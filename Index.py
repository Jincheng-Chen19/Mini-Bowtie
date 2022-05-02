import numpy as np
from tqdm import tqdm


class Index:
    def __init__(self):
        """
        the class for index
        idea partially from https://blog.csdn.net/stormlovetao/article/details/7048481

        :attribute message: The basic information of the reference
        :attribute B: save the last column of the BWT
        :attribute S: save the order before sorted for every 36 rows
        :attribute C: save the number of the nucleotides smaller than current nucleotide
        :attribute start: save the line whose previous order is 0
        :attribute cp: save the checkpoint for every 448 rows
        """
        self.message = ""
        self.B = ""
        self.S = []
        self.C = {}
        self.start = int
        self.cp = [[0], [0], [0], [0]]

    def readFa(self, file):
        """
        Read in the fasta file, and return the genomic sequence.

        :param file: the path to the fasta file
        :return: genomic sequence with a "$" at the end.
        """
        fa = open(file, "r")
        seq = ""
        for line in fa:
            if line.startswith(">"):
                self.message += line[:-1]
            else:
                seq += line[:-1]
        seq += "$"  # add the special symbol at the end of the sequence
        return seq

    def BWT(self, file):
        """
        Do the Burrows-Wheeler transform, and save the last column, the order for every 36 rows,
        and the number of the nucleotides smaller than the current nucleotide.

        :param file: the path to the fasta file
        """
        sequence = self.readFa(file)  # Get the sequence
        size = len(sequence)
        # sort the BW matrix by the forward column. This BW matrix is left rotated.
        # this code cited https://stackoverflow.com/questions/59913437/create-a-bowtie-pattern-with-python-3-0
        sort = sorted(range(size), key=lambda i: sequence[i:])
        self.start = sort.index(0)  # record the row for the start row before rotation
        for i in range(0, len(sort)):
            self.B += sequence[sort[i] - 1]  # record the last column
            if i % 36 == 0:
                self.S.append(sort[i])  # record the order every 36 rows
        self.C["A"] = 0
        self.C["C"] = self.B.count("A")
        self.C["G"] = self.B.count("C") + self.C["C"]
        self.C["T"] = self.B.count("G") + self.C["G"]

    def makecp(self, bp: str, pos: int):
        """
        Calculate the checkpoint.
        The checkpoint is that before this row, how many current nucleotides in the last column.

        :param bp: the current nucleotide.
        :param pos: count the number before this row.
        :return: the number of current nucleotides before pos.
        """
        n = "ACGT"
        l = len(self.cp[3]) - 1
        c = self.cp[n.find(bp)][l]
        for i in range(l * 448, pos):
            if self.B[i] == bp:
                c += 1
        return c

    def setcp(self):
        """
        Save the checkpoint for every 448 rows.
        """
        for i in range(448, len(self.B), 448):
            self.cp[0].append(self.makecp("A", i))
            self.cp[1].append(self.makecp("C", i))
            self.cp[2].append(self.makecp("G", i))
            self.cp[3].append(self.makecp("T", i))

    def Occ(self, bp: str, pos: int):
        """
        Count before this row, how many current nucleotides in the last column with the help of checkpoint list.

        :param bp: the current nucleotide.
        :param pos: count the number before this row.
        :return: the count.
        """
        n = "ACGT"
        point = pos // 448
        c = self.cp[n.find(bp)][point]
        for i in range(point * 448, pos):
            if self.B[i] == bp:
                c += 1
        return c

    def LF(self, pos):
        """
        find the position of the previous nucleotide in last column.

        :param pos: the current row
        :return: the row for the previous nucleotide
        """
        bp = self.B[pos]
        return self.C[bp] + self.Occ(bp, pos) + 1

    def deBWT(self):
        """
        Use the LF function to recover the genomic sequence

        :return: the genomic sequence
        """
        pos = 0  # From the last nucleotide
        ref = ""
        while self.B[pos] != "$":
            ref = self.B[pos] + ref # insert the
            pos = self.LF(pos)  # get the position for the previous nucleotide
        return ref

    def match(self, seed):
        """
        Find the BW ranges of the seed.

        :param seed: extracted from the read
        :return: the BW range of the seed
        """
        n = "ACGT"
        bp = seed[-1]
        spos = self.C[bp] + 1  # lower limit
        if bp == "T":
            epos = len(self.B)
        else:
            bp = n[n.index(bp) + 1]
            epos = self.C[bp] + 1  # upper limit
        for i in range(1, len(seed)):
            bp = seed[-(i + 1)]
            spos = self.C[bp] + self.Occ(bp, spos) + 1
            epos = self.C[bp] + self.Occ(bp, epos) + 1
            if spos >= epos or epos > len(self.B):
                return None, None
        return spos, epos

    def findpos(self, pos: int):
        """
        Find the position in the genomic sequence according to the row in BW matrix.

        :param pos: the row in BW matrix.
        :return: the position in the genomic sequence
        """
        t = 0
        while pos % 36 != 0:  # save the position every 36 rows
            if pos == self.start:
                return 0 + t  # if get the starting row, return the starting position
            pos = self.LF(pos)  # find the previous nucleotide
            t += 1  # record how many steps go previous and plus this to get the actual match position
        point = pos // 36  # get the index
        return self.S[point] + t

    def extension(self, start: int, end: int):
        """
        Based on the starting and the end position, return the reference sequence between it, and the actual starting position.

        :param start: starting position in the genomic sequence
        :param end: end position in the genomic sequence
        :return: the sequence and the actual starting position
        """
        end += 10  # more nucleotides at end for SW, for deletion situation
        # find the proper end position and get the end row
        while end not in self.S:
            if end > len(self.B) - 1:
                end = len(self.B) - 1
                break
            end += 1
        pos = self.S.index(end) * 36
        start -= 10  # more nucleotides at start for SW, for deletion situation
        # find the proper starting position and get the starting row
        while start not in self.S:
            if start <= 0:
                start = 0
                break
            start -= 1
        if start == 0:
            poss = self.start
        else:
            poss = self.S.index(start) * 36
        range = ""
        #  recover the sequence from the end row to the starting row to get the extended reference sequence
        while not pos == poss:
            range = self.B[pos] + range
            pos = self.LF(pos)
        return range, start

