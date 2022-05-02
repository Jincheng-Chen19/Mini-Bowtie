import numpy as np
from Index import Index
import simplesam
from simplesam import Writer, OrderedDict
from tqdm import tqdm


class Align:
    def __init__(self):
        """
        The class for alignment

        :attribute index: get index for alignment
        :attribute message: record the message for each read
        :attribute reads: record sequence for each read
        :attribute quality: record quality for each read
        :attribute seedF: consider the reads are from the reference sequence and extract the seeds,
        :attribute seedR: consider the reads are from the reverse complement sequence and extract the seeds
        :attribute start: after filtrating and save the starting position of the extended ref sequence
        :attribute EC:  after filtrating and save the sequence of the extended ref sequence
        :attribute flag: a column in sam file, recorded the match situation
                         0:  match to forward reference strand. 4 no any match. 16 match to reverse reference strand
        :attribute startpos: save the actual starting position for each reads
        :attribute mark: save the mark got from SW matrix for each reads
        :attribute CIGAR: for the sam file, record the number of match, mismatch, insertion, and deletion
        """
        self.index = Index()
        self.message = []
        self.reads = []
        self.quality = []
        self.seedF = []
        self.seedR = []
        self.start = []
        self.EC = []
        self.flag = []
        self.startpos = []
        self.mark = []
        self.CIGAR = []

    def setIndex(self, file):
        """
        Set the index based on the class Index

        :param file: the path to the FASTA file containing genomic sequence
        """
        self.index.BWT(file)
        self.index.setcp()

    def readFq(self, file, reads: int):
        """
        Get the message, sequence, and quality of the reads from fastq file

        :param file: the path to the fastq file
        :param reads: optional. default is None, which will read in all the reads.
                      if give this parameter, it will only read in and align this number of reads at started.
        """
        fq = open(file, "r")
        all = []
        for line in fq:
            all.append(line[:-1])  # read in all the messages
        # default, read in all the reads.
        if reads is None:
            for i in range(0, len(all), 4):
                self.message.append(all[i])
                self.reads.append(all[i + 1])
                self.quality.append(all[i + 3])
                self.flag.append([])
        # according to the parameter "reads" to read in
        else:
            for i in range(0, 4 * reads, 4):
                self.message.append(all[i])
                self.reads.append(all[i+1])
                self.quality.append(all[i+3])
                self.flag.append([])

    def getseed(self, read: str, seed: int, gap: int):
        """
        Get all the seeds for one read.

        :param read: the sequence of this read
        :param seed: the length of each seed, default 16.
        :param gap: the gap between each seed, default 10.
        :return: FL: the seed when considering the read is aligned to the reference sequence
                 RL：the seed when considering the read is aligned to the reverse complement reference sequence
        """
        FL, RL = [], []
        # get the reverse complement sequence for the read
        rr = read[::-1].lower()
        rr = rr.replace("a", "T")
        rr = rr.replace("c", "G")
        rr = rr.replace("g", "C")
        rr = rr.replace("t", "A")
        # according to the parameter to get seeds
        for i in range(0, len(read), gap):
            if i + seed < len(read):
                FL.append(read[i: i + seed])
                RL.append(rr[i: i + seed])
        return FL, RL

    def seeds(self, seed: int, gap: int):
        """
        With the function getseed, get and save all the seeds for all reads.

        :param seed: the length of each seed, default 16.
        :param gap: the gap between each seed, default 10.
        """
        for i in self.reads:
            self.seedF.append(self.getseed(i, seed, gap)[0])
            self.seedR.append(self.getseed(i, seed, gap)[1])

    def match(self, seed: list):
        """
        Get the BW ranges for a list of seeds

        :param seed: a list of seeds
        :return: a list containing all the BW ranges for these seeds
        """
        match = []
        for i in range(0, len(seed)):
            match.append(self.index.match(seed[i]))
        return match

    def findpos(self, seed: list):
        """
        Accordign to the BW ranges, find the positions in the reference genome as extension candidates.

        :param seed: a list of seeds
        :return: a list containing all the positions for each seed
        """
        match = self.match(seed)
        pos = []
        for i in range(0, len(match)):
            p = []
            # when a seed has the BW range, use function in class Index to find the positions
            if match[i][0] is not None:
                for j in range(match[i][0], match[i][1]):
                    p.append(self.index.findpos(j))
                pos.append(p)
        return pos

    def extension(self, seed: list, length: int, gap: int):
        """
        Based on all these positions, filtrate to get proper positions,
        and then extend these position in the reference sequence to make the extended sequence.

        :param seed: a list of seeds
        :param length: the length of the read that these seeds were extracted from it
        :param gap: the gap between each seed, default 10.
        :return: ext: all the extended reference sequence,
                 startpos: the starting position of the extension reference sequence
        """
        pos = self.findpos(seed)
        ext = []
        startpos = []
        start = []
        # based on the matched positions of each seed, get the probable starting positions for reads.
        for i in range(0, len(pos)):
            if pos[i]:
                for j in pos[i]:
                    s = j - i * gap
                    start.append(s)
        # Filtrate the starting positions.
        startunique = []
        while start:
            s = start.pop(0)
            add = False
            # considering the mismatch, 10bp around one position is considered as the same starting position.
            for i in range(s - 5, s + 5):
                while i in start:
                    start.remove(i)
                    add = True
            if add:
                startunique.append(s)
        # after filtrating, extend the position to get extenend reference sequence and its starting position
        for start in startunique:
            end = start + length
            e, s = self.index.extension(start, end)
            # make sure the extenend reference sequence is longer than reads.
            if len(e) < length:
                continue
            ext.append(e)
            startpos.append(s)
        return ext, startpos

    def allextension(self, gap: int):
        """
        Use function extension to get extended reference sequence and starting position for all the reads.

        :param gap: the gap between each seed, default 10.
        """
        for i in tqdm(range(0, len(self.seedF))): # tqdm package learnt from https://www.jb51.net/article/203570.htm
            exF, startF = self.extension(self.seedF[i], len(self.reads[i]), gap)
            exR, startR = self.extension(self.seedR[i], len(self.reads[i]), gap)
            self.EC.append(exF)
            self.EC.append(exR)
            self.start.append(startF)
            self.start.append(startR)

    def SWmark(self, ref: str, seq: str, match: float, mismatch: float, indel: float):
        """
        Make the Smith-Waterman mark matrix for extended reference sequence and sequence of the read.

        :param ref: the extended reference sequence
        :param seq: the sequence of the read
        :param match: how many score will add if match
        :param mismatch: how many score will reduce if mismatch
        :param indel: how many score will reduce if insertion and if deletion
        :return: the whole mark matrix will be return
        """
        ref = " " + ref
        seq = " " + seq
        x, y = len(ref), len(seq)
        ScoreMat = np.zeros(shape=(x, y))
        for i in range(1, x):
            for j in range(i, y):
                m = int
                # if match
                if ref[i] == seq[j]:
                    m = ScoreMat[i - 1, j - 1] + match
                # if mismatch
                if ref[i] != seq[j]:
                    m = ScoreMat[i - 1, j - 1] - mismatch
                xgap = ScoreMat[i - 1, j] - indel  # for insertion
                ygap = ScoreMat[i, j - 1] - indel  # for deletion
                ScoreMat[i, j] = max(0, m, xgap, ygap)  # find the highest score
        return ScoreMat

    def traceback(self, ScoreMat, max_x, max_y,  ref: str, seq: str, match: float, mismatch: float, indel: float):
        """
        with the mark matrix, and the position of the max, traceback to get the match situation,
        return the start position and the match situation.

        :param ScoreMat: the SW mark matrix
        :param max_x: the row of the max
        :param max_y: the column of the max
        :param ref: the extended reference sequence
        :param seq: the sequence of the read
        :param match: how many score will add if match
        :param mismatch: how many score will reduce if mismatch
        :param indel: how many score will reduce if insertion and if deletion
        :return: CIGAR: record the match situation, four number refer to the number of match, mismatch,
                        insertion, and deletion, respectively
                 x: refer to after traceback, the start position of the match.
        """
        ref = " " + ref
        seq = " " + seq
        S = np.max(ScoreMat)
        x = max_x
        y = max_y
        CIGAR = [0, 0, 0, 0]
        while S > 0:
            # for the match situation
            if ScoreMat[x, y] == ScoreMat[x - 1, y - 1] + match and ref[x] == seq[y]:
                x -= 1
                y -= 1
                S = ScoreMat[x, y]
                CIGAR[0] += 1
            # for the mismatch situation
            elif ScoreMat[x, y] == ScoreMat[x - 1, y - 1] - mismatch and ref[x] != seq[y]:
                x -= 1
                y -= 1
                S = ScoreMat[x, y]
                CIGAR[1] += 1
            # for the insertion situation
            elif ScoreMat[x, y] == ScoreMat[x - 1, y] - indel:
                x -= 1
                S = ScoreMat[x, y]
                CIGAR[2] += 1
            # for the deletion situation
            elif ScoreMat[x, y] == ScoreMat[x, y-1] - indel:
                y -= 1
                S = ScoreMat[x, y]
                CIGAR[3] += 1
        return CIGAR, x

    def markres(self, lowerlimit: float, match=1, mismatch=-1/3, indel=1):
        """
        Use Smith-Waterman algorithm to calculate the mark of the matches,
        and further find the best matches for all the reads.

        :param lowerlimit: the lowest mark that was considered as matched
        :param match: how many score will add if match
        :param mismatch: how many score will reduce if mismatch
        :param indel: how many score will reduce if insertion and if deletion
        """
        for i in tqdm(range(0, len(self.reads))):
            flag = []
            matrix = []
            startpos = []
            markF = []
            markR = []
            CIGAR = []
            maxF, maxR = 0, 0
            rr = ""
            # for the extended sequence from foward seeds
            if self.EC[2 * i]:
                # get the highest mark for forward sequence
                for j in range(0, len(self.EC[2 * i])):
                    ScoreMat = self.SWmark(self.EC[2 * i][j], self.reads[i], match, mismatch, indel)
                    markF.append(np.max(ScoreMat))
                    matrix.append(ScoreMat)
                maxF = max(markF)
            # for the extended sequence from reverse complement seeds
            if self.EC[2 * i + 1]:
                # reverse complement sequence of the read
                rr = self.reads[i][::-1].lower()
                rr = rr.replace("a", "T")
                rr = rr.replace("c", "G")
                rr = rr.replace("g", "C")
                rr = rr.replace("t", "A")
                # get the highest mark for reversed sequence
                for j in range(0, len(self.EC[2 * i + 1])):
                    ScoreMat = self.SWmark(self.EC[2 * i + 1][j], rr, match, mismatch, indel)
                    markR.append(np.max(ScoreMat))
                    matrix.append(ScoreMat)
                maxR = max(markR)
            # if all the max mark lower than the lowest score
            if maxF < lowerlimit and maxR < lowerlimit:
                self.mark.append(0)
                self.startpos.append([])
                self.flag[i] = [4]
                self.CIGAR.append([[0, 0, 0, 0]])
                continue
            # if forward matches have higher score
            if maxF > maxR:
                F = np.where(markF == np.max(markF))[0]  # find the highest scores
                self.mark.append(np.max(markF))
                for j in F:
                    ScoreMat = matrix[j]
                    m = np.where(ScoreMat == np.max(ScoreMat))  # get the position in SW matrix for highest scores
                    for k in range(0, len(m[0])):
                        # traceback to get the match situation and the starting position in extended strain
                        c, s = self.traceback(ScoreMat, m[0][k], m[1][k], self.EC[2 * i][j], self.reads[i], match, mismatch, indel)
                        flag.append(0)
                        # add the starting position of the extended sequence and the starting position of the match together
                        # to get the real start position
                        startpos.append(self.start[2 * i][j] + s + 1)
                        CIGAR.append(c)
            # if reverse complement matches have higher score
            if maxR > maxF:
                R = np.where(markR == np.max(markR))[0]  # find the highest scores
                self.mark.append(np.max(markR))
                for j in R:
                    ScoreMat = matrix[len(markF) + j]
                    m = np.where(ScoreMat == np.max(ScoreMat))  # get the position in SW matrix for highest scores
                    for k in range(0, len(m[0])):
                        # traceback to get the match situation and the starting position
                        c, s = self.traceback(ScoreMat, m[0][k], m[1][k], self.EC[2 * i + 1][j], rr, match, mismatch, indel)
                        flag.append(16)
                        startpos.append(self.start[2 * i + 1][j] + s + 1)
                        CIGAR.append(c)
            # both situation have highest score
            if maxF == maxR:
                F = np.where(markF == np.max(markF))[0]
                self.mark.append(np.max(markF))
                for j in F:
                    ScoreMat = matrix[j]
                    m = np.where(ScoreMat == np.max(ScoreMat))
                    for k in range(0, len(m[0])):
                        c, s = self.traceback(ScoreMat, m[0][k], m[1][k], self.EC[2 * i + 1][j], self.reads[i], match, mismatch, indel)
                        flag.append(0)
                        startpos.append(self.start[2 * i][j] + s + 1)
                        CIGAR.append(c)
                R = np.where(markR == np.max(markR))[0]
                self.mark.append(np.max(markR))
                for j in R:
                    ScoreMat = matrix[len(markF) + j]
                    m = np.where(ScoreMat == np.max(ScoreMat))
                    for k in range(0, len(m[0])):
                        c, s = self.traceback(ScoreMat, m[0][k], m[1][k], self.EC[2 * i + 1][j], rr, match, mismatch, indel)
                        flag.append(16)
                        startpos.append(self.start[2 * i + 1][j] + s + 1)
                        CIGAR.append(c)
            startpos = list(set(startpos))
            self.startpos.append(startpos)
            self.flag[i] = flag
            self.CIGAR.append(CIGAR)

    def write(self, samfile, name, sequence, flag, startpos, mapping_quality, quality, cigar):
        """
        Using the package simplesam to write for every read into a .sam file

        :param samfile: the name of the Writer
        :param name: the name of the read
        :param sequence: the sequence of the read
        :param flag: the flag number
        :param startpos: the start position of match in reference genome
        :param mapping_quality: the mapping_quality got by  SW matrix
        :param quality: the quality given by fastq file
        :param cigar: mapping situation
        :return:
        """
        line = simplesam.Sam()
        line.qname = name[1:]
        line.flag = flag
        if flag == 0:
            line.seq = sequence
        # if match to reverse reference strand, write in the reverse complement sequence of the reads
        if flag == 16:
            sequence = sequence[::-1].lower()
            sequence = sequence.replace("a", "T")
            sequence = sequence.replace("c", "G")
            sequence = sequence.replace("g", "C")
            sequence = sequence.replace("t", "A")
            line.seq = sequence
        line.rname = "chr1"
        line.mapq = mapping_quality
        line.pos = startpos
        line.qual = quality
        # write the CIGAR in right format
        c = [0, 0, 0, 0, 0]
        c[0] = cigar[0] + cigar[1]
        c[1] = cigar[3]
        c[2] = cigar[2]
        c[3] = cigar[0]
        c[4] = cigar[1]
        line.cigar = str(c[0]) + "M" + str(c[1]) + "I" + str(c[2]) + "D" + str(c[3]) + "=" + str(c[4]) + "X"
        samfile.write(line)

    def writesamfile(self, outputfile):
        """
        Write a samfile for the whole alignment

        :param outputfile: the path of the output file
        """
        out_file = open(outputfile, 'w')
        samfile = Writer(out_file)
        for i in tqdm(range(0, len(self.reads))):
            for j in range(0, len(self.startpos[i])):
                self.write(samfile, self.message[i], self.reads[i], self.flag[i][j], self.startpos[i][j], self.mark[i], self.quality[i], self.CIGAR[i][j])
        samfile.close()
