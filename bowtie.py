import numpy as np
from Index import Index
from align import Align
import simplesam
from simplesam import Writer, OrderedDict
import sys
from optparse import OptionParser
from tqdm import tqdm


def Bowtie(reads, FASTA, FASTQ, outputfile, seed=16, gap=10, match=1, mismatch=1 / 3, indel=1, lowerlimit=50):
    """
    The final function use class Align to achieve the alignment.

    :param reads: optional. default is None, which will read in all the reads.
                  if give this parameter, it will only read in and align this number of reads at started.
    :param FASTA: the path to the fasta file
    :param FASTQ: the path to the fastq file
    :param outputfile: the path to the output file
    :param seed: the length of each seed, default 16.
    :param gap: the gap between each seed, default 10.
    :param match: how many score will add if match
    :param mismatch: how many score will reduce if mismatch
    :param indel: how many score will reduce if insertion and if deletion
    :param lowerlimit: the lowest mark that was considered as matched
    """
    bowtie = Align()
    bowtie.setIndex(FASTA)
    bowtie.readFq(FASTQ, reads)
    bowtie.seeds(seed, gap)
    bowtie.allextension(gap)
    bowtie.markres(lowerlimit, match, mismatch, indel)
    bowtie.writesamfile(outputfile)


# use OptionParser to make the code command-line executable
parser = OptionParser()
parser.add_option("-r", "--reads", action="store", type="int", dest="reads", default=None,
                  help="align the number of reads at started")
parser.add_option("-R", "--FASTA", action="store", type="string", dest="FASTA", help="the path to the FASTA file")
parser.add_option("-U", "--FASTQ", action="store", type="string", dest="FASTQ", help="the path to the FASTQ file")
parser.add_option("-o", "--outputfile", action="store", type="string", dest="outputfile",
                  help="the path to the output file")
parser.add_option("-s", "--seed", action="store", type="int", dest="seed", default=16, help="the length of each seed")
parser.add_option("-g", "--gap", action="store", type="int", dest="gap", default=10, help="the gap between seeds")
parser.add_option("-m", "--match", action="store", type="float", dest="match", default=1,
                  help="add this score if match")
parser.add_option("-n", "--mismatch", action="store", type="float", dest="mismatch", default=1 / 3,
                  help="reduce this score if mismatch")
parser.add_option("-i", "--indel", action="store", type="float", dest="indel", default=1,
                  help="reduce this score if insertion and if deletion")
parser.add_option("-l", "--lowerlimit", action="store", type="float", dest="lowerlimit", default=20,
                  help="the lowest score considered as match")
options, args = parser.parse_args()

# The code to achieve the alignment
Bowtie(reads=options.reads, FASTA=options.FASTA, FASTQ=options.FASTQ, outputfile=options.outputfile, seed=options.seed,
       gap=options.gap, match=options.match, mismatch=options.mismatch, indel=options.indel,
       lowerlimit=options.lowerlimit)
