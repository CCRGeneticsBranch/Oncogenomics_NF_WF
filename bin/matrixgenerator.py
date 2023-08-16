#!/usr/bin/python

import sys


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import os
import numpy

#with open(sys.argv[1], "r") as infile:
#    in_name = os.path.basename(infile.name)

#abso = os.path.abspath(str(in_name))
#inpath = os.path.dirname(abso)



#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', rsync=False, bash=True)




matrices = matGen.SigProfilerMatrixGeneratorFunc(sys.argv[1], "GRCh37", sys.argv[2], plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=True, seqInfo=True, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc(in_name, "GRCh37", str(inpath), plot=True, exome=True, bed_file=None, chrom_based=False, tsb_stat=True, seqInfo=True, cushion=100)


