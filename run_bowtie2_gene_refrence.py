#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#
#        Moore Laboratory
#        UMASS Medical School / HHMI
#        RNA Therapeutics Institute
#        Albert Sherman Center, ASC4-1009
#        368 Plantation Street
#        Worcester, MA 01605
#        USA
#
#################################################################

import os, argparse
from bal.reference.bowtie2_gene_reference import Bowtie2GeneReference
from bal.settings import *

from bal.reference.prepare import get_executables

###################################################################

###################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Makes Bowtie2 Reference for eac fasta file in the input directory
   ''')
   parser.add_argument("-i" ,
                       help = "Input Directory. Contains fasta files" ,
                       required = True ,
                       metavar = "input_directory" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_directory",
                       type = str)
   arguments = parser.parse_args()

   executables = get_executables(os.path.dirname(os.path.realpath(__file__)))

   this_job = Bowtie2GeneReference(
                 name             = "Bowtie2_Gene_Reference",
                 input_files      = [ arguments.i] ,
                 output_directory = arguments.o ,
                 executable       = executables['bowtie2-build']
                  )

   this_job.run()

############################################################################

if __name__ == "__main__":
    main()