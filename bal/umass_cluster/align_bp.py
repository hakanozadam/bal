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
####################################################################

import argparse
import os
import operator
import sys
import prepare

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':

    sys.path.append( bal_dir  )
    from bal.core.bowtie2 import Bowtie2
    from bal.settings import *
else:
    exit(1)

####################################################################


####################################################################


####################################################################

def main():

   executables = prepare.get_executables(bal_dir)

   parser = argparse.ArgumentParser(description=
   '''
   Align the given fastq file(s) against the provided bowtie2 refrence
   ''')
   parser.add_argument("-x" ,
                       help = "Bowtie2 reference base" ,
                       required = True ,
                       metavar = "genome_fasta_file" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_directory" ,
                       type = str)
   parser.add_argument("-a" ,
                       help = "bowtie2 arguments" ,
                       required = False ,
                       default = " ",
                       metavar = "bowtie2_Arguments" ,
                       type = str)
   parser.add_argument("-p" ,
                       help = "bowtie2 threads" ,
                       required = False ,
                       default = 1,
                       metavar = "bowtie2_threads" ,
                       type = int)
   parser.add_argument("-U" ,
                       help = "Input Fastq File" ,
                       required = True ,
                       metavar = "input_fastq_file" ,
                       type = str)
   arguments = parser.parse_args()

   bt2_job = Bowtie2(name="bt2_bp", input_files=[],
                     output_directory = os.path.abspath(arguments.o),
                     executable = executables['bowtie2'] ,
                     executable_arguments = settings["align_bp_bt2_arguments"] + " " + arguments.a ,\
                     genome_reference = arguments.x ,
                     single_end_fastq = arguments.U,
                     paired_end_fastqs = list(),
                     p = arguments.p)
   bt2_job.run()

#####################################################################

if __name__ == '__main__':
    main()