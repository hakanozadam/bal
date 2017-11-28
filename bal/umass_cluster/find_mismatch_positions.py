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
import sys
import pysam

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':
    sys.path.append( bal_dir  )
    from bal.core.find_mismatch_positions import FindMismatchPositions
    from bal.settings import *
else:
    exit(1)

####################################################################

####################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Extract the branchpoints in bed format from the given sam file(s).
   The sam file reference must be formatted properly for this.
   The input sam file comes from aligning reads to the reference branchpoints.
   ''')

   parser.add_argument("-f" ,
                    metavar = 'fasta_file' ,
                    help = "Reference sequence fasta file" ,
                    required = True ,
                    type = str)
   parser.add_argument("-s" ,
                    metavar = 'sam_file' ,
                    help = "Alignment file in sam format" ,
                    required = True ,
                    type = str)
   parser.add_argument("-o" ,
                    metavar = 'output_directory' ,
                    help = "Output directory" ,
                    required = True ,
                    type = str)

   arguments = parser.parse_args()

   this_job = FindMismatchPositions(name             = 'Find_Mismatch_Positions',
                                    fasta_file       = arguments.f,
                                    sam_file         = arguments.s,
                                    output_directory = arguments.o)

   this_job.run()

#####################################################################

if __name__ == '__main__':
    main()