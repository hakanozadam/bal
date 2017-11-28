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
    from bal.core.bp_bed_from_sam import BpBedFromSam
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

   parser.add_argument("-1" ,
                    metavar = 'mate_1_file' ,
                    help = "First sam file. Must be provided" ,
                    required = True ,
                    dest = 'mate_1',
                    type = str)
   parser.add_argument("-2" ,
                    metavar = 'mate_2_file' ,
                    help = "If the data are paired end, mate_2 file is provided in this argument" ,
                    required = False ,
                    dest = 'mate_2' ,
                    type = str)
   parser.add_argument("-o" ,
                    metavar = 'output_bed' ,
                    help = "Output bed file " ,
                    required = False ,
                    dest = 'output' ,
                    type = str)

   arguments = parser.parse_args()

   input_files = [arguments.mate_1]
   if arguments.mate_2:
       input_files.append(arguments.mate_2)

   this_job = BpBedFromSam(name             = 'Bed_From_Sam_Ref',
                           input_files      = input_files,
                           output_directory = arguments.output)

   this_job.run()

#####################################################################

if __name__ == '__main__':
    main()