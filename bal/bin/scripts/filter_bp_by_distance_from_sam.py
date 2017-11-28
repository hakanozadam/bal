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
from collections import defaultdict

import pysam

#####################################################################

def handle_arguments():
   parser = argparse.ArgumentParser(description=
   '''
   only picks those reads mapping to a branchpoint whose distance is less than or equal to d
   The input and output files are both in sam format.
   ''')
   parser.add_argument("-i" ,
                       help = "Input sam file" ,
                       required = True ,
                       metavar = "input_sam_file" ,
                       type = str)
   parser.add_argument("-d" ,
                       help = "Distance Threshold" ,
                       required = True ,
                       metavar = "radius" ,
                       type = int)
   parser.add_argument("-o" ,
                       help = "Output file" ,
                       required = True ,
                       metavar = "output_fastq_file" ,
                       type = str)

   arguments = parser.parse_args()
   return arguments

#####################################################################

def main():
    arguments = handle_arguments()
    distance_index = 5

    input_file     = pysam.AlignmentFile(arguments.i, 'r')
    output_sam_stream = pysam.AlignmentFile(arguments.o, "wh", header = input_file.header)

    for read in input_file.fetch():
        aligned_ref_label    = input_file.getrname(read.reference_id)
        ref_contents = aligned_ref_label.split("__")
        if int(ref_contents[distance_index]) <= int(arguments.d):
            output_sam_stream.write(read)

    input_file.close()
    output_sam_stream.close()


#####################################################################ß

if __name__ == '__main__':
    main()
else:
    exit(1)

####################################################################
