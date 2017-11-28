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
   This script filters out alignments based on a minumum mapping quality score,
   alignment score, read length etc.
   ''')
   parser.add_argument("-i" ,
                       help = "Input sam file" ,
                       required = True ,
                       metavar = "input_sam_file" ,
                       type = str)
   parser.add_argument("--dist-min" ,
                       help = "Minimum Distance Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "dist_min" ,
                       dest = "dist_min",
                       type = int)
   parser.add_argument("--dist-max" ,
                       help = "Maximum Distance Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "dist_max" ,
                       dest = "dist_max",
                       type = int)
   parser.add_argument("--mapq-min" ,
                       help = "Minimum Mapping Quality Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "mapq_min" ,
                       dest = "mapq_min",
                       type = int)
   parser.add_argument("--mapq-max" ,
                       help = "Maximum Mapping Quality Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "mapq_max" ,
                       dest = "mapq_max",
                       type = int)
   parser.add_argument("--length-min" ,
                       help = "Minimum Length Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "length_min" ,
                       dest = "length_min",
                       type = int)
   parser.add_argument("--length-max" ,
                       help = "Maximum Length Threshold. Not implemented yet" ,
                       required = False ,
                       metavar = "length_max" ,
                       dest = "length_max",
                       type = int)
   parser.add_argument("-o" ,
                       help = "Output file" ,
                       required = True ,
                       metavar = "output_sam_file" ,
                       type = str)

   arguments = parser.parse_args()
   return arguments

#####################################################################

def main():
    arguments = handle_arguments()
    distance_index = 5

    base_attributes = ['dist', 'length', 'mapq']
    min_attributes  = map( lambda x: x + '_min', base_attributes )
    max_attributes  = map( lambda x: x + '_max', base_attributes )
    min_filters = dict()

    input_file     = pysam.AlignmentFile(arguments.i, 'r')
    output_sam_stream = pysam.AlignmentFile(arguments.o, "wh", header = input_file.header)

    for read in input_file.fetch():
        aligned_ref_label    = input_file.getrname(read.reference_id)
        ref_contents = aligned_ref_label.split("__")
        this_distance = int(ref_contents[distance_index])
        #print(dir(read))
        if arguments.mapq_min != None and read.mapping_quality < arguments.mapq_min:
            continue
        if arguments.mapq_max != None and read.mapping_quality > arguments.mapq_max:
            continue
        if arguments.dist_min != None and this_distance < arguments.dist_min:
            continue
        if arguments.dist_max != None and this_distance > arguments.dist_max:
            continue
        if arguments.length_min != None and read.inferred_length < arguments.length_min:
            continue
        if arguments.length_max != None and read.inferred_length > arguments.length_max:
            continue
        output_sam_stream.write(read)

    input_file.close()
    output_sam_stream.close()


#####################################################################ß

if __name__ == '__main__':
    main()
else:
    exit(1)

####################################################################
