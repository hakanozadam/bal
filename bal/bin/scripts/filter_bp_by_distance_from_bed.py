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

#####################################################################

def handle_arguments():
   parser = argparse.ArgumentParser(description=
   '''
   Given a bed file of branchpoints, it filters out those
   ''')
   parser.add_argument("-i" ,
                       help = "Input bed file" ,
                       required = True ,
                       metavar = "input_bed_file" ,
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
   parser.add_argument("-g" ,
                       help = "Output distances strictly GREATER than d (defualt is False)."
                              "That is, in the absence of this option, it reports BP's"
                              "with distance d or less." ,
                       action = "store_true",
                       dest = "greater_than")

   arguments = parser.parse_args()
   return arguments

#####################################################################

def main():
    arguments = handle_arguments()
    distance_index = 5

    with open(arguments.i, 'r') as input_stream,\
        open(arguments.o, 'w') as output_stream:

        if arguments.greater_than:
            print("Reporting branchpoints with distance strictly greater than " + str( arguments.d ))
            for line in input_stream:
                contents = line.strip().split('__')
                if len(contents) < 6:
                    continue
                if int(contents[distance_index]) > int(arguments.d):
                    print(line.strip(), file = output_stream)
        else:
            print("Reporting branchpoints with distance less than " + str( arguments.d ))
            for line in input_stream:
                contents = line.strip().split('__')
                if len(contents) < 6:
                    continue
                if int(contents[distance_index]) <= int(arguments.d):
                    print(line.strip(), file = output_stream)

#####################################################################ß

if __name__ == '__main__':
    main()
else:
    exit(1)

####################################################################
