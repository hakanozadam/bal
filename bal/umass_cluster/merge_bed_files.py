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

####################################################################


####################################################################

def merge_bed_files(input_list, output_bed):
    output_bed_entries = dict()
    for file in input_list:
        with open(file, 'r') as input_stream:
            for line in input_stream:
                contents = line.strip().split()
                contents[4] = int(contents[4])
                print(contents)
                if output_bed_entries.get(contents[3]):
                    output_bed_entries[contents[3]][4] = \
                        output_bed_entries[contents[3]][4] +\
                        contents[4]
                else:
                    output_bed_entries[contents[3]] = contents

    sorted_entries = sorted(output_bed_entries, key = operator.itemgetter(0))
    print('---sorted_entries---\n', sorted_entries)
    print('---output_entries---\n', output_bed_entries)

    with open(output_bed, 'w') as output_stream:
        for entry in sorted_entries:
           output_bed_entries[entry][4] = str(output_bed_entries[entry][4])
           print("\t".join(output_bed_entries[entry]) , file = output_stream)

####################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Merge the bed files given as a list in the input file and write the result
   to the output given in the -o argument. The input list should be one file per line.
   ''')
   parser.add_argument("-i" ,
                       help = "Input bed file list One file per line." ,
                       required = True ,
                       metavar = "input_file_list" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Merged bed file" ,
                       required = True ,
                       metavar = "output_file" ,
                       type = str)
   arguments = parser.parse_args()

   input_file_list = list()

   with open(arguments.i, 'r') as input_list_stream:
       for line in input_list_stream:
           file_path = line.strip()
           if not os.path.isfile(file_path):
               raise(FileNotFoundError("The bed file\n", file_path,\
                                       "doesn't exist!") )
           input_file_list.append(file_path)

   merge_bed_files(input_file_list, arguments.o)


#####################################################################

if __name__ == '__main__':
    main()