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
import re

####################################################################

####################################################################

def get_file_name(file_path):
    file_name, extension = "", ""
    file_base = os.path.basename(file_path)
    paired_case = re.compile(r"(?P<head>.*)(?P<tail>(([._][1,2])\.f(ast)?q))$" , flags = re.IGNORECASE)
    single_case = re.compile(r"(?P<head>.*)(?P<tail>(\.f(ast)?q))$" , flags = re.IGNORECASE)
    paired_result = paired_case.search(file_base)
    single_result = single_case.search(file_base)
    if paired_result:
        file_name, extension = paired_result.group("head"), paired_result.group("tail")
    elif single_result:
        file_name, extension = single_result.group("head"), single_result.group("tail")
    else:
        print("Problem with the file name ", file_path)
        exit(1)
    return [file_name, extension]

####################################################################

def main():
    ''' Partition the given fastq file into little pieces'''

    parser = argparse.ArgumentParser(description=
    '''
    BAL Pipeline For Umass Cluster

    Partition the given fastq file into smaller pieces.
    The number of reads per piece is defined by the user.
    ''')

    parser.add_argument("-i" ,
                    metavar  = 'input_fastq_file' ,
                    help     = "Input fastq file." ,
                    required = True ,
                    type     = str)
    parser.add_argument("-n" ,
                    metavar  = 'reads_per_file' ,
                    help     = "Number of reads per file" ,
                    required = True ,
                    type     = int)
    parser.add_argument("-o" ,
                    metavar  = 'Output_directory' ,
                    help     = "Output directory" ,
                    required = True ,
                    type     = str)
    arguments = parser.parse_args()

    if arguments.n < 1:
        arguments.n = 2

    if not os.path.isfile(arguments.i):
        print("Couldn't find the input file", arguments.i)

    os.makedirs(arguments.o, exist_ok=True)

    file_counter = 1
    read_counter = 0

    file_name, extension = get_file_name(arguments.i)

    current_output_file = os.path.join(arguments.o, file_name + "-" + str(file_counter) + extension)
    current_output_stream = open(current_output_file, "w")

    with open(arguments.i, "r") as input_stream:
        line_mode = 0
        for line in input_stream:
            if line_mode == 0:
                if read_counter >= arguments.n:
                    current_output_stream.close()
                    file_counter += 1
                    current_output_file = os.path.join(arguments.o, file_name + "-" + str(file_counter) + extension)
                    current_output_stream = open(current_output_file, "w")
                    read_counter = 0
                read_counter += 1
            print(line, file = current_output_stream, end="")
            line_mode += 1
            if line_mode == 4:
                line_mode = 0

    current_output_stream.close()


#####################################################################

if __name__ == '__main__':
    main()