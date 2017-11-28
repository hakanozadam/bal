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

#####################################################################

def handle_arguments():
   parser = argparse.ArgumentParser(description=
   '''
   convert a given fasta file to a fastq file where the qualities will all be set to I
   ''')
   parser.add_argument("-i" ,
                       help = "Input fasta file" ,
                       required = True ,
                       metavar = "input_fasta_file" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output fastq file" ,
                       required = True ,
                       metavar = "output_fastq_file" ,
                       type = str)

   arguments = parser.parse_args()
   return arguments

#####################################################################

def main():
    arguments = handle_arguments()
    with FastaFile(arguments.i) as fasta_input_stream,\
        open(arguments.o, 'w') as fastq_output_file:

        for fasta_entry in fasta_input_stream:
            new_fastq_entry = FastqEntry(header = "Read_" + fasta_entry.header,
                                         sequence = fasta_entry.sequence,
                                         plus = '+',
                                         quality = "I"*len(fasta_entry.sequence))
            print(new_fastq_entry, file = fastq_output_file)

#####################################################################ß

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':

    sys.path.append( bal_dir  )
    from genomic_io.fasta import FastaFile, FastaEntry
    from genomic_io.fastq import FastqFile, FastqEntry
    
    main()
else:
    exit(1)

####################################################################
