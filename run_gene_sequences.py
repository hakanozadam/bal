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
from bal.reference.extract_gene_sequences import ExtractGeneSequences

###################################################################

###################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Extracts Gene Sequences and Coordinates
   ''')
   parser.add_argument("-g" ,
                       help = "Gtf File" ,
                       required = True ,
                       metavar = "Gtf_File" ,
                       type = str)
   parser.add_argument("-f" ,
                       help = "Fasta File" ,
                       required = True ,
                       metavar = "Fasta_File" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_directory",
                       type = str)
   arguments = parser.parse_args()

   this_job = ExtractGeneSequences(
                 name = "Extract_Gene_Sequences",
                 input_files  = [arguments.g, arguments.f] ,
                 output_directory = arguments.o
                  )

   this_job.run()

############################################################################

if __name__ == "__main__":
    main()