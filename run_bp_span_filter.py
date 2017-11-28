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
from bal.genomic_io.functions import bp_coverage_filter_sam

###################################################################

###################################################################

def main():
   parser = argparse.ArgumentParser(description=
   '''
   Filter a given sam file using bp_span_coverage_filter
   ''')
   parser.add_argument("-i" ,
                       help = "Input sam file" ,
                       required = True ,
                       metavar = "input_sam_file" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output sam file" ,
                       required = True ,
                       metavar = "output_sam_file",
                       type = str)
   arguments = parser.parse_args()

   bp_coverage_filter_sam(arguments.i, arguments.o)

############################################################################

if __name__ == "__main__":
    main()