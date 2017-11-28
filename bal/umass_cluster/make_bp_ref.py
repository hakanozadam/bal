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
import sys
import prepare

script_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)) )
bal_dir = os.path.split( os.path.split( script_directory)[0] )[0]

if __name__ == '__main__':

    sys.path.append( bal_dir  )
    from bal.core.make_bp_reference import BpReference
else:
    exit(1)

####################################################################


####################################################################


####################################################################

def main():
   executables = prepare.get_executables(bal_dir)

   parser = argparse.ArgumentParser(description=
   '''
   Given a genome fasta file and branchpoint annotation in bed format, make bowtie2 branchpoint references.
   ''')
   parser.add_argument("-f" ,
                       help = "Genome fasta file" ,
                       required = True ,
                       metavar = "genome_fasta_file" ,
                       type = str)
   parser.add_argument("-b" ,
                       help = "Branchpoint bed file" ,
                       required = True ,
                       metavar = "brancpoint_bed_file" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output_directory" ,
                       type = str)
   parser.add_argument("-n" ,
                       help = "number of nucleotides on either side of the branchpoint" ,
                       required = False ,
                       default = 100,
                       metavar = "number_of_nucs" ,
                       type = int)
   arguments = parser.parse_args()

   make_ref_job = BpReference(name = "make_bp_ref", input_files = [arguments.f, arguments.b],
                              output_directory = arguments.o, executable = executables['bowtie2-build'],
                              executable_arguments = '',
                              number_of_nucleotides = arguments.n)
   make_ref_job.run()

#####################################################################

if __name__ == '__main__':
    main()