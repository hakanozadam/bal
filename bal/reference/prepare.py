#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#        Rachel Brown
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

import argparse
import os
from shutil import which
from sys import platform as _os

#################################################################

def get_commandline_arguments():
    ''' Parse and return the command line arguments'''

    parser = argparse.ArgumentParser(description=
    '''
    BAL Reference Prepare

    This script creates bowtie2 and HISAT references for BAL.
    In order to prepare HISAT and bowtie2 reference,
    BAL needs a whole genome reference in, fasta format, with exon annotation in a GTF file.
    BAL locally aligns the reads against the first N nucleotide of the introns.
    By default, N = 20 but this can be modified in the N parameter.
    ''')

    parser.add_argument("-g" ,
                    metavar  = 'gtf file' ,
                    help     = "GTF file annotating the exons in the genome of interest." ,
                    required = True ,
                    type     = str)
    parser.add_argument("-f" ,
                    metavar = 'Genomic Fasta File' ,
                    help = "The fasta file that contains the genomic sequence" ,
                    required = True ,
                    type = str)
    parser.add_argument("-N" ,
                    metavar = 'Number of five prime intron nucleotides' ,
                    help = "This is the number of five prime nucleotides in the intron where the reaqds are going to "
                           "be locally aligned against." ,
                    required = False ,
                    default = 20,
                    type = int)
    parser.add_argument("-o" ,
                    metavar = 'Output Directory' ,
                    help = "Output directory" ,
                    required = True ,
                    type = str)
    return parser.parse_args()

#################################################################################

def check_HISAT_files(ref_base):
   ''' TODO: Check for the existence of other files as well'''
   result = list()
   suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
   for suffix in suffixes:
      if (not os.path.isfile(ref_base + suffix) ) and\
         (not os.path.isfile(ref_base + suffix + "l")):
          result.append("Couldn't find the HISAT reference: " + ref_base + suffix + " or " +
          ref_base + suffix + "l")

   return result

#################################################################################

def process_commandline_arguments(cmd_args):
   ''' Check if the input files exist or not and do some consistency checks '''
   error_messages = list()
   if not os.path.isfile(cmd_args.f):
      error_messages.append("Couldn't find the fasta file " + cmd_args.f)
   if not os.path.isfile(cmd_args.g):
      error_messages.append("Couldn't find the gtf file " + cmd_args.g)

   if error_messages:
      print("Error!\nThe following error(s) occurred:")
      for error in enumerate(error_messages):
         print("{n}) {e}".format(n = error[0] +  1, e = error[1]))
      exit(1)

   return cmd_args

##################################################################################
def get_arguments():
   return process_commandline_arguments(get_commandline_arguments())

###################################################################################



###################################################################################

def get_executables(bin_directory):
   ''' Check the existence of executables: hisat, bowtie2
   Put their paths in a dictionary and return it'''

   #check the os and define bin variables for executables accordingly
   
   if _os == "linux" or _os == "linux2":
       hisat_relative_path    = 'bal/bin/hisat/linux_x86_64'
       bowtie2_relative_path  = 'bal/bin/bowtie2/linux_x86_64'
       bowtie2_build_relative_path  = 'bal/bin/bowtie2/linux_x86_64/bowtie2-build'

       
   elif _os == "darwin":
       hisat_relative_path    = 'bal/bin/hisat/mac_os_x_x86_64'
       bowtie2_relative_path  = 'bal/bin/bowtie2/mac_os_x_x86_64'
       bowtie2_build_relative_path  = 'bal/bin/bowtie2/mac_os_x_x86_64/bowtie2-build'
       print(bowtie2_build_relative_path)

   executables = dict()
   error_messages = list()

   executables['hisat']   = os.path.join(bin_directory, hisat_relative_path, 'hisat')
   executables['hisat-build']   = os.path.join(bin_directory, hisat_relative_path, 'hisat-build')
   executables['hisat_extract_splice_sites'] = os.path.join(bin_directory, hisat_relative_path,\
                                                            'extract_splice_sites.py')

   executables['bowtie2'] = os.path.join(bin_directory, bowtie2_relative_path,'bowtie2')
   executables['bowtie2-build'] = os.path.join(bin_directory, bowtie2_build_relative_path)

   for executable, path in executables.items():
      if not which(path):
         error_messages.append("Couldn't find the {executable} executable at {path}"\
            .format(executable = executable, path = path))

   if(error_messages):
      print('The following executable(s) are missing. If you have the files in the indicated path,'
            'make sure that the files are executable.')
      print("\n".join(error_messages))
      exit(1)

   return executables
