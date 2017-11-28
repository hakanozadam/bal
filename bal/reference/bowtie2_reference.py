
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

from shutil import which
import os
import subprocess
import re
from ..core.step import Step
from ..core.exceptions import *

#################################################################
#################################################################

class Bowtie2Reference(Step):

   ###############################################################
   def __init__(self , name , input_files , output_directory , executable , executable_arguments,
                genome_reference_name ):
       self.pre_input_files = input_files
       self.input_files = []
       super().__init__(name , self.input_files , output_directory , executable, executable_arguments)
       self.executable_arguments  = executable_arguments
       self.genome_reference_name = genome_reference_name
       self.fasta_file = input_files[0]
       self.reference_base = os.path.join(self.output_directory , self.genome_reference_name)

   ################################################################

   def prepare(self):
      if not os.path.isfile(self.fasta_file):
          raise(StepError("Couldn't find the fasta file", self.fasta_file))

      self.command = " ".join( (self.executable , self.executable_arguments ,
                              self.fasta_file , self.reference_base) )

   ###################################################################

   def post_run(self):

       missing_references = list()
       suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
       error_messages = list()

       for suffix in suffixes:
          if (not os.path.isfile(self.reference_base + suffix) ) :
              missing_references.append("Couldn't find the bowtie2 reference: " + self.reference_base + suffix)
       if len(missing_references) > 0:
           error_messages.append("Couldn't find the following bowtie2 reference(s):\n" +\
                  "\n".join(missing_references))

       if len(error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
       else:
           subprocess.call('touch ' + self.success_file , shell=True )
       self.error_messages = error_messages

   #####################################################################

