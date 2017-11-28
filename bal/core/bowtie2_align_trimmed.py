
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
from glob import glob

from .step import Step
from .exceptions import *
from ..settings import *

#################################################################
#################################################################

class Bowtie2AlignTrimmed(Step):

   ###############################################################
   def __init__(self, name, input_files, output_directory, bt2_reference_directory ,
                fastq_directory, executable , executable_arguments ,
                p = 10):
      super().__init__(name, [], output_directory, executable, executable_arguments)

      self.bt2_reference_directory = os.path.abspath(bt2_reference_directory)
      if not os.path.isdir(self.bt2_reference_directory):
          raise(FileNotFoundError("The reference directory does not exist: " +\
                                  self.bt2_reference_directory))
      self.fastq_directory         = os.path.abspath(fastq_directory)
      if not os.path.isdir(self.fastq_directory):
          raise( FileNotFoundError("The fastq directory doesn't exist!: " +\
              self.fastq_directory) )

      self.executable_arguments = executable_arguments
      self.genes    = list()
      self.commands = list()
      self.command  = ""

      self.alignment_directory = os.path.join(self.output_directory , "alignments")
      os.makedirs(self.alignment_directory, exist_ok=True)

      self.sam_files = list()

      if p < 1 :
          self.p = 1
      else:
          self.p = p

   ################################################################

   ################################################################

   def _check_reference(self, ref_base):
       missing_references = list()
       suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
       for suffix in suffixes:
           this_ref_file = ref_base + suffix
           if not os.path.isfile(this_ref_file):
               missing_references.append(this_ref_file)
       if len(missing_references) > 0:
           error_message = "Missing bowtie2 reference files: " +\
               "\n".join(missing_references)
           raise(FileNotFoundError(error_message))

   ################################################################

   def prepare(self):

      fastq_files = glob(self.fastq_directory + "/*.fastq")
      for file in fastq_files:
          this_gene = os.path.splitext( os.path.basename(file) )[0]
          this_ref_base = os.path.join(self.bt2_reference_directory, this_gene , "gene")
          self._check_reference(this_ref_base)
          output_directory = os.path.join(self.alignment_directory, this_gene)
          os.makedirs(output_directory, exist_ok=True)
          sam_file = os.path.join(output_directory, "alignments.sam")
          self.sam_files.append(sam_file)
          unaligned_fastq = os.path.join(output_directory, "unaligned.fastq")
          command = "{executable} {arguments} -x {index} -S {sam_output} -U {input_fastq} " \
                    "--un {unaligned_fastq} ".format( executable      = self.executable,
                                                      arguments       = self.executable_arguments,
                                                      index           = this_ref_base,
                                                      sam_output      = sam_file,
                                                      input_fastq     = file,
                                                      unaligned_fastq = unaligned_fastq)
          self.commands.append(command)


   ####################################################################

   def _module_run(self):
       self.run_parallel_subprocess(self.p)


   ####################################################################

   def post_run(self):
       missing_files = list()
       self.error_messages = list()

       for file in self.sam_files:
           if not os.path.isfile(file):
               missing_files.append(file)

       if len(missing_files) > 0 :
           self.error_messages.append("There are missing alignment files:\n" +
                                       "\n".join(missing_files) )
           subprocess.call('touch ' + self.failure_file , shell=True )
           self.returncode = 1
       else:
           subprocess.call('touch ' + self.success_file , shell=True )
           self.returncode = 0

   #####################################################################