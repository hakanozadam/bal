
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

from .step import Step
from .exceptions import *

#################################################################
#################################################################

class Bowtie2(Step):

   ###############################################################
   def __init__(self, name, input_files, output_directory, executable , executable_arguments ,\
                  genome_reference , single_end_fastq, paired_end_fastqs, p = 1):
      super().__init__(name, input_files, output_directory, executable, executable_arguments)

      self.executable_arguments = executable_arguments
      self.single_end_fastq     = single_end_fastq
      self.paired_end_fastqs    = paired_end_fastqs

      self.sam_output = os.path.join(self.output_directory, "alignment.sam")
      self.unaligned_fastq_single_end = os.path.join(self.output_directory, "unaligned.fastq")
      self.unaligned_fastq_mate_1 = os.path.join(self.output_directory, "unaligned.1.fastq")
      self.unaligned_fastq_mate_2 = os.path.join(self.output_directory, "unaligned.2.fastq")

      self.pre_genome_reference = genome_reference
      if p < 1 :
          self.p = 1
      else:
          self.p = p
      self.p = str(self.p)

   ################################################################

   @property
   def genome_reference(self):
       return self._genome_reference

   ################################################################

   @genome_reference.setter
   def genome_reference(self, ref_base):
       missing_references = list()
       suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
       for suffix in suffixes:
          if (not os.path.isfile(ref_base + suffix) ) and\
             (not os.path.isfile(ref_base + suffix + "l")):
              missing_references.append("Couldn't find the bowtie2 reference: " + ref_base + suffix + " or " +
              ref_base + suffix + "l")
       if len(missing_references) > 0:
           raise StepError("Couldn't find the following bowtie2 reference(s):\n" +\
                  "\n".join(missing_references))
       self._genome_reference = ref_base


   ################################################################

   def prepare(self):
      self.genome_reference = self.pre_genome_reference

      if self.single_end_fastq:
          if not os.path.isfile(self.single_end_fastq):
              raise(StepError("Couldn't find the input fastq file " + self.single_end_fastq))
      elif self.paired_end_fastqs[0] and self.paired_end_fastqs[1]:
          if not os.path.isfile(self.paired_end_fastqs[0]):
              raise(StepError("Couldn't find the mate 1 file " + self.paired_end_fastqs[0]))
          if not os.path.isfile(self.paired_end_fastqs[1]):
              raise(StepError("Couldn't find the mate 2 file " + self.paired_end_fastqs[1]))
      else:
          raise(StepError("No input files have been specified properly.\n"
                          "self.single_end_fastq = {single}\n"
                          "self.paired_end_fastqs[0] = {pair_0}\n"
                          "self.paired_end_fastqs[1] = {pair_1}\n"\
                          .format( single = self.single_end_fastq,
                                   pair_0 = self.paired_end_fastqs[0] ,
                                   pair_1 = self.paired_end_fastqs[1])))

      self.command = self.executable + " -p " + self.p + " " + self.executable_arguments +\
                    " -x {index} -S {sam_output} " \
                          " --no-unal "\
                    .format(index = self.genome_reference, sam_output = self.sam_output)

      if self.single_end_fastq:
          self.command += " -U {fastq_file} --un {unaligned_fastq_single_end}"\
                          .format(fastq_file = self.single_end_fastq ,
                                  unaligned_fastq_single_end = self.unaligned_fastq_single_end)
      else:
          self.command += " -1 {mate_1} -2 {mate_2} --un-conc {unaligned_paired_base} "\
                          .format(mate_1                = self.paired_end_fastqs[0] ,
                                  mate_2                = self.paired_end_fastqs[1] ,
                                  unaligned_paired_base = self.unaligned_fastq_single_end)

   ###################################################################

   def post_run(self):

       error_messages = list()
       if self.returncode:
           error_messages.append("The return code is {} ".format(self.returncode) )

       if self.single_end_fastq:
           if not os.path.isfile(self.unaligned_fastq_single_end):
               error_messages.append("Couldn't find the unaligned fastq file " +\
                   self.unaligned_fastq_single_end)
           if not os.path.isfile(self.sam_output):
               error_messages.append("Couldn't find the sam file " +\
                   self.sam_output)
       else:
           if (not os.path.isfile(self.unaligned_fastq_mate_1)) or\
               (not os.path.isfile(self.unaligned_fastq_mate_2) ):
               error_messages.append("Couldn't find the unaligned paired end fastq files: " +
                " , " + (self.unaligned_fastq_mate_1 , self.unaligned_fastq_mate_2))

       if len(error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
       else:
           subprocess.call('touch ' + self.success_file , shell=True )
           self._parse_hisat_summary()

       self.error_messages = error_messages

   #####################################################################

   #TODO: Write unittests for this function!
   def _parse_hisat_summary(self):
       if self.single_end_fastq:
           token_specs = [
               ("TOTAL_READS"       , "TOTAL_PERCENTAGE"    , "were unpaired"),
               ("UNALIGNED"         , "UNALIGNED_PERCENTAGE", "aligned 0 times"),
               ("UNIQUE_ALIGNMENTS" , "UNIQUE_PERCENTAGE"   , "aligned exactly 1 time"),
               ("MULTI_ALIGNMENTS"  , "MULTI_PERCENTAGE"    , "aligned >1 times"  )
           ]
       else:
           token_specs = [
               ("TOTAL_READS" ,        "TOTAL_PERCENTAGE" , "were paired"),
               ("ZERO_CONC_COUNT" ,    "ZERO_CONC_PERCENT" , "aligned concordantly 0 times"),
               ("UNIQUE_CONC_COUNT" ,  "UNIQUE_CONC_PERCENT" , "aligned concordantly exactly 1 time"),
               ("MULTI_CONC_COUNT" ,   "MULTI_CONC_PERCENT" , "aligned concordantly >1 times"),
               ("ALIGNED_DISC_COUNT" , "ALIGNED_DISC_PERCENT" , "aligned discordantly 1 time"),
               ("ZERO_MATE_COUNT" , "ZERO_MATE_PERCENT" , "aligned 0 times"),
               ("UNIQUE_MATE_COUNT" , "UNIQUE_MATE_PERCENT" , "aligned exactly 1 time"),
               ("MULTI_MATE_COUNT" , "MULTI_MATE_PERCENT" , "aligned >1 times")
           ]


       re_string = ''
       overall_string = r"(?P<OVERALL_PERCENT>.*)%\s+overall.*"
       for triple in token_specs:
           re_string += r'\s+(?P<' + triple[0] +\
                        r'>\d+)(?=.*\s+' + triple[2] + r'.*)|' +\
                        r'\((?P<' + triple[1] + r'>.*)%\)(?=\s+' + triple[2] + r'.*)|'

       re_string += overall_string + r"|(?P<SKIP>[\n]+)"
       parser_re = re.compile(re_string)
       results = dict()

       for mo in parser_re.finditer(self.std_err):
           kind = mo.lastgroup
           value = mo.group(kind)
           if kind != "SKIP":
               results[kind] = value

       self.alignment_statistics = results
       return results