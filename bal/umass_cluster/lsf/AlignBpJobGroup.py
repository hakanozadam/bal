#!/bin/env python3

# AUTHOR:
#        Hakan Ozadam
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
from time import gmtime, strftime , localtime
import re
import shutil
from collections import namedtuple, defaultdict
import glob
import pysam


###############################################################################
import sys
if __name__ == '__main__':
    script_directory = os.path.dirname(os.path.realpath(__file__))
    sys.path.append( os.path.split( script_directory)[0]  )
    import JobGroup , Job
else:
    from . import JobGroup, Job


#################################################################################

def get_directories(maindir):
    contents = ( filter(lambda x: os.path.isdir(os.path.join(maindir, x)), os.listdir(maindir)) )
    return map( lambda y: os.path.join(maindir, y) , contents )

#################################################################################

class AlignBpJobGroupException(Exception):
   pass

#################################################################################

class AlignBpJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory ,
                  arguments, reference,
                  alignment_mode, threads,
                  name = "Align_Bp" , run_time = 4,
                  executable = "", memory = 8192   ,
                  rna_strandness = "F",
                  **kwargs ):

      self.input_extensions = ["fastq" , "fq"]

      super().__init__( name = name, input_directory = input_directory ,
                        output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        file_selection_mode = "directories",
                        **kwargs)

      if arguments:
          self.arguments       = arguments
      else:
          self.arguments = " "

      self.reference       = os.path.join(os.path.abspath(reference) , "bp")
      self.alignment_mode  = alignment_mode
      self.run_threads     = threads
      self.run_time        = run_time * 60 # Convert hours to minutes

      if self.run_time < 230:
          self.run_time = 230
      if self.run_threads < 1:
          self.run_threads = 1

      self.memory = 8092 if memory else memory

      self.exclude_system_files()
      self.jobs = list()

      self.run_success_file_name = ".success"
      self.sam_files_1 = defaultdict(list)
      self.sam_files_2 = defaultdict(list)

      self.individual_sam_files_directory = os.path.join(self.output_directory, "individual_sam_files")
      os.makedirs(self.individual_sam_files_directory, exist_ok = True)

      if rna_strandness not in ('R', 'F'):
         raise(Exception('Wrong rna_strandness. It can be either R or F'))
      self.rna_strandness = rna_strandness

      self.libraries_directory = os.path.join(self.output_directory, 'libraries')
      os.makedirs(self.libraries_directory, exist_ok=True)

      forbidden_aligner_parameter = re.compile("(?P<param>\s+-p\s+)")
      forbidden_parameter_match   = forbidden_aligner_parameter.search(self.arguments)
      if forbidden_parameter_match:
           raise(Bowtie2JobGroupException("You can not specify -p in the aligner parameters."))

   ###########################################################################

   def prepare_jobs(self):

      self.jobs = list()

      mate_1_option = " --norc "
      mate_2_option = " --nofw "

      if self.rna_strandness == 'R':
         mate_1_option = " --nofw "
         mate_2_option = " --norc "


      for directory in self.input_files:
          subdirectories = get_directories(directory)
          library_name = os.path.basename(directory)

          for subdir in subdirectories:
              arranged_input_files = self.arrange_input_files(
                  os.path.join(subdir, "hisat")  )
              file = arranged_input_files[0]
              input_fastq              = file.strand_1_file
              output_directory_base    = os.path.join(self.libraries_directory ,
                                                      os.path.basename(directory),
                                                      os.path.basename(subdir) )

              if self.alignment_mode == "single":
                  output_directory_1 = os.path.join(output_directory_base, 'align_bp_single_end')
                  success_file_1 = os.path.join(output_directory_1, self.run_success_file_name)
                  if os.path.isfile(success_file_1):
                      continue
              else:
                  output_directory_1 = os.path.join(output_directory_base, 'mate_1')
                  output_directory_2 = os.path.join(output_directory_base, 'mate_2')
                  os.makedirs(output_directory_2, exist_ok=True )
                  self.sam_files_2[library_name].append(os.path.join(output_directory_2, "alignment.sam"))
                  success_file_1 = os.path.join(output_directory_1, self.run_success_file_name)
                  success_file_2 = os.path.join(output_directory_2, self.run_success_file_name)
                  if os.path.isfile(success_file_1) and os.path.isfile(success_file_2):
                      continue
                  summary_file_2             = os.path.join(output_directory_2,
                                                            "%s.pipeline_summary"%file.name)


              self.sam_files_1[library_name].append(os.path.join(output_directory_1, "alignment.sam"))
              os.makedirs(output_directory_1, exist_ok=True )

              summary_file_1             = os.path.join(output_directory_1,
                                                        "%s.pipeline_summary"%file.name)
              command_1 = "{executable} -a \\\" {arguments} \\\" -x {reference} -o {output_directory} " \
                          "-p {threads} -U {input_fastq} &> {summary}".format(
                  executable = self.executable,
                  arguments = self.arguments + mate_1_option,
                  reference = self.reference,
                  output_directory = output_directory_1,
                  threads = self.run_threads,
                  input_fastq = input_fastq,
                  summary = summary_file_1
              )

              job_1 = Job.Job(command          = command_1 ,
                              output_directory = self.cluster_out_directory ,
                              job_name         = file.name,
                              time_limit       = self.run_time ,
                              cores            = self.run_threads ,
                              memory           = self.memory)

              self.jobs.append(job_1)

              if self.alignment_mode == "paired":
                  input_fastq_2              = file.strand_2_file

                  command_2 = "{executable} -a \\\" {arguments}\\\" -x {reference} -o {output_directory} " \
                              "-p {threads} -U {input_fastq} &> {summary}".format(
                      executable = self.executable,
                      arguments = self.arguments + mate_2_option,
                      reference = self.reference,
                      output_directory = output_directory_2,
                      threads = self.run_threads,
                      input_fastq = input_fastq_2,
                      summary = summary_file_2
                  )

                  job_2 = Job.Job(command          = command_2 ,
                                  output_directory = self.cluster_out_directory ,
                                  job_name         = file.name + "_2",
                                  time_limit       = self.run_time ,
                                  cores            = self.run_threads ,
                                  memory           = self.memory)

                  self.jobs.append(job_2)

   ##########################################################################

   def _merge_sam_helper(self, output_sam_file, sam_files):
       template_file     = pysam.AlignmentFile(list(self.sam_files_1.values())[0][0], 'r')
       output_sam_stream = pysam.AlignmentFile(output_sam_file, "wh", header = template_file.header)

       for file in sam_files:
           samfile = pysam.AlignmentFile(file, "r")

           for read in samfile.fetch():
               if not read.is_unmapped:
                   output_sam_stream.write(read)
           samfile.close()

       output_sam_stream.close()
       template_file.close()

   ##########################################################################

   def merge_sam_files(self):
       all_sam_files_1 = list()
       all_sam_files_2 = list()

       if self.alignment_mode == "single":
           for lib_name, sam_files in self.sam_files_1.items():
               this_output_directory = os.path.join(self.individual_sam_files_directory, lib_name)
               os.makedirs(this_output_directory, exist_ok = True)
               output_sam_file = os.path.join(this_output_directory , "alignment.sam")
               self._merge_sam_helper(output_sam_file, sam_files)
               all_sam_files_1 += sam_files
           output_sam_file = os.path.join(self.output_directory, "alignment.sam")
           self._merge_sam_helper(output_sam_file, all_sam_files_1)
       else:
           for lib_name, sam_files in self.sam_files_1.items():
               this_output_directory = os.path.join(self.individual_sam_files_directory, lib_name)
               os.makedirs(this_output_directory, exist_ok = True)
               output_sam_file_1 = os.path.join(this_output_directory , "alignment_1.sam")
               self._merge_sam_helper(output_sam_file_1, sam_files)
               all_sam_files_1 += sam_files

           for lib_name, sam_files in self.sam_files_2.items():
               this_output_directory = os.path.join(self.individual_sam_files_directory, lib_name)
               os.makedirs(this_output_directory, exist_ok = True)
               output_sam_file_2 = os.path.join(this_output_directory , "alignment_2.sam")
               self._merge_sam_helper(output_sam_file_2, sam_files)
               all_sam_files_2 += sam_files

           output_sam_file_1 = os.path.join(self.output_directory, "alignment_1.sam")
           output_sam_file_2 = os.path.join(self.output_directory, "alignment_2.sam")
           self._merge_sam_helper(output_sam_file_1, all_sam_files_1)
           self._merge_sam_helper(output_sam_file_2, all_sam_files_2)


   ##########################################################################

   def post_run(self):

       self.missing_files = list()
       self.error_message = ''

       all_sam_files = list()
       for lib_sam_files in self.sam_files_1.values():
          all_sam_files += lib_sam_files
       for lib_sam_files in self.sam_files_2.values():
          all_sam_files += lib_sam_files

       for sam_file in all_sam_files:
           if not os.path.isfile(sam_file):
               self.missing_files.append(sam_file)

       if len(self.missing_files) > 0 :
           self.error_message = "The following alignment files do not exist.\n" +\
                                "\n".join(self.missing_files)

       self.merge_sam_files()

       super().post_run()
      
   ###########################################################################

   def _get_file_2(self, lib_name, directory ):
       dir_files = glob.glob( directory + "/*" )
       ext2 = re.compile(r"(?P<ext2>.*)[._]2\.f(ast)?q$" , flags = re.IGNORECASE)

       for file in dir_files:
           file_base = os.path.basename(file)
           this_search = ext2.search(file_base)
           if this_search and lib_name == this_search.group("ext2"):
               return file

       return ""

###############################################################################

   def arrange_input_files(self, input_directory):

       input_files = list()
       for extension in self.input_extensions:
            input_files += glob.glob(input_directory + "/*." + extension)

       Library = namedtuple("Library" , ["name" , "strand_1_file" , "strand_2_file"])
       ext1 = re.compile(r"(?P<ext1>.*)[._]1\.f(ast)?q$" , flags = re.IGNORECASE)
       arranged_input_files = list()

       for file in sorted(input_files):
           file_name = os.path.splitext(os.path.basename(file))[0]

           if self.alignment_mode == "single":
               arranged_input_files.append( Library(file_name, file, "") )
           else:
               name_with_extension = os.path.basename(file)
               ext1_search = ext1.search(name_with_extension)
               if ext1_search:
                   lib_name = ext1_search.group("ext1")
                   file2 = self._get_file_2(lib_name, input_directory)
                   if not file2:
                       raise(BalJobGroupException("Couldn't find the strand 2 file for %s"%file))
                   arranged_input_files.append(Library(lib_name, file, file2) )

       return(arranged_input_files)

###############################################################################
