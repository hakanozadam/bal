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
from collections import namedtuple
import glob


###############################################################################
import sys
if __name__ == '__main__':
    script_directory = os.path.dirname(os.path.realpath(__file__))
    sys.path.append( os.path.split( script_directory)[0]  )
    import JobGroup , Job
else:
    from . import JobGroup, Job

class BalJobGroupException(Exception):
   pass

class BalJobGroup(JobGroup.JobGroup):

   def __init__(self , input_directory, output_directory ,
                  arguments, reference,
                  alignment_mode, threads,
                  name = "BAL" , run_time = 1,
                  executable = "", memory = 8192   ,
                  file_selection_mode = "directories" , **kwargs ):

      self.input_extensions = ["fastq" , "fq"]


      super().__init__( name = name, input_directory = input_directory , output_directory = output_directory ,
                        executable = executable, input_extensions = self.input_extensions ,
                        file_selection_mode = "directories",
                        **kwargs)

      if arguments:
          self.arguments       = arguments
      else:
          self.arguments = " "
      self.reference       = os.path.abspath(reference)
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

      self.bp_candidates_bed_file = os.path.join(self.output_directory, 'bp_candidates.bed')

      forbidden_aligner_parameter = re.compile("(?P<param>\s+-p\s+)")
      forbidden_parameter_match   = forbidden_aligner_parameter.search(self.arguments)
      if forbidden_parameter_match:
           raise(Bowtie2JobGroupException("You can not specify -p in the aligner parameters."))

   ###########################################################################

   def prepare_jobs(self):

      self.jobs = list()

      if self.alignment_mode == "single":
          for directory in self.input_files:
             arranged_input_files = self.arrange_input_files(directory)

             for file in arranged_input_files:

                 input_fastq              = file.strand_1_file
                 output_directory    = os.path.join(self.output_directory ,os.path.basename(directory), file.name)
                 success_file = os.path.join(output_directory, self.run_success_file_name)
                 if os.path.isfile(success_file):
                     continue
                 os.makedirs(output_directory, exist_ok=True )
                 summary_file             = os.path.join(output_directory,
                                                         "%s.summary"%file.name)
                 command = "{executable} -p {threads} {arguments} -x {reference} -o {output_directory} "\
                           " -U {input} &> {summary}".format(executable = self.executable ,\
                                                threads          = self.run_threads ,\
                                                arguments        = self.arguments ,\
                                                output_directory = output_directory,\
                                                reference        = self.reference ,\
                                                input            = input_fastq ,\
                                                summary          = summary_file )

                 job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = file.name,
                                       time_limit = self.run_time , cores = self.run_threads , memory = self.memory)
                 self.jobs.append(job)
      else:
          for directory in self.input_files:
              arranged_input_files = self.arrange_input_files(directory)
              for file in arranged_input_files:
                 input_1_fastq              = file.strand_1_file
                 input_2_fastq              = file.strand_2_file
                 output_directory      = os.path.join(self.output_directory , os.path.basename(directory) , file.name)
                 success_file = os.path.join(output_directory, self.run_success_file_name)
                 if os.path.isfile(success_file):
                     continue
                 os.makedirs(output_directory, exist_ok=True )
                 summary_file             = os.path.join(output_directory ,
                                                         "%s.summary"%file.name)
                 command = "{executable} -p {threads} {arguments} -x {reference} -o {output_directory} "\
                           " -1 {input_1} -2 {input_2}  &> {summary}".format(executable = self.executable ,\
                                                threads          = self.run_threads ,\
                                                arguments        = self.arguments ,\
                                                reference        = self.reference ,\
                                                output_directory = output_directory,\
                                                input_1          = input_1_fastq ,\
                                                input_2          = input_2_fastq ,\
                                                summary          = summary_file )

                 job = Job.Job(command = command ,
                                       output_directory = self.cluster_out_directory ,
                                       job_name = file.name,
                                       time_limit = self.run_time , cores = self.run_threads , memory = self.memory)
                 self.jobs.append(job)

   ###########################################################################

   ##########################################################################
   def post_run(self):
       self.error_message = ""
       missing_files = list()
       self.sam_directory = self.output_directory
       os.makedirs(self.sam_directory, exist_ok = True)

       for directory in self.input_files:
           arranged_input_files = self.arrange_input_files(directory)
           lib_bp_bed_files_list = os.path.join(self.output_directory ,os.path.basename(directory),
                                           "candidate_bp_files_list.txt")
           with open(lib_bp_bed_files_list, 'w') as lib_bp_files_stream:
              for file in arranged_input_files:
                  success_file = os.path.join(self.output_directory ,os.path.basename(directory), file.name,
                                                '.success')
                  bp_bed_file = os.path.join(self.output_directory ,os.path.basename(directory), file.name,
                                                'candidate_branchpoints', 'bp_candidates.bed')
                  if not os.path.isfile(success_file):
                        missing_files.append(success_file)
                  if not os.path.isfile(bp_bed_file):
                        missing_files.append(bp_bed_file)
                  print(bp_bed_file, file = lib_bp_files_stream)

       if len(missing_files) > 0 :
           self.error_message += "Error: Some jobs have failed because the following files are missing: %s"\
                                 " "%"\n".join(missing_files)

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

def main():
   parser = argparse.ArgumentParser(description=
   '''
   BAL Pipeline For Umass Cluster
   ''')
   parser.add_argument("-i" ,
                       help = "Input directory.." ,
                       required = True ,
                       metavar = "input directory" ,
                       type = str)
   parser.add_argument("-o" ,
                       help = "Output directory" ,
                       required = True ,
                       metavar = "output directory",
                       type = str)
   parser.add_argument("-a" ,
                       help = "Bal arguments" ,
                       required = False ,
                       metavar = "Arguments",
                       type = str)
   parser.add_argument("-x" ,
                       help = "Bal reference" ,
                       required = True ,
                       metavar = "Arguments",
                       type = str)
   parser.add_argument("-m" ,
                       help = "Alignment Mode" ,
                       required = True ,
                       metavar = "Mode",
                       choices = ['single' , 'paired'],
                       type = str)
   parser.add_argument("-p" ,
                       help = "Number of Threads" ,
                       required = False ,
                       default = 1,
                       metavar = "threads",
                       type = int)
   parser.add_argument("-t" ,
                       help = "Run time in HOURS" ,
                       required = False ,
                       metavar = "run_time",
                       default = 1,
                       type = int)
   parser.add_argument("--memory" ,
                       help = "Required Memory in MB" ,
                       required = False ,
                       metavar = "memory",
                       default = 8192,
                       type = int)
   args = parser.parse_args(JobGroup.correct_command_line_arguments("-a"))

   if not args.t:
       run_time = 1
   else:
       run_time = args.t

   jobGroup = BalJobGroup(
                 input_directory  = args.i ,
                 output_directory = args.o ,
                 arguments        = args.a ,
                 reference        = args.x ,
                 alignment_mode   = args.m ,
                 threads          = args.p ,
                 run_time         = args.t ,
                 memory           = args.memory
                  )

   jobGroup.run_in_main()
   ######################################################################################

if __name__ == '__main__':
   main()
