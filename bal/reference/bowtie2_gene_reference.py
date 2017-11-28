
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

from shutil import which
import os
from subprocess import Popen, PIPE, call
import re
import glob
from ..core.step import Step
from ..core.exceptions import *
from ..settings import *

#################################################################
#################################################################

class Bowtie2GeneReference(Step):

   ###############################################################
   def __init__(self , name , input_files , output_directory , executable ,
                executable_arguments = "",
                genome_reference_name = "gene" , processes = 10 ):
       self.pre_input_files = input_files
       self.input_files = []
       super().__init__(name , self.input_files , output_directory , executable, executable_arguments)
       self.executable_arguments  = executable_arguments
       self.genome_reference_name = genome_reference_name
       self.fasta_directory       = os.path.abspath(input_files[0])


       self.bowtie2_ref_directory = os.path.join(self.output_directory, settings['bowtie2_gene_sub_directory'])
       os.makedirs(self.bowtie2_ref_directory, exist_ok=True)

       # Number of processes run in parallel
       self.processes = processes


   ################################################################

   def _executable_run(self):

        if not os.path.isdir(self.fasta_directory):
          raise FileNotFoundError("Could not find the directory " + self.fasta_directory )
        self.fasta_files = glob.glob(self.fasta_directory + "/*.fa")
        gene_with_extension = map(os.path.basename, self.fasta_files)
        self.genes = [ os.path.splitext(f)[0] for f in gene_with_extension ]

        i=0
        commands = list()
        active_processes = list()


        # Prepare the commands to be run
        for g in self.genes:
            reference_dir = os.path.join(self.bowtie2_ref_directory, g)
            os.makedirs(reference_dir, exist_ok=True)
            fasta_file = os.path.join(self.fasta_directory, g + ".fa")
            if not os.path.isfile(fasta_file):
                raise(FileNotFoundError("Could not find the fasta file " + fasta_file))
            reference_base = os.path.join(reference_dir, "gene")

            this_command = " ".join( (self.executable, fasta_file, reference_base) )
            commands.append(this_command)

        # Start running them in parallel so that there is at most (and preferably exactly, mostly)
        # self.processes many processes at a time.

        if len(self.genes) <= self.processes:
            for c in commands:
                active_processes.append(Popen([c + self.executable_arguments], stdout=PIPE,
                                                  stderr=PIPE, close_fds=True,
                                                  shell = True))
            while active_processes:
                for p in active_processes:
                    if p.poll() is not None:
                        p.stdout.close()
                        active_processes.remove(p)
        else:
            # First initialize the processes
            command_counter = 0
            for c in commands[0:self.processes]:
                active_processes.append(Popen([c + self.executable_arguments], stdout=PIPE,
                                                  stderr=PIPE, close_fds=True,
                                                  shell = True))

            for c in commands[self.processes:]:
                process_removed = False
                while not process_removed:
                    for p in active_processes:
                        if p.poll() is not None:
                            p.stdout.close()
                            process_removed = True
                            active_processes.remove(p)
                            active_processes.append(Popen([c + self.executable_arguments], stdout=PIPE,
                                                  stderr=PIPE, close_fds=True,
                                                  shell = True))
                            break

            # wait for the last active processes to finish
            while active_processes:
                for p in active_processes:
                    if p.poll() is not None:
                        p.stdout.close()
                        active_processes.remove(p)

        self.std_out = ""
        self.std_err = ""
        self.returncode = 0

        suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')

        # Now check if the refrence files have been created
        missing_references = list()

        for g in self.genes:
            reference_dir = os.path.join(self.bowtie2_ref_directory, g)
            reference_base = os.path.join(reference_dir, "gene")
            for suffix in suffixes:
                if (not os.path.isfile(reference_base + suffix) ) :
                    missing_references.append("Couldn't find the bowtie2 reference: " +\
                                              reference_base + suffix)

        if len(missing_references) > 0:
            error_message = g + " : The following bt2 references are missing " +\
                "\n".join(missing_references)
            raise(FileNotFoundError(error_message))

   ###################################################################

   def prepare(self):
      # to trick the parent class to run _executable_run function
      self.command = "dummy"

   ###################################################################

   def post_run(self):
       call('touch ' + self.success_file , shell=True )
