
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
from .bowtie2_reference import Bowtie2Reference
from .prepare import get_executables
from ..settings import *

#################################################################
#################################################################

class HISATReference(Bowtie2Reference):

    def __init__(self, known_splicesite_script, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.known_splicesite_script = known_splicesite_script

    def prepare(self):
        super().prepare()
        self.gtf_file                = self.pre_input_files[1]
        if not os.path.isfile(self.gtf_file):
            raise(StepError("Couldn't find the gtf file: " + self.gtf_file))
        self.known_splicesite_infile = os.path.join(self.output_directory,
                                                    settings['hisat_known_splice_sites_file'])
        self.command =  self.known_splicesite_script + " " + \
                           self.gtf_file + " > " + self.known_splicesite_infile + " ; " +\
                           self.command
   ###################################################################

    def post_run(self):
       missing_references = list()
       suffixes = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2',  '.rev.2.bt2')
       error_messages = list()

       if not os.path.isfile(self.known_splicesite_infile):
           error_messages.append("Couldn't find the known splice site file ",
                                 self.known_splicesite_infile)

       for suffix in suffixes:
          if (not os.path.isfile(self.reference_base + suffix) ) and\
             (not os.path.isfile(self.reference_base + suffix + "l")):
              missing_references.append("Couldn't find the hisat reference: " + self.reference_base + suffix)

       if len(missing_references) > 0:
           error_messages.append("Couldn't find the following hisat reference(s):\n" +\
                  "\n".join(missing_references))

       if len(error_messages) > 0:
           subprocess.call('touch ' + self.failure_file , shell=True )
       else:
           subprocess.call('touch ' + self.success_file , shell=True )
       self.error_messages = error_messages

   #####################################################################

