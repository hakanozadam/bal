
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

from collections import OrderedDict

#################################################################
#################################################################

class SequentialEngine(OrderedDict):
    def run(self):
        for step, job in self.items():
            if job.did_success:
                print(step, " already run. Skipping...")
                continue
            else:
                print('Running', step)
                job.run()
                if job.did_fail:
                    print(step, 'failed')
                    exit(1)

