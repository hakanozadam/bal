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

from sys import stdin
### Fastq Reading Functions   ###############

class FastqEntry:
    def __init__(self , header , sequence , plus , quality):
        self.header   = header
        self.sequence = sequence
        self.plus     = plus
        self.quality  = quality

    def reverse_comp_sequence(self):
        complement_table = {"A" : "T" , "C" : "G" , "G" : "C", "T" : "A", "N" : "N" ,
                            "a" : "T" , "c" : "G" , "g" : "C", "t" : "C", "n" : "N"}
        raw_result = self.sequence[::-1]
        result = list()
        for nucleotide in raw_result:
            result.append(complement_table[nucleotide])

        return ''.join(result)

    def reverse_complement(self):
        self.sequence = self.reverse_comp_sequence()
        self.quality  = self.quality[::-1]

    def __str__(self ):
        return("\n".join( ('@' + self.header , self.sequence , self.plus , self.quality ) ) )

class FastqFile:
    def __init__(self , file):
        if(file):
            self.f = open(file , "r")
        else:
            #self.f = "a"
            self.f = stdin

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __del__(self):
        self.f.close()
         
    def __getitem__(self , index ):
        header   = self.f.readline().rstrip()
        if(header == ""):
            self.f.close()
            raise(IndexError)
        if header[0] != '@':
            raise(IOError('The first character of the fastq file must me a @'))
        sequence = self.f.readline().rstrip()
        plus     = self.f.readline().rstrip()
        quality  = self.f.readline().rstrip()
        return(FastqEntry( header[1:] , sequence , plus , quality  ))


