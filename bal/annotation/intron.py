#!/bin/env python3
#
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

from collections import namedtuple, defaultdict, OrderedDict

from ..genomic_io.gtf import GtfEntry, GtfFile
from ..genomic_io.fasta import FastaFile, FastaEntry
from ..annotation.gtf_tools import find_transcripts_with_max_exons, sort_gtf_file
from ..settings import *

#################################################################

def get_intron_boundries_from_exon_boundires(exon_boundries):
    '''
    The input is an iterable of pairs where the first element of the pair
    is the exon start location and the second eleement of the pair is exon end location
    This function outputs a list of pairs where the first element of the pair is
    an intron start and the second element is an intron end

    The boundries of exons and introns are one-based and both inclusive
    '''
    intron_boundries = list()
    if len(exon_boundries) < 2:
        return list()
    intron_start = exon_boundries[0][1] + 1
    for exon in exon_boundries[1:]:
        intron_end = exon[0] - 1
        intron_boundries.append( (intron_start, intron_end) )
        intron_start = exon[1] + 1

    return intron_boundries

#################################################################

#################################################################

def get_intron_sequences(gtf_file, fasta_file, five_prime_slice_len,
                         intron_sequence_file, five_prime_slice_file):
    '''
    Given
    i)   gtf annotation in gtf file
    ii)  genome sequence in the fasta file
    iii) the number of nucleotides to be taken from the five prime end of the intron
          ( five_prime_slice_len )

    it creates two files in fasta format
    1) intron_sequence_file:
    Each entry comes from an intron
    The naming convention is
       chromosome_gene_intronIndex_strand_starPosition_endPosition
          chromosome, gene and strand are the obvious genomic features of the intron
          strand       : + , - : the strand of the gene the intron belongs to
          startPosition: The start position of the intron on the chromosome, 1-based, Inclusive
          endPosition  : The end position of the intron on the chromosome, 1-based, inclusive

    '''

    selected_transcripts = find_transcripts_with_max_exons(gtf_file)


    ## We want to reduce the memory footprint of our script.
    # So we'll read the fasta file chromosome-by-chromosome
    # and process the transcripts chromosome-by-chromosome
    # So this way, there will be only one chromosome sequence kept in memory at a time
    # To reduce the cost of search of transcripts,
    # we initially partition them into chromosomes in a dictionary
    # Also note that we want to find the intron intervals using other (more sophisticated)
    # methods in the future. Right now we are  using the introns of the transcripts where
    # those trasncripts have the max number of introns within the gene

    transcripts_by_chromosomes = defaultdict(list)
    for transcript in selected_transcripts:
        transcripts_by_chromosomes[transcript.chr].append(transcript)

    with FastaFile(fasta_file) as fasta_input,\
         open(intron_sequence_file, "w") as intron_seq_output,\
         open(five_prime_slice_file, "w") as five_prime_output:

        for fasta_entry in fasta_input:
            try:
                chr_transcripts = transcripts_by_chromosomes[fasta_entry.header]
            except KeyError:
                raise KeyError("There is a problem with the annotation. "
                      "Couldn't find a fasta entry with the header " + fasta_entry +\
                      ". It seems like chromosome names in your gtf file and fasta file are not compatible.")


            for transcript in chr_transcripts:
                if len(transcript.exons) < 2:
                    continue
                intron_boundries = get_intron_boundries_from_exon_boundires(transcript.exons)

                intron_counter = 1

                for intron_end_points in intron_boundries:
                    intron_header = "{chr}_{gene}_{intron_index}_{strand}_{startPos}_{endPos}".format(\
                                     chr = fasta_entry.header,
                                     gene = transcript.gene ,
                                     intron_index = intron_counter,
                                     strand = transcript.strand ,
                                     startPos = intron_end_points[0],
                                     endPos = intron_end_points[1])

                    intron_sequence = fasta_entry.sequence[intron_end_points[0] - 1 : (intron_end_points[1] ) ]

                    # If the intron sequence is too short, discard it
                    if len(intron_sequence) < 2 * five_prime_slice_len:
                        continue

                    # Fix for bowtie2:
                    # If the intron sequence contains many N's discard it.

                    intron_entry = FastaEntry(header = intron_header, sequence = intron_sequence)

                    if transcript.strand == '-':
                       intron_entry.reverse_complement()
                    five_prime_entry = FastaEntry(header = intron_header ,
                                                  sequence = intron_entry.sequence[0:five_prime_slice_len ]  )
                    N_counter = 0
                    for nucleotide in five_prime_entry.sequence:
                        if nucleotide not in ['A', 'C', 'G', 'T', 'a' , 'c', 'g', 't']:
                            N_counter += 1
                    if N_counter >= settings['max_N_in_intron']:
                        continue

                    print( str(intron_entry), file=intron_seq_output)
                    print( str(five_prime_entry), file=five_prime_output)

                    intron_counter += 1
    return selected_transcripts

#################################################################

def get_intron_five_prime_sequences(gtf_file, fasta_file, five_prime_slice_len,
                                    five_prime_slice_file, intron_bed_file):
    field_separator = "__"
    gtf_entry_list       = sort_gtf_file(gtf_file)



    # The keys are of the form
    # chr strand gene five_prime_location
    # the keys are lists that hold the 3 prime end of the corresponding intron starting at the specific location
    intron_list = defaultdict(list)

    intron_start        = gtf_entry_list[0].end
    intron_end          = -1
    previous_gene       = (gtf_entry_list[0]).attribute_contents["gene_id"]
    previous_transcript = (gtf_entry_list[0]).attribute_contents["transcript_id"]
    current_gene        = previous_gene
    current_transcript  = previous_transcript
    current_chr         = (gtf_entry_list[0]).seqname
    current_strand      = (gtf_entry_list[0]).strand

    introns_grouped_by_chr        = defaultdict( list )
    # if the intron is added to the dictionary introns_grouped_by_chr, then
    # we put the entry in the dictionary introns_in_grouped_by_chr_dic
    # where the key is the intron identifier and the value is 1
    # meaning that the intron is recorded.
    # We use this dict to make sure that we don't record introns in the
    # list inside the dictionary introns_grouped_by_chr
    introns_in_grouped_by_chr_dic = dict()

    for entry in gtf_entry_list[1:]:
        current_gene       = entry.attribute_contents["gene_id"]
        current_transcript = entry.attribute_contents["transcript_id"]
        current_chr        = entry.seqname
        current_strand     = entry.strand

        if current_transcript == previous_transcript:
            # We keep adding introns from the same transcript
            # first add the new intron
            intron_end = entry.start

            # We need to avoid the order nucleotides coming from the exons
            # So we add or subtract 1 from the boundary values as necessary
            if current_strand == '+':
                this_entry_five_prime_location  = intron_start + 1
                this_entry_three_prime_location = intron_end   - 1
            else:
                this_entry_five_prime_location  = intron_end   - 1
                this_entry_three_prime_location = intron_start + 1

            this_intron_identifier = field_separator.join((current_chr, current_strand, current_gene,
                                                           str(this_entry_five_prime_location) ))

            if not (this_entry_three_prime_location in intron_list[this_intron_identifier] ):
                intron_list[this_intron_identifier].append(this_entry_three_prime_location)

            if len( intron_list[this_intron_identifier] ) == 1 and \
                    introns_in_grouped_by_chr_dic.get(this_intron_identifier) != 1:
                introns_grouped_by_chr[current_chr].append(this_intron_identifier)
                introns_in_grouped_by_chr_dic[this_intron_identifier] = 1

            # Then update variables for the next possible intron
            intron_start = entry.end
            intron_end   = -1

        else:
            # We are on a new transcript
            intron_start        = entry.end
            intron_end          = -1
            previous_transcript = current_transcript
            previous_gene       = current_gene


    # for key, value in intron_list.items():
    #     print(key, " --> ", value)
    #     five_p_loc = key.split(field_separator)[3]
    #     print(five_p_loc)

    # Now that we have all the intron coordinates, we can pick the sequences from the fasta file
    with FastaFile(fasta_file) as fasta_input,\
        open(five_prime_slice_file , 'w') as fasta_output,\
        open(intron_bed_file, 'w') as bed_output:

        for fasta_entry in fasta_input:
            print(fasta_entry.header)
            if introns_grouped_by_chr.get(fasta_entry.header) == None:
                print("WARNING: Couldn't find the chromosome", fasta_entry.header,
                        "in the gtf file", gtf_file)
                continue

            for intron in introns_grouped_by_chr[fasta_entry.header]:
                contents         = intron.split(field_separator)
                five_prime_loc   = int(contents[3])
                strand           = contents[1]
                three_prime_ends = intron_list[intron]
                three_prime_ends.sort()

                # Note that gtf coordinates are 1 based and python string indeces are 0-based
                # So we needed the -1 adjustments
                if strand == "+":
                    entry_start = five_prime_loc - 1
                    entry_end   = five_prime_loc + five_prime_slice_len - 1
                elif strand == "-":
                    entry_start = five_prime_loc - five_prime_slice_len
                    entry_end = five_prime_loc
                    three_prime_ends.reverse()
                else:
                    print("WARNING: strand of the intron", intron, "is not + or - . It is ", strand)
                    continue

                intron_header     = intron + field_separator + field_separator.join( map(str, three_prime_ends) )
                intron_sequence   = fasta_entry.sequence[entry_start : entry_end]
                this_intron_entry = FastaEntry(intron_header, intron_sequence)

                N_counter = 0
                for nucleotide in intron_sequence:
                    if nucleotide not in ['A', 'C', 'G', 'T', 'a' , 'c', 'g', 't']:
                        N_counter += 1
                if N_counter >= settings['max_N_in_intron']:
                    continue

                if strand == "-":
                    this_intron_entry.reverse_complement()
                    bed_file_entry = "\t".join((fasta_entry.header, str(entry_start ), str(entry_end ),
                                            intron_header, "0", strand, this_intron_entry.sequence ))
                else:
                    bed_file_entry = "\t".join((fasta_entry.header, str(entry_start ), str(entry_end ),
                                            intron_header, "0", strand, this_intron_entry.sequence ))

                print(this_intron_entry, file = fasta_output)
                print(bed_file_entry, file    = bed_output)

##################################################################


def main():
    pass

if __name__ == "__main__":
    main()