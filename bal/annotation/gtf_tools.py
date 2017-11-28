#!/bin/env python3

# AUTHORS:
#        Hakan Ozadam
#        Rachel Brown
#        Moore Laboratory
#        UMASS Medical School / HHMI
#        RNA Therapeutics Institute
#        Albert Sherman Center, ASC4-1009
#        368 Plantation Street
#        Worcester, MA 01605
#        USA
#
#################################################################

from ..genomic_io.gtf import GtfEntry, GtfFile
from collections import namedtuple, defaultdict
from functools import partial

#################################################################
#################################################################
# TODO
# Fix the transcript annotation end issue
# It shouldn't effect any results now though
transcript_annotation = namedtuple('transcript_annotation',['name', 'chr', 'gene', 'start', 'end' , 'strand', 'exons'])

##################################################################
def sort_gtf_file(gtf_file):
    '''
    Sorts the given gtf file according to seqname, gene_id,stranscript_id, start and end
    So we ignore fields not coming from transcripts.
    For this, we take only the features named "exon" and 'cds'
    :param gtf_file:
    :return: a sorted list of gtf entries
    '''
    gtf_entry_list = list()
    with GtfFile(gtf_file) as gtf_input:
        for gtf_entry in gtf_input:
            if gtf_entry.feature.lower() != "exon":
                continue
            gtf_entry_list.append(gtf_entry)

    gtf_entry_list.sort( key =  lambda x : (x.seqname , x.attribute_contents['gene_id'],\
                                            x.attribute_contents['transcript_id'], x.start, x.end) )
    return gtf_entry_list

####################################################################
def initialize_transcript_annotation(gtf_entry):
    this_annotation = transcript_annotation(name = gtf_entry.attribute_contents['transcript_id'],
                                            chr  = gtf_entry.seqname ,
                                            gene = gtf_entry.attribute_contents['gene_id'],
                                            start = gtf_entry.start, end = gtf_entry.end,
                                            strand = gtf_entry.strand,
                                            exons = [ (gtf_entry.start , gtf_entry.end) ])
    return this_annotation

#####################################################################

def find_transcripts_with_max_exons(gtf_file):
    '''
    For each gene, finds the transcript having maximum number of exons.
    Note that this will also gives us the transcript having maximum number of introns
    in that gene.
    The results are reported in a named tuple called transcript

    We don't assume that the gtf file is sorted. So we sort it ourselves.
    We assume that each transcript and gene has a unique name

    If you decide not to use ssort_gtf_file function, make sure that entries are coming from exons!
    '''

    gtf_entry_list       = sort_gtf_file(gtf_file)
    selected_transcripts = list()

    first_entry  = gtf_entry_list[0]
    current_gene = first_entry.attribute_contents['gene_id']
    previous_gene = current_gene
    current_transcript = first_entry.attribute_contents['transcript_id']
    previous_transcript = current_transcript
    current_transcript_annotation = initialize_transcript_annotation(first_entry)
    # make a deepcopy or initialize like this dont use current_.. = previous_...
    previous_transcript_annotation = initialize_transcript_annotation(first_entry)
    current_transcript_exon_count  = 1
    previous_transcript_exon_count = 1

    for gtf_entry in gtf_entry_list[1:]:
        current_gene = gtf_entry.attribute_contents['gene_id']
        if current_gene != previous_gene:
            if current_transcript_exon_count > previous_transcript_exon_count:
                selected_transcripts.append(current_transcript_annotation)
            else:
                selected_transcripts.append(previous_transcript_annotation)

            previous_gene       = current_gene
            current_transcript  = gtf_entry.attribute_contents['transcript_id']
            previous_transcript = current_transcript
            current_transcript_annotation  = initialize_transcript_annotation(gtf_entry)
            previous_transcript_annotation = initialize_transcript_annotation(gtf_entry)
            current_transcript_exon_count  = 1
            previous_transcript_exon_count = 1
            continue

        if current_transcript != gtf_entry.attribute_contents['transcript_id']:
            if current_transcript_exon_count > previous_transcript_exon_count:
                previous_transcript = current_transcript
                previous_transcript_annotation = current_transcript_annotation
                previous_transcript_exon_count = current_transcript_exon_count

            current_transcript = gtf_entry.attribute_contents['transcript_id']
            current_transcript_annotation  = initialize_transcript_annotation(gtf_entry)
            current_transcript_exon_count  = 1
        else:
            current_transcript_exon_count += 1
            current_transcript_annotation.exons.append( (gtf_entry.start, gtf_entry.end) )

    if current_transcript_exon_count > previous_transcript_exon_count:
        selected_transcripts.append(current_transcript_annotation)
    else:
        selected_transcripts.append(previous_transcript_annotation)

    return(selected_transcripts)

#######################################################################

#######################################################################

def main():
    pass

if __name__ == "__main__":
    main()