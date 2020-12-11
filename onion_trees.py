import pandas as pd
import re
import numpy as np
from Bio import SeqIO, AlignIO, Phylo, Align
from itertools import groupby
import more_itertools as mit
from bjorn_support import map_gene_to_pos

def identify_deletions(input_filepath: str, patient_zero: str, min_del_len: int=2,
                       start_pos: int=265, end_pos: int=29674) -> pd.DataFrame:
    """Identify deletions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_del_len: minimum length of deletions to be identified"""
    # read MSA file
    consensus_data = AlignIO.read(input_filepath, 'fasta')
    # prcess MSA to remove insertions and fix position coordinate systems
    seqs, ref_seq = process_cns_seqs(consensus_data, patient_zero, start_pos, end_pos)
    # load into dataframe
    seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), columns=['sequence'])
                .reset_index().rename(columns={'index': 'idx'}))
    seqsdf['seq_len'] = seqsdf['sequence'].str.len()
    seqsdf['del_positions'] = seqsdf['sequence'].apply(find_deletions)
    # sequences with one or more deletions
    del_seqs = seqsdf.loc[seqsdf['del_positions'].str.len() > 0]
    del_seqs = del_seqs.explode('del_positions')
    # compute length of each deletion
    del_seqs['del_len'] = del_seqs['del_positions'].apply(len)
    # only consider deletions longer than 2nts
    del_seqs = del_seqs[del_seqs['del_len'] > min_del_len]
    # fetch coordinates of each deletion
    del_seqs['relative_coords'] = del_seqs['del_positions'].apply(get_coords)
    # group sample by the deletion they share
    del_seqs = (del_seqs.groupby(['relative_coords', 'del_len'])
                        .agg(samples=('idx', 'unique'),       # list of sample IDs with the deletion
                             num_samples=('idx', 'nunique'))  # num of samples with the deletion
                        .reset_index()
                        .sort_values('num_samples'))
    del_seqs['type'] = 'deletion'
    # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
    del_seqs['absolute_coords'] = del_seqs['relative_coords'].apply(adjust_coords, args=(start_pos,))
    # record the 5 nts before each deletion (based on reference seq)
    del_seqs['prev_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
    # record the 5 nts after each deletion (based on reference seq)
    del_seqs['next_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
    # approximate the gene where each deletion was identified
    del_seqs['gene'] = del_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0])).apply(map_gene_to_pos)
    return del_seqs


def identify_insertions(input_filepath: str, patient_zero: str, min_ins_len: int=2,
                       start_pos: int=265, end_pos: int=29674) -> pd.DataFrame:
    """Identify insertions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_ins_len: minimum length of insertions to be identified"""
    # read MSA file
    consensus_data = AlignIO.read(input_filepath, 'fasta')
    # load into dataframe
    ref_seq = get_seq(consensus_data, patient_zero)[start_pos:end_pos]
    insert_positions = identify_insertion_positions(ref_seq)
    seqs = get_seqs(consensus_data)
    seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), columns=['sequence'])
                .reset_index().rename(columns={'index': 'idx'}))
    seqsdf['seq_len'] = seqsdf['sequence'].str.len()
    seqsdf['ins_positions'] = seqsdf['sequence'].apply(find_insertions, args=(insert_positions,))
    # sequences with one or more deletions
    ins_seqs = seqsdf.loc[seqsdf['ins_positions'].str.len() > 0]
    ins_seqs = ins_seqs.explode('ins_positions')
    # compute length of each deletion
    ins_seqs['ins_len'] = ins_seqs['ins_positions'].apply(len)
    # only consider deletions longer than 2nts
    ins_seqs = ins_seqs[ins_seqs['ins_len'] > min_ins_len]
    # fetch coordinates of each deletion
    ins_seqs['relative_coords'] = ins_seqs['ins_positions'].apply(get_coords)
    # group sample by the deletion they share
    ins_seqs = (ins_seqs.groupby(['relative_coords', 'ins_len'])
                        .agg(samples=('idx', 'unique'),       # list of sample IDs with the deletion
                             num_samples=('idx', 'nunique'))  # num of samples with the deletion
                        .reset_index()
                        .sort_values('num_samples'))
    ins_seqs['type'] = 'insertion'
    # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
    ins_seqs['absolute_coords'] = ins_seqs['relative_coords'].apply(adjust_coords, args=(start_pos,))
    # record the 5 nts before each deletion (based on reference seq)
    ins_seqs['prev_5nts'] = ins_seqs['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
    # record the 5 nts after each deletion (based on reference seq)
    ins_seqs['next_5nts'] = ins_seqs['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
    # approximate the gene where each deletion was identified
    ins_seqs['gene'] = ins_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0])).apply(map_gene_to_pos)
    return ins_seqs

def process_cns_seqs(cns_data: Align.MultipleSeqAlignment, patient_zero: str,
                     start_pos: int, end_pos: int) -> (dict, str):
    """Process aligned consensus sequences to prepare them for identifying deletions. 
    The reference sequence is used to identify insertion positions, 
    which are then removed and position numbers are updated."""
    # sequence for patient zero (before removing pseudo deletions)
    ref_seq = get_seq(cns_data, patient_zero)
    # identify insertions (by identifyin 'fake' deletions in the aligned reference sequence)
    insertion_positions = identify_insertion_positions(ref_seq)
    # remove insertions from each sequence to consolidate correct nt positions
    for rec in cns_data:
        rec.seq = remove_insertions(str(rec.seq), insertion_positions)
    # sanity check: ensure that there are no "fake" deletions in reference sequence
    ref_seq = get_seq(cns_data, patient_zero)
    assert not identify_insertion_positions(ref_seq)
    # grab sequences from MSA
    seqs = get_seqs(cns_data, start_pos, end_pos)
    return seqs, ref_seq


# support functions
def get_seqs(bio_seqs: Align.MultipleSeqAlignment, min_pos: int=265, max_pos: int=29674) -> dict:
    """Parse aligned sequences from Bio.Align.MultipleSeqAlignment to a dict object.
    The keys are sample names and values are their consensus sequences. 
    Each sequence is trimmed from both ends using `min_pos` and `max_pos`"""
    seqs = {}
    for row in bio_seqs:
        sample_name = str(row.id)
        s = str(row.seq)
        seqs[sample_name] = s[min_pos:max_pos]
    return seqs


def find_deletions(x):
    del_positions = [m.start() for m in re.finditer('-', x)]
    deletions = [list(deletion) for deletion in mit.consecutive_groups(del_positions)]
    return deletions


def find_insertions(x, insert_positions: list):
    ins_positions = [m for m in insert_positions if x[m]!='-' and x[m]!='n']
    insertions = [list(insert) for insert in mit.consecutive_groups(ins_positions)]
    return insertions


def get_seq(all_seqs: Align.MultipleSeqAlignment, sample_name: str) -> str:
    """Fetches the aligned sequence of a specific sample name"""
    for rec in all_seqs:
        if rec.name == sample_name:
            seq = rec.seq
            break
    return str(seq)


def identify_insertion_positions(ref_seq: str) -> list:
    """helper function to identify positions where '-' was found in a sequence"""
    return [m.start() for m in re.finditer('-', str(ref_seq))]


def remove_insertions(seq: str, positions: list) -> str:
    for i, pos in enumerate(positions):
        seq = seq[:pos-i] + seq[pos+1-i:]
    return seq


def get_coords(x):
    """helper function to get deletion coordinates (start:end) from a Pandas Series containing deletion positions"""
    min_pos = np.min(x)
    max_pos = np.max(x)
    return f'{min_pos}:{max_pos}'
    
    
def adjust_coords(x, start_pos: int):
    """helper function to adjust deletion coordinates by adding the number of nucleotides that were trimmed (265nts)"""
    start = int(x.split(':')[0])
    end = int(x.split(':')[1])
    return f'{start+start_pos}:{end+start_pos}'


def find_del_positions(x):
    return [m.start() for m in re.finditer('-', x)]
    
    
def find_deletions_old(x):
    del_positions = [m.start() for m in re.finditer('-', x)]
    return [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(del_positions), 
                                                            lambda x: x[0]-x[1])]

def cross_join(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """helper function to perform a cross-join between two dataframes
    Useful for computing pairwise relationships...etc."""
    df1 = df1.assign(key=0)
    df2 = df2.assign(key=0)
    return pd.merge(df1, df2, on='key').drop(columns='key')


def is_deletion_common(x):
    return x['del_positions_x']==x['del_positions_y']