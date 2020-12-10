import pandas as pd
import re
import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from itertools import groupby
import more_itertools as mit

def identify_deletions(input_filepath: str, patient_zero: str, 
                       min_del_len: int=2) -> pd.DataFrame:
    # read MSA file
    consensus_data = AlignIO.read(input_filepath, 'fasta')
    # prcess MSA to remove insertions and fix position coordinate systems
    seqs, ref_seq = process_cns_seqs(consensus_data, patient_zero)
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
    del_seqs['del_coords'] = del_seqs['del_positions'].apply(get_deletion_coord)
    # group sample by the deletion they share
    del_seqs = (del_seqs.groupby(['del_coords', 'del_len'])
                        .agg(samples=('idx', 'unique'),       # list of sample IDs with the deletion
                             num_samples=('idx', 'nunique'))  # num of samples with the deletion
                        .reset_index()
                        .sort_values('num_samples'))
    return del_seqs

def process_cns_seqs(cns_data, patient_zero) -> (dict, str):
    """Process aligned consensus sequences to prepare them for identifying deletions. The reference sequence is used to identify insertion positions, which are then removed and position numbers are updated."""
    # sequence for patient zero (before removing pseudo deletions)
    ref_seq = get_seq(cns_data, patient_zero)
    # identify insertions
    insertion_positions = identify_insertion_positions(ref_seq)
    # remove insertions from each sequence to consolidate correct nt positions
    for rec in cns_data:
        rec.seq = remove_insertions(str(rec.seq), insertion_positions)
    # sanity check: ensure that there are no "fake" deletions in reference sequence
    ref_seq = get_seq(cns_data, patient_zero)
    assert not identify_insertion_positions(ref_seq)
    # grab sequences from MSA
    seqs = get_seqs(cns_data)
    return seqs, ref_seq


# support functions
def get_seqs(bio_seqs, min_pos: int=265, max_pos: int=29674) -> dict:
    """Parse aligned sequences from Bio.Align.MultipleSeqAlignment to a dict object.
    The keys are sample names and values are their consensus sequences. 
    Each sequence is trimmed from both ends using `min_pos` and `max_pos`"""
    seqs = {}
    for row in bio_seqs:
        sample_name = str(row.id)
        s = str(row.seq)
        seqs[sample_name] = s[min_pos:max_pos]
    return seqs

def find_del_positions(x):
    return [m.start() for m in re.finditer('-', x)]

def find_deletions(x):
    del_positions = [m.start() for m in re.finditer('-', x)]
    deletions = [list(deletion) for deletion in mit.consecutive_groups(del_positions)]
    return deletions

def get_seq(all_seqs, sample_name: str) -> str:
    for rec in all_seqs:
        if rec.name == sample_name:
            seq = rec.seq
            break
    return str(seq)

def identify_insertion_positions(ref_seq: str) -> list:
    return [m.start() for m in re.finditer('-', str(ref_seq))]

def remove_insertions(seq: str, positions: list) -> str:
    for i, pos in enumerate(positions):
        seq = seq[:pos-i] + seq[pos+1-i:]
    return seq

def get_deletion_coord(x):
    min_pos = np.min(x)
    max_pos = np.max(x)
    return f'{min_pos}:{max_pos}'
    
def adjust_coords(x):
    start = int(x.split(':')[0])
    end = int(x.split(':')[1])
    return f'{start+265}:{end+265}'
    
def find_deletions_old(x):
    del_positions = [m.start() for m in re.finditer('-', x)]
    return [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(del_positions), 
                                                            lambda x: x[0]-x[1])]

def cross_join(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    df1 = df1.assign(key=0)
    df2 = df2.assign(key=0)
    return pd.merge(df1, df2, on='key').drop(columns='key')

def is_deletion_common(x):
    return x['del_positions_x']==x['del_positions_y']
