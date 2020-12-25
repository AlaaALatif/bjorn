import math
import re
import numpy as np
import pandas as pd
import more_itertools as mit
from Bio import SeqIO, AlignIO, Phylo, Align
from bjorn_support import map_gene_to_pos



GENE2POS = {
            '5UTR': {'start': 0, 'end': 265},
            'ORF1ab': {'start': 265, 'end': 21555},
            'S': {'start': 21562, 'end': 25384},
            'ORF3a': {'start': 25392, 'end': 26220},
            'E': {'start': 26244, 'end': 26472},
            'M': {'start': 26522, 'end': 27191},
            'ORF6': {'start': 27201, 'end': 27387},
            'ORF7a': {'start': 27393, 'end': 27759},
            'ORF7b': {'start': 27755, 'end': 27887},
            'ORF8': {'start': 27893, 'end': 28259},
            'N': {'start': 28273, 'end': 29533},
            'ORF10': {'start': 29557, 'end': 29674},
            '3UTR': {'start': 29674, 'end': 29902}
           }


def identify_replacements(input_fasta, 
                          meta_fp,
                          patient_zero: str='NC_045512.2', 
                          gene2pos: dict=GENE2POS):
    print(f"Loading Alignment file at: {input_fasta}")
    cns = AlignIO.read(input_fasta, 'fasta')
    print(f"Initial cleaning...")
    seqs, ref_seq = process_cns_seqs(cns, patient_zero,
                                     start_pos=0, end_pos=30000)
#     ref_seq = get_seq(cns, patient_zero)
#     seqs = get_seqs(cns, 0, 30000)
    print(f"Creating a dataframe...")
    seqsdf = (pd.DataFrame(index=seqs.keys(), 
                           data=seqs.values(), 
                           columns=['sequence'])
                .reset_index()
                .rename(columns={'index': 'idx'}))
    print(f"Identifying mutations...")
    # for each sample, identify list of substitutions (position:alt)
    seqsdf['replacements'] = seqsdf['sequence'].apply(find_replacements, args=(ref_seq,))
    # wide-to-long data manipulation
    seqsdf = seqsdf.explode('replacements')
    # initialize position column
    seqsdf['pos'] = -1
    # populate position column
    seqsdf.loc[~seqsdf['replacements'].isna(), 'pos'] = (seqsdf.loc[~seqsdf['replacements'].isna(), 'replacements']
       .apply(lambda x: int(x.split(':')[0])))
    # filter out non-substitutions
    seqsdf = seqsdf.loc[seqsdf['pos']!=-1]
    print(f"Mapping Genes to mutations...")
    # identify gene of each substitution
    seqsdf['gene'] = seqsdf['pos'].apply(map_gene_to_pos)
    seqsdf = seqsdf.loc[~seqsdf['gene'].isna()]
    # filter our substitutions in non-gene positions
    seqsdf = seqsdf.loc[seqsdf['gene']!='nan']
    print(f"Compute codon numbers...")
    # compute codon number of each substitution
    seqsdf['codon_num'] = seqsdf.apply(compute_codon_num, args=(GENE2POS,), axis=1)
    print(f"Fetch reference codon...")
    # fetch the reference codon for each substitution
    seqsdf['ref_codon'] = seqsdf.apply(get_ref_codon, args=(ref_seq, GENE2POS), axis=1)
    print(f"Fetch alternative codon...")
    # fetch the alternative codon for each substitution
    seqsdf['alt_codon'] = seqsdf.apply(get_alt_codon, args=(GENE2POS,), axis=1)
    print(f"Map amino acids...")
    # fetch the reference and alternative amino acids
    seqsdf['ref_aa'] = seqsdf['ref_codon'].apply(get_aa)
    seqsdf['alt_aa'] = seqsdf['alt_codon'].apply(get_aa)
    # filter out substitutions with non-amino acid alternates (bad consensus calls)
    seqsdf = seqsdf.loc[seqsdf['alt_aa']!='nan']
    print(f"Fuse with metadata...")
    # load and join metadata
    meta = pd.read_csv(meta_fp)
    seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
    # clean and process sample collection dates
    seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                   & (seqsdf['collection_date']!='1900-01-00')]
    seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
    seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
    # aggregate on each substitutions, compute number of samples and other attributes
    subs = (seqsdf.groupby(['gene', 'pos', 'ref_aa', 'codon_num', 'alt_aa'])
    .agg(
     num_samples=('ID', 'nunique'),
     first_detected=('date', 'min'),
     last_detected=('date', 'max'),
#      locations=('location', uniq_locs),
     location_counts=('location', lambda x: np.unique(x, return_counts=True)),
     samples=('ID', 'unique')
    )
    .reset_index())
    subs['locations'] = subs['location_counts'].apply(lambda x: list(x[0]))
    subs['location_counts'] = subs['location_counts'].apply(lambda x: list(x[1]))
    # 1-based nucleotide position coordinate system
    subs['pos'] = subs['pos'] + 1
    return subs


def find_replacements(x, ref):
    return [f'{i}:{n}' for i, n in enumerate(x) 
            if n!=ref[i] and n!='-' and n!='n']


def compute_codon_num(x, gene2pos: dict):
    pos = x['pos']
    ref_pos = gene2pos[x['gene']]['start']
    return math.ceil((pos - ref_pos + 1) / 3)


def get_ref_codon(x, ref_seq, gene2pos: dict):
    ref_pos = gene2pos[x['gene']]['start']
    codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
    return ref_seq[codon_start: codon_start+3].upper()


def get_alt_codon(x, gene2pos: dict):
    ref_pos = gene2pos[x['gene']]['start']
    codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
    return x['sequence'][codon_start: codon_start+3].upper()


def get_aa(codon: str):
    CODON2AA = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    return CODON2AA.get(codon, 'nan')


def uniq_locs(x):
    return list(np.unique(x))


def loc_counts(x):
    _, counts = np.unique(x, return_counts=True)
    return list(counts)


def identify_deletions(input_filepath: str, 
                       meta_fp: str,
                       patient_zero: str='NC_045512.2', 
                       min_del_len: int=2,
                       start_pos: int=265, 
                       end_pos: int=29674) -> pd.DataFrame:
    """Identify deletions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_del_len: minimum length of deletions to be identified"""
    # read MSA file
    consensus_data = AlignIO.read(input_filepath, 'fasta')
    # prcess MSA to remove insertions and fix position coordinate systems
    seqs, ref_seq = process_cns_seqs(consensus_data, patient_zero, start_pos, end_pos)
    # load into dataframe
    seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), 
                           columns=['sequence'])
                .reset_index().rename(columns={'index': 'idx'}))
    # load and join metadata
    meta = pd.read_csv(meta_fp)
    seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
    # clean and process sample collection dates
    seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                   & (seqsdf['collection_date']!='1900-01-00')]
    seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
    seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
    # compute length of each sequence
    seqsdf['seq_len'] = seqsdf['sequence'].str.len()
    # identify deletion positions
    seqsdf['del_positions'] = seqsdf['sequence'].apply(find_deletions)
    # sequences with one or more deletions
    del_seqs = seqsdf.loc[seqsdf['del_positions'].str.len() > 0]
    del_seqs = del_seqs.explode('del_positions')
    # compute length of each deletion
    del_seqs['del_len'] = del_seqs['del_positions'].apply(len)
    # only consider deletions longer than 2nts
    del_seqs = del_seqs[del_seqs['del_len'] >= min_del_len]
    # fetch coordinates of each deletion
    del_seqs['relative_coords'] = del_seqs['del_positions'].apply(get_indel_coords)
    # group sample by the deletion they share
    del_seqs = (del_seqs.groupby(['relative_coords', 'del_len'])
                        .agg(samples=('idx', 'unique'),
                             num_samples=('idx', 'nunique'),
                             first_detected=('date', 'min'),
                             last_detected=('date', 'max'),
#                              locations=('location', uniq_locs),
                             location_counts=('location', lambda x: np.unique(x, return_counts=True)))
                        .reset_index()
                        .sort_values('num_samples'))
    del_seqs['locations'] = del_seqs['location_counts'].apply(lambda x: list(x[0]))
    del_seqs['location_counts'] = del_seqs['location_counts'].apply(lambda x: list(x[1]))
    del_seqs['type'] = 'deletion'
    # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
    del_seqs['absolute_coords'] = del_seqs['relative_coords'].apply(adjust_coords, args=(start_pos+1,))
    del_seqs['pos'] = del_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0]))
    # approximate the gene where each deletion was identified
    del_seqs['gene'] = del_seqs['pos'].apply(map_gene_to_pos)
    del_seqs = del_seqs.loc[~del_seqs['gene'].isna()]
    # filter our substitutions in non-gene positions
    del_seqs = del_seqs.loc[del_seqs['gene']!='nan']
    # compute codon number of each substitution
    del_seqs['codon_num'] = del_seqs.apply(compute_codon_num, args=(GENE2POS,), axis=1)
    # fetch the reference codon for each substitution
    del_seqs['ref_codon'] = del_seqs.apply(get_ref_codon, args=(ref_seq, GENE2POS), axis=1)
    # fetch the reference and alternative amino acids
    del_seqs['ref_aa'] = del_seqs['ref_codon'].apply(get_aa)
    # record the 5 nts before each deletion (based on reference seq)
    del_seqs['prev_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
    # record the 5 nts after each deletion (based on reference seq)
    del_seqs['next_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
    return del_seqs[['type', 'gene', 'absolute_coords', 'del_len', 'pos', 
                     'ref_aa', 'codon_num', 'num_samples',
                     'first_detected', 'last_detected', 'locations',
                     'location_counts', 'samples',
                     'ref_codon', 'prev_5nts', 'next_5nts'
                     ]]


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
    if insert_positions:
        seqs = get_seqs(consensus_data)
        seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), columns=['sequence'])
                    .reset_index().rename(columns={'index': 'idx'}))
        seqsdf['seq_len'] = seqsdf['sequence'].str.len()
        seqsdf['ins_positions'] = seqsdf['sequence'].apply(find_insertions, args=(insert_positions,))
#             # approximate the gene where each deletion was identified
#         del_seqs['gene'] = del_seqs['pos'].apply(map_gene_to_pos)
#         del_seqs = del_seqs.loc[~del_seqs['gene'].isna()]
#         # filter our substitutions in non-gene positions
#         del_seqs = del_seqs.loc[del_seqs['gene']!='nan']
#         # compute codon number of each substitution
#         del_seqs['codon_num'] = del_seqs.apply(compute_codon_num, args=(GENE2POS,), axis=1)
#         # fetch the reference codon for each substitution
#         del_seqs['ref_codon'] = del_seqs.apply(get_ref_codon, args=(ref_seq, GENE2POS), axis=1)
#         # fetch the reference and alternative amino acids
#         del_seqs['ref_aa'] = del_seqs['ref_codon'].apply(get_aa)
#         # record the 5 nts before each deletion (based on reference seq)
#         del_seqs['prev_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
#         # record the 5 nts after each deletion (based on reference seq)
#         del_seqs['next_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
        # sequences with one or more deletions
        ins_seqs = seqsdf.loc[seqsdf['ins_positions'].str.len() > 0]
        ins_seqs = ins_seqs.explode('ins_positions')
        # compute length of each deletion
        ins_seqs['ins_len'] = ins_seqs['ins_positions'].apply(len)
        # only consider deletions longer than 2nts
        ins_seqs = ins_seqs[ins_seqs['ins_len'] >= min_ins_len]
        # fetch coordinates of each deletion
        ins_seqs['relative_coords'] = ins_seqs['ins_positions'].apply(get_indel_coords)
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
#         # approximate the gene where each deletion was identified
#         ins_seqs['gene'] = ins_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0])).apply(map_gene_to_pos)
        return ins_seqs#[['type', 'gene', 'absolute_coords', 'del_len', 'pos', 
#                      'ref_aa', 'codon_num', 'num_samples',
#                      'first_detected', 'last_detected',
#                      'locations', 'location_counts', 'samples',
#                      'ref_codon', 'prev_5nts', 'next_5nts'
#                      ]]
    else:
        return pd.DataFrame()


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


def get_indel_coords(x):
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