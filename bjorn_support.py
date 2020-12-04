import glob
import subprocess
import path
import numpy as np
import pandas as pd

def get_variant_filepaths(sample_ids: list, analysis_path: str='/home/gk/analysis') -> dict:
    """Takes list of sample IDs and returns filepaths of variant data for corresponding samples"""
    variant_paths = {}
    for s_id in sample_ids:
        f = glob.glob(f"{analysis_path}/**/variants/illumina/*{s_id}*.tsv")
        variant_paths[s_id] = f[0]
    return variant_paths


def find_loc(f: str):
    """helper function to fetch location from filepath"""
    return f.split('/')[-1].split('_')[0].split('-')[-1]


def get_variant_data(variant_filepaths: dict):
    """Takes dict of variant filepaths and loads all variant data into dataframe"""
    df = (pd.concat((pd.read_csv(f, sep='\t')
                     .assign(sample=s_id, location=find_loc(f)) for s_id, f in variant_filepaths.items())))
    return df


def concat_fasta(in_dir, out_dir):
    """Concatenate fasta sequences into single fasta file"""
    cat_cmd = f"cat {in_dir}/*.fa > {out_dir}.fa"
    return subprocess.check_call(cat_cmd, shell=True)

def align_fasta(fasta_filepath):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file"""
    msa_cmd = f"mafft --auto --thread {num_cpus} {fasta_filepath}.fa > {fasta_filepath}_aligned.fa"
    return subprocess.check_call(msa_cmd, shell=True)