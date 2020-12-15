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
    cat_cmd = f"cat {in_dir}/*.fa* > {out_dir}.fa"
    subprocess.check_call(cat_cmd, shell=True)
    return f"{out_dir}.fa"


def align_fasta(fasta_filepath, num_cpus=8):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file using mafft"""
    out_filepath = fasta_filepath.split('.')[0] + '_aligned.fa'
    msa_cmd = f"mafft --auto --thread {num_cpus} {fasta_filepath} > {out_filepath}"
    subprocess.check_call(msa_cmd, shell=True)
    return out_filepath


def align_fasta_reference(fasta_filepath, num_cpus=8, ref_fp: str=''):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file using mafft"""
    out_filepath = fasta_filepath.split('.')[0] + '_aligned.fa'
    msa_cmd = f"mafft --auto --thread {num_cpus} --keeplength --addfragments {fasta_filepath} {ref_fp} > {out_filepath}"
    subprocess.check_call(msa_cmd, shell=True)
    return out_filepath


def add_gene_column(df: pd.DataFrame) -> pd.DataFrame:
    """Takes dataframe containing intra-host variant information and adds a column indicating the gene for each mutation
     Uses the GFF_FEATURE to infer gene, and when GFF_FEATURE is missing (i.e. for indels), it uses nucleotide position 
     to infer gene. The `gff2gene` mapping was obtained from the file /home/al/data/hcov19/gff/NC_045512.2.gff3 and using 
     the create_gff2gene_mapping() function defined below."""
    gff2gene = {
        'cds-YP_009724389.1': 'ORF1ab',
        'cds-YP_009725295.1': 'ORF1ab',
        'cds-YP_009724390.1': 'S',
        'cds-YP_009724396.1': 'ORF8',
        'cds-YP_009724397.2': 'N',
        'cds-YP_009724395.1': 'ORF7a',
        'cds-YP_009724391.1': 'ORF3a',
        'cds-YP_009724393.1': 'M',
        'cds-YP_009724394.1': 'ORF6',
        'cds-YP_009725255.1': 'ORF10',
        'cds-YP_009725318.1': 'ORF7b',
        'cds-YP_009724392.1': 'E'
    }
    # infer the gene from GFF_FEATURE
    try:
        df['gene'] = df['GFF_FEATURE'].apply(lambda x: gff2gene.get(x, 'nan'))
    except:
        raise KeyError('GFF_FEATURE column not found in the input dataframe.')
    # infer the gene from position when GFF_FEATURE is missing
    try:
        df.loc[df['gene']=='nan', 'gene'] = df.loc[df['gene']=='nan', 'POS'].apply(map_gene_to_pos)
    except:
        raise KeyError('POS column not found in the input dataframe.')
    return df
    

def map_gene_to_pos(x):
    """helper function to infer the gene based on nucleotide position of SARS-CoV-2 genome"""
    pos = x
    if pos >= 269 and pos <= 21555:
        return 'ORF1ab'
    elif pos >= 21564 and pos <= 25382:
        return 'S'
    elif pos >= 25410 and pos <= 26214:
        return 'ORF3a'
    elif pos >= 26247 and pos <= 26471:
        return 'E'
    elif pos >= 26523 and pos <= 27187:
        return 'M'
    elif pos >= 27203 and pos <= 27382:
        return 'ORF6'
    elif pos >= 27398 and pos <= 27754:
        return 'ORF7a'
    elif pos >= 27757 and pos <= 27887:
        return 'ORF7b'
    elif pos >= 27896 and pos <= 28257:
        return 'ORF8'
    elif pos >= 28289 and pos <= 29528:
        return 'N'
    elif pos >= 29564 and pos <= 29670:
        return 'ORF10'
    return 'nan'


def create_gff2gene_mapping(variant_data: pd.DataFrame, gff_filepath: str) -> dict:
    fn = gffutils.example_filename(gff_filepath)
    db = gffutils.create_db(fn, dbfn='test.db', merge_strategy='merge', force=True)
    all_gffs = variant_data['GFF_FEATURE'].dropna().unique()
    gff2gene = {}
    for gff in all_gffs:
        if gff:
            gff2gene[gff] = db[gff]['gene'][0]
    return gff2gene