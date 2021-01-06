from path import Path
import pandas as pd
import bjorn_support as bs
import mutations as bm
import subprocess
import shlex


date = '01-01-2021'
num_cpus = 25
in_alab_seqs = Path('/home/al/code/HCoV-19-Genomics/consensus_sequences/')
in_alab_meta = Path('/home/al/code/HCoV-19-Genomics/metadata.csv')
gisaid_seqs = Path('/home/al/analysis/gisaid/sequences_2021-01-01_08-20.fasta')
gisaid_meta = Path('/home/al/analysis/gisaid/metadata_2021-01-01_08-12.tsv')
ref_fp = Path('/home/al/data/hcov19/NC045512.fasta')
patient_zero = 'NC_045512.2'
out_dir = Path(f'/home/al/analysis/alab_mutations_{date}')
if not Path.isdir(out_dir):
    Path.mkdir(out_dir)
    print(f"Created results directory: {out_dir}")
else:
    print(f"Results directory {out_dir} already exists...Continuing...")
fa_fp = out_dir/'seqs.fa'
if not Path.isfile(fa_fp):
    fa_fp = bs.concat_fasta(in_alab_seqs, out_dir/'seqs')
print(f"Concatenated all sequences and wrote to {fa_fp}")
msa_fp = Path(fa_fp.split('.')[0] + '_aligned.fa')
if not Path.isfile(msa_fp):
    print(f"Aligning sequences with reference...")
    msa_fp = bs.align_fasta_reference(fa_fp, msa_fp, ref_fp=ref_fp, num_cpus=num_cpus)
print(f"Multiple sequence alignment of A-lab samples with reference saved in {msa_fp}")
msa2_fp = Path(fa_fp.split('.')[0] + '_aligned_absolute.fa')
if not Path.isfile(msa2_fp):
    print(f"Aligning sequences without reference...")
    msa2_fp = bs.align_fasta(fa_fp, msa2_fp, num_cpus=num_cpus)
print(f"Multiple sequence alignment of A-lab samples without reference saved in {msa2_fp}")
tree_fp = msa_fp + '.treefile'
if not Path.isfile(tree_fp):
    print(f"Computing phylogenetic tree...")
    tree_fp = bs.compute_tree(msa_fp, num_cpus=num_cpus)
print(f"Phylogenetic tree of A-lab samples saved in {tree_fp}")
subs_long_fp = out_dir/f'alab_substitutions_long_{date}.csv'
subs_long, _ = bm.identify_replacements_per_sample(msa_fp, in_alab_meta, bm.GENE2POS, data_src='alab')
subs_long.to_csv(subs_long_fp, index=False)
subs_wide = bm.identify_replacements(msa_fp, in_alab_meta)
subs_wide_fp = out_dir/f'alab_substitutions_wide_{date}.csv'
subs_wide.sort_values('num_samples', ascending=False).to_csv(subs_wide_fp, index=False)
print(f"Substitution-based mutations of A-lab samples saved in {subs_wide_fp}")
dels_long_fp = out_dir/f'alab_deletions_{date}.csv'
dels_long, _ = bm.identify_deletions_per_sample(msa_fp, in_alab_meta, patient_zero, bm.GENE2POS)
dels_long.to_csv(dels_long_fp, index=False)
dels_wide = bm.identify_deletions(msa_fp, in_alab_meta)
dels_wide_fp = out_dir/'deletions_aggregated.csv'
dels_wide.sort_values('num_samples', ascending=False).to_csv(dels_wide_fp, index=False)
print(f"Deletion-based mutations of A-lab samples saved in {dels_wide_fp}")
print(f"Aligning GISAID Sequences from {gisaid_seqs}...")
gisaid_msa_fp = Path(gisaid_seqs.split('.')[0] + '_aligned.fa')
if not Path.isfile(gisaid_msa_fp):
    gisaid_msa_fp = bs.align_fasta_reference(gisaid_seqs, num_cpus=25, ref_fp=ref_fp)
print(f"Multiple sequence alignment of GISAID Sequences saved in {gisaid_msa_fp}")
print("Analyzing Mutations...")
gisaid_subs_wide_fp = out_dir/f'gisaid_substitutions_wide_{date}.csv'
if not Path.isfile(gisaid_subs_wide_fp):
    print("Identifying substitution-based mutations - wide (aggregated)...")
    gisaid_subs = bm.identify_replacements(gisaid_msa_fp, gisaid_meta, data_src='gisaid')
    gisaid_subs.to_csv(gisaid_subs_wide_fp, index=False)
gisaid_subs_long_fp = out_dir/f'gisaid_substitutions_long_{date}.csv'
if not Path.isfile(gisaid_subs_long_fp):
    print("Identifying substitution-based mutations - long...")
    gisaid_subs_long, _ = bm.identify_replacements_per_sample(gisaid_msa_fp, gisaid_meta, bm.GENE2POS, data_src='gisaid')
    gisaid_subs_long.to_csv(gisaid_subs_long_fp, index=False)

gisaid_dels_wide_fp = out_dir/f'gisaid_deletions_wide_{date}.csv'
if not Path.isfile(gisaid_dels_wide_fp):
    print("Identifying deletion-based mutations - wide (aggregated)...")
    gisaid_dels = bm.identify_deletions(gisaid_msa_fp, gisaid_meta, data_src='gisaid')
    gisaid_dels.to_csv(gisaid_dels_wide_fp, index=False)
gisaid_dels_long_fp = out_dir/f'gisaid_deletions_long_{date}.csv'
if not Path.isfile(gisaid_subs_long_fp):
    print("Identifying deletion-based mutations - long...")
    gisaid_dels_long, _ = bm.identify_deletions_per_sample(gisaid_msa_fp, gisaid_meta, patient_zero, bm.GENE2POS, data_src='gisaid')
    gisaid_dels_long.to_csv(gisaid_dels_long_fp, index=False)
print(f"Deletion-based mutations of GISAID data saved in {gisaid_subs_fp}")