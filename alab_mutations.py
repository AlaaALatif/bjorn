from path import Path
import pandas as pd
import bjorn_support as bs
import mutations as bm
import subprocess
import shlex

in_alab_seqs = Path('/home/al/code/HCoV-19-Genomics/consensus_sequences/')
in_alab_meta = Path('/home/al/code/HCoV-19-Genomics/metadata.csv')
in_gisaid_seqs = Path('/home/al/analysis/gisaid/sequences_2021-01-01_08-20.fasta')
in_gisaid_meta = Path('/home/al/analysis/gisaid/metadata_2021-01-01_08-12.tsv')
ref_fp = Path('/home/al/data/hcov19/NC045512.fasta')
patient_zero = 'NC_045512.2'
out_dir = Path('/home/al/analysis/alab_mutations_01-01-2020')
Path.mkdir(out_dir)
print(f"Created results directory: {out_dir}")
fa_fp = bs.concat_fasta(in_alab_seqs, out_dir/'seqs')
print(f"Concatenated sequences and wrote to {fa_fp}")
print(f"Aligning sequences...")
msa_fp = bs.align_fasta_reference(fa_fp, num_cpus=25, ref_fp=ref_fp)
print(f"Multiple sequence alignment saved in {msa_fp}. Bye")
subs_long, _ = bm.identify_replacements_per_sample(msa_fp, in_alab_meta, bm.GENE2POS)
subs_long.to_csv(out_dir/'substitutions.csv', index=False)
subs_wide = bm.identify_replacements(msa_fp, in_alab_meta)
subs_fp = out_dir/'substitutions_aggregated.csv'
subs_wide.sort_values('num_samples', ascending=False).to_csv(subs_fp, index=False)
print(f"Substitution-based mutations saved in {subs_fp}")
dels_long, _ = bm.identify_deletions_per_sample(msa_fp, in_alab_meta, patient_zero, bm.GENE2POS)
dels_long.to_csv(out_dir/'deletions.csv', index=False)
dels_wide = bm.identify_deletions(msa_fp, in_alab_meta)
dels_fp = out_dir/'deletions_aggregated.csv'
dels_wide.sort_values('num_samples', ascending=False).to_csv(dels_fp, index=False)
print(f"Deletion-based mutations saved in {subs_fp}")
print(f"Aligning GISAID Sequences from {in_gisaid_seqs}...")
gisaid_msa_fp = bs.align_fasta_reference(in_gisaid_seqs, num_cpus=25, ref_fp=ref_fp)
print(f"Multiple sequence alignment of GISAID Sequences saved in {gisaid_msa_fp}. Bye")
