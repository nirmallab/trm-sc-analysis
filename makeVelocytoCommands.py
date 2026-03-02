# %%
from pathlib import Path
# %%
sample = 'BRI-2937'
num = 'BRI-2937'

baseCmd = f'docker run --rm \
  -v /home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell:/data \
  -v /data/trm-sc-analysis/data/external:/ref \
  velocyto:latest run \
    -b /data/{sample}/per_sample_outs/{num}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
    /data/{sample}/per_sample_outs/{num}/count/sample_alignments.bam \
    /ref/Mus_musculus.GRCm38.93.gtf'


# %%
samples = ['BRI-2937', 'BRI-2939', 'BRI-2941']
experimentDir = Path('/home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell')
allCmds = []
for sample in samples:
    outsDir = experimentDir / sample / 'per_sample_outs'
    for numDir in outsDir.iterdir():
        num = numDir.stem
        baseCmd = f'docker run --rm \
  -v /home/tyj566/mnt/h/lsp-analysis/tyler/jasonSingleCell:/data \
  -v /data/trm-sc-analysis/data/external:/ref \
  velocyto:latest run \
    -b /data/{sample}/per_sample_outs/{num}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
    /data/{sample}/per_sample_outs/{num}/count/sample_alignments.bam \
    /ref/Mus_musculus.GRCm38.93.gtf;'

        allCmds.append(baseCmd)

print('\n'.join(allCmds))