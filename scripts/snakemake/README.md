```bash
# dryrun
bash /gpfs/projects/bsc83/utils/wrappers/Slurm-wrapper-Ruben-unstable_PauModified.sh \
  -j main_snakemake \
  -u 1 \
  -q debug \
  -t 2:00:00 \
  -c 'bash submit_dryrun.sh'

# run
bash /gpfs/projects/bsc83/utils/wrappers/Slurm-wrapper-Ruben-unstable_PauModified.sh \
  -j main_snakemake \
  -u 1 \
  -q bsc_ls \
  -t 2:00:00 \
  -c 'bash submit.sh'
```

```bash
# module load intel mkl impi gcc/7.2.0
# module load miniconda3/py39_4.10.3 && source activate snakemake
conda activate pt_snakemake
snakemake \
  -s snakemake/Snakefile \
  -j 10 \
  --latency-wait 120 \
  --cluster "sbatch \
  -q gp_bscls \
  -c 1  \
  -A bsc83 \
  --time=2:00:00"  
```
