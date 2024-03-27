```bash
bash /gpfs/projects/bsc83/utils/wrappers/Slurm-wrapper-Ruben-unstable_PauModified.sh \
  -j main_snakemake \
  -u 1 \
  -q debug \
  -t 2:00:00 \
  -c 'bash submit.sh'
```

```bash
module load intel mkl impi gcc/7.2.0
module load miniconda3/py39_4.10.3 && source activate snakemake
snakemake \
-s Snakefile \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q debug \
  -c {resources.threads}  \
  --time=2:00:00" \
  -n
```
