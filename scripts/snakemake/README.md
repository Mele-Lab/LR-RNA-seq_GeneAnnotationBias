```bash
snakemake \
-s Snakefile.smk \
-j 100 \
--latency-wait 120 \
--cluster "sbatch \
  -q debug \
  -c {resources.threads}  \
  --time=2:00:00" \
  -n
```
