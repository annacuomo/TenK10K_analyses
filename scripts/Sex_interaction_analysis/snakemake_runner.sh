SNAKEFILE="./snakemake_sex_interaction_singlecells.smk"
CONFIG="cluster.json"
LOG="/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/logs"

snakemake \
  --snakefile $SNAKEFILE \
  --keep-going \
  --rerun-incomplete \
  --jobs 20 \
  --cluster \
    "qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -pe smp 20 \
    -l tmp_requested=20G \
    -l mem_requested=20G \
    -e $LOG \
    -o $LOG" \
  > $LOG/snake_`date +%Y-%m-%d.%H:%M:%S`.log 
