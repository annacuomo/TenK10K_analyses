SNAKEFILE="./snakemake_sex_interaction_singlecells.smk"
CONFIG="cluster.json"
LOG="/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/logs"

nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --jobs 20 \
    --cluster \
      "qsub -S /bin/bash \
      -q short.q \
      -r yes \
      -pe smp 20 \  # number of threads
      -l tmp_requested=20G \
      -l mem_requested=20G \
      -o $LOG \
      -j y \  # this is to join error and log
      -V"
  > $LOG/snake_`date +%Y-%m-%d.%H:%M:%S`.log 
