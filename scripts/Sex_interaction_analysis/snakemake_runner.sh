SNAKEFILE="./snakemake_sex_interaction_singlecells.smk"
CONFIG="./cluster.json"
LOG="/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/logs"

nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --latency-wait 30 \
    --cluster-config $CONFIG \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --jobs 20 \
    --cluster \
      "qsub -S /bin/bash \
      -q {queue} \
      -r yes \
      -pe smp {threads} \  # number of threads
      -l tmp_requested={memory}G \
      -l mem_requested={memory}G \
      -o $LOG \
      -j y \  # this is to join error and log
      -V" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
