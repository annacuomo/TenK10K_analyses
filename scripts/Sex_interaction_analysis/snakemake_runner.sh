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
      -q {cluster.queue} \
      -r yes \
      -pe smp {cluster.threads} \  
      -l tmp_requested={cluster.memory}G \
      -l mem_requested={cluster.memory}G \
      -o $LOG \
      -j y \  
      -V" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
