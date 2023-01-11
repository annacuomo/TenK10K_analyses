SNAKEFILE="./snakemake_sex_interaction_singlecells.smk"
CONFIG="cluster.json"
LOG="/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/Monocytes/logs"

nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --configfile $CONFIG \
    --rerun-incomplete \
    --jobs 20 \
    --restart-times 2 \
    --keep-going \
    --cluster \
        "qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -pe smp {cluster.threads} \
        -l tmp_requested={cluster.tmp_memory} \
        -l mem_requested={cluster.memory} \
        -e $LOG \
        -o $LOG \
        -j y \
        -V" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
