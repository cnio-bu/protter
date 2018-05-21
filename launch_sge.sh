read -p "Is the conda environment activated? " -n 1 -r
echo  
if [[ $REPLY =~ ^[Yy]$ ]]
then
    LOGDIR=log/cluster
    mkdir -p $LOGDIR
    snakemake -s Snakefile --drmaa-log-dir $LOGDIR --drmaa ' -cwd -pe smp {cluster.threads} -q {cluster.queue}' --cores 10 --cluster-config cluster_sge.yaml
fi
