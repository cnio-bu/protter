read -p "Is the conda environment activated? " -n 1 -r
echo  
if [[ $REPLY =~ ^[Yy]$ ]]
then
    LOGDIR=log/cluster
    mkdir -p $LOGDIR
    snakemake -s Snakefile --drmaa-log-dir $LOGDIR --drmaa ' -cwd -pe smp {cluster.threads} -l h_vmem={cluster.mem}M -q {cluster.queue}' --cores 20 --cluster-config cluster_sge.yaml $1
fi
