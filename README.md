# **Pipelines**

## *Collection of pipelines used in RNA/DNA Genomics Analysis*

# *0. Prerequisites and General Info.*

1. Download and install Nextflow: 
```

# Download using curl and run bash script
curl -s https://get.nextflow.io | bash

# Move nextflow to directory witihin the scope of $PATH variable
# To check directories in the path 
echo $PATH

# To move nextflow to the bin
mv nextflow /usr/bin
```

2. Install Singularity : https://docs.sylabs.io/guides/3.0/user-guide/installation.html#

OR 

Docker : https://docs.docker.com/engine/install/

**IF using Docker over Singularity you must go into the nextflow.config file and specify docker.enabled = true**
Otherwise the pipeline will try to use singularity container. As this pipeline was designed to run on a HPC 
it was by default designed using singularity but can run with docker just as eaisly.

3. The HPC management system I used when designing was SLURM, this can be changed howver by manipulating the executor of each process. 
For instance in the nextflow.config file 'slurm' can be changed to 'sge' to use SGE executor.

A detailed list of executors can be found here : https://www.nextflow.io/docs/latest/executor.html

4. For HPC users: must specify -profile cluster for jobs to be submitted to the cluster, some processes will still be executed locally.
You must manually specify the cluster options in the nextflow.config file. 

# *1. Usage*

*git clone this repo*
```
git clone https://github.com/benpastore/nextflow_pipelines.git
```

There are currently two pipelines set up 
1. CHIP-sequencing end-to-end analysis
2. mRNA seuqencing analysis

To run a given workflow use associated shell script in the directory with the main.nf for each workflow.

For example:
```
cd nextflow_pipelines

sh run.sh workflows/chip/run.sh -profile cluster -design design.tsv --single_end false 
```

Any params found in the nextflow.config file can be overridden by specifing them using '--' in the command line. 

For example if one wants to override single_end mode default (which is false) to true one could do:
```
cd nextflow_pipelines

sh run.sh workflows/chip/run.sh -profile cluster -design design.tsv --single_end true 
```

If one wanted to resume a failed/stalled pipeline one can include the -resume directive
```
cd nextflow_pipelines

sh run.sh workflows/chip/run.sh -profile cluster -design design.tsv --single_end false -resume
```



# *2. Directory Structure*

When in the nextflow_pipelines directory the structure is as follows: 

1. Modules
In the modules directory (nextflow_pipelines/modules) each software (i.e. bowtie, trim_galore...) has its own directory associated with a min.nf
file. 

```
./modules 
  |
  bowtie
    |
    main.nf --> Associated processes
      
2. Workflows 
Each workflow has its own names directory under workflows/. Each workflow has an associated run.sh script to run the associated pipeline.

./workflows
  |
  chip
    |
    run.sh
    main.nf
    nextflow.config
```

