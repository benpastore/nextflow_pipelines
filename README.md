# **Pipelines Galore**

## *Collection of pipelines used in RNA/DNA Genomics Analysis*

# *0. Prerequisites*

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


