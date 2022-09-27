#!/bin/bash

which nextflow || curl -s https://get.nextflow.io | bash

singularity pull rnaseq.sif docker://benpasto/rnaseq:latest