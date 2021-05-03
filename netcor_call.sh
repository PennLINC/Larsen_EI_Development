#!/bin/bash
#$ -o /cbica/projects/alpraz_EI/output/job_output/
#$ -e /cbica/projects/alpraz_EI/output/job_output/
atlas=$1
prefix=$2
inset=$3

3dNetCorr -in_rois "$atlas" -prefix "$prefix" -inset "$inset"