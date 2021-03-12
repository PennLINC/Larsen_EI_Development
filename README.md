# Guide the E:I development project  
The goal of this project is to measure the development of the E:I ratio using fMRI. First, we use a pharmalogical imaging dataset with the GABAergic benzodiazepine alprazolam to empirically generate a multivariate model for the effect of shifts in the E:I ratio on patterns of fMRI connectivity. We then apply this model to a large developmental and assess how model-predicted E:I ratio changes with age. This is a guide to the scripts used to execute the primary  analyses.  

## File overview
01: This downloads the output of all the flywheel processing (including fmriprep and XCP processing). Those scripts can be found in xcp_fw directory.  
02: This generates connectivity matrices for all the existing parcellations we have. (Atlases were generated with `createCoverageAtlash.sh`).  
03: These do the same as the above scripts but for all the PNC data.  
04: This trains the classifier models and then applies them to the PNC data.  
05: This evaluates the the results of the PNC model-generated labels and classification distances. It assesses age relationships as well as clinical relationships.  


## Notes on the steps for replication 
All the scripts should be run from CUBIC. The directory there are stored in is `/cbica/projects/alpraz_EI/scripts`.  
Many of the scripts pull input data from `/cbica/projects/alpraz_EI/input` or from `/cbica/projects/alpraz_EI/data`.  
The scripts are numbered according to the order they should run.  

1. This step downloads preprocessed data from flywheel.  

2. This script launches `qsub` jobs that generate connectivity matrices for all the pharmacological sessions and for all the atlases. 
The script uses AFNI's `3dNetCorr` function. For convenience reasons with qsub, the actual call to `3dNetCorr` is in the script `netcor_call.sh`.  
For replication, you can probably just generate a matrix for one or two subjects using the schaefer400 atlas. For example, you could run the following:  

3. The `03` scripts do the same as step 1 and 2 for the PNC data.  

4. This script does the SVM classification. My final analyses use 3 different classification models:  
- a classifier with all regions. 
- a classifier with transmodal areas only. 
- a classifier with unimodal areas only   

This should be run 3 times with different options each time.  
First, set `subdivide = FALSE` and `subdivision = "all"`.  
Second, set `subdivide = TRUE` and `subdivision = "transmodal25"`.  
Last, set `subdivide = TRUE` and `subdivision = "unimodal25"`.    
More information about this is included in the script comments.  

There is functionality to run permutation tests for each of these.  
`perm_test="permute_off"` for no permutations.   

`perm_test="permute_on"` for permutations.  

5. Finally, this Rmd script takes the output from step `4` and does some analyses.  
- It relies on the script `Alpraz_viz_functions.R` to make most of the figures/visualizations. 
  - That script contains all of the functions to visualize classifier results.
- It also fits the GAM models for the age effects. 

