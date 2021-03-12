# Guide the Alpraz E:I balance project

## File overview
01: This downloads the output of all the flywheel processing (including fmriprep and XCP processing). Those scripts can be found in xcp_fw directory.  
02: This generates connectivity matrices for all the existing parcellations we have. (Atlases were generated with `createCoverageAtlash.sh`).  
03: These do the same as the above scripts but for all the PNC data.  
04: This does the trains the classifier models and then applies them to the PNC data.  
05: This evaluates the the results of the PNC model-generated labels and classification distances. It assesses age relationships as well as clinical relationships.  


## Notes on the steps for replication 
All the scripts should be run from CUBIC. The directory there are stored in is `/cbica/projects/alpraz_EI/scripts`.  
Many of the scripts pull input data from `/cbica/projects/alpraz_EI/input` or from `/cbica/projects/alpraz_EI/data`.  
The scripts are numbered according to the order they should run.  

1. This step downloads data from flywheel and does not need to be replicated.  

2. This script launches a ton of `qsub` jobs that generate connectivity matrices for all the alpraz sessions and for all the atlases. 
The script uses AFNI's `3dNetCorr` function. For convenience reasons with qsub, the actual call to `3dNetCorr` is in the script `netcor_call.sh`.  
For replication, you can probably just generate a matrix for one or two subjects using the schaefer400 atlas. For example, you could run the following:  
``` bash
atlas=schaefer400x7_aal_threshold_0.95.nii.gz
file_out=output_name
input_file=/cbica/projects/alpraz_EI/data/TASK_GSR/xcpengine/sub-012969/ses-001561/task-emotionid/space-MNI152NLin2009cAsym/task/stats/sub-012969_ses-001561_task-emotionid_space-MNI152NLin2009cAsym_res4d.nii.gz

3dNetCorr -in_rois "$atlas" -prefix $file_out -inset $input_file
```

Alternatively, you can use your preferred method of generating connectivity matrices and compare outputs. For comparison, the outputs from my script are stored like:
`/cbica/projects/alpraz_EI/data/TASK_GSR/12969/1561/schaefer400x7_aal_CC_GSR_000.netcc` for the above example subject.  

3. The `03` scripts do the same as step 1 and 2 for the PNC data. Not necesarry to run these.  

4. This script does the SVM classification. My final analyses use 3 different classification models:  
- a classifier with all regions. 
- a classifier with transmodal areas only. 
- a classifier with unimodal areas only   

This should be run 3 times with different options each time.  
First, set `subdivide = FALSE` and `subdivision = "all"`.  
Second, set `subdivide = TRUE` and `subdivision = "transmodal25"`.  
Last, set `subdivide = TRUE` and `subdivision = "unimodal25"`.    
More information about this is included in the script comments.  

There is functionality to run permutation tests for each of these, but I have it turned off by default because it takes a while. I don't think the permutations are necessary to replicate.  
`perm_test="permute_off"` for no permutations.   

`perm_test="permute_on"` for permutations.  

The important thing to replicate here is that you get significant classifier accuracy. The accuracies should be in the mid to high 60s.  There is code in the script for the next step that will print the classification accuracy and significance.  

5. Finally, this Rmd script takes the output from step `4` and does some analyses.  
- It relies on the script `Alpraz_viz_functions.R` to make most of the figures/visualizations. 
  - That script is basically a bunch of functions. It isn't documented right now, sorry, but it is just making images.
- The important things to replicate in script `05` are the GAM models that look at how classification distance is related to `age`. This is the central result of the paper. You can see in the script that there is a function that runs GAM models and prints the plots. This should be done for the Rmd sections that run the function for `Schaefer 400`. 

## For the blind replication
The most important set of analyses to replicate with your own scripts would be:
1. Running the classification model
2. Applying the model to the PNC data
3. Looking at the relationship between `age` and `distance` from the hyperplane.

### 1. Running the model(s)  
You can use the `data.frame`s that I have saved already. They are here:
Schaefer400: `/cbica/projects/alpraz_EI/input/CorMats/schaefer400x7_aal_GSR_all.rds` (all the pairwise connections).  
Schaefer400 transmodal: `/cbica/projects/alpraz_EI/input/CorMats/schaefer400x7_aal_GSR_transmodal25.rds` (all the connections to transmodal regions).  
Schaefer400 unimodal: `/cbica/projects/alpraz_EI/input/CorMats/schaefer400x7_aal_GSR_unimodal25.rds` (all the connections to unimodal regions).  

These can be loaded in `R` with 
``` R
df <- readRDS('/cbica/projects/alpraz_EI/input/CorMats/schaefer400x7_aal_GSR_transmodal25.rds')
```
for example.  
If you prefer to use Matlab, you can just open that file in `R` and then save as a csv with 
``` R
write.table(df,"output_filename.csv", row.names=F,sep=",")
```

The dataframe is one row per observation. The columns are as follows:
The first 2 are: "subid" (subject ID), "sesid" (the session ID, 2 per subject).  
The third column is the `y` variable, `drug`. This is the label for each session (`0` = `drug`; `1` = `placebo`).  
All the remaining columns are the features. The labels are meaningless for these. There is one column per pairwise connection.  
So, for the classifier, `y` is column 3 and `x` is column 4:end.  

### 2. Apply to the PNC data.
There are `rds` files saved for the PNC data that match the format of the alpraz files. They can be found here:  
`/cbica/projects/alpraz_EI/input/CorMats/PNC/schaefer400x7_aal_all.rds`  
`/cbica/projects/alpraz_EI/input/CorMats/PNC/schaefer400x7_aal_transmodal25.rds`  
`/cbica/projects/alpraz_EI/input/CorMats/PNC/schaefer400x7_aal_unimodal25.rds`  
The only difference is these don't have a column for `drug` because there is no drug manipulation here. (Column 1 and 2 are `subid` and `sesid`, the rest are features).  


### 3. Run age analyses.  
See the top section of `05_PNC_prediction_analysis.Rmd` to see how to load the demographic data, identify exclusions, and get things merged together.  Hopefully this is clear from the code.  
The critical model to run is:
`Distance ~s(age,k=4)+ sex + FD+ idemoBehAllEmoNrCount`  
`Distance` is the classification distance for each subject.  
`age` is the age of the subj.  
`sex` is subj sex.  
`FD` is mean FD.  
`idemoBehAllEmoNrCount` is the number of skipped responses from the `idemo` task. This is a proxy for alertness of the subject.  

This could be run, for example, as: 
``` R
model <- gam(Distance ~s(age,k=4)+ oSex+FD+ idemoBehAllEmoNrCount, data = my_data_frame, subset = exclusions==0)  
summary(model)
```
Again, see the script to get the exclusions.  
