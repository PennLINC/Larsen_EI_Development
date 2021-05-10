#!/bin/bash

## Extract the GABRA5 values for each parcel from different atlases

for atlas in schaefer200x7_aal_threshold_0.95 schaefer400x7_aal_threshold_0.95 glasser360_aal_threshold_0.95 gordon333_threshold_0.95; do
	3dROIstats -1DRformat -nobriklab -mask ../${atlas}.nii.gz -nomeanout -nzmean GABRA5_2mm_MNI.nii.gz[0]  > GABRA5_${atlas}_vals.txt
	Rscript clean_ROIstats.R GABRA5_${atlas}_vals.txt
done
