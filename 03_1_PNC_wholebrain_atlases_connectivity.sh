#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/alpraz_EI/output/job_output/
#$ -e /cbica/projects/alpraz_EI/output/job_output/
#$ -l h_vmem=50.5G,s_vmem=50.3G
#$ -V
#$ -cwd

sublist="/cbica/projects/alpraz_EI/input/n1601_sublist.csv" 

datadir=/cbica/projects/alpraz_EI/input/CorMats/PNC/
mkdir -p $datadir
atlasdir=/cbica/projects/alpraz_EI/input/atlases/

for atlas in $(ls $atlasdir/schaefer*wholebrain_2mm.nii.gz); do
	fullname=$(basename $atlas)
	atlasName=${fullname%_2mm*}
	while IFS=, read subid sesid other; do
		outDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/${subid}/${sesid}"
		dataDir="/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/xcpengine/sub-${subid}/ses-PNC1/task-idemo/space-MNI152NLin2009cAsym/task/stats/"
		dataFile="${dataDir}/sub-${subid}_ses-PNC1_task-idemo_space-MNI152NLin2009cAsym_res4d.nii.gz"

		#Make the connectivity matrix
		if [ -f "${outDir}/${atlasName}_CC_000.netcc" ]; then
			continue
		else
			thisName="sub${subid}_${sesid}_${atlasName}"
			call="qsub -l short -N "$thisName" netcor_call.sh "$atlas" "${outDir}/${atlasName}_CC" "${dataFile}""
			# echo $subid $sesid $atlasName
			eval $call
			# 3dNetCorr -in_rois "$atlas" -prefix "${outDir}/${atlasName}_CC_${datavers}" -inset "${dataDir}/stats/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_res4d.nii.gz"
		fi
		wait		
	done <$sublist
done


