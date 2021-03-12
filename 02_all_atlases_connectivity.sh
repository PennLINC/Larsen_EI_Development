#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/alpraz_EI/output/job_output/
#$ -e /cbica/projects/alpraz_EI/output/job_output/
#$ -l h_vmem=50.5G,s_vmem=50.3G
#$ -V
#$ -cwd

sublist="/cbica/projects/alpraz_EI/scripts/alpraz_sublist.csv" 

datadir=/cbica/projects/alpraz_EI/input/CorMats/
mkdir -p $datadir
atlasdir=/cbica/projects/alpraz_EI/input/atlases/

for atlas in $(ls $atlasdir/*_threshold_0.95.nii.gz); do
	fullname=$(basename $atlas)
	atlasName=${fullname%_thresh*}
	while IFS=, read subid sesid other; do
		if [ $sesid == "na" ];
		then continue
		else

		for datavers in "rest_noGSR" "rest_GSR" "GSR" "noGSR"; do
			if [ $datavers = "GSR" ]; then
				outDir="/cbica/projects/alpraz_EI/data/TASK_GSR/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/TASK_GSR/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/task/stats/"
				dataFile="${dataDir}/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_res4d.nii.gz"
			elif [ $datavers = "rest_GSR" ]; then
				outDir="/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS_GSR/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS_GSR/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/regress/"
				dataFile="${dataDir}/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_residualised.nii.gz"
			elif [ $datavers = "rest_noGSR" ]; then
				outDir="/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/regress/"
				dataFile="${dataDir}/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_residualised.nii.gz"
			else
				outDir="/cbica/projects/alpraz_EI/data/TASK/${subid}/${sesid}"
				dataDir="/cbica/projects/alpraz_EI/data/TASK/xcpengine/sub-0${subid}/ses-00${sesid}/task-emotionid/space-MNI152NLin2009cAsym/task/stats/"
				dataFile="${dataDir}/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_res4d.nii.gz"
			fi
			#Make the connectivity matrix
			if [ -f "${outDir}/${atlasName}_CC_${datavers}_000.netcc" ]; then
				continue
			else
				thisName="${datavers}_${subid}_${sesid}_${atlasName}"
				call="qsub -l short -N "$thisName" netcor_call.sh "$atlas" "${outDir}/${atlasName}_CC_${datavers}" "${dataFile}""
				# echo $subid $sesid $atlasName
				eval $call
				# 3dNetCorr -in_rois "$atlas" -prefix "${outDir}/${atlasName}_CC_${datavers}" -inset "${dataDir}/stats/sub-0${subid}_ses-00${sesid}_task-emotionid_space-MNI152NLin2009cAsym_res4d.nii.gz"
			fi
		done
		wait
	fi

		
	done <$sublist
done

