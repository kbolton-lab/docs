function W12_FINAL_PROJ {
	# array2=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/10/" | grep -v "crai" | cut -d'_' -f1 | sort | head -n100 | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/log))
	# for eid in ${array2[@]}; do W10_FINAL $eid 10; done
	# BQSR2
	# 1=eid
	# 2=project

	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	subfolder=${1:0:2}
	project=$2
	echo $1 $subfolder $wf_folder $project
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $1 $subfolder $wf_folder $project

	# dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f $wf_folder/input_json/$subfolder/$1.json -y --priority low
	echo dx run workflow-G4xy25QJQ2884J0qJ2204fJ6 -f $wf_folder/input_json/$subfolder/$1.json -y --priority low
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12
}

function W12_COMBINE {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	subfolder=${1:0:2}
	project=$2
	echo $1 $subfolder $wf_folder $project
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $1 $subfolder $wf_folder $project

	echo dx run workflow-G4xzYxjJQ28FFXY59fv1ZxK7 -f $wf_folder/input_json/$subfolder/$1.json -y --priority low
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12
}

function TRANSPLANT {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	project=$KELLY
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt

	i=0
	for sample in $(grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/Transplant/eids.txt | tail -n+3); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 1200
		fi
	done
}