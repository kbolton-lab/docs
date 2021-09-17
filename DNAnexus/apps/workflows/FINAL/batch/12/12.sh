function __W12_FINAL {
	# array2=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/10/" | grep -v "crai" | cut -d'_' -f1 | sort | head -n100 | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/log))
	# for eid in ${array2[@]}; do W10_FINAL $eid 10; done
	# BQSR2
	# 1=eid
	# 2=folder_prefix=$3 # i.e. 10-60
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	2=${1:0:2}
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $1 $2 $wf_folder

	# mutect _scatter
	#dx run workflow-G4qkQJQJ6XG0P2Bz4Qj3zxk2 -f $wf_folder/input_json/$2/$1.json -y --priority low
	# mutect_single
	dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f $wf_folder/input_json/$2/$1.json -y --priority low
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12
}

# dx run workflow-G4v1kZjJQ286B56406bBFxvQ -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/12/1208696.json -y --priority low

# kelly
function W12_FINAL_K {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt
	sample=$1
	subfolder=${sample:0:2}
	project=$2
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $sample $subfolder $wf_folder $project
	echo dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY
	echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
}

# brian 
function W12_FINAL_B {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt
	sample=$1
	subfolder=${sample:0:2}
	project=$2
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $sample $subfolder $wf_folder $project
	
	echo dx run workflow-G4x8k0QJ6XG88z719jyVVX77 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $BRIAN
	echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
}