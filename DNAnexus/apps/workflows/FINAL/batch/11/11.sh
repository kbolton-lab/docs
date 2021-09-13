function W11_FINAL {
	# array2=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/10/" | grep -v "crai" | cut -d'_' -f1 | sort | head -n100 | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/log))
	# for eid in ${array2[@]}; do W10_FINAL $eid 10; done
	# BQSR2
	# 1=eid
	# 2=folder_prefix=$3 # i.e. 10-60
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $1 $2 $wf_folder

	dx run workflow-G4gzV7QJ6XG2Q4yX4K4g4FPQ -f $wf_folder/input_json/$2/$1.json -y --name FINAL_11
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log11
}
#1100024