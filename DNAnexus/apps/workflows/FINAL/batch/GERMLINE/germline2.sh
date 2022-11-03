

workflow=workflow-G775ggjJQ288g5QY22PYQQx4
function GERMLINE {
	i=1
	batch=$1
	workflow=$2
	project=$3
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400

	batch_total=$(wc -l $batch | awk '{print $1}')
	for sample in $(cat $batch); do
		subfolder=${sample:0:2}
		cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram | jq -r '.[].describe.id')
		cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram.crai | jq -r '.[].describe.id')
		echo $sample
		echo cram=$cram
		echo cram_index=$cram_index
		
		if [[ $i -lt 81 ]]; then
			printf "$sample	${sample}_23153_0_0.cram\t${sample}_23153_0_0.cram.crai\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80
		elif [[ $i -lt 161 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160
		elif [[ $i -lt 241 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240
		elif [[ $i -lt 321 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320
		else
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400
		fi
		
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}

		if [[ $batch_total -lt 80 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 80  ]]; 
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 160 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 160 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 240 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 240 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 320 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 320 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 400 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 400 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		fi
		i=$((i+1))
	done
}
export -f GERMLINE 


dx ls project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.tsv" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/brian.complete.txt
for batch in aa ac ae ag ai ak am ao aq as au aw; do
    grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/brian.complete.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch} > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{a,c,e,g,i,k,m,o,q,s,u,w}_remaining
for batch in aw; do
    GERMLINE /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining $workflow $BRIAN
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{a,c,e,g,i,k,m,o,q,s,u,w}_remaining

--instance-type final_annotation_and_declutter=mem1_ssd1_v2_x4
--instance-type stage-G4KqPz0J6XGJZGB842qJVYQK=mem1_ssd1_v2_x4
workflow=workflow-G77BY70J6XG09X1BPy08GzFz

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_aq_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_au_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_aw_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type stage-G4KqPz0J6XGJZGB842qJVYQK=mem1_ssd1_v2_x4


 dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_aq_remaining.batch.header.801  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'


dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.tsv" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/kelly.complete.txt
for batch in ab ad af ah aj al an ap ar at av; do
    grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/kelly.complete.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch} > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{b,d,f,h,j,l,n,p,r,t,v}_remaining
for batch in ab ad af ah aj al an ap ar at av; do
    GERMLINE /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining $workflow $KELLY
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{b,d,f,h,j,l,n,p,r,t,v}_remaining

dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ap_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'











for file in $(ls *.final.tsv | head -n1 ); do
	echo $file
	head -n1 $file | TRANSPOSE | sha256sum
done

rm change.cols.txt
touch change.cols.txt
for file in $(ls *.final.tsv); do
	sha=$(head -n1 $file | sha256sum | awk '{print $1}')
	if [[ $sha != "31fb67749b92af0e715a62b2d6fc586c44cb68b4841ad3dc70c4648d83bf770b" ]]; then
		echo $file
		echo $file >> change.cols.txt
	fi
done

echo {40..50}0000

for file in $(cat change.cols.txt | tail -n+11); do
	eid=$(echo $file | cut -d_ -f1)
	./cols.R ${eid}_23153_0_0.final.tsv 1000051_23153_0_0.columns.txt ${eid}_23153_0_0.final.updated.tsv
	sha=$(head -n1 ${eid}_23153_0_0.final.updated.tsv | sha256sum | awk '{print $1}')
	if [[ $sha == "31fb67749b92af0e715a62b2d6fc586c44cb68b4841ad3dc70c4648d83bf770b" ]]; then
		printf "$file\tGOOD\n"
		mv ${eid}_23153_0_0.final.updated.tsv ${eid}_23153_0_0.final.tsv
	else
		printf "$file\tBAD\n"
	fi
	echo
done

for file in $(ls *.final.tsv); do
	ncols=$(head -n1 $file | TRANSPOSE | wc -l)
	if [[ $ncols -ne 277 ]]; then
		echo $file
	fi
done


for i in {1..28}0; do
	for file in 1000051_23153_0_0.final.tsv 1004632_23153_0_0.final.tsv; do
		printf "$file\t$i\t"
		head -n1 $file | TRANSPOSE | head -n $i | sha256sum
	done
	echo
done

for file in 1000051_23153_0_0.final.tsv 1004632_23153_0_0.final.tsv; do
	head -n1 $file | TRANSPOSE | head -n 52 | tail -n1
done


(samtools view DNA01917.kraken_filtered.ends_trimmed.rehead.bam -H; samtools view DNA01917.kraken_filtered.ends_trimmed.rehead.bam -h | gawk -F'\t' 'gsub(/[0-9]{1,3}H/,"",$6) {print}' OFS='\t') | samtools view -@ 16 -b -o DNA01917.kraken_filtered.ends_trimmed.rehead.cigar_edit.bam && samtools index DNA01917.kraken_filtered.ends_trimmed.rehead.cigar_edit.bam

(samtools view DNA01917.kraken_filtered.ends_trimmed.rehead.11000001_11500000.bam -H; samtools view DNA01917.kraken_filtered.ends_trimmed.rehead.11000001_11500000.bam -h | head -n419314 | tail -n1) | samtools view -@ 16 -h


cat hotspots.germline.DV.tsv | while read key eid index; do
	GERMLINE $key $eid $index
done

cat germline.ch_pd2.DV.tsv | while read key eid index; do
	GERMLINE $key $eid $index
done > germline.ch_pd2.GERMLINE.tsv

BRIAN:
Transplant in: 
	$BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/13/
	$BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/
Lung in:
	$BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/

KELLY:
	$KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/
	$KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant/
	$KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/

i=1
for sample in 1002297 1240236 1300667 4137322; do


i=1
for sample in $(cat germline_eids_already_ran.uniq.txt); do
	subfolder=${sample:0:2}
	printf "$i\t$subfolder\t$sample\t"
	if [[ $subfolder == 12 ]]; then
		file=$(dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12 --json --name ${sample}_23153_0_0.final.tsv --all-projects | jq -cr '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12") | .describe.id')
		printf "12\t"
	elif [[ $subfolder == 13 ]]; then
		file=$(dx find data --class file --path $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/13 --json --name ${sample}_23153_0_0.final.tsv | jq -cr '.[].describe.id')
		printf "13\t"
	else
		file=$(dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung --json --name ${sample}_23153_0_0.final.tsv --all-projects | jq -cr '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung") | .describe.id')
		tl="Lung"
		if [[ $file == "" ]]; then
			file=$(dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant --json --name ${sample}_23153_0_0.final.tsv --all-projects | jq -cr '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant") | .describe.id')
			tl="Transplant"
		fi
		printf "$tl\t"
	fi
	printf "$file\n"
	i=$((i+1))
done

/storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/filtered/
1323241_23153_0_0.final.tsv

dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung --json --name ${sample}_23153_0_0.final.tsv --all-projects

dx find data --class file --json --name ${sample}_23153_0_0.final.tsv --all-projects --norecurse 

dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/${sample}_23153_0_0.final.tsv --json --all-projects


cat test.json | jq -c 'select( .[].describe.folder == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/")'

cat test.json | jq -c '.[].describe.folder == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung"'

cat test.json | jq -c '.[].describe.folder == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung"'

dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung --json --name ${sample}_23153_0_0.final.tsv --all-projects > test.json

cat test.json | jq -c '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung")'


for sample in 1277250 1323241 1379716; do
	subfolder=${sample:0:2}
	printf "$i\t$subfolder\t$sample\t"

	file=$(dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung --json --name ${sample}_23153_0_0.final.tsv --all-projects | jq -cr '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung") | .describe.id')
	if [[ $file == "" ]]; then
		file=$(dx find data --class file --path /CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant --json --name ${sample}_23153_0_0.final.tsv --all-projects | jq -cr '.[] | select((.describe.folder | tostring) == "/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant") | .describe.id')
	fi
	# fi
	printf "$file\n"
	i=$((i+1))
done

12 length = 38997831
13 length = 38997831
Tr length = 38997831
Lu length = 39652397

./alreadyRan.R \
	-i /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.final.tsv \
	--output /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.updated.tsv \
	--target-length 38997831 \
	--col-order /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/1000051_23153_0_0.columns.txt \
	-T ~/Bolton/data/gene_census_TSG.txt \
	--oncoKB-curated ~/Bolton/data/all_curated_genes_v2.0.tsv \
	-p ~/Bolton/data/pd_table_kbreview_bick_trunc2.txt \
	--pan-myeloid ~/Bolton/data/panmyeloid_variant_counts.vep.annotated.vcf.tsv \
	--bolton-bick-vars ~/Bolton/data/bick.bolton.vars3.txt \
	--mut2-bick ~/Bolton/data/topmed.n2.mutation.c.p.txt \
	--mut2-kelly ~/Bolton/data/kelly.n2.mutation.c.p.txt \
	--matches2 ~/Bolton/data/matches.2.c.p.txt

./alreadyRan.R \
	-i /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.final.tsv \
	--output /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.updated.tsv \
	--target-length 38997831 \
	--col-order /scratch1/fs1/bolton/UKBB/germline/already/1000051_23153_0_0.columns.txt \
	-T /scratch1/fs1/bolton/UKBB/germline/already/gene_census_TSG.txt \
	--oncoKB-curated /scratch1/fs1/bolton/UKBB/germline/already/all_curated_genes_v2.0.tsv \
	-p /scratch1/fs1/bolton/UKBB/germline/already/pd_table_kbreview_bick_trunc2.txt \
	--pan-myeloid /scratch1/fs1/bolton/UKBB/germline/already/panmyeloid_variant_counts.vep.annotated.vcf.tsv \
	--bolton-bick-vars /scratch1/fs1/bolton/UKBB/germline/already/bick.bolton.vars3.txt \
	--mut2-bick /scratch1/fs1/bolton/UKBB/germline/already/topmed.n2.mutation.c.p.txt \
	--mut2-kelly /scratch1/fs1/bolton/UKBB/germline/already/kelly.n2.mutation.c.p.txt \
	--matches2 /scratch1/fs1/bolton/UKBB/germline/already/matches.2.c.p.txt

./alreadyRan.R \
	-i $input \
	--output $output \
	--target-length $target_length \
	--col-order $col_order \
	-T $T \
	--oncoKB-curated $oncoKB_curated \
	-p $p \
	--pan-myeloid $pan_myeloid \
	--bolton-bick-vars $bolton_bick_vars \
	--mut2-bick $mut2_bick \
	--mut2-kelly $mut2_kelly \
	--matches2 $matches2

for eid in $(cat failed.txt); do
	bsub -n8 -oo long.log -g /bwileytest3 -G compute-timley -q general -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a "docker(kboltonlab/annotate_wes_ch:3.0)" /bin/bash -c "/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/germline/NSAMPS/run.sh \
	/scratch1/fs1/bolton/UKBB/lung/pass/pon_pass_myeloid/${eid}_23153_0_0.pon.pass.myeloid.tsv \
	/scratch1/fs1/bolton/UKBB/lung/pass/newest/${eid}_23153_0_0.newest.tsv \
	39652397 \
	/scratch1/fs1/bolton/UKBB/germline/already/gene_census_TSG.txt \
	/scratch1/fs1/bolton/UKBB/germline/already/all_curated_genes_v2.0.tsv /scratch1/fs1/bolton/UKBB/germline/already/pd_table_kbreview_bick_trunc2.txt \
	/scratch1/fs1/bolton/UKBB/germline/already/panmyeloid_variant_counts.vep.annotated.vcf.tsv /scratch1/fs1/bolton/UKBB/germline/already/bick.bolton.vars3.txt \
	/scratch1/fs1/bolton/UKBB/germline/already/topmed.n2.mutation.c.p.txt \
	/scratch1/fs1/bolton/UKBB/germline/already/kelly.n2.mutation.c.p.txt \
	/scratch1/fs1/bolton/UKBB/germline/already/matches.2.c.p.txt"
done

(head -n1 ../newest/${eid}_23153_0_0.newest.tsv; awk -F '\t' '$2=="102248627" {print}' ../newest/${eid}_23153_0_0.newest.tsv) | TRANSPOSE

if [[ $project == "Lung" ]]; then
	length=39652397
else
	length=38997831
fi
printf "$project\t$length\n"

 

cat germline.already.ran.download2.txt | tail -n+2 | while read index sub eid project file; do
	if [[ $project == "Lung" ]]; then
		length=39652397
	else
		length=38997831
	fi
	printf "$project\t$length\n"
done
	bsub -oo $eid.log -g /bwileytest3 -G compute-timley -q general -M 64G -R 'select[mem>64G] span[hosts=1] rusage[mem=64G]' -a "docker(kboltonlab/annotate_wes_ch:3.0)" ./run.sh /scratch1/fs1/bolton/UKBB/germline/already/inputs/${eid}_23153_0_0.final.tsv /scratch1/fs1/bolton/UKBB/germline/already/outputs/${eid}_23153_0_0.updated.tsv $length /scratch1/fs1/bolton/UKBB/germline/already/1000051_23153_0_0.columns.txt /scratch1/fs1/bolton/UKBB/germline/already/gene_census_TSG.txt /scratch1/fs1/bolton/UKBB/germline/already/all_curated_genes_v2.0.tsv /scratch1/fs1/bolton/UKBB/germline/already/pd_table_kbreview_bick_trunc2.txt /scratch1/fs1/bolton/UKBB/germline/already/panmyeloid_variant_counts.vep.annotated.vcf.tsv /scratch1/fs1/bolton/UKBB/germline/already/bick.bolton.vars3.txt /scratch1/fs1/bolton/UKBB/germline/already/topmed.n2.mutation.c.p.txt /scratch1/fs1/bolton/UKBB/germline/already/kelly.n2.mutation.c.p.txt /scratch1/fs1/bolton/UKBB/germline/already/matches.2.c.p.txt
done

head -n1 sample_files/4898186_23153_0_0.final.tsv | awk -F'\t' '{$0=gensub(/\s*\S+/,"",242)}1' OFS='\t' | awk -F'\t' '$278="eid" {$0=gensub(/\s*\S+/,"",267)}1' OFS='\t'

for sample in $(cat 65.txt | head -n1); do
	sample=${sample}_23153_0_0
	tail -n+2 $sample.final.tsv | head 
	head -n2 $sample.final.tsv | tail -n1 | awk -F'\t' '{$0=gensub(/\s*\S+/,"",242)}1' OFS='\t' | awk -F'\t' -v sample=$sample -F'\t' '$278=sample {$0=gensub(/\s*\S+/,"",267)}1' OFS='\t' | TRANSPOSE

head -n10 $sample.final.tsv | awk -F'\t' '{$0=gensub(/\s*\S+/,"",268)}1' OFS='\t' | awk -F'\t' -v sample=$sample -F'\t' '{$0=gensub(/\s*\S+/,"",242)}1' OFS='\t' | awk -F'\t' '{print NF}'

head -n10 $sample.final.tsv | awk -F'\t' '{$0=gensub(/\s*\S+/,"",3)}1' OFS='\t' | awk -F'\t' '{print NF}'

sed -E 's/(\s+)?\S+//268g' $sample.final.tsv | sed -E 's/(\s+)?\S+//242g' > $sample.final.new.tsv
head -n2 $sample.final.new.tsv | tail -n1 | tr '\t' '\n' | wc -l

awk -F'\t' '{$242=$268=""}1' OFS='\t' $sample.final.tsv | head -n3 | TRANSPOSE
head -n2 $sample.final.tsv | tail -n1 | awk -F'\t' '{$242=$268="";gsub("\t+","\t",$0)}1' OFS='\t'


head -n1 $sample.final.tsv > 
cut -f1-241,243-268,270-278 -d$'\t' $sample.final.tsv

head -n1 4898186_23153_0_0.final.tsv | awk -F'\t' '{$0=gensub(/\s*\S+/,"",269)}1' OFS='\t' | awk -F'\t' '$278="eid" {$0=gensub(/\s*\S+/,"",242)}1' OFS='\t' | TRANSPOSE | sha256sum


head -n1 combined.pon.germline2.tsv | TRANSPOSE | sha256sum

head -n1 $sample.final.tsv | cut -f1-241,243-268,270-278 -d$'\t' | awk -F'\t' '$277="eid" {print}' OFS='\t' | TRANSPOSE | sha256sum
head -n1 $sample.final.tsv | cut -f1-241,243-268,270-278 -d$'\t' | awk -F'\t' '$277="eid" {print}' OFS='\t' > new.header.65.txt
for file in $(ls *.tsv); do
	head -n1 $file | sha256sum
done

for sample in $(cat 65.txt); do
	sample=${sample}_23153_0_0
	(cat new.header.65.txt; tail -n+2 $sample.final.tsv | cut -f1-241,243-268,270-278 -d$'\t' | awk -F'\t' -v eid=$sample '$277=eid {print}' OFS='\t') > $sample.final.new.tsv
done



mutect_passed combined.pon.germline2.tsv germline2.ch_pd2.pon.nsamples.removed.NOT.stream_gene_variant.tsv # 1657 lines