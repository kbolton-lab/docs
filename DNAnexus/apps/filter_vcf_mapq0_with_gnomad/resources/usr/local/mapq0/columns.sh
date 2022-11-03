10
6 4b452252b34db4c00c26c4a5c88bbbe988f6ee2a337dd6295c25026608080eaa  - # 271
130 60b363fd16bd03badd2929d935da5ab363319aabdeae43a18f661688174de757  - # 273
3648 c6c626d02949e74ecd6a81cd9a370b2849f40bc359cb2601e7d795cc6f1d1455  - # 273
150 dfe50f87e8967acb557fa0229ba120241ca7d3ec42b249abe7b2e48d020a1b25  - # 271

for file in $(grep -v c6c626d02949e74ecd6a81cd9a370b2849f40bc359cb2601e7d795cc6f1d1455 out.txt | awk '{print $3}' | cut -d. -f1 | head); do
   ./columns.R $file.final.tsv columns.txt $file.final.updated.tsv
done

for sample in $(ls *.updated* | cut -d. -f1); do
   one=$(awk -F'\t' '{print $1,$2,$3,$4}' $sample.final.tsv | sha256sum | awk '{print $1}')
   two=$(awk -F'\t' '{print $1,$2,$3,$4}' $sample.final.updated.tsv | sha256sum | awk '{print $1}')
   if [[ one==2 ]]; then 
      echo $sample.final.tsv: yes
      mv $sample.final.updated.tsv $sample.final.tsv
   else
      echo no
   fi
done

11
8 4b452252b34db4c00c26c4a5c88bbbe988f6ee2a337dd6295c25026608080eaa  -
135 60b363fd16bd03badd2929d935da5ab363319aabdeae43a18f661688174de757  -
3662 c6c626d02949e74ecd6a81cd9a370b2849f40bc359cb2601e7d795cc6f1d1455  -
188 dfe50f87e8967acb557fa0229ba120241ca7d3ec42b249abe7b2e48d020a1b25  -

awk '$2 != 273 && $1 !="c6c626d02949e74ecd6a81cd9a370b2849f40bc359cb2601e7d795cc6f1d1455"' out.txt


for file in $(ls *.tsv); do 
   wc=$(head -n1 $file | tr '\t' '\n' | wc -l)
   hash=$(head -n1 $file | sha256sum | awk '{print $1}')
   printf "$hash\t$wc\t$file\n"
done > out.txt




#######
#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Missing", call.=FALSE)
}

input = read.table(args[1], header=T, quote="", comment.char="", sep="\t")
columns = read.table(args[2], header=T)$V1
if (is.null(input$Mutect2_ROQ)) {
   input$Mutect2_ROQ = NA
}
if (is.null(input$Mutect2_samtools_DP)) {
   input$Mutect2_samtools_DP = NA
}
input = input[,columns]
write.table(input, args[3], sep="\t", quote=F, row.names=F)