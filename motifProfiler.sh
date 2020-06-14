#/usr/bin/env bash

while getopts i:o:t:a:h:d:b:w:c:g: option
do
case "${option}"
in
o) OUTDIR=${OPTARG};;
i) InputFile=${OPTARG};;
t) threshold=${OPTARG};;
a) GENOMEFA=${OPTARG};;
h) HOCOMOCOthresholds=${OPTARG};;
d) PWMDIR=${OPTARG};;
b) HOCOMOCOtoTF=${OPTARG};;
w) BIGWIGFILE=${OPTARG};;
c) CoreCoordinates=${OPTARG};;
g) GencodeFile=${OPTARG};;


esac
done

display_usage() { 
echo -e "\nUsage: pipeline.sh -i InputDir -o OutputDir -t threshold -a genome.fa -h HOCOMOCOthresholds -d folderWithPWMs -b HOCOMOCOfullannotations -w bigwigFile -c coreCoordinateFile -g geneCodeGeneAnnotationFile\n"
echo -e "\nExample: pipeline.sh -i ATAC.XXX.rep1_peaks.narrowPeak -o /home/user/outputDir -t 4.99 -a hg19.fa -h HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt -d HOCOMOCOv11_core_pwm -b HOCOMOCOv11_full_annotation_HUMAN_mono.tsv -w ATAC.XXX.rep1_peaks.bw -c core.txt -g gencode.v19.annotation.gtf\n"
echo "You can get all necessary HOCOMOCO files (PWMs and thresholds) here: "
echo "PWMs: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwm_HUMAN_mono.tar.gz "
echo "Thresholds: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt "
echo "HOCOMOCO to TF: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv " 
} 

#get path to the directory with .py code files
SCRIPTPATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ ! -f "$HOCOMOCOtoTF" ]
then
echo -e "\nERROR: You have to set the file with HOCOMOCO full annotation (e.g. hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv)"
display_usage
exit 1
fi

if [ ! -f "$HOCOMOCOthresholds" ]
then
echo -e "\nERROR: You have to set the file with HOCOMOCO thresholds (e.g. http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt)"
display_usage
exit 1
fi

if [ ! -d "$PWMDIR" ]
then
echo -e "\nERROR: You have to set the directory with HOCOMOCO PWM (e.g. get it here: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwm_HUMAN_mono.tar.gz)"
display_usage
exit 1
fi

if [ ! -f "$GENOMEFA" ]
then
echo -e "\nERROR: You have to set the .fasta file with the genomic sequences"
display_usage
exit 1
fi

if [ ! -f "$InputFile" ]
then
echo -e "\nERROR: You have to set the .bed file with ATAC-seq peaks (e.g. narrowPeak file from the HMCan output)"
display_usage
exit 1
fi

if [ -z "$threshold" ]
then
echo -e "\nERROR: You have to set the score threshold for the .bed file with ATAC-seq peaks. Set 0 if you do not wish to perform filtering by score"
display_usage
exit 1
fi

if [ -z "$OUTDIR" ]
then
echo -e "\nERROR: You have to set the output directory"
display_usage
exit 1
fi

echo "TFhubsFinder -- input parameters:"
echo "Input file: $InputFile"
echo "Threshold: $threshold"
echo "Output Directory: $OUTDIR"
echo "FASTA file with the genome: $GENOMEFA"
echo "File with HOCOMOCO thresholds: $HOCOMOCOthresholds"
echo "File with full HOCOMOCO annotations: $HOCOMOCOtoTF"

echo "..Will create the output directory if it does not exist"

mkdir $OUTDIR || true

rm -rf $OUTDIR/FASTA/ || true

mkdir $OUTDIR/FASTA/ || true

rm -rf $OUTDIR/TFBS/ || true

mkdir $OUTDIR/TFBS/ || true

output_fasta=$OUTDIR/FASTA/
TFBS_output=$OUTDIR/TFBS/


file=$(basename $InputFile)
bed6_output=$OUTDIR/${file}.bed6.bed
fasta_output=$output_fasta/${file}.fa
expressedGenes=$OUTDIR/peaksInExpressedGenes.bed


echo "1. filtering $file with threshold $threshold..."
cat $InputFile | awk "BEGIN {OFS=\"\t\"} {if (\$5 > $threshold) {print \$1,\$2,\$3,\$4,\$4,\$10}}" > $bed6_output # Keeping just 5 columns; Change if needed more


echo "1.2. determining promoter regions and finding expressed TFs based on ATAC-seq peaks..."
#create a temporary file with promoter region coordinates (+- 750 bp around TSS)
tail -n +8 $GencodeFile | awk '
BEGIN {OFS="\t"}{
	if ($7 == "+") {print $1,$3,$4,$4+1,$18}

	else if ($7 == "-") {print $1,$3,$5-1,$5,$18}
}' | grep -wE "(transcript)" | awk 'BEGIN {OFS="\t"}{print $1,$2,$3,$4,$5}' | awk 'BEGIN {OFS="\t"}{print $1,$3-750,$4+750,"TSS",$5}' | sed "s/[\";]//g" > $OUTDIR/promoter_regions.bed

cat $bed6_output $OUTDIR/promoter_regions.bed | sort -k1,1 -k2,2n | python3 $SCRIPTPATH/intersect.py TSS > $expressedGenes
# output: chrom, start, end, peak, gene


echo "2. making fasta files..."
echo $GENOMEFA
bedtools getfasta -fi $GENOMEFA -bed $bed6_output -name -fo $fasta_output
sed -i 's/:.*$//g' $fasta_output
echo "done!"
  
echo "3. looking for motif hits..."
date

python3 $SCRIPTPATH/TFBS_finder.py $HOCOMOCOthresholds $PWMDIR $fasta_output regions $HOCOMOCOtoTF
    
echo "done!"
date
	
echo "4. calculating motif statistics..."
echo "python3 $SCRIPTPATH/Crossref_v3.py $bed6_output $TFBS_output $CoreCoordinates $HOCOMOCOthresholds $HOCOMOCOtoTF $BIGWIGFILE"

python3 $SCRIPTPATH/Crossref_v3.py $bed6_output $TFBS_output $CoreCoordinates $HOCOMOCOthresholds $HOCOMOCOtoTF $BIGWIGFILE $expressedGenes

echo "done!"
date
