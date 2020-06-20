
cd /home/boeva/NAS_public/data/projects/EnhancerPromoterInteractions/TunaProject

narrowPeakFile=testData/ATAC.GM12878.50Kcells.rep1_peaks.narrowPeak 
ATACseqThreshold=4.99


ATACseqThreshold=10
bigwigFile=testData/ATAC.GM12878.50Kcells.rep1.wig.bw
outputDir=testOutput
genomeFa=../../../annotations/Human/hg19/hg19.fa
hocomocoThresholds=annotationFiles/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt
hocomocoPWMdir=annotationFiles/HOCOMOCOv11_core_pwm
hocomocoPWMdescriptionFile=annotationFiles/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv
coreCoordinateFile=annotationFiles/core.txt

. code/ValentinaCode/motifProfiler.sh -i $narrowPeakFile -t $ATACseqThreshold -w $bigwigFile -o $outputDir  -a $genomeFa -h $hocomocoThresholds -d $hocomocoPWMdir -b $hocomocoPWMdescriptionFile -c $coreCoordinateFile -g /home/boeva/NAS_public/data/annotations/Human/hg19/gencode.v19.annotation.gtf
