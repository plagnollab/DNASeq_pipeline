java=/share/apps/jdk1.7.0_25/bin/java
vcftools=/cluster/project8/vyp/vincent/Software/vcftools_0.1.11/bin/vcftools
#GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
GATK=/cluster/project8/vyp/vincent/Software/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
bundle=/scratch2/vyp-scratch2/reference_datasets/GATK_bundle
target=/cluster/project8/vyp/exome_sequencing_multisamples/target_region/data/merged_exome_target_cleaned.bed
Rfilter=/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/trio/filter_de_novo_file.R
Rannot=/cluster/project8/vyp/vincent/Software/pipeline/GATK_v2/trio/final_annotate_file.R

fasta=/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta

exec=no
Rbin=/share/apps/R-3.0.2/bin/R
script=default
#argument needed: code, ped file, oFolder
##returns 4 different files per trio

until [ -z "$1" ]; do
        # use a case statement to test vars. we always test $1 and shift at the end of the for block.
    case $1 in
	--script)
	    shift
	    script=$1;;
	--makeVCF)
	    shift
	    makeVCF=$1;;
	--BAMlist)
	    shift
	    BAMlist=$1;;
        --output)
	    shift
	    output=$1;;
	--pedfile)
	    shift
	    pedFile=$1;;
	-* )
            echo "Unrecognized option: $1"
            exit 1;;
    esac
    shift
    if [ "$#" = "0" ]; then break; fi
done

for folder in temp cluster cluster/submission cluster/out cluster/error ${oFolder}; do
    if [ ! -e $folder ]; then mkdir $folder; fi
done


echo "Running awk script"
awk '{print $2}' $pedFile > ${output}_samples.tab
echo "Running awk script"

if [[ "$script" == "default" ]]; then 
    code=`basename $output`
    script=cluster/submission/${code}.sh
fi

echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -pe smp 1
#$ -l tmem=10G,h_vmem=10G
#$ -l h_rt=72:0:0
#$ -V
#$ -R y
" >  $script


if [[ "$makeVCF" == "yes" ]]; then
    
    if [ ! -e $BAMlist ]; then echo "The file $BAMlist does not exist"; exit; fi
    tmpDir=/scratch0/test

    echo "
if [ ! -e $tmpDir ]; then mkdir $tmpDir; fi

$java -Djava.io.tmpdir=$tmpDir -Xmx8g -jar $GATK -R $fasta -T HaplotypeCaller -I $BAMlist \
 --dbsnp ${bundle}/dbsnp_137.b37.vcf \
     -stand_call_conf 30.0 \
     -stand_emit_conf 10.0 \
     -L $target \
     --activeRegionExtension 100 \
     -o ${output}.vcf

rm -rf $tmpDir
" >> $script

fi


echo "
echo \"Run GATK to compute positions that seem Mendel incompatible\"
$java -Xmx2g -jar $GATK -R $fasta -T PhaseByTransmission --DeNovoPrior 0.0005 -V ${output}.vcf -ped $pedFile -o ${output}_phased.vcf --MendelianViolationsFile ${output}_noMendel.tab

echo \"Filter the GATK output file to extract the proper de novo variants\"    
$Rbin CMD BATCH --no-save --no-restore --input.file=${output}_noMendel.tab --output.file=${output}_deNovo.tab ${Rfilter} ${output}.Rlog 

echo \"Now extract from the main VCF the positions of interest\" 
$vcftools --positions ${output}_deNovo.tab.positions --vcf ${output}.vcf --out ${output}_deNovo --recode

echo \"Now annotate the variants using ANNOVAR\"
/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/convert2annovar.pl --allallele -format vcf4 --includeinfo ${output}_deNovo.recode.vcf > ${output}_db

/cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/summarize_annovar_VP.pl -ver1000g 1000g2012apr -verdbsnp 137 -veresp 6500si -alltranscript -buildver hg19 --genetype ensgene --remove ${output}_db /cluster/project8/vyp/vincent/Software_heavy/annovar_Feb2013/humandb_hg19/

echo \"Final annotation of de novo variants\"
$Rbin CMD BATCH --no-save --no-restore --input.file=${output}_db.genome_summary.csv --output.file=${output}_final_de_novo_annotated.csv ${Rannot} ${output}_annot.Rlog 

" >> $script



echo "Submission script in $script"

