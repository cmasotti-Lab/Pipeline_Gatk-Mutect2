SCRATCH60="/home/scratch60/vlira_07Fev2023/"
SCRATCH45="/home/scratch45/vlira_27jan/"

SAMPLES_DIR="$SCRATCH60/preprocessing_FINAL_result/"
SAMPLES_FILE=$(find "$SAMPLES_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam')
#SAMPLES_FILE=$(find "/home/scratch45/vlira_17dez/DownSample_ROP/Result3_DownSample-ROP/DownsampleSam/" -maxdepth 1 -mindepth 1  -name 'ROP-65-*.dedup.tags.bqsr.bam')

REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/"
GATK="$SCRATCH60/tools/gatk-4.3.0.0/./gatk"
TARGET="$SCRATCH60/references/xgen-exome-research-panel-v2-targets-hg38.bed"
PON="$SCRATCH45/PanelOfNormals/PoN.Mutect2_COVID/PoN.100COVID.vcf.gz"
GNOMAD="$SCRATCH60/references/af-only-gnomad.hg38.vcf.gz"
ANNOVAR="/home/scratch60/vlira_07Fev2023/tools/annovar/table_annovar.pl"
ANNOVAR_DB="/home/scratch60/vlira_07Fev2023/humandb/"
CROSS_REFERENCE="$SCRATCH60/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

OUTPUT_DIR="Result_Mutect2.PoN.100COVID"
mkdir $OUTPUT_DIR

find "$SAMPLES_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam' > $OUTPUT_DIR/samples.list
#find "/home/scratch45/vlira_17dez/DownSample_ROP/Result3_DownSample-ROP/DownsampleSam/" -maxdepth 1 -mindepth 1  -name 'ROP-65-*.dedup.tags.bqsr.bam' > $OUTPUT_DIR/samples.list

TIME_FILE="$OUTPUT_DIR/$OUTPUT_DIR.log"

export OUTPUT_DIR
export TIME_FILE
export PON
export GNOMAD
export TARGET
export GATK
export REF_FASTA
export ANNOVAR
export ANNOVAR_DB
export CROSS_REFERENCE



stage_Mutect2 (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

    $GATK --java-options "-Xmx20G" Mutect2 \
          -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
          -L $TARGET \
          -I $SAMPLE \
          --germline-resource $GNOMAD \
          --panel-of-normals $PON   \
          --f1r2-tar-gz $OUTPUT_DIR/Mutect2/$NAME.f1r2.tar.gz \
          -bamout $OUTPUT_DIR/Mutect2/$NAME.bamout.bam \ 
          -O $OUTPUT_DIR/Mutect2/$NAME.unfiltered.vcf.gz 2> $OUTPUT_DIR/Mutect2/$NAME.unfiltered.log
}
export -f stage_Mutect2

stage_LearnReadOrientationModel (){
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_LearnReadOrientationModel<<<<<<" >> $TIME_FILE
  date >> $TIME_FILE

  local SAMPLE_READORIENTATION=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.f1r2.tar.gz')
  local SAMPLE_F1R2=$(echo $SAMPLE_READORIENTATION| sed 's/\s/ -I  /g')

  $GATK --java-options "-Xmx20G" LearnReadOrientationModel \
        -I $SAMPLE_F1R2 \
        -O "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz 2> $OUTPUT_DIR/LearnReadOrientationModel/read-orientation-model.tar.gz.log
} 


stage_GetPileupSummaries (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando GetPileupSummaries para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK --java-options "-Xmx100G" GetPileupSummaries \
        -I $SAMPLE  \
        -V $GNOMAD \
        -L $GNOMAD \
        -O $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.table 2> $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.log
}
export -f stage_GetPileupSummaries

stage_CalculateContamination (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando CalculateContamination para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE
 
  $GATK CalculateContamination \
        -I $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.table \
        -tumor-segmentation $OUTPUT_DIR/CalculateContamination/$NAME.segments.table \
        -O $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.table 2> $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.log
}
export -f stage_CalculateContamination


stage_FilterMutectCalls (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE
 
  $GATK FilterMutectCalls \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -V $OUTPUT_DIR/Mutect2/$NAME.unfiltered.vcf.gz \
        --tumor-segmentation $OUTPUT_DIR/CalculateContamination/$NAME.segments.table \
        --contamination-table $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.table \
        --stats $OUTPUT_DIR/Mutect2/$NAME.unfiltered.vcf.gz.stats \
        --ob-priors "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        -O $OUTPUT_DIR/FilterMutectCalls/$NAME.filtered.vcf.gz  2> $OUTPUT_DIR/FilterMutectCalls/$NAME.filtered.vcf.log
}
export -f stage_FilterMutectCalls



left_normalization () {
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando Normalization para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

    bcftools norm -m-both -O z -o $OUTPUT_DIR/left_normalization/$NAME.norm_Step1.vcf.gz $OUTPUT_DIR/FilterMutectCalls/$NAME.filtered.vcf.gz 2> $OUTPUT_DIR/left_normalization/$NAME.norm_Step1.log
    bcftools norm -O z -f $REF_FASTA/Homo_sapiens_assembly38.fasta -o $OUTPUT_DIR/left_normalization/$NAME.norm_Step2.vcf.gz $OUTPUT_DIR/left_normalization/$NAME.norm_Step1.vcf.gz 2> $OUTPUT_DIR/left_normalization/$NAME.norm_Step2.log
    bcftools index $OUTPUT_DIR/left_normalization/$NAME.norm_Step2.vcf.gz

}
export -f left_normalization


annotation (){

  # Merge all genotyped samples to get a single VCF/BCF using bcftools merge
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando Vcftools para todas juntar os vcf das amostras<<<<<<" >> $TIME_FILE
  date >> $TIME_FILE
  
  local SAMPLE_CALL_GENO=$(find "$OUTPUT_DIR"/left_normalization/ -maxdepth 1 -mindepth 1  -name '*.norm_Step2.vcf.gz')

  #bcftools merge -m id -O z -o $OUTPUT_DIR/annotation/mutect.merged.vcf $SAMPLE_CALL_GENO 2> $OUTPUT_DIR/annotation/mutect.merged.log
  #bcftools index "$OUTPUT_DIR/annotation/mutect.merged.vcf"

  #vcf-merge $SAMPLE_CALL_GENO > $OUTPUT_DIR/annotation/mutect.merged2.vcf 2> $OUTPUT_DIR/annotation/mutect.merged2.log
  #bcftools index "$OUTPUT_DIR/annotation/mutect.merged2.vcf"


  bcftools norm -m-both -O z -o $OUTPUT_DIR/annotation/mutect.merged.norm_Step1.vcf.gz $OUTPUT_DIR/annotation/mutect.merged.vcf 2> $OUTPUT_DIR/annotation/mutect.merged.norm_Step1.log
  bcftools norm -O z -f $REF_FASTA/Homo_sapiens_assembly38.fasta -o $OUTPUT_DIR/annotation/mutect.merged.norm_Step2.vcf.gz $OUTPUT_DIR/annotation/mutect.merged.norm_Step1.vcf.gz 2> $OUTPUT_DIR/annotation/mutect.merged.norm_Step2.log
  bcftools index $OUTPUT_DIR/annotation/mutect.merged.norm_Step2.vcf.gz


  # # annovar
   echo "" >> $TIME_FILE
   echo ">>>>>> Executando annovar para todas juntar os vcf das amostras<<<<<<" >> $TIME_FILE
  
   $ANNOVAR  --vcfinput $OUTPUT_DIR/annotation/mutect.merged.norm_Step2.vcf.gz $ANNOVAR_DB -buildver hg38 --remove \
   --protocol refGene,avsnp147,gnomad_exome,abraom,cosmic95,icgc28,dbnsfp42a  \
   --operation gx,f,f,f,f,f,f --arg '-splicing 5',,,,,, --polish \
   --xreffile $CROSS_REFERENCE --otherinfo --thread 15 --outfile $OUTPUT_DIR/annotation/annovar.norm 2> $OUTPUT_DIR/annotation/annovar.norm.log

   sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/annovar.norm.hg38_multianno.vcf| sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/annovar.norm.hg38_multianno.correct.vcf 

  date >> $TIME_FILE
  

  date >> $TIME_FILE
  echo ">>>>>> Executando SnpSift para todas juntar os vcf das amostras<<<<<<" >> $TIME_FILE
  
  java -jar -Xmx50G /home/scratch60/vlira_07Fev2023/tools/snpEff/SnpSift.jar extractFields  "$OUTPUT_DIR"/annotation/annovar.norm.hg38_multianno.correct.vcf \
    -e . "CHROM" "POS" "ID" "REF" "ALT" "FILTER" "AC" "AN" "DP" "Func.refGene" "Gene.refGene" "GeneDetail.refGene" "ExonicFunc.refGene" "AAChange.refGene" \
    "COSMIC_Census_Gene.refGene" "Role_in_Cancer.refGene" "Translocation_Partner.refGene" "Therapeutic_Agents.refGene" "Cancer_Syndromes.refGene" \
    "panel.refGene" "gnomAD_exome_ALL" "abraom_freq" "cosmic95" "ICGC_Id" "GEN[*].GT"  > "$OUTPUT_DIR"/Final.mutect2.txt 2> "$OUTPUT_DIR"/Final.mutect2.log

  echo "" >> $TIME_FILE

}
export -f annotation



echo "                                                     >>>>>> Starting Pipeline  to Run GATK-MUTECT2 <<<<<<" >> $TIME_FILE
date >> $TIME_FILE

#mkdir $OUTPUT_DIR/Mutect2/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P10 bash -c 'stage_Mutect2  "$@"' 'stage_Mutect2'

#mkdir $OUTPUT_DIR/LearnReadOrientationModel/
#stage_LearnReadOrientationModel

#mkdir $OUTPUT_DIR/GetPileupSummaries/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P2 bash -c 'stage_GetPileupSummaries  "$@"' 'stage_GetPileupSummaries'

#mkdir $OUTPUT_DIR/CalculateContamination/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P5 bash -c 'stage_CalculateContamination  "$@"' 'stage_CalculateContamination'

#mkdir $OUTPUT_DIR/FilterMutectCalls/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P5 bash -c 'stage_FilterMutectCalls  "$@"' 'stage_FilterMutectCalls'

#mkdir $OUTPUT_DIR/left_normalization/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P15 bash -c 'left_normalization  "$@"' 'left_normalization'

mkdir $OUTPUT_DIR/annotation/
annotation



echo "" >> $TIME_FILE
echo ">>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE
