#Pipeline atualizado em 14/06/2024
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado
#2. Alterado o diretorio onde estava o PON na scratch45 para /home/users/vlira/PanelOfNormals/PoN.100COVID.vcf.gz
#3. Agora o script recebe 2 agumentos: 1- lista de amostras; 2- diretorio SCRATCH
#4. Utualizado os databases: hg38_cosmic98_coding,hg38_avsnp150,hg38_clinvar_20220320, hg38_gnomad40_exome

export WD="/home/scratch90/vlira_13may2024/"

export MEM=200
export JOBS=5

export DATA=$(date "+%F")
export DATA="2024-06-14"  # EDITE AQUI SE QUISER USAR UMA PASTA DE UMA DATA ESPECIFICA
export OUTPUT_DIR=${WD}"/Result_Mutect2.TOY.${DATA}"

export REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/"
export GATK="$WD/tools/gatk-4.6.0.0/./gatk"
export PICARD="java -jar ${WD}/tools/picard-3.2.0/picard.jar"

export TARGET="$WD/references/xgen-exome-research-panel-v2-targets-hg38.bed"
export PON="/home/users/vlira/PanelOfNornal/PoN.100COVID.vcf.gz"
export GNOMAD="$WD/references/af-only-gnomad.hg38.vcf.gz"
#export GNOMAD="$WD/references/af-only-gnomad.SABE1171.Abraom.hg38.vcf.gz"
export ANNOVAR="$WD/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$WD/humandb/"
export CROSS_REFERENCE="$WD/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

export INDEL_KNOWN="/home/projects2/LIDO/molPathol/oncoseek/nextseq/references/Mills_and_1000G_gold_standard.indels.hg38.vcf"
export DBSNP="$WD/references/dbsnp_151.hg38.vcf.gz"
export GNOMAD2="$WD/references/af-only-gnomad.SABE1171.Abraom.hg38.new.vcf"

export TIME_FILE="$OUTPUT_DIR.log"

mkdir $OUTPUT_DIR
#$(find "$WD/preprocessing_FINAL_result/" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam.cram') > $OUTPUT_DIR/samples.list.cram

ls $WD/preprocessing_FINAL_result/*.dedup.tags.bqsr.bam.cram > $OUTPUT_DIR/samples.cram
export SAMPLE_LIST="$OUTPUT_DIR/samples.cram"

ls $WD/preprocessing_FINAL_result/*.dedup.tags.bqsr.bam > $OUTPUT_DIR/samples.list.bam
export SAMPLE_LIST_BAM="$OUTPUT_DIR/samples.list.bam"

ls $WD/preprocessing_TOY_result/*.bam > $OUTPUT_DIR/samples.TOY.bam
export SAMPLE_TOY_BAM="$OUTPUT_DIR/samples.TOY.bam"


#######################################
## PRE-PROCECING
#######################################
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

#1. Raw Unmapped Reads
#2. Map to Reference
#3. Raw Mapped Reads
#4. MarkDuplicates

#5. Recalibreta Base Qality Score
   #  SetNmMdAndUqTags

stage_SetNmMdAndUqTags (){
  local SAMPLE=$1
  # SAMPLE="TOY-56-DownsampleSam.bam"
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_SetNmMdAndUqTags para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  ${PICARD} SetNmMdAndUqTags \
    R=$REF_FASTA/Homo_sapiens_assembly38.fasta \
    I=$SAMPLE \
    CREATE_INDEX=true \
    O=$WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bam" \
    2> $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bam.log"
}
export -f  stage_SetNmMdAndUqTags

echo "            >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " > $TIME_FILE
mkdir $OUTPUT_DIR/preprocessing_TOY_result/
xargs -a ${SAMPLE_TOY_BAM} -t -n1 -P${JOBS} bash -c 'stage_SetNmMdAndUqTags  "$@"' 'stage_SetNmMdAndUqTags'

#   "GATK BaseRecalibrator|PrintReads"

stage_base_recalibration(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_base_recalibration para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  $GATK --java-options "-Xmx${MEM}G" BaseRecalibrator \
    -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
    -I $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bam" \
    -known-sites ${INDEL_KNOWN} \
    -known-sites ${DBSNP} \
    -known-sites ${GNOMAD} \
    -L ${TARGET} \
    -O $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.table" \
    > $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.table.log" \
    2> $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.table.log2"

  $GATK --java-options "-Xmx${MEM}G" ApplyBQSR \
    -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
    -I $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bam" \
    --bqsr-recal-file $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.table"\
    -O $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.bam" \
    > $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.bam.log" \
    2> $WD/preprocessing_TOY_result/"${NAME%.*}.dedup.tags.bqsr.bam.log2"
}
export -f  stage_base_recalibration

mkdir $OUTPUT_DIR/preprocessing_TOY_result/
xargs -a ${SAMPLE_TOY_BAM} -t -n1 -P${JOBS} bash -c 'stage_base_recalibration  "$@"' 'stage_base_recalibration'



stage_Mutect2 (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

    $GATK --java-options "-Xmx${MEM}G" Mutect2 \
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
  echo ">>>>>> STAGE_LearnReadOrientationModel<<<<<<  $(date) " >> $TIME_FILE

  local SAMPLE_READORIENTATION=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.f1r2.tar.gz')
  local SAMPLE_F1R2=$(echo $SAMPLE_READORIENTATION| sed 's/\s/ -I  /g')

  $GATK --java-options  "-Xmx${MEM}G"  LearnReadOrientationModel \
        -I $SAMPLE_F1R2 \
        -O "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz 2> $OUTPUT_DIR/LearnReadOrientationModel/read-orientation-model.tar.gz.log
} 
export -f stage_LearnReadOrientationModel

stage_GetPileupSummaries (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando GetPileupSummaries para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  $GATK --java-options -Xmx500G GetPileupSummaries \
        -I $SAMPLE  \
        -V $GNOMAD \
        -L $GNOMAD \
        -O $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.table 2> $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.log
}
export -f stage_GetPileupSummaries

for i in `cat ${SAMPLE_LIST_BAM}`; do
 stage_GetPileupSummaries $i
done

################################## TESTE ##################################



$GATK --java-options -Xmx500G  GetPileupSummaries \
      -I /home/scratch90/rtorreglosa_25jun2024/preprocessing_READ_result/ROP-98-ExC85-xgenV2_S66.dedup.tags.bam \
      -V $GNOMAD \
      -L $GNOMAD \
      -O $OUTPUT_DIR/GetPileupSummaries.erro/TESTE3.getpileupsummaries.table 2> $OUTPUT_DIR/GetPileupSummaries.erro/TESTE3.getpileupsummaries.log

$GATK --java-options -Xmx500G  GetPileupSummaries \
      -I /home/scratch90/rtorreglosa_25jun2024/preprocessing_READ_result/ROP-98-ExC85-xgenV2_S66.dedup.tags.bqsr.bam \
      -V $GNOMAD \
      -L $GNOMAD \
      -O $OUTPUT_DIR/GetPileupSummaries.erro/TESTE7.getpileupsummaries.table 2> $OUTPUT_DIR/GetPileupSummaries.erro/TESTE7.getpileupsummaries.log


$GATK --java-options -Xmx500G  GetPileupSummaries \
      -I /home/scratch90/vlira_13may2024/preprocessing_FINAL_result/ROP-98-ExC85-xgenV2_S66.dedup.tags.bqsr.bam \
      -V $GNOMAD \
      -L $GNOMAD \
      -O $OUTPUT_DIR/GetPileupSummaries.erro/TESTE4.getpileupsummaries.table 2> $OUTPUT_DIR/GetPileupSummaries.erro/TESTE4.getpileupsummaries.log


####################################################################
# 
#cp -r MUTECT2_ROP/Result_Mutect2.PoN.100COVID/GetPileupSummaries/ Result_Mutect2.ROP.2024-06-14/

stage_CalculateContamination (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando CalculateContamination para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  $GATK --java-options "-Xmx${MEM}G" CalculateContamination \
        -I $OUTPUT_DIR/GetPileupSummaries/$NAME.getpileupsummaries.table \
        -tumor-segmentation $OUTPUT_DIR/CalculateContamination/$NAME.segments.table \
        -O $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.table 2> $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.log
}
export -f stage_CalculateContamination


stage_FilterMutectCalls (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE
 
  $GATK --java-options "-Xmx${MEM}G" FilterMutectCalls \
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
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  # # annovar
   echo "" >> $TIME_FILE
   echo ">>>>>> Executando annovar para amostra $NAME: <<<<<< $(date)" >> $TIME_FILE
  
   $ANNOVAR  --vcfinput $OUTPUT_DIR/FILTER_PASS/$NAME.vcf.gz $ANNOVAR_DB -buildver hg38 --remove \
   --protocol refGene,avsnp150,gnomad41_exome_filt,abraom,cosmic99,icgc28,dbnsfp42a_filt,clinvar_20220320  \
   --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
   --xreffile $CROSS_REFERENCE --otherinfo --thread 5 \
   --outfile $OUTPUT_DIR/annotation/$NAME > $OUTPUT_DIR/annotation/$NAME.log 2> $OUTPUT_DIR/annotation/$NAME.log2

   sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/$NAME.hg38_multianno.vcf| sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/$NAME.hg38_multianno.correct.vcf 
}
export -f annotation



annotationXargs (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  # # annovar
   echo "" >> $TIME_FILE
   echo ">>>>>> Executando annovar para amostra $NAME: <<<<<< $(date)" >> $TIME_FILE
  
   $ANNOVAR  --vcfinput $OUTPUT_DIR/FILTER_PASS/$NAME.vcf.gz $ANNOVAR_DB -buildver hg38 --remove \
   --protocol refGene,avsnp150,gnomad41_exome_filt,abraom,cosmic99,icgc28,dbnsfp42a_filt,clinvar_20220320  \
   --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
   --xreffile $CROSS_REFERENCE --otherinfo --thread 5 \
   --outfile $OUTPUT_DIR/annotationXargs/$NAME > $OUTPUT_DIR/annotationXargs/$NAME.log 2> $OUTPUT_DIR/annotationXargs/$NAME.log2

   sed 's/\\x3b/;/g' $OUTPUT_DIR/annotationXargs/$NAME.hg38_multianno.vcf| sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotationXargs/$NAME.hg38_multianno.correct.vcf 
}
export -f annotationXargs






echo "            >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " >> $TIME_FILE

mkdir $OUTPUT_DIR/Mutect2/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'stage_Mutect2  "$@"' 'stage_Mutect2'

mkdir $OUTPUT_DIR/LearnReadOrientationModel/
stage_LearnReadOrientationModel

mkdir $OUTPUT_DIR/GetPileupSummaries/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P1 bash -c 'stage_GetPileupSummaries "$@"' 'stage_GetPileupSummaries'


 mkdir $OUTPUT_DIR/CalculateContamination/
 xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'stage_CalculateContamination  "$@"' 'stage_CalculateContamination'

 mkdir $OUTPUT_DIR/FilterMutectCalls/
 xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'stage_FilterMutectCalls  "$@"' 'stage_FilterMutectCalls'

 mkdir $OUTPUT_DIR/left_normalization/
 xargs -a ${SAMPLE_LIST_BAM}  -t -n1 -P${JOBS} bash -c 'left_normalization  "$@"' 'left_normalization'

# É PRECISSO REFINIR A ORDEM DAS ETAPAS

selectPASS(){
    local SAMPLE=$1
    NAME="${SAMPLE##*/}"
#    echo $NAME
#    echo "> Executando para Amostra: "$NAME" <<< $(date) " >> $LOG_FILE
    bcftools view -f PASS -O z $OUTPUT_DIR/left_normalization/$NAME.norm_Step2.vcf.gz > $OUTPUT_DIR/FILTER_PASS/$NAME.vcf.gz
    bcftools index $OUTPUT_DIR/FILTER_PASS/$NAME.vcf.gz 
}
export -f selectPASS

mkdir -p $OUTPUT_DIR/FILTER_PASS
#echo ">>>>>> selectPASS : <<< $(date) " >> $LOG_FILE
xargs -a $OUTPUT_DIR/samples.list.bam -t -n1 -P5 bash -c 'selectPASS  "$@"' 'selectPASS'


#zgrep -v "#" FILTER_PASS/ROP-98-ExC85-xgenV2_S66.dedup.tags.bqsr.bam.gz| cut -f 10|  cut -d ":" -f 1| sort| uniq -c



mkdir $OUTPUT_DIR/annotationXargs/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P5 bash -c 'annotationXargs "$@"' 'annotationXargs'

# USEI O LOOP FOR ABAIXO PARA ANNOTARA AS VARIANTES
for i in `cat ${SAMPLE_LIST_BAM}`; do
 annotation $i
done

SAMPLES=$(find $OUTPUT_DIR/annotation/ -maxdepth 1 -mindepth 1  -name '*.hg38_multianno.txt')
for SAMPLE in $SAMPLES; do
  NAME="${SAMPLE##*/}"
  awk -v N="$NAME" 'BEGIN { FS="\t"; OFS="\t" } $45 == "PASS" { print N, $0 }' "$SAMPLE" >> all.ROP_annotated.PASS.tsv
done


awk '{if ($46 == "PASS") print $_ }' all.ROP_annotated.txt | cut -f "49" | cut -d ":" -f 1| sort| uniq -c
awk '{if ($46 == "PASS") print $_ }' all.ROP_annotated.txt > all.ROP_annotated.PASS.tsv
cut -f "49" all.ROP_annotated.PASS.txt | cut -d ":" -f 1| sort| uniq -c

cut -f 1-18,36,38,42,45,48,83,111,148,152-249 $OUTPUT_DIR/annotation/annovar.norm.hg38_multianno.txt > $OUTPUT_DIR/all.snv.annotated.tsv


echo "" >> $TIME_FILE
echo "${TIME} >>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE
