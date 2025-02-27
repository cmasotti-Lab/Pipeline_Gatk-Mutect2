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
export OUTPUT_DIR=${WD}"/Result_Mutect2.GATK4.6.${DATA}"

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

export TIME_FILE="$OUTPUT_DIR.log"

mkdir $OUTPUT_DIR
#$(find "$WD/preprocessing_FINAL_result/" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam.cram') > $OUTPUT_DIR/samples.list.cram

ls $WD/preprocessing_FINAL_result/*.dedup.tags.bqsr.bam.cram > $OUTPUT_DIR/samples.cram
export SAMPLE_LIST="$OUTPUT_DIR/samples.cram"

ls $WD/preprocessing_FINAL_result/*.dedup.bam > $OUTPUT_DIR/samples.list.bam
export SAMPLE_LIST_BAM="$OUTPUT_DIR/samples.list.bam"

ls $WD/preprocessing_READ_result/*.dedup.tags.bqsr.bam > $OUTPUT_DIR/samples.list.bam
export SAMPLE_LIST_BAM="$OUTPUT_DIR/samples.list.bam"

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
  echo "" >> $WD/preprocessing_READ_result/preprocessing.log
  echo ">>>>>> Executando stage_SetNmMdAndUqTags para Amostra: "$NAME" <<<  $(date) " >> $WD/preprocessing_READ_result/preprocessing.log

  ${PICARD} SetNmMdAndUqTags \
    R=$REF_FASTA/Homo_sapiens_assembly38.fasta \
    I=$SAMPLE \
    CREATE_INDEX=true \
    O=$WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bam" \
    2> $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bam.log"
}
export -f  stage_SetNmMdAndUqTags

    # "GATK BaseRecalibrator|PrintReads"
stage_base_recalibration(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $WD/preprocessing_READ_result/preprocessing.log
  echo ">>>>>> Executando stage_base_recalibration para Amostra: "$NAME" <<<  $(date) " >> $WD/preprocessing_READ_result/preprocessing.log

  $GATK --java-options "-Xmx${MEM}G" BaseRecalibrator \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -I $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bam" \
        -known-sites ${INDEL_KNOWN} \
        -known-sites ${DBSNP} \
        -known-sites ${GNOMAD} \
        -L ${TARGET} \
        -O $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.table" \
        > $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.table.log" \
        2> $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.table.log2"

  $GATK --java-options "-Xmx${MEM}G" ApplyBQSR \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -I $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bam" \
        --bqsr-recal-file $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.table"\
        -O $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.bam" \
        > $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.bam.log" \
        2> $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.bam.log2"
}
export -f  stage_base_recalibration


echo " >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " > $WD/preprocessing_READ_result/preprocessing.log

xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'stage_SetNmMdAndUqTags  "$@"' 'stage_SetNmMdAndUqTags'

xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'stage_base_recalibration  "$@"' 'stage_base_recalibration'



stage_Mutect2 (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  $GATK --java-options "-Xmx${MEM}G" Mutect2 \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -L $TARGET \
        -I $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.bam" \
        --germline-resource $GNOMAD \
        --panel-of-normals $PON   \
        --f1r2-tar-gz $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.f1r2.tar.gz" \
        -bamout $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.bamout.bam" \
        -O $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.unfiltered.vcf.gz" 2> $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.unfiltered.log"
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

  $GATK --java-options -Xmx300G GetPileupSummaries \
        -I $WD/preprocessing_READ_result/"${NAME%.dedup*}.dedup.tags.bqsr.bam" \
        -V $GNOMAD \
        -L $GNOMAD \
        -O $OUTPUT_DIR/GetPileupSummaries/"${NAME%.dedup*}.getpileupsummaries.table" 2> $OUTPUT_DIR/GetPileupSummaries/"${NAME%.dedup*}".getpileupsummaries.log
}
export -f stage_GetPileupSummaries


stage_CalculateContamination (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando CalculateContamination para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE

  $GATK --java-options "-Xmx${MEM}G" CalculateContamination \
        -I $OUTPUT_DIR/GetPileupSummaries/"${NAME%.dedup*}.getpileupsummaries.table" \
        -tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${NAME%.dedup*}.segments.table" \
        -O $OUTPUT_DIR/CalculateContamination/"${NAME%.dedup*}.calculatecontamination.table" 2> $OUTPUT_DIR/CalculateContamination/"${NAME%.dedup*}.calculatecontamination.log"
}
export -f stage_CalculateContamination


stage_FilterMutectCalls (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$NAME" <<<  $(date) " >> $TIME_FILE
 
  $GATK --java-options "-Xmx${MEM}G" FilterMutectCalls \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -V $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.unfiltered.vcf.gz" \
        --tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${NAME%.dedup*}.segments.table" \
        --contamination-table $OUTPUT_DIR/CalculateContamination/"${NAME%.dedup*}.calculatecontamination.table" \
        --stats $OUTPUT_DIR/Mutect2/"${NAME%.dedup*}.unfiltered.vcf.gz.stats" \
        --ob-priors "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        -O $OUTPUT_DIR/FilterMutectCalls/"${NAME%.dedup*}.filtered.vcf.gz" 2> $OUTPUT_DIR/FilterMutectCalls/"${NAME%.dedup*}.filtered.vcf.log"
}
export -f stage_FilterMutectCalls


left_normalization () {
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando Normalization para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE
  bcftools norm -m-both -O z -o $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step1.vcf.gz" \
    $OUTPUT_DIR/FilterMutectCalls/"${NAME%.dedup*}.filtered.vcf.gz" 2> $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step1.log"
  bcftools norm -O z -f $REF_FASTA/Homo_sapiens_assembly38.fasta -o $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step2.vcf.gz" \
    $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step1.vcf.gz" 2> $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step2.log"
  bcftools index $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step2.vcf.gz"
}
export -f left_normalization


selectPASS(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "> Executando selectPASS para Amostra: "$NAME" <<< $(date) " >> $TIME_FILE
  bcftools view -f PASS -O z $OUTPUT_DIR/left_normalization/"${NAME%.dedup*}.norm_Step2.vcf.gz" > $OUTPUT_DIR/FILTER_PASS/"${NAME%.dedup*}.pass.vcf.gz"
  bcftools index $OUTPUT_DIR/FILTER_PASS/"${NAME%.dedup*}.pass.vcf.gz" 
}
export -f selectPASS


selectVar(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "> Executando selectVar para Amostra: "$NAME" <<< $(date) " >> $TIME_FILE
  zgrep "#" $OUTPUT_DIR/FILTER_PASS/"${NAME%.dedup*}.pass.vcf.gz" > $OUTPUT_DIR/FILTER_VAR/"${NAME%.dedup*}.pass.var.vcf"
  zgrep -v "#" $OUTPUT_DIR/FILTER_PASS/"${NAME%.dedup*}.pass.vcf.gz" |  grep -v -E "0/0:|0/\.|\.\/\." >> $OUTPUT_DIR/FILTER_VAR/"${NAME%.dedup*}.pass.var.vcf"
  bgzip $OUTPUT_DIR/FILTER_VAR/"${NAME%.dedup*}.pass.var.vcf" 
  bcftools index $OUTPUT_DIR/FILTER_VAR/"${NAME%.dedup*}.pass.var.vcf.gz" 
}
export -f selectVar


annotation (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  # # annovar
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando annovar para amostra $NAME: <<<<<< $(date)" >> $TIME_FILE
  
  $ANNOVAR  \
    --vcfinput $OUTPUT_DIR/FILTER_VAR/"${NAME%.dedup*}.pass.var.vcf.gz"  $ANNOVAR_DB -buildver hg38 --remove \
    --protocol refGene,avsnp150,gnomad41_exome_filt,abraom,cosmic99,icgc28,dbnsfp42a_filt,clinvar_20220320  \
    --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
    --xreffile $CROSS_REFERENCE --otherinfo --thread 10 \
    --outfile $OUTPUT_DIR/annotation/${NAME%.dedup*} > $OUTPUT_DIR/annotation/${NAME%.dedup*}.log 2> $OUTPUT_DIR/annotation/${NAME%.dedup*}.log2
  sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/${NAME%.dedup*}.hg38_multianno.vcf | sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/${NAME%.dedup*}.hg38_multianno.correct.vcf 
}
export -f annotation





echo "            >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " >> $TIME_FILE

mkdir $OUTPUT_DIR/Mutect2/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P5 bash -c 'stage_Mutect2  "$@"' 'stage_Mutect2' &

mkdir $OUTPUT_DIR/LearnReadOrientationModel/
stage_LearnReadOrientationModel 

mkdir $OUTPUT_DIR/GetPileupSummaries/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P1 bash -c 'stage_GetPileupSummaries "$@"' 'stage_GetPileupSummaries' &

mkdir $OUTPUT_DIR/CalculateContamination/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'stage_CalculateContamination  "$@"' 'stage_CalculateContamination'

mkdir $OUTPUT_DIR/FilterMutectCalls/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'stage_FilterMutectCalls  "$@"' 'stage_FilterMutectCalls'

mkdir $OUTPUT_DIR/left_normalization/
xargs -a ${SAMPLE_LIST_BAM}  -t -n1 -P10 bash -c 'left_normalization  "$@"' 'left_normalization'

mkdir -p $OUTPUT_DIR/FILTER_PASS
echo ">>>>>> selectPASS : <<< $(date) " >> $TIME_FILE
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'selectPASS  "$@"' 'selectPASS'

mkdir -p $OUTPUT_DIR/FILTER_VAR 
echo ">>>>>> selectVAR : <<< $(date) " >> $TIME_FILE
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'selectVar  "$@"' 'selectVar'

mkdir $OUTPUT_DIR/annotation/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'annotation "$@"' 'annotation'

# USEI O LOOP FOR ABAIXO PARA ANNOTARA AS VARIANTES
for i in `cat ${SAMPLE_LIST_BAM}`; do
 annotation $i
done


# PASSO  - Concatenar tabelas do annovar

SAMPLES=$(find $OUTPUT_DIR/annotation -maxdepth 1 -mindepth 1  -name '*.hg38_multianno.txt')
for SAMPLE in $SAMPLES; do
  NAME="${SAMPLE##*/}"
  awk -OFS="\t" -v N=${NAME%.hg38*} '{print N,$_ }' $SAMPLE >> $OUTPUT_DIR/Final_GATK.4.6_annotated.txt
done

sed -i 's/\s/\t/' $OUTPUT_DIR/Final_GATK.4.6_annotated.txt

echo "" >> $TIME_FILE
echo "${TIME} >>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE


mkdir $OUTPUT_DIR/vcfs_toMutalisk
SAMPLES=$(find $OUTPUT_DIR/mutations_to_filter_fromVCF/synonymous -maxdepth 1 -mindepth 1  -name '*.synonymous.bed')

for SAMPLE in ${SAMPLES}; do
  NAME="${SAMPLE##*/}"
  vcftools --gzvcf $OUTPUT_DIR/FILTER_VAR/"${NAME%.synonymous.bed*}.pass.var.vcf.gz" \
   --positions ${SAMPLE} \
   --out $OUTPUT_DIR/vcfs_toMutalisk/"${NAME%.synonymous.bed*}" --recode --keep-INFO-all \
   2> $OUTPUT_DIR/vcfs_toMutalisk/"${NAME%.synonymous.bed*}.log2"
done

   bedtools intersect -wa \
    -a $OUTPUT_DIR/FILTER_VAR/"${NAME%.synonymous.bed*}.pass.var.vcf.gz" \
    -b  ${SAMPLE} \
    > $OUTPUT_DIR/Mutation_toMutalisk/$NAME.bedtools
    
/home/scratch90/vlira_13may2024//Result_Mutect2.GATK4.6.2024-06-14/FILTER_VAR/ROP-97-ExC85-xgenV2_S65.pass.var.vcf.gz


http://mutalisk.org/result.php?rid=Snc2VuzCmJ

#synonymous
http://mutalisk.org/result.php?rid=YJqifIJ7TW

# PCAWG - SigProfiler 
http://mutalisk.org/result.php?rid=dn13TvlmTd