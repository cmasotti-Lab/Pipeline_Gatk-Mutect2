#!/bin/bash 

export WD="/home/venus/mar/vlira"

export CRAM_DIR="$WD/samples/"
export CRAM_FILES=$(find "$CRAM_DIR" -maxdepth 1 -mindepth 1  -name '*.cram' )
export REF_FASTA="$WD/reference/"
export GATK="$WD/gatk-4.3.0.0/./gatk"
#export TARGET="$WD/reference/xgen-exome-research-panel-v2-targets-hg38.bed"   # usado para criar _2024/
export TARGET="$WD/reference/baits_optimized_hg38.bed"

export OUTPUT_DIR="RESULT_PON-GATK4.3_Mutect2.ToPureCN"
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/Mutect2/
mkdir $OUTPUT_DIR/GenomicsDBImport/
mkdir $OUTPUT_DIR/PanelOfNormals/

export TIME_FILE="$OUTPUT_DIR/$OUTPUT_DIR.log"


export CRAM_LIST=$1



echo "                                            >>>>>> Starting Pipeline  to create PON Mutect2 <<<<<<" >> $TIME_FILE
date >> $TIME_FILE


STAGE_Mutect2(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_Mutect2 <<<" >> $TIME_FILE
  echo ">>>>>> Executando para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK Mutect2  \
       -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
       -I $SAMPLE \
       --genotype-germline-sites true --genotype-pon-sites true \
       -max-mnp-distance 0 \
       -O $OUTPUT_DIR/Mutect2/$NAME.vcf.gz 2> $OUTPUT_DIR/Mutect2/$NAME.log
}
export -f STAGE_Mutect2

STAGE_GenomicsDB(){
 local SAMPLES_VCF=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.vcf.gz')
 local SAMPLES=$(echo $SAMPLES_VCF| sed 's/\s/ -V  /g')

  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_GenomicsDB <<<" >> $TIME_FILE
  echo ">>>>>> Montando PoN para as Amostras: "$SAMPLES_VCF" <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK GenomicsDBImport \
       -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
       -L $TARGET \
       --merge-input-intervals true \
       --genomicsdb-workspace-path $OUTPUT_DIR/pon_db \
       -V $SAMPLES   2> $OUTPUT_DIR/GenomicsDBImport/GenomicsDBImport.log
}
export -f STAGE_GenomicsDB

STAGE_CreateSomaticPanelOfNormals(){
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_CreateSomaticPanelOfNormals <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK  --java-options "-Xmx6500m" CreateSomaticPanelOfNormals\
      -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
      --germline-resource $REF_FASTA/af-only-gnomad.hg38.vcf.gz \
      -V gendb://$OUTPUT_DIR/pon_db \
      -O $OUTPUT_DIR/PanelOfNormals/pon.vcf.gz 2> $OUTPUT_DIR/PanelOfNormals/my.pon.log
}
export -f STAGE_CreateSomaticPanelOfNormals


DepthOfCoverage(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> DepthOfCoverage <<<" >> $TIME_FILE
  echo ">>>>>> Executando para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK DepthOfCoverage  \
       -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
       -I $SAMPLE \
       -L $TARGET \
       --output-format TABLE \
       --interval-merging-rule OVERLAPPING_ONLY \
       -O /home/venus/mar/vlira/DepthOfCoverage/$NAME 2> /home/venus/mar/vlira/DepthOfCoverage/$NAME.log
}
export -f DepthOfCoverage



#for run in $CRAM_FILES; do 
#  STAGE_Mutect2 "$run"  >> $TIME_FILE
#done
# xargs -a /home/venus/mar/vlira/samples.55.list -t -n1 -P2 bash -c 'STAGE_Mutect2  "$@"' 'STAGE_Mutect2'

STAGE_GenomicsDB

STAGE_CreateSomaticPanelOfNormals

mkdir $OUTPUT_DIR/DepthOfCoverage
xargs -a /home/venus/mar/vlira/samples.list -t -n1 -P3 bash -c 'DepthOfCoverage  "$@"' 'DepthOfCoverage'


echo "" >> $TIME_FILE
echo ">>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE



# mkdir $OUTPUT_DIR/DepthOfCoverage
#  $GATK DepthOfCoverage  \
#        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
#        -I /home/venus/mar/vlira/samples/C029228-ExC89-xgenV2.hg38.final.cram \
#        -L $TARGET \
#        --output-format TABLE \
#        --interval-merging-rule OVERLAPPING_ONLY \
#        -O /home/venus/mar/vlira/DepthOfCoverage/gakt-4.2.6.1 2> /home/venus/mar/vlira/DepthOfCoverage/gakt-4.2.6.1.log


# $GATK DepthOfCoverage  \
#        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
#        -I /home/venus/mar/vlira/samples/C029228-ExC89-xgenV2.hg38.final.cram \
#        -L /home/venus/mar/vlira/baits_hg38_intervals.txt \
#        --output-format TABLE \
#        --interval-merging-rule OVERLAPPING_ONLY \
#        -O /home/venus/mar/vlira/DepthOfCoverage/gakt-4.2.6.1.baits 2> /home/venus/mar/vlira/DepthOfCoverage/gakt-4.2.6.1.baits.log &


# export WD="/home/scratch60/vlira_18jan2024/"
# ${WD}/tools/gatk-4.3.0.0/./gatk  DepthOfCoverage  \
#        -R /home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/Homo_sapiens_assembly38.fasta \
#        -I ${WD}/PureCN_Docker/NORMAL/C028869-ExC69-xgenV1.hg38.final.cram \
#        -L ${WD}/PureCN_Docker/reference_files/xgen-exome-research-panel-v2-targets-hg38.bed \
#        --output-format TABLE \
#        --interval-merging-rule OVERLAPPING_ONLY \
#        -O ${WD}/PureCN_Docker/NORMAL/DepthOfCoverage2/C028869-ExC69-xgenV1.hg38.final.cram.depth > ${WD}/PureCN_Docker/NORMAL/DepthOfCoverage2/C028869-ExC69-xgenV1.hg38.final.cram.depth.erro 2> ${WD}/PureCN_Docker/NORMAL/DepthOfCoverage2/C028869-ExC69-xgenV1.hg38.final.cram.depth.log &
