WD="/home/venus/mar/vlira"

CRAM_DIR="$WD/samples/"
CRAM_FILES=$(find "$CRAM_DIR" -maxdepth 1 -mindepth 1  -name '*.cram')
REF_FASTA="$WD/reference/"
GATK="$WD/gatk-4.3.0.0/./gatk"
TARGET="$WD/reference/xgen-exome-research-panel-v2-targets-hg38.autossome.bed"

OUTPUT_DIR="RESULT_PON-GATK4.3_Mutect2"
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/Mutect2/
mkdir $OUTPUT_DIR/GenomicsDBImport/
mkdir $OUTPUT_DIR/PanelOfNormals/

TIME_FILE="$OUTPUT_DIR/$OUTPUT_DIR.log"




echo "                                                     >>>>>> Starting Pipeline  to create PON Mutect2 <<<<<<" >> $TIME_FILE
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
       -max-mnp-distance 0 \
       -O $OUTPUT_DIR/Mutect2/$NAME.vcf.gz 2> $OUTPUT_DIR/Mutect2/$NAME.log
}


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


#for run in $CRAM_FILES; do 
#  STAGE_Mutect2 "$run"  >> $TIME_FILE
#done

STAGE_GenomicsDB

STAGE_CreateSomaticPanelOfNormals

echo "" >> $TIME_FILE
echo ">>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE
