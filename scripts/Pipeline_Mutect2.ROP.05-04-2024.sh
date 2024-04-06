#Pipeline atualizado em 05/04/2024
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado
#2. Alterado o diretorio onde estava o PON na scratch45 para /home/users/vlira/PanelOfNormals/PoN.100COVID.vcf.gz
#3. Agora o script recebe 2 agumentos: 1- lista de amostras; 2- diretorio SCRATCH
#4. Utualizado os databases: hg38_cosmic98_coding,hg38_avsnp150,hg38_clinvar_20220320, hg38_gnomad40_exome

export SCRATCH60="/home/scratch60/vlira_15mar2024/"

# Verifica se dois argumentos foram passados
if [ "$#" -eq 2 ]; then
    export SAMPLE_LIST=${1}
    export SCRATCH="${2}/vlira_15mar2024/"
else
    # Valores definidos por padrão
    export SAMPLE_LIST="${SCRATCH60}/samples.lists/samples.list"
    export SCRATCH=${SCRATCH60}
fi

export MEM=200
export JOBS=5

export DATA=$(date "+%F")
export DATA="2024-04-05"  # EDITE AQUI SE QUISER USAR UMA PASTA DE UMA DATA ESPECIFICA
export OUTPUT_DIR=${SCRATCH}"/Result_Mutect2.ROP.toPureCN.${DATA}"

export REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/"
export GATK="$SCRATCH60/tools/gatk-4.3.0.0/./gatk"

export TARGET="$SCRATCH60/references/xgen-exome-research-panel-v2-targets-hg38.bed"
export PON="/home/users/vlira/PanelOfNornal/PON_ToPureCN/PON_Mutect2/pon.vcf.gz"
export GNOMAD="$SCRATCH60/references/af-only-gnomad.hg38.vcf.gz"
#export GNOMAD="$SCRATCH60/references/af-only-gnomad.SABE1171.Abraom.hg38.vcf.gz"
export ANNOVAR="$SCRATCH60/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$SCRATCH60/humandb/"
export CROSS_REFERENCE="$SCRATCH60/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

export TIME_FILE="$OUTPUT_DIR.log"

mkdir $OUTPUT_DIR


stage_Mutect2 (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

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
  echo ">>>>>> STAGE_LearnReadOrientationModel<<<<<<" >> $TIME_FILE
  date >> $TIME_FILE

  local SAMPLE_READORIENTATION=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.f1r2.tar.gz')
  local SAMPLE_F1R2=$(echo $SAMPLE_READORIENTATION| sed 's/\s/ -I  /g')

  $GATK --java-options  "-Xmx${MEM}G"  LearnReadOrientationModel \
        -I $SAMPLE_F1R2 \
        -O "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz 2> $OUTPUT_DIR/LearnReadOrientationModel/read-orientation-model.tar.gz.log
} 


stage_GetPileupSummaries (){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> Executando GetPileupSummaries para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK --java-options  "-Xmx${MEM}G"  GetPileupSummaries \
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
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE
 
  $GATK --java-options "-Xmx${MEM}G" FilterMutectCalls \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -V $OUTPUT_DIR/Mutect2/$NAME.unfiltered.vcf.gz \
        --tumor-segmentation $OUTPUT_DIR/CalculateContamination/$NAME.segments.table \
        --contamination-table $OUTPUT_DIR/CalculateContamination/$NAME.calculatecontamination.table \
        --stats $OUTPUT_DIR/Mutect2/$NAME.unfiltered.vcf.gz.stats \
        --ob-priors "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        -O $OUTPUT_DIR/FilterMutectCalls/$NAME.filtered.vcf  2> $OUTPUT_DIR/FilterMutectCalls/$NAME.filtered.vcf.log
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
   --protocol refGene,avsnp150,gnomad40_exome,abraom,cosmic98_coding,icgc28,dbnsfp42a,clinvar  \
   --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
   --xreffile $CROSS_REFERENCE --otherinfo --thread 5 --outfile $OUTPUT_DIR/annotation/annovar.norm 2> $OUTPUT_DIR/annotation/annovar.norm.log

   sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/annovar.norm.hg38_multianno.vcf| sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/annovar.norm.hg38_multianno.correct.vcf 

   date >> $TIME_FILE
  

  date >> $TIME_FILE
  echo ">>>>>> Executando SnpSift para todas juntar os vcf das amostras<<<<<<" >> $TIME_FILE
  
  # java -jar -Xmx50G $SCRATCH/tools/snpEff/SnpSift.jar extractFields  "$OUTPUT_DIR"/annotation/annovar.norm.hg38_multianno.correct.vcf \
  #   -e . "CHROM" "POS" "ID" "REF" "ALT" "FILTER" "AC" "AN" "DP" "Func.refGene" "Gene.refGene" "GeneDetail.refGene" "ExonicFunc.refGene" "AAChange.refGene" \
  #   "COSMIC_Census_Gene.refGene" "Role_in_Cancer.refGene" "Translocation_Partner.refGene" "Therapeutic_Agents.refGene" "Cancer_Syndromes.refGene" \
  #   "panel.refGene" "gnomAD_exome_ALL" "abraom_freq" "cosmic95" "ICGC_Id" "GEN[*].GT"  > "$OUTPUT_DIR"/Final.mutect2.txt 2> "$OUTPUT_DIR"/Final.mutect2.log

  echo "" >> $TIME_FILE

}
export -f annotation




echo "                       >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<<" >> $TIME_FILE
date >> $TIME_FILE

#mkdir $OUTPUT_DIR/Mutect2/
#xargs -a ${SAMPLE_LIST} -t -n1 -P${JOBS} bash -c 'stage_Mutect2  "$@"' 'stage_Mutect2'

#mkdir $OUTPUT_DIR/LearnReadOrientationModel/
#stage_LearnReadOrientationModel

#mkdir $OUTPUT_DIR/GetPileupSummaries/
#xargs -a ${SAMPLE_LIST} -t -n1 -P${JOBS} bash -c 'stage_GetPileupSummaries  "$@"' 'stage_GetPileupSummaries'

#mkdir $OUTPUT_DIR/CalculateContamination/
#xargs -a ${SAMPLE_LIST} -t -n1 -P${JOBS} bash -c 'stage_CalculateContamination  "$@"' 'stage_CalculateContamination'

#mkdir $OUTPUT_DIR/FilterMutectCalls/
#xargs -a ${SAMPLE_LIST} -t -n1 -P${JOBS} bash -c 'stage_FilterMutectCalls  "$@"' 'stage_FilterMutectCalls'

# mkdir $OUTPUT_DIR/left_normalization/
# xargs -a ${SAMPLE_LIST} -t -n1 -P${JOBS} bash -c 'left_normalization  "$@"' 'left_normalization'

mkdir $OUTPUT_DIR/annotation/
annotation


echo "" >> $TIME_FILE
echo "${TIME} >>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE
