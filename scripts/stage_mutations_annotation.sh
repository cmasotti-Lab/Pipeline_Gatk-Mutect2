export SCRATCH60="/home/scratch60/vlira_15mar2024/"
export SCRATCH="/home/scratch60/vlira_15mar2024/"
export DATA="2024-04-05"  # EDITE AQUI SE QUISER USAR UMA PASTA DE UMA DATA ESPECIFICA
export OUTPUT_DIR=${SCRATCH}"/Result_Mutect2.ROP.toPureCN.${DATA}"

export ANNOVAR="$SCRATCH60/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$SCRATCH60/humandb/"
export CROSS_REFERENCE="$SCRATCH60/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"


stage_mutations_annotation (){
	local SAMPLE=$1
	# /home/scratch60/vlira_15mar2024/preprocessing_FINAL_result/ROP-100-ExC85-xgenV2_S67.dedup.tags.bqsr.bam
	local NAME="${SAMPLE##*/}"
	#ROP-1-ExC85-xgenV2_S7.dedup.tags.bqsr.bam
	local SAMPLE_ID=$(echo $NAME | cut -d'.' -f1)

	vcftools --gzvcf $OUTPUT_DIR/left_normalization/${SAMPLE_ID}.dedup.tags.bqsr.bam.norm_Step2.vcf.gz --bed $OUTPUT_DIR/mutations_to_filter_fromVCF/${SAMPLE_ID}.bed --recode --stdout | gzip > $OUTPUT_DIR/mutations_to_filter_fromVCF/${SAMPLE_ID}.recode.vcf.gz

	mkdir $OUTPUT_DIR/mutations_annotation
	$ANNOVAR  --vcfinput $OUTPUT_DIR/mutations_to_filter_fromVCF/${SAMPLE_ID}.recode.vcf.gz \
	   $ANNOVAR_DB -buildver hg38 --remove \
	   --protocol refGene,dbnsfp42a,clinvar_20220320  \
	   --operation gx,f,f --arg '-splicing 5',, --polish \
	   --xreffile $CROSS_REFERENCE --otherinfo --thread 5 --outfile $OUTPUT_DIR/mutations_annotation/${SAMPLE_ID} 2> $OUTPUT_DIR/mutations_annotation/${SAMPLE_ID}.log

	sed 's/\\x3b/;/g' $OUTPUT_DIR/mutations_annotation/${SAMPLE_ID}.hg38_multianno.vcf | sed 's/\\x3d/=/g' > $OUTPUT_DIR/mutations_annotation/${SAMPLE_ID}.hg38_multianno.correct.vcf 
}
export -f stage_mutations_annotation


find "$OUTPUT_DIR/mutations_to_filter_fromVCF/" -maxdepth 1 -mindepth 1  -name '*.bed' > $OUTPUT_DIR/mutations_to_filter_fromVCF/samples.list

stage_mutations_annotation
xargs -a $OUTPUT_DIR/mutations_to_filter_fromVCF/samples.list -t -n1 -P15 bash -c 'stage_mutations_annotation  "$@"' 'stage_mutations_annotation'
