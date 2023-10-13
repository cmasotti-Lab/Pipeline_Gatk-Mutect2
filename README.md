# Pipeline_Gatk-Mutect2

A Pipeline foi desenvolvida para identificação de variantes somaticas em dados de Exoma tumor-only.
Tem como base o GATK Mutect2.
Utilizamos a ferramenta Mutect2 com uma abordagem de tumor-only para identificar mutações somáticas. 
O Mutect2 utiliza um painel de normais (PoN) para excluir as variantes germinativas raras, então utilizamos o PoN com 100 indivíduos do projeto Covid-Brasil. 



Referencias:
https://gatk.broadinstitute.org/hc/en-us/articles/360047232772--Notebook-Intro-to-using-Mutect2-for-somatic-data
https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136
https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf
https://gatk.broadinstitute.org/hc/en-us/articles/360035890491?id=11127
https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow

## O pipeline esta divido em duas etapa
  - Criação do Painel de Normal.
  - Identificação das variantes somática.

