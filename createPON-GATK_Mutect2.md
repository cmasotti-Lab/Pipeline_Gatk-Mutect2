# Pipeline_Gatk-Mutect2

A Pipeline foi desenvolvida para identificação de variantes somaticas em dados de Exoma tumor-only.

As amostras utilizadas são do projeto de cancer de reto localmente avançado.

O pipeline tem como base a abordagem do GATK Mutect2 para identificaça mutações somáticas em amostras tumorais em a presença do controle normal.

O Mutect2 utiliza um painel de normais (PoN) com não relacionadas para saparar as variantes germinativas das somaticas.

Para o desenvolvimento do PoN utilizamos amostras de 100 indivíduos (não-cancer) do projeto SARS-CoV2-Brasil (SECOLIN, et al. 2021).  

### Documentação do GATK Mutect2:

https://gatk.broadinstitute.org/hc/en-us/articles/360047232772--Notebook-Intro-to-using-Mutect2-for-somatic-data
https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136
https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf
https://gatk.broadinstitute.org/hc/en-us/articles/360035890491?id=11127
https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow

## O pipeline esta divido em duas etapa
  - Criação do Painel de Normal.
  
  - Identificação das variantes somática.

