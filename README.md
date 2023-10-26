# Pipeline_Gatk-Mutect2

![Workflow](https://github.com/cmasotti-Lab/Pipeline_Gatk-Mutect2/assets/11162991/74970b4a-c88d-4d17-8957-6b3824d61d9f)
[Detailed Extended Pipelines](https://drive.google.com/file/d/100eEe_oiofVWKpySCTJfEtS2OvwRBM9m/view?usp=sharing)

This pipeline has been developed for the identification of somatic variants in tumor-only exome data.



The samples used are from the locally advanced rectal cancer project.

The pipeline is based on the GATK Mutect2 approach for identifying somatic mutations in tumor samples in the presence of normal controls.

Mutect2 uses a Panel of Normals (PoN) with unrelated samples to separate germline variants from somatic ones.

For the development of the PoN, we used samples from 100 non-cancer individuals from the SARS-CoV2-Brasil project (SECOLIN, et al. 2021).

### GATK Mutect2 Documentation:

- [Intro to using Mutect2 for somatic data](https://gatk.broadinstitute.org/hc/en-us/articles/360047232772--Notebook-Intro-to-using-Mutect2-for-somatic-data)
- [How to Call somatic mutations using GATK4 Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
- [Mutect2 Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035889791?id=11136)
- [Mutect2 PDF Manual](https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf)
- [Tumor-Only Variant Calling Workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035890491?id=11127)
- [DNA Sequencing Variant Calling Pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow)

## Pipeline Stages

This pipeline is divided into two stages:

- Creation of the Panel of Normals (PoN).
- Identification of somatic variants.

[Detailed Extended Pipelines](https://drive.google.com/file/d/100eEe_oiofVWKpySCTJfEtS2OvwRBM9m/view?usp=sharing)
