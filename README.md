# Mouse Gut Microbiome Metranscriptomics

## Files overview

**rnaseq_0922.smk**

A workflow for preprocessing fastq files with fastp, bowtie2 alignment, and generating featurecounts. Information on how to obtain the reference genomes and fastq files along with further details on the workflow can be found in the Rmd files. 

**stm_dge.Rmd**

Markdown for differential gene expression of *Salmonella enterica subsp enterica* serovar Tyhphimurium using DeSeq2. 

**btheta_dge.Rmd**

Markdown for differential gene expression of *Bacteroides thetaiotamicron* using DeSeq2. 

**stm_1110.txt**

Featurecounts table for S. Tm. 

**bth_1110.txt**

Featurecounts table for B. thetaiotamicron. 

**exp_stm.csv**

Sample names and conditions for S. Tm differential gene expression. 

**exp_bth.csv**

Sample names and conditions for B. thetaiotamicron differential gene expression. 

**prelim_DGE.Rmd** 

A scratch work file. Ignore. 

The remaining csv files are raw data of the DeSeq2 results, results sorted by log2 fold change and filted by p-values <0.05, and raw data for MA plots, volcano plots, and heatmaps.
