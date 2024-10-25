# Welcome to SAMAR-lite !

SAMAR-lite lets you perform ultra rapid differential gene expression analysis in organisms when you don't have a reference transcriptome sequence.

For non-model organisms, conventional differential expression analysis pipelines require de-novo transcriptome assembly, which is computationally expensive. 
SAMAR-lite provides a shortcut in which RNA-seq reads are directly quasi-mapped to a reference proteome database containing protein sequences which would have been otherwise used for annotation and functional inference in the assembly-based approach). 
The use of quasi-mapping instead of full alignments makes it faster than [SAMAR](https://github.com/bioinfodlsu/samar), but less accurate.

# 1. Installation
This pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
All other dependencies are handled automatically by Snakemake.

## 1.1. Install Conda 
Conda can be downloaded and installed from [Miniforge](https://github.com/conda-forge/miniforge).
Once installation is complete, you can test your Conda installation by running:
```
$ conda list
```

## 1.2. Install Snakemake
Snakemake recommends installation via Conda:
```
$ conda create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```
## 1.3. Download/clone this repository
```
$ git clone https://github.com/bioinfodlsu/samar_lite.git
$ cd samar_lite
```
## 1.4. Check if SAMAR_lite works
A toy example is provided in the `test_data` folder. 
In this example, there are 6 RNA-seq samples, of which 3 are the "control" group replicates and 3 are the "treated" group replicates. 
The paths to the files and information about the experiment design is specified in config_test.yaml. 
From the top-level directory of SAMAR_lite, run:
```
$ snakemake -p --configfile testdata/config_test.yaml --use-conda --cores all 
```
Since this is the first-ever run, it will take a while (5-10 minutes) for Conda to download and install all dependencies.
Once the run is complete, you can check if some outputs are generated in the newly created `test_output` folder.
The contents of the output folder are described in the section below.

# 2. User Guide
 
### 2.1. Input
The pipeline requires, at the very least: (1) RNA-seq fastq reads, (2) a reference protein database in fasta format, (3) a k-mer size, and (4) a coverage threshold.
These and other input parameters are specified via a YAML-format config file -- a template config.yaml is provided in the config folder within the workflow folder.
 
### 2.2. Output
The output of the pipeline is stored inside the folder specified in the out_dir entry of the config file. 
Inside it, the mappings can be found in the mappings folder, the count data in the counts folder, and the results of differential expression analysis in the DEanalysis folder. 
 
In the DEanalysis folder in test_output, you will find the following files: 
 
- DE_count_summary.txt, which contains the count of up- and downregulated genes 
- Test_result-Group_treated_vs_control.csv, which contains the table of the results of hypothesis testing for each gene in the reference. 
- DESeq2_fit.RDS which can be loaded using R, if further analysis is required.
 
### 2.3. Config Parameters
A template configuration file is provided as a YAML file: `workfow/config/config.yaml`.
The parameters to be provided by the user are:
 
| Keyword       |   Possible values         | Default  |  Description  |
| ------------- |------------------------| ------ |  ------------|
| reference | path/to/reference/file.fasta  | - | Reference protein database in fasta format |
| k-size | positive integer value | 7 | Size of the k-mers |
| threshold | positive float value | 40.0 | Minimum mapping coverage |
| reads: sample_i | path/to/query/file_i.fastq | - | RNA-seq fastq reads |
| out_dir | path/to/output/dir | - | Directory where outputs will be stored where, the mappings can be found in the mapping folder, the count data in the counts folder, and the results of differential expression analysis in the DEanalysis folder |
| deanalysis | yes, no | yes |  If set to "no", the pipeline does not proceed to DE anlaysis and halts after counting. If keyword is not provided, value defaults to "yes" |
|factor| [factor1, factor2, .. ] | - | Name of factors in the study |
| sample_info: sample_i: |  [level of factor1, level of factor2, ..] | -| Describes the factor levels which the sample corresponds to |
 
## 2.4 Running Pipeline 
After fully constructing the config.yaml file and with the snakemake conda environment you created earlier activated, you can call the pipeline from the top-level directory of SAMAR_lite:
```
$ snakemake -p --configfile workflow/config/config.yaml --use-conda --cores all 
```
 


# 3 Publication
Santiago, K.C.L., Shrestha, A.M.S. DNA-protein quasi-mapping for rapid differential gene expression analysis in non-model organisms. BMC Bioinformatics 25 (Suppl 2), 335 (2024). [Link](https://doi.org/10.1186/s12859-024-05924-1)

This work was presented at [GIW XXXI/ISCB-Asia V, Dec 12 – 14, 2022, Tainan, Taiwan](https://www.iscb.org/giw-iscb-asia2022).


# 4 Contact
For comments, suggestions, issues, please send an email to anish.shrestha--atmark--dlsu.edu.ph
