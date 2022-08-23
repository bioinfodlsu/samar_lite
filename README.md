Welcome to SAMAR-lite !

SAMAR-lite lets you perform quick-and-dirty differential gene expression analysis in organisms when you don't have a reference transcriptome sequence.

For non-model organisms, conventional DE analysis pipelines requires de-novo transcriptome assembly, which requires massive computational resources. 
This pipeline provides a shortcut in which RNA-seq reads are directly mapped to the reference proteome which would have been otherwise used for annotation and functional inference in the assembly-based approach. 
It is faster than [SAMAR](https://bitbucket.org/project_samar/samar/), but less accurate.

# 1. Installation
The read mapping part of SAMAR-lite is written in Rust.
## 1.1. Install  Rust 
### 1.1.1 Using Conda
This is probably the easiest way to install Rust. 

```
$ conda create --name rust
$ conda activate rust
$ conda install -c conda-forge rust
```
Test Rust by running:
```
$ cargo
```

### 1.1.2 As recommended by Rust:
Or you can install Rust directly:
```
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
After completing the rust installation, test Rust by running:
```
$ cargo
```

## 1.2. Download/clone this repository

```
$ git clone https://bitbucket.org/project_samar/samar_lite.git
$ cd samar_lite
```

## 1.3. Compiling the Binaries
In the `samar-lite` folder, run the `make` command.
```
$ make
```
This will automatically compile the binaries of the `reference` and `alignr` folders.
After compilation, the binaries will be placed inside the `bin` folder.

# 2. Running Project
From the `samar-lite` folder, go to the `bin` folder.
```
$ cd bin
```

## 2.1. Reference Generation
First the reference will be generated. 
Note that the input files for the reference section is: 
- protein reference fasta file

### 2.1.1. Running the Reference Binary 
To generate the reference, run this command: 
```
$ ./ref-align <path-to-reference>.fasta <path-to-output>.json <kmer size>
```
A json file in the specified path with be generated using the reference and kmer size. 

## 2.2. Pseudo Alignment Generation
Next, using the previously generated reference JSON, the Pseudo Alignments will be generated.
Note that the input files for the pseudo alignment section are: 
- reference JSON file
- interleaved DNA fastq query file 

### 2.1.1. Running the Reference Binary 
To generate the Pseudo Alignments, run this command: 
```
$ ./alignr <path-to-reference>.json <path-to-query>.fastq <threshold of coverage> > <path-to-desired-output>.txt
```
A text file in the specified path with be generated using the reference and threshold of coverage. 

# 3. Snakemake Pipeline

## 3.1. Install Snakemake
Snakemake recommends installation via Conda:
```
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```

## 3.2 Config File 
 
### 3.2.1. Input
The pipeline requires, at the very least: (1) RNA-seq fastq **interleaved** read, (2) a reference protein database in fasta format, (3) a k-mer size, and (4) a coverage threshold.
These and other input parameters are specified via a YAML-format config file -- a template config.yaml is provided in the config folder within the workflow folder.
 
### 3.2.2. Output
The output of the pipeline is stored inside the folder specified in the out_dir entry of the config file. 
Inside it, the mappings can be found in the mappings folder, the count data in the counts folder, and the results of differential expression analysis in the DEanalysis folder. 
 
In the DEanalysis folder in test_output, you will find the following files: 
 
- DE_count_summary.txt, which contains the count of up- and downregulated genes 
- Test_result-Group_treated_vs_control.csv, which contains the table of the results of hypothesis testing for each gene in the reference. 
- DESeq2_fit.RDS which can be loaded using R, if further analysis is required.
 
### 3.2.3. Config Parameters
The parameters in the YAML config file can be seen below:
 
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
 
## 3.3 Running Pipeline 
After constructing a config.yaml file and with the snakemake conda environment you created earlier activated, you can call the pipeline from the top-level directory of SAMAR_lite:
```
$ snakemake --configfile <workflow/config/config.yaml> --use-conda --cores all 
```
 
### 3.3.1 Example 
A toy example is provided in the test_data folder. 
In this example, there are 6 RNA-seq samples, of which 3 are the "control" group replicates and 3 are the "treated" group replicates. 
The paths to the files and information about the experiment design is specified in config_test.yaml. 
From the top directory of the pipeline and with the snakemake conda environment activated, run:
```
$ snakemake --configfile testdata/config.yaml --use-conda --cores all 
```
Once the run is complete, you can check the output in the folder test_output which is created in the top-level directory of SAMAR-lite.
