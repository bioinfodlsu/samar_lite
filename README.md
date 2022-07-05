# 1. Installation

## 1.1. Install  Rust 
### 1.1.1 Using Conda
This is probably the most painless way to install Rust. 

```
conda create --name rust
conda activate rust
conda install -c conda-forge rust
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

