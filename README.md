
# 1. Software Requirement
C++ Build Tools and Rust are required for this project. All other dependencies will be handled by rust. 

## 1.1. Install C++ Build Tools (for Rust)

Download the C++ build tools [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/).

## 1.2. Install Rust 

Download Rust [here](https://www.rust-lang.org/).
After completing the rust installation, test Rust by running:
```
$ cargo
```

## 1.3. Download/clone this repository

```
$ git clone https://kylesans@bitbucket.org/project_samar/samar_lite.git
$ cd samar_lite
```

# 2. Running Project
The project is separated into 2 sections: A reference section and a pseudo alignment section. 

## 2.1. Reference Generation
First the reference will be generated. 
Note that the input files for the reference section are: 
- protein reference fasta file

### 2.1.1. Reference folder
The reference section of samar_lite can be found in the reference folder: 

```
$ cd reference
```

### 2.1.2. Compiling the Binary
Once in the reference folder, the reference section can now be compiled into a binary file. 
To compile the reference section use:
```
$ cargo build --release
```
A new folder will appear called "target". Further in that folder is a new folder called "release". 
The binary of the reference section will be found in the release folder: 
```
$ cd target/release 
```

### 2.1.3. Running the Reference Binary 
Once inside the release folder, run this command: 
```
$ ./ref-align <path-to-reference>.fasta <path-to-output>.json <kmer size>
```
A json file in the specified path with be generated using the reference and kmer size. 

## 2.2. Pseudo Alignment Generation
Next, using the previously generated reference JSON, the Pseudo Alignments will be generated.
Note that the input files for the pseudo alignment section are: 
- reference JSON file
- interleaved DNA fastq query file 

### 2.2.1. Alignment folder
The pseudo alignment section of the project can be found in the alignment folder: 
```
$ cd alignr
```
### 2.1.2. Compiling the Binary
Similar to the Reference section, the compilation of the binary is the same.

Once in the reference folder, the reference section can now be compiled into a binary file. 
To compile the reference section use:
```
$ cargo build --release
```
A new folder will appear called "target". Further in that folder is a new folder called "release". 
The binary of the reference section will be found in the release folder: 
```
$ cd target/release 
```
### 2.1.3. Running the Reference Binary 
Once inside the release folder, run this command: 
```
$ ./alignr <path-to-reference>.json <path-to-query>.fastq <threshold of coverage> > <path-to-desired-output>.txt
```
A text file in the specified path with be generated using the reference and threshold of coverage. 

