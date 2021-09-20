## Compiling to Binary 

In the root directory of the project, open the terminal and input the commnad:

```
$ cargo build
```

After running the command above, a new folder called debug will be created containing the binary called "new-p-align". 

## Running the binary

From the root directory, navigate to the debug folder, open the terminal and input the command:

```
./new-p-align <path-to-reference>.fasta <path-to-query>.fastq
```

Both paths require the file extension. 

It is assumed that the input query sequence is already interleaved. 
