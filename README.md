# minimizers
Compute set of minimizers for DNA or Protein sequences.

The implementation is based on the article from [Roberts et al.](https://academic.oup.com/bioinformatics/article/20/18/3363/202143) that first described how to determine the set of minimizers for biological sequences.

## Usage

```
Purpose
-------

Compute the set of minimizers for DNA or Protein sequences.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE         Path to a FASTA file with DNA or protein sequences.
  -o OUTPUT_FILE        Path to output to save results.
  --of {binary,csv}     Output file format. `binary` will use pickle to save
                        results into a binary file. `csv` will save results to
                        a CSV file with sequence identifiers in the first
                        field followed by a minimizer per field.
  --k K_VALUE           Value for the size of the kmers.
  --w WINDOW_SIZE       Window size value. Minimizers will be determined based
                        on groups of `w` adjacent kmers.
  --p                   If the start position of the kmers should be returned.
  --m {default,skipper}
                        Defines function that will be used to determine
                        minimizers.
  --t THREADS           Number of CPU cores/threads used to parallelize
                        minimizer computation.
```
