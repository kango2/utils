# Hard-mask specific sequences in a fasta file with Ns

Use cases: sex chromosome masking

```
masksequences.pl input.fasta output.fasta maskIDs
```

# Split large chromosomes (>512Mb) and corresponding annotations into smaller chunks

Use cases: 
1. `samtools` is unable to work with very large chromsosomes 

```
splitchrann.pl input.fasta output.fasta inputannotation outputannotations annotationtype maxchrsize
```
