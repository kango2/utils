# Wrap sequence into 60 character strings in Fasta file

Fasta files can sometimes contain the whole sequence in a single line. `samtools faidx` and `subread-index` utilities don't like this for good reasons. When you want to peek into the file, it will print the entire sequence and fill up the terminal. Also, opossum chromosome is >750Mb in size and Indian Muntjac, type of deer, [has only three chromosomes](https://doi.org/10.1038/s42003-020-1096-9) for a mammalian sized genome. Following Perl oneliner will read in your Fasta file and write to STDOUT sequences in 60 character long strings.

```
perl -lne 'if ($_=~/(>\S+)/){ $header=$1; $seq{$header}=""; \
push(@headers,$header)}else{$sequence{$header}.=$_} \
END { foreach $h (@headers) { \
print "$h"; for ($i=0;$i<length($sequence{$h});$i+=60){ \
print substr($sequence{$h},$i,60)}}}' \
input.fasta >output.fasta
```
