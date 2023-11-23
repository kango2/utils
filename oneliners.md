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

# Create `bed` file with fixed size windows from SAM/BAM/CRAM file

Often times, we need to calculate coverage/depth etc from SAM/BAM?CRAM files over a sliding window of fixed size. We can create a `bed` file from header information in SAM/BAM files to quickly generate such file for use.

```
export WINDOW=1000000
samtools view -H INPUT.bam | perl -lne 'if ($_=~/SN:(\S+)\tLN:(\d+)/){ $c=$1;$l=$2; for ($i=0;$i<$l;$i+=$ENV{"WINDOW"}) { print "$c\t$i\t". ((($i+$ENV{"WINDOW"}) > $l) ? $l : ($i+$ENV{"WINDOW"}))  }} '
```
Alternately, one can use `.faidx` file as well to achieve the same outcome.

```
export WINDOW=1000000
cat genome.fasta.faidx | perl -lne '@a=split("\t",$_); $c=$a[0];$l=$a[1]; for ($i=0;$i<$l;$i+=$ENV{"WINDOW"}) { print "$c\t$i\t". ((($i+$ENV{"WINDOW"}) > $l) ? $l : ($i+$ENV{"WINDOW"}))  }} '
```

# Get summary statistics from a fasta file

```
inputfile=/path/to/input.fasta
echo ${inputfile},$(cat ${inputfile} | perl -lne 'if($_=~/>/){push(@l,length($s));$s="";}else{$s.=$_;$total+=length($_);$gc+=$_=~tr/[GC]//;$ncount+=$_=~tr/N//;}END{push(@l,length($s));shift@l;@l=sort{$b<=>$a}@l;$n9p=0;$n5p=0;foreach(@l){$n+=$_;if($n>=$total*0.5 && $n5p==0){$n5p=$_}elsif($n>=$total*0.9 && $n9p==0){$n9p=$_}}; print "$total,".scalar(@l).",".(sprintf "%.02f", (($gc*100)/$total)).",$ncount,".(sprintf "%.02f", $total/scalar(@l)).",".($l[int(scalar(@l)/2)]).",$l[0],$l[$#l],$n5p,$n9p"}')
```

For fastq files, we can convert them to fasta on the fly and use them as follows

```
inputfile=/path/to/input.fasta
echo ${inputfile},$(cat ${inputfile} | sed -n '1~4s/^@/>/p;2~4p' | perl -lne 'if($_=~/>/){push(@l,length($s));$s="";}else{$s.=$_;$total+=length($_);$gc+=$_=~tr/[GC]//;$ncount+=$_=~tr/N//;}END{push(@l,length($s));shift@l;@l=sort{$b<=>$a}@l;$n9p=0;$n5p=0;foreach(@l){$n+=$_;if($n>=$total*0.5 && $n5p==0){$n5p=$_}elsif($n>=$total*0.9 && $n9p==0){$n9p=$_}}; print "$total,".scalar(@l).",".(sprintf "%.02f", (($gc*100)/$total)).",$ncount,".(sprintf "%.02f", $total/scalar(@l)).",".($l[int(scalar(@l)/2)]).",$l[0],$l[$#l],$n5p,$n9p"}')
```

NOTES:
1. If the file is a zipped file, then use `zcat` instead of `cat`
2. Output comma separated columns: `inputfile,# of bases,# of reads,GC,Ncount,meanLen,medianLen,longestLen,shortestLen,N50,N90

