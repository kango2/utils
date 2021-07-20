library(GenomicRanges)
library(Biostrings)
library(tidyverse)
windowwidth <- 1000000
fastafile <- "/pat/to/fasta/file.fasta.gz"
fasta <- readDNAStringSet(fastafile)
names(fasta) <- str_match(names(fasta), pattern = "^\\S+")
granges <- GenomicRanges::GRanges(seqnames = str_match(names(fasta), pattern = "^\\S+"), ranges = IRanges(start = 1, width = width(fasta)))
swindows <- unlist(GenomicRanges::slidingWindows(granges, width = windowwidth, step = windowwidth))
gc <- data.frame(alphabetFrequency(fasta[swindows]))
gc <- gc %>% mutate(gc = (G + C) / (A + C + G + T) )
gc <- bind_cols(data.frame(swindows) %>% select(chr = seqnames, start, end), select(gc,gc))
