library(tidyverse)
library(plotly)
library(ggplotlyExtra)
a <- read_delim("~/Downloads/vsGCF_004115215.2_mOrnAna1.pri.v4_genomic.target.chainpairs.tab", delim = "\t", col_names = F)
cnames <- read_delim("~/Downloads/platypuschrnames.txt", delim="\t", col_names = F)
clen <- read_delim("~/Downloads/platypussize.txt", delim = "\t", col_names = F)
colnames(cnames) <- c("cname", "cacc", "nothing", "contigid", "nothng")
colnames(clen) <- c("contigid", "len")

cinfo <- left_join(cnames,clen)
cinfo <- cinfo %>% rowid_to_column()
cinfo <- cinfo %>% mutate(offset = lag(len), offset = replace_na(offset, 0), offset = cumsum(offset))
cinfo <- mutate(cinfo, cname = str_replace(cname, "Chromosome ", "chr"))
colnames(a) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")
a <- left_join(a, dplyr::select(cinfo, tacc = contigid, tchr = cname, trank = rowid, toffset=offset), by = "tacc")
a <- left_join(a, dplyr::select(cinfo, qacc = contigid, qchr = cname, qrank = rowid, qoffset=offset), by = "qacc")
a <- mutate(a, tlen=tend-tstart, qlen=qend-qstart)

minalnlen <- 1000
g <-   drop_na(a) %>% dplyr::filter(qlen>minalnlen | tlen>minalnlen) %>%
  mutate(adjstart = if_else(strand=="+", qstart, qend), adjend=if_else(strand=="+",qend,qstart)) %>% 
  ggplot() + 
  geom_segment(aes(x=tstart+toffset, y=adjstart+qoffset, xend=tend+toffset, yend=adjend+qoffset,color=strand)) + 
  theme_bw() +
  scale_x_continuous(breaks = drop_na(a) %>% arrange(trank) %>% pull(toffset) %>% unique(), minor_breaks = NULL, labels = drop_na(a) %>% arrange(trank) %>% pull(tchr) %>% unique()) + 
  scale_y_continuous(breaks = drop_na(a) %>% arrange(qrank) %>% pull(qoffset) %>% unique(), minor_breaks = NULL, labels = drop_na(a) %>% arrange(qrank) %>% pull(qchr) %>% unique())

ggplotly(g)


b <- drop_na(a) %>% group_by(tchr,qchr) %>% summarise(tsum = sum(tlen), qsum = sum(qlen))
b <- left_join(b, select(cinfo, tchr = cname, tlen = len))
b <- left_join(b, select(cinfo, qchr = cname, qlen = len))
pdf("~/Downloads/platypusselfaln.pdf", 10, 10)
b %>% ggplot(aes(x=qchr,y=tchr, fill=qsum/qlen)) + geom_tile() + scale_fill_viridis_c() + theme_bw() + theme(text = element_text(size=12), axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

for (i in cinfo$cname) {
  cband <- cinfo %>% mutate(chromStart=0, name="gneg", score=1,strand="+") %>% select(chrom=cname,chromStart, chromEnd=len,name,score,strand)
  cband <- data.frame(cband)
  minalnlen <- 10000
  print(i)
  pdf(paste("~/Downloads/",i,".pdf",sep = ""), height = 5, width = 5)
  circos.initializeWithIdeogram(cband)
  circos.genomicLink(data.frame(filter(c,tchr==i & tlen>minalnlen) %>% select(tchr,tstart,tend)), data.frame(filter(c,tchr==i & tlen>minalnlen) %>% select(qchr,qstart,qend)), col = "green")
  dev.off()
}
