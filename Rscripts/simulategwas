# Load the necessary library
library(ggplot2)
library(ggmanh)


# Create synthetic data for humans
chrLen <- c(248387328,242696752,201105948,193574945,182045439,172126628,160567428,146259331,150617247,134758134,135127769,133324548,113566686,101161492,99753195,96330374,84276897,80542538,61707364,66210255,45090682,51324926,154259566,62460029)
chrNames <- paste("chr", c(1:22, "X", "Y"), sep = "")

chrInfo <- data.frame(chrNames, chrLen)

gwas_pval <- 5e-8
prs_pval <- 1e-3


n <- 100000
gwas_data <- sample(chrInfo$chrNames, n, replace = T, prob = chrInfo$chrLen/sum(chrInfo$chrLen))
gwas_data <- left_join(tibble(chrNames = gwas_data), chrInfo) %>% arrange(chrNames)
gwas_data <- gwas_data %>% rowwise() %>% mutate(Position = sample(chrLen, 1))
gwas_data$p_values <- runif(n, min = prs_pval, max = 1)
##beta distribution of p-values
#gwas_data$p_values <- rbeta(n, shape1 = 5, shape2 = 1)^7
gwas_data <- data.frame(gwas_data)


gwas_significant_n <- 5
gwas_sig <- sample(n, gwas_significant_n)
gwas_sig <- gwas_data[gwas_sig,]
gwas_sig$p_values <- runif(gwas_significant_n, min = 1e-10, max = 1e-8)
tmp <- list()  
for (i in 1:nrow(gwas_sig)) {
  neighbours <- 40
  pos <- sample((gwas_sig[i,]$Position - 10000):(gwas_sig[i,]$Position + 10000), neighbours)
  pval <- rbeta(neighbours, shape1 = 2, shape2 = 1)^3
  t <- data.frame(chrNames = gwas_sig[i,]$chrNames, chrLen = gwas_sig[i,]$chrLen, Position = pos, p_values = pval)
  tmp[[i]] <- t
}
tmp <- bind_rows(tmp)
gwas_sig <- bind_rows(gwas_sig, tmp)



prs_significant_n <- 10
prs_sig <- sample(n, prs_significant_n)
prs_sig <- gwas_data[prs_sig,]
prs_sig$p_values <- runif(prs_significant_n, min = gwas_pval, max = prs_pval)
tmp <- list()  
for (i in 1:nrow(prs_sig)) {
  neighbours <- 1000
  pos <- sample((prs_sig[i,]$Position - 10000):(prs_sig[i,]$Position + 10000), neighbours)
  pval <- rbeta(neighbours, shape1 = 5, shape2 = 1)^7
  t <- data.frame(chrNames = prs_sig[i,]$chrNames, chrLen = prs_sig[i,]$chrLen, Position = pos, p_values = pval)
  tmp[[i]] <- t
}
tmp <- bind_rows(tmp)
prs_sig <- bind_rows(prs_sig, tmp)

gwas_data <- bind_rows(gwas_data, gwas_sig)
gwas_data <- bind_rows(gwas_data, prs_sig)

gwas_data$chrNames <- factor(gwas_data$chrNames, pull(chrInfo, chrNames))

m <- manhattan_plot(x = gwas_data, pval.colname = "p_values", chr.colname = "chrNames", pos.colname = "Position", 
               plot.title = "GWAS Example", 
               y.label = "-log10(P)",
               signif = c(5e-08, 1e-03), point.size = 2 
               )
m + geom_hline(yintercept = 3, color = "red") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
