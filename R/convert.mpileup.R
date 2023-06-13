#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# read and convert mpileups.csv to likelihoods of each of the 10 genotypes
#--------------------------------------------------------------------------------------------------------------------------------
source('functions.R')

d1 <- read.csv('../data/AADRv50.PigmAlleles.mpileups.tsv',sep='\t')
d2 <- read.csv('../data/Lehti.PigmAlleles.mpileups.tsv',sep='\t')

x1 <- convert.mpileups(d1)
x2 <- convert.mpileups(d2)

summary <- rbind(x1$summary,x2$summary)
meta <- rbind(x1$meta, x2$meta)

write.csv(summary, '../data/summary.csv')
write.csv(meta, '../data/meta.csv', row.names=F)
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
dim(d1)

install.packages("gtools")
