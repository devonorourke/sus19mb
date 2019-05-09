# Script to make the files needed for emperor 
# from a phyloseq object

library(phyloseq)
library(tidyverse)


# Number of principle components to use
numPC <- 5

# Read phyloseq object
ps <- readRDS('../Downloads/Phyloseq_filtered.rds')


# Do MDS calculations
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")


a <- out.bc.log$vectors
a <- as.data.frame(a)
a <- a[,1:numPC]
a$pc <- row.names(a)
a <- as_tibble(a)
a <- select(a, pc, everything())
#names(a)
names(a)[1] <- 'pc vector number'


# Write PCA data to file
write.table(a, 'emperor/pca.dat', quote=F, sep="\t", row.names = F)


# Prep eigenvalue and percent explained data
b <- as.data.frame(out.bc.log$values)

eig <- b$Eigenvalues
eig <- eig[1:numPC]
eig <- c('eigvals',eig)

varexp <- b$Relative_eig
varexp <- 100 * varexp
varexp <- varexp[1:numPC]
varexp <- c('% variation explained', varexp)



# Add that to PCA data file
cat('\n\n', 
    paste0(eig, collapse='\t'), '\n',
    paste0(varexp, collapse='\t'),
    file='emperor/pca.dat','\n', append=T)

# Do meta data file

meta <- sample_data(ps)
meta <- as_tibble(meta)
meta <- rename(meta, '#SampleID'=SampleID)
meta <- dplyr::select(meta, `#SampleID`, everything())

write_tsv(meta, 'emperor/meta.txt')	


# On the command line,
# conda activate emperor
# make_emperor.py -i pca.dat -m meta.dat -o outdir
# conda deactivate
