# Interactive trees

library(phyloseq)
library(tidyverse)
library(plotly)
library(data.table)
library(scales)

map <- read_csv('bact_alldata_mapfile.csv')
tax <- read_csv('bact_alldata_taxatable_wTax.csv')


otu <- as.data.frame(select(tax, -X1,-taxonomy))
row.names(otu) <- tax$X1
otu <- otu_table(otu, taxa_are_rows = T)

taxonomy <- data.frame(taxonomy=tax$taxonomy)
row.names(taxonomy)<- tax$X1

# Split taxonomy into separate columns
taxonomy <- data.frame(separate(taxonomy, col=taxonomy, into= c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                sep = ";"))

# Change spaces to NA for missing data
taxonomy <- apply(taxonomy, 2, function(x) gsub("^$|^ $", NA, x))

taxonomy <- as.matrix(taxonomy)
#head(taxonomy)

taxonomy <- tax_table(taxonomy)

names(map)[1] <- 'Sample'
row.names(map) <- map$Sample
map[7:9] <- lapply(map[7:9] , factor)

meta <- sample_data(map)


input_tree <- read_tree('Data_Files/pruned_tree.tre')

ps <- phyloseq(otu, taxonomy, meta, input_tree)


#ps_H <- readRDS('../../Downloads/Phyloseq_filtered.rds')


#ps_alpha <- subset_taxa(ps, Class=='c__Alphaproteobacteria')
ps_ricket <- subset_taxa(ps, Order=='o__Rickettsiales')



Lowest_taxa <- do.call(paste, c(as.data.frame(tax_table(ps_ricket))[1:7], sep=";"))

Lowest_taxa <- gsub("(;NA)*$", "", Lowest_taxa, perl=T)
Lowest_taxa <- gsub("(;.__)*$", "", Lowest_taxa, perl=T)
Lowest_taxa <- gsub("^.*;", "", Lowest_taxa, perl=T)

colnames(tax_table(ps_ricket))[1] <- 'Lowest_taxa'
tax_table(ps_ricket)[,'Lowest_taxa'] <- Lowest_taxa


pt <- interactive_plot_tree(ps_ricket, ladderize="left", color="Vertposition", size='abundance', nodelabf=nodeplotblank, tooltip = 'Lowest_taxa')  + 
    scale_color_manual(values=c("black","saddlebrown","tan2","green1","green4"))
  
guides(color = guide_legend(override.aes = list(size=4))) 
  scale_color_brewer(palette = 'Dark2') 



plotly_tree<-ggplotly(p=pt, tooltip = 'Lowest_taxa') 

plotly_tree

htmlwidgets::saveWidget(plotly_tree, 'tree_with_hover.html')



