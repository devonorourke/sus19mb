library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(gplots)


#ps <- readRDS('../Downloads/Phyloseq_filtered.rds')

map <- read_csv('bact_alldata_mapfile.csv')
tax <- read_csv('bact_alldata_taxatable_wTax.csv')


otu <- as.data.frame(select(tax, -X1,-taxonomy))
row.names(otu) <- tax$X1
#head(otu)
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

head(map)

meta <- sample_data(map)

ps <- phyloseq(otu, taxonomy, meta)

sample_data(ps)



otut <- as.data.frame(otu_table(ps, matrix))
otut <- t(otut)

nmds_out <- metaMDS(otut, distance='bray', k=2, trymax=100 ,maxit=1000)



NMDS <- data.frame(x=nmds_out$point[,1],y=nmds_out$point[,2])

NMDS <- cbind(NMDS, meta)
head(NMDS)


# Note: "names" in the geom_point() command is not recognized by the ggplot2 package (Warning: Ignoring unknown aesthetics: names). 
# This is actually recognized by the plotly package and is the way to define information that is displayed 
# by hovering over points in the graph. 

p <- ggplot(data = NMDS, aes(x, y)) + 
  geom_point(size = 2.8, alpha=0.8, aes(color = Vertposition, shape = litterorother, names = Sample))+
#  geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
#               arrow = arrow(length = unit(0.2, "cm")),colour="grey") + 
#  geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species), size=3.8)+
  coord_fixed()+
#  scale_color_gradient2(low = "yellow", high = "red", mid = "orange",  midpoint = 9.32)+
#  guides(colour=guide_colourbar(title="Alpha Acid %"))+
  guides(shape=guide_legend(title=NULL))+
  theme(legend.position = "bottom")+
  ylab("NMDS Axis 2")+
  xlab("NMDS Axis 1")+
  theme_minimal()+
  theme(axis.text=element_text(family="A"))+
  theme(axis.title.y = element_text(family="A"))+
  theme(legend.text=element_text(face = "italic", family="A", size =12), legend.title= element_text(family="A", size=15))+
  theme(text = element_text(family="A"))

p

#Interactive plot with Plotly Package 
plotlyNMDS<-ggplotly(p=p, tooltip = "Sample")

plotlyNMDS












