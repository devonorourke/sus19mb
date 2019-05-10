library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(gplots)
library(vegan)

# Get data
otu <- read_tsv('Data_Files/otu.tsv')
meta <- read_tsv('Data_Files/metadata.tsv')


# The otu table needs to be transposed first
otu <- t(otu)

# The meta data has some columns that should be factors
meta[7:9] <- lapply(meta[7:9] , factor)


#Do your ordination
nmds_out <- metaMDS(otu, distance='bray', k=2, trymax=100 ,maxit=1000)


# Get the coordinates of the points 
NMDS <- data.frame(x=nmds_out$point[,1],y=nmds_out$point[,2])

# Add meta data to the data frame
NMDS <- cbind(NMDS, meta)
head(NMDS)


# Create the plot using ggplot.

# Note: "names" in the geom_point() command is not recognized by the ggplot2 package (Warning: Ignoring unknown aesthetics: names). 
# This is actually recognized by the plotly package and is the way to define information that is displayed 
# by hovering over points in the graph. 

p <- ggplot(data = NMDS, aes(x, y)) + 
  geom_point(size = 2.8, alpha=0.8, aes(color = Vertposition, shape = litterorother, names = Sample))+
  coord_fixed()+
  guides(shape=guide_legend(title=NULL))+
  theme(legend.position = "bottom")+
  ylab("NMDS Axis 2")+
  xlab("NMDS Axis 1")+
  theme_minimal() 


# Interactive plot with Plotly Package 
# The tooltip needs to match what you put in as names above in the ggplot
plotlyNMDS<-ggplotly(p=p, tooltip = "Sample")

# Show the plot
plotlyNMDS

# You can save the plot as html
htmlwidgets::saveWidget(plotlyNMDS, "test.html")














