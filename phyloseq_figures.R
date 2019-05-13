# Some pics
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(gplots)


ps <- readRDS('../Downloads/Phyloseq_filtered.rds')

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



row.names(map) <- map$Sample
map[7:9] <- lapply(map[7:9] , factor)

head(map)

meta <- sample_data(map)

ps <- phyloseq(otu, taxonomy, meta)

sample_data(ps)

#sample_data(ps)$Treatment <- gsub("BpO", "B+O", sample_data(ps)$Treatment)
#saveRDS(ps, 'dada2/Phyloseq_filtered.rds')


#png('pics/alpha_by_day.png')

#p <- plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"), color="Day") + theme_bw() + geom_point(size=2)

p <- plot_richness(ps, x="Vertposition", measures=c("Shannon", "Simpson")) + theme_bw() + geom_point(size=2)
p + geom_boxplot(data=p$data, aes(x=Vertposition, y=value, color=NULL), alpha = 0.1) + facet_grid( variable ~ litterorother, scales='free')

#dev.off()



#png('pics/tree.png', width=1200, height = 1200)
plot_tree(ps, ladderize="left", color="Vertposition", size='abundance') + coord_polar(theta="y") + 
  guides(color = guide_legend(override.aes = list(size=5))) + scale_color_brewer(palette = 'Dark2') 
#dev.off()

ps1 <- subset_samples(ps, Day == 1)
png('pics/tree_Day1.png', width=1200, height = 1200)
plot_tree(ps1, ladderize="left", color="Treatment", size='abundance') + coord_polar(theta="y") + 
  guides(color = guide_legend(override.aes = list(size=5))) + scale_color_brewer(palette = 'Dark2') 
dev.off()

ps4 <- subset_samples(ps, Day == 4)
png('pics/tree_Day4.png', width=1200, height = 1200)
plot_tree(ps4, ladderize="left", color="Treatment", size='abundance', label.tips = 'taxa_names') + coord_polar(theta="y") + 
  guides(color = guide_legend(override.aes = list(size=5))) + scale_color_brewer(palette = 'Dark2') 
dev.off()



# Barplot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~Vertposition, scales="free_x")+geom_bar(stat="identity")
plot_bar(ps.top20, x="Sample", fill="Phylum") + facet_wrap(~Vertposition, scales="free_x")+geom_bar(stat="identity")


# Modifying plot_bar code
mdf = psmelt(ps.top20)

gp <- colorRampPalette(brewer.pal(11, "Spectral"))
myPal <- gp(20)
myPal <- c("#ff6287",	 "#feade5",
           "#984f00",	 "#9b9fff",
           "#4f3c00",	 "#f2bd71",
           "#00315e",	 "#91d78c",
           "#6a7a00",	 "#78003f",
           "#ff7e67",	 "#b0c4ff",
           "#018a86",	 "#206dff",
           "#ae0dd2",	 "#7c000e",
           "#2a175e",	 "#ff42db",
           "#00961e",	 "#00b48d")


png('pics/Family_bar.png')
p = ggplot(mdf, aes_string(x = 'SampleID', y = 'Abundance', fill = "Family", color='Family'))
p = p + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values=myPal) + scale_color_manual(values=myPal)
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
p <- p + facet_wrap(  'Treatment', scales='free_x')
p	
dev.off()


png('pics/Phyla_bar.png')
p = ggplot(mdf, aes_string(x = 'SampleID', y = 'Abundance', fill = "Phylum", color='Phylum'))
p = p + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values=myPal) + scale_color_manual(values=myPal)
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
p <- p + facet_wrap(  'Treatment', scales='free_x')
p	
dev.off()

p = ggplot(mdf, aes_string(x = 'SampleID', y = 'Abundance', fill = "Phylum", color='Phylum'))
p = p + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values=myPal) + scale_color_manual(values=myPal)
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
p <- p + facet_wrap(  'Treatment')
p	+ coord_polar('y', start=0)



# Making Pie
psd <- subset_samples(ps, Day == 4)
psd <- prune_taxa(taxa_sums(psd) > 0, psd)
#psd <- prune_taxa(!is.na(taxa_sums(psd)), psd)


#psd <- merge_samples(psd, 'Treatment')
#sample_data(psd)$Treatment <- sample_names(psd)

top20 <- names(sort(taxa_sums(psd), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(psd, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

#pd <- merge_samples(ps.top20, 'Treatment')
#sample_data(pd)$Treatment <- sample_names(pd)



forPie = psmelt(ps.top20)

fp <- forPie %>% group_by(Treatment, Phylum, Family) %>% 
  summarize(Abundance=sum(Abundance)/n()) %>% 
  #		mutate(Abundance=Abundance / sum(Abundance))
  print()	

png('pics/Day4_Family_pie.png', width=1400, height=1100)		
p = ggplot(fp, aes(x='', y = Abundance, fill = Family))
p = p + geom_bar(stat = "identity", position = "stack") 
p = p + scale_fill_manual(values=myPal) + scale_color_manual(values=myPal)

#p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
p <- p + facet_wrap(  'Treatment')
p <- p + theme(axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               panel.border = element_blank(),
               panel.grid=element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_blank())
p	+ coord_polar('y', start=0)

dev.off()

#gp.ch = subset_taxa(ps, Phylum == "p__Firmicutes")
#plot_bar(gp.ch, fill="as.character(OTU)") + facet_wrap(~OTU, nrow=ntaxa(gp.ch))
#plot_bar(gp.ch, fill="Genus") + facet_wrap(~OTU, nrow=ntaxa(gp.ch))


# Looking at f__Lachnospiraceae

psd <- subset_samples(ps, Day == 4)
psd <- subset_taxa(psd, Family == 'f__Lachnospiraceae')

tax_table(psd)


p = ggplot(mdf, aes_string(x = 'SampleID', y = 'Abundance', fill = "Family", color='Family'))
p = p + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values=myPal) + scale_color_manual(values=myPal)
p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0)) 
p <- p + facet_wrap(  'Treatment', scales='free_x')
p	

plot_bar(psd, x="SampleID", fill="Species") + facet_wrap(~Treatment, scales="free_x")+geom_bar(stat="identity")



# Abundance Heatmap

#ps <- readRDS('dada2/Phyloseq_filtered.rds')

prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

# Are any phyla minimally represented?
#plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# or

myPrev <- prevdf %>% group_by(Phylum) %>% summarize(n=n(),sum=sum(Prevalence),mean=mean(Prevalence),max=max(Prevalence))
myPrevF <- prevdf %>% group_by(Family) %>% summarize(n=n(),sum=sum(Prevalence),mean=mean(Prevalence),max=max(Prevalence))
myPrev %>% print(n=100)


t <- data.frame(tax_table(ps))
head(t)

t %>% dplyr::filter(Family == 'f__mitochondria')


otu <- as.matrix(otu_table(ps))
otu <- log10((otu + 1))
bw <- colorRampPalette(c('white','blue'))

dim(otu)


metaColor <- sample_data(ps)
#rowColor <- rep(c('red','green'),14)

rownames(otu)
metaColor$SampleID

png('pics/Abundance_heatmap.png', height=1000, width=1400)
heatmap.2(otu, trace="none", col = bw, margin=c(6, 6), density.info = 'none',
          RowSideColors = metaColor$Vertposition, labCol=F, key.title=NA, key.xlab='Log Abundance', cexRow = 1.7)
dev.off()




png('pics/Abundance_heatmap_legend.png')
plot(NA,NA,xlim=c(0,150), ylim=c(0,150),
     xlab=NA, ylab=NA, axes=F)

rect(20,10,60,30,col='khaki')
rect(80,10,120,30,col='gold')
rect(20,40,60,60,col='coral')
rect(80,40,120,60,col='firebrick')
rect(20,70,60,90,col='deepskyblue')
rect(80,70,120,90,col='dodgerblue4')
rect(20,100,60,120,col='chartreuse')
rect(80,100,120,120,col='chartreuse4')

text(40,130,'Day 1')
text(100,130,'Day 4')
text(10,20, "CB+O3")
text(10,50, 'CB')
text(10,80, 'O3')
text(10,110, 'Air')

dev.off()
