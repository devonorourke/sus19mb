#Load McToolsR
library("mctoolsr", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")


#Load .txt file of biom, but I think you can just load the biom file
otu_table_fp1 <- "bact_mctoolsr.txt"
#Load mapping file
mapping_fp1 <- "bact_alldata_mapfile.txt"
fulltable1 <- load_taxa_table(otu_table_fp1, mapping_fp1) 

View(fulltable1)

sort(colSums(fulltable1$data_loaded))

#You can filter our data or just run this code and calculate distance matrix and ordination
Bac = filter_data(fulltable1)
#Create bray Curtis Distance Matrix
dmbac = calc_dm(Bac$data_loaded)
#Create NMDS Scores
ordbac <- calc_ordination(dmbac, 'nmds')
plot_ordination(Bac, ordbac, 'Vertposition', hulls = F) #hulls draw lines around groups, change to T if you want lines

tax_table_imp<-read.csv("Load a CSV with the Bacterial Taxa you want to create vectors of.csv")
fit1<-envfit(ordbac, tax_table_imp, perm = 999, na.rm = TRUE)
##save this data
Bact_vectors<-data.frame((fit1$vectors)$arrows, (fit1$vectors)$r, (fit1$vectors)$pvals)
colnames(Bact_vectors)<-c("NMS1","NMS2","Rsq","pval")
write.csv(Bact_vectors,"Bact_NMS_vector_results.csv", row.names=TRUE)
NMS_coordinates<-scores(ordbac,display="sites")


#create vectors of the major orders I analyzed previously
#add the proportion of total sequences of each vector to the "for_plotting" object below
for_ploting<-as.data.frame(cbind(NMS_coordinates,tax_table_imp))
str(tax_table_imp)

##COLORFUL
#now plot these data
png(file="FILE NAME YOU WANT TO SAVE.png", width = 6000, height = 6000, res = 1200)
par(mar=c(4,4,1,1))
plot(for_ploting$MDS2 ~ for_ploting$MDS1,
     xlab = "NMS1",
     ylab = "NMS2",
     font=2,
     font.lab=2,
     cex.axis=1,
     pch = c(8, 20, 17, 9)[as.factor(for_ploting$FunctionalGroup)], cex=.8, #change Functional Group to whatever variable you are interested in
     col =c("black","saddlebrown","green4","blue")[as.factor(for_ploting$FunctionalGroup)],  # different 'pch' types 
     data = for_ploting)
#This put 95% confidence ellipses around your centroid on the figure
#you need lty and col to be equal to the number of variables you have
#I had 4 functional groups so I have 4 colors
ordiellipse(ordbac, group=for_ploting$FunctinalGroup,kind = "se", 
            conf=0.95, lwd=1.9, lty= c(1,5,2,4), col =c("black","saddlebrown","green4","blue"))

#This plots vectors of the significant vectors you are interested in
#p.max means only variables with p-values less than 0.05 will be plotted
#Change chem1 to fit1 to look at bacterial or other groups you want
plot(chem1, p.max = 0.05, col = "black", cex = .6)
#Create a legend
#Change legend to your variables
legend(
  x ="bottomright",
  legend = c("All","Control","Grass","Mixed"), # for readability of legend
  pch = c(8, 20, 17, 9),
  col =c("black","saddlebrown","green4","blue"),
  cex = .80 # scale the legend to look attractively sized
)

dev.off()


#Perform multivariate analysis on your data
######Adonis
adonis(dmbac ~ Fertilizer, data=for_ploting, permutations=9999)


#Load a file with environmental variables
#Make sure that the samples all need to be in the same order
chem1 <- read.csv("Environmental Variable File .csv")
View(chem1)
#Load NMDS1 and NMDS2 File
nmdssec <- read.csv("nmdsbac.csv")
View(nmdssec)
#Load mapping file
supermapfile <- read.csv("nmdsmap.csv")
View(supermapfile)
library("vegan", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

chem = envfit(nmdssec, chem1, na.rm = TRUE, perm = 999)
chem

