#Sources for this script:
#https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#content
#http://joey711.github.io/phyloseq/import-data#import_functions
#https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
#http://joey711.github.io/phyloseq/plot_ordination-examples.html
#https://rstudio-pubs-static.s3.amazonaws.com/713954_d40760746cd3402fb0d1012ff67e5ab9.html
#https://rstudio-pubs-static.s3.amazonaws.com/330760_8bba830836324bf6b100d4e76f49e3d2.html#alpha-diversity

setwd("C:/Users/User/Downloads/Circadiano_micro/Analises_Microbiota_Katia_R")

library("readxl")			# For importing data from Excel file
library(ggplot2)			# For data visualization
library(vegan)				# Ecological data analysis
library(RColorBrewer)		# Color palettes
library(reshape2)			# Data transformation
library(scales)				# Scaling data
library(data.table)			# Fast data manipulation
library(microbiome)			# Microbiome analysis
library(dplyr)				# Data manipulation
library(phyloseq)			# Handling microbiome data
library(DT)					#for interactive tables
library(microbiomeutilities) # Additional utilities for microbiome analysis
library(mirlyn)				# Data rarefaction
library(tibble)				# Handling data frames efficiently
library(MetaCycle)			# Searches for oscillation in the dataset following a given period
library(GUniFrac)			# Data rarefaction
#library(devtools)			# Used to install pairwiseAdonis from github
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") #downloading pairwiseAdonis
library("pairwiseAdonis")	# For pairwise PERMANOVA tests
#BiocManager::install("DESeq2")
library(DESeq2)				# For differencial abundance analysis
library(pals)				# Color palettes that have different colors next to each other rather than a continuous spectrum
library(Polychrome)			# Color palettes that have different colors next to each other rather than a continuous spectrum



##### Preparing phyloseq object #####

# Load ASV (Amplicon Sequence Variant) table from Excel
otu_mat<- read_excel("asv_otu5_paraR_coral.xlsx")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")     # Set ASV column as row names
otu_mat <- as.matrix(otu_mat)    # Convert to matrix format
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)  # Create phyloseq OTU table

#Loading taxonomy spreadsheet
tax_mat<- read_excel("asv_tax5_paraR.xlsx")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
tax_mat <- as.matrix(tax_mat)
TAX = tax_table(tax_mat)  # Create phyloseq taxonomy table

#Loading metadata spreadsheet
samples_df <- read_excel("metadados_paraR_coral.xlsx")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
labels <- row.names(samples_df)
samples_df <- cbind (labels, samples_df) #para facilitar juntar a diversidade por iteracao depois de fazer mirl
samples = sample_data(samples_df)   # Convert metadata to phyloseq format
samples_df$Time <- as.integer(samples_df$Time)
samples_df$Depth_CORRIGIDA<-as.integer(samples_df$Depth_CORRIGIDA)
samples_df$Salinity_CORRIGIDA<-as.integer(samples_df$Salinity_CORRIGIDA)
samples_df$Temperature<-as.integer(samples_df$Temperature)

micro_coral <- phyloseq(OTU, TAX, samples)
sample_variables(micro_coral)
print(micro_coral) # confirm whether number of taxa is the same as the number of ASVs in the spreadsheets, if number of samples is the same in metadata and ASV spreadsheets, if taxonomic ranks are 6 in taxonomy spreadsheet and if there are 8 variables in the metadata spreadsheet


# Check if there is any ASV which is absent in all samples
any(taxa_sums(micro_coral) == 0)

#if true, remove them
micro_coral_1 <- prune_taxa(taxa_sums(micro_coral) > 0, micro_coral)
any(taxa_sums(micro_coral_1) == 0) #this line has to return FALSE

rank_names(micro_coral_1) # check the taxonomic rank information

# Remove unwanted taxa (e.g., Chloroplasts, Mitochondria, Archaea)
micro_coral_2 <- subset_taxa(micro_coral_1,Order!="Chloroplast")
micro_coral_3 <- subset_taxa(micro_coral_2,Family!="Mitochondria")
micro_coral_5 <- subset_taxa(micro_coral_3,Kingdom!="Archaea")
micro_coral_4 <- subset_taxa(micro_coral_5,Kingdom!="NA")

any(tax_table(micro_coral_4) == 'Archaea') 
any(tax_table(micro_coral_4) == 'Chloroplast')
any(tax_table(micro_coral_4) == 'Mitochondria')
any(tax_table(micro_coral_4)[,1] == 'NA')
#To confirm there are no more items belonging to the omitted categories (it has to return FALSE)


# Export processed OTU and taxonomy tables to CSV files
otu.df <- as.data.frame(otu_table(micro_coral_4))  # make a dataframe for ASV information.
tax.df <- as.data.frame(tax_table(micro_coral_4))  
write.csv(otu.df, "./otu.df_phyloseq_4.csv", row.names=TRUE)
write.csv(tax.df, "./tax.df_phyloseq_4.csv", row.names=TRUE)

ntaxa(micro_coral_4) #Counting the amount of ASVs in the dataset


##### Filtering and rarefaction #####

# Remove ASVs with less than 5 counts across all samples
reads_per_OTU4 <- taxa_sums(micro_coral_4)
print(length(reads_per_OTU4[reads_per_OTU4 < 5]))

micro_coral_4.f <- filter_taxa(micro_coral_4, function(x) sum(x > 5) > 0, TRUE)

ntaxa(micro_coral_4.f)
#[1] 25833

# Exporting ASV and taxonomy tables to excel
otu.df.4.f <- as.data.frame(otu_table(micro_coral_4.f))  # make a dataframe for OTU information.
tax.df.4.f <- as.data.frame(tax_table(micro_coral_4.f))  
write.csv(otu.df.4.f, "./otu.df_phyloseq_4.f.csv", row.names=TRUE)
write.csv(tax.df.4.f, "./tax.df_phyloseq_4.f.csv", row.names=TRUE)

# Rarefaction analysis
Rarefy_whole_rep_circmicro.f <- rarefy_whole_rep(micro_coral_4.f, rep = 100)

# Generate rarefaction curve
Rarecurve_circmicro.f <- rarecurve(Rarefy_whole_rep_circmicro.f, sample = "Sample")


# Generate rarefaction curves for each species
#Mde
test_M_decactis.f <- Rarefy_whole_rep_circmicro.f [Rarefy_whole_rep_circmicro.f$Species=="M_decactis",]
Rarecurve_Mde.f <- rarecurve(test_M_decactis.f, sample = "Sample")
Rarecurve_Mde.f #exported as pdf in A4 landscape

#Mhi
test_M_hispida.f <- Rarefy_whole_rep_circmicro.f [Rarefy_whole_rep_circmicro.f$Species=="M_hispida",]
Rarecurve_Mhi.f <- rarecurve(test_M_hispida.f, sample = "Sample")
Rarecurve_Mhi.f #exported as pdf in A4 landscape

#Tc
test_T_coccinea.f <- Rarefy_whole_rep_circmicro.f [Rarefy_whole_rep_circmicro.f$Species=="T_coccinea",]
Rarecurve_Tco.f <- rarecurve(test_T_coccinea.f, sample = "Sample")
Rarecurve_Tco.f #exported as pdf in A4 landscape

#Tt
test_T_tagusensis.f <- Rarefy_whole_rep_circmicro.f [Rarefy_whole_rep_circmicro.f$Species=="T_tagusensis",]
Rarecurve_Tta.f <- rarecurve(test_T_tagusensis.f, sample = "Sample")
Rarecurve_Tta.f #exported as pdf in A4 landscape


## Removing samples that had few reads or did not reach plateau 

sort(sample_sums(micro_coral_4.f)) #prints total read count per sample (library size)

micro_coral_4.f.sub <-subset_samples(micro_coral_4.f, #phyloseq object
                                    sample_names(micro_coral_4.f) !="Mde45" & sample_names(micro_coral_4.f) !="Mde33") #removing Mde45 e Mde33 because few reads (even though Mde33 reached plateau)

sort(sample_sums(micro_coral_4.f.sub)) 
ntaxa(micro_coral_4.f.sub)

## Rarefying dataset
# Disclaimer: Rarefaction may lead to differences in abundances, causing different results. Use alphadiv_mirl_4.f.sub.csv to get the same result.

mirl_micro_coral_4.f.sub <- mirl(micro_coral_4.f.sub, libsize=20371, rep = 100, set.seed = 120, replace=FALSE)

#ASVs with no counts in any samples will not be trimmed in the end, analysis will be performed withou replacement (i.e. default)

##### Alpha diversity #####

## Calculate alpha diversity metrics using mirlyn package
# Generates dataframe of alpha-diversity metric from mirl_object
alphadiv_df_mirl_4.f.sub <- alphadivDF(mirl_micro_coral_4.f.sub)

# Exporting rarefied dataset to excel file
write.csv(alphadiv_df_mirl_4.f.sub, "alphadiv_mirl_4.f.sub.csv")
head(alphadiv_df_mirl_4.f.sub)


##Plotting alpha diversity as a box plot (not included in the manuscript)
p <- alphawichVis(alphadiv_df_mirl_4.f.sub, xvar = "Species", colorvar = "Light")
p <- p + geom_point(aes(color = Light)) + geom_boxplot(aes(fill = "Light"), alpha = 0) + scale_fill_manual(values = c("#CBD588", "#5F7FC7")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#  scale_x_discrete(expand = c(2, 4))
print(p)
ggsave("./mirl_micro_coral_4_f_sub_Alpha_Diversity.pdf", height = 6, width = 7)


## Plotting alpha diversity as box plot with one data point per sample rather than one point per mirl function replicate

#Calculating average diversity between the 100 replicates belonging to each sample
library(dplyr)

#	Because there were no sample label but it is possible to know (printing mirl_micro_coral_4.f.sub) that there are 69 samples. Also,
#	each sample has 100 replicates and in the exported excel sheet (a data frame in R), it is possible to see that mirl samples are 
#	in the same order. Hence, the solution was including an index column in this data frame. This column contains numbers from 1 to 69 until
#	the 69th sample and, on the 70th line (which corresponds to the first sample again), it goes back to 1. Each sample will have an index
#	between 1 and 69 and diversity values will be grouped by sample.

replicate_index = c(1:69)
alphadiv_df_mirl_4.f.sub.i <- cbind (replicate_index, alphadiv_df_mirl_4.f.sub)
write.csv(alphadiv_df_mirl_4.f.sub.i, "alphadiv_mirl_4.f.sub.indexado.csv") #file with indexes, each sample has the same number


# Average alpha diversity by sample including other columns from the metadata table
av.a.diver_mirl_4.f.sub_Sa <- alphadiv_df_mirl_4.f.sub.i %>% 
  group_by(replicate_index, Species, Time, Day, Light, Symbiont, Temperature, Depth_CORRIGIDA, Salinity_CORRIGIDA) %>%
  dplyr::summarise(SD = sd(DiversityIndex), mean = mean(DiversityIndex)) %>%
  as.data.frame()

# Exporting file including alpha diversity averages per sample
write.csv(av.a.diver_mirl_4.f.sub_Sa, "alphadiv_mirl_4.f.sub_bySample.csv") 


# Plotting box plot
ggplot(av.a.diver_mirl_4.f.sub_Sa, aes(x = Species, y = mean, color = Light)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 3) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.shape = NA) +
  labs(title = "Shannon diversity by species and light",
       x = "Species",
       y = "Shannon",
       color = "Light") +
  theme_minimal()
ggsave("./mirl_micro_coral_4_f_sub_Alpha_Diversity_1dotPerSample.pdf", height = 6, width = 7)


##### Normality tests (performed with original data and log transformed data) #####

##Histogram

hist(av.a.diver_mirl_4.f.sub_Sa$mean, col='#f0868c', main='Shannon\n1 sample = 100 mirl iterations averaged') 
hist(log(av.a.diver_mirl_4.f.sub_Sa$mean), col='#f0868c', main='Shannon\n1 sample = 100 mirl iterations averaged')

##Shapiro-Wilk test

shapiro.test(av.a.diver_mirl_4.f.sub_Sa$mean)
#W = 0.95949, p-value = 0.025
#not normally distributed

column_values_Shannon_log <- log(av.a.diver_mirl_4.f.sub_Sa$mean)
shapiro.test(column_values_Shannon_log)
#W = 0.95984, p-value = 0.02613
#not normally distributed


##### Kruskal-Wallis test for alpha diversity within each species comparing day and night #####

## Comparing species
kruskal.test(av.a.diver_mirl_4.f.sub_Sa$mean ~ sample_data(av.a.diver_mirl_4.f.sub_Sa)$Species)
# Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 26.055, df = 3, p-value = 9.286e-06
pairwise.wilcox.test(av.a.diver_mirl_4.f.sub_Sa$mean, sample_data(av.a.diver_mirl_4.f.sub_Sa)$Species, p.adj = "bonf")
# Pairwise comparisons using Wilcoxon rank sum exact test 
# data:  av.a.diver_mirl_4.f.sub_Sa$mean and sample_data(av.a.diver_mirl_4.f.sub_Sa)$Species 
#             M_decactis M_hispida T_coccinea
#  M_hispida  0.02514    -          -         
#  T_coccinea 0.00012   1.00000     -         
#  T_tagusensis 5e-06   0.07555     0.70921   
#  P value adjustment method: bonferroni 


## Comparing day and night in each species
# Kruskal-Wallis test for overall comparison among species and treatments
kruskal.test(av.a.diver_mirl_4.f.sub_Sa$mean ~ interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa)$Species, sample_data(av.a.diver_mirl_4.f.sub_Sa)$Light))
#  Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 28.787, df = 7, p-value = 0.0001582


# Pairwise Wilcoxon tests for comparisons among treatments within each species
pairwise.wilcox.test(av.a.diver_mirl_4.f.sub_Sa$mean, interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa)$Species, sample_data(av.a.diver_mirl_4.f.sub_Sa)$Light), p.adj = "bonf")
#                M_decactis.Day M_hispida.Day T_coccinea.Day T_tagusensis.Day M_decactis.Night M_hispida.Night T_coccinea.Night
#  M_hispida.Day      0.5787         -             -              -                -                -               -               
#  T_coccinea.Day     0.0691         1.0000        -              -                -                -               -               
#  T_tagusensis.Day   0.0161         1.0000        1.0000         -                -                -               -               
#  M_decactis.Night   1.0000         1.0000        0.0438         0.0092           -                -               -               
#  M_hispida.Night    1.0000         1.0000        1.0000         0.5252           1.0000           -               -               
#  T_coccinea.Night   0.1037         1.0000        1.0000         0.2972           0.4262           1.0000          -               
#  T_tagusensis.Night 0.0438         1.0000        1.0000         1.0000           0.0691           1.0000          1.0000          
#Significantly different pairs (<0.05): Tta.Day x Mde.Night; Tco.Day x Mde.Night; Tta.Night x Mde.Day; Tta.Day x Mde.Day --- none in the same species

## Subsetting dataset in order to compare day vs night and time in each species
av.a.diver_mirl_4.f.sub_Sa_Mde <- subset(av.a.diver_mirl_4.f.sub_Sa, Species == 'M_decactis')
av.a.diver_mirl_4.f.sub_Sa_Mhi <- subset(av.a.diver_mirl_4.f.sub_Sa, Species == 'M_hispida')
av.a.diver_mirl_4.f.sub_Sa_Tco <- subset(av.a.diver_mirl_4.f.sub_Sa, Species == 'T_coccinea')
av.a.diver_mirl_4.f.sub_Sa_Tta <- subset(av.a.diver_mirl_4.f.sub_Sa, Species == 'T_tagusensis')


# M_decactis
# Kruskal-Wallis test for overall comparison among species and treatments
kruskal.test(av.a.diver_mirl_4.f.sub_Sa_Mde$mean ~ interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa_Mde)$Light))
#  Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 0.89338, df = 1, p-value = 0.3446


# M_hispida
# Kruskal-Wallis test for overall comparison among species and treatments
kruskal.test(av.a.diver_mirl_4.f.sub_Sa_Mhi$mean ~ interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa_Mhi)$Light))
#  Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 0.0092593, df = 1, p-value = 0.9233


# T_tagusensis
# Kruskal-Wallis test for overall comparison among species and treatments
kruskal.test(av.a.diver_mirl_4.f.sub_Sa_Tta$mean ~ interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa_Tta)$Light))
#  Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 1.0312, df = 1, p-value = 0.3099


# T_coccinea
# Kruskal-Wallis test for overall comparison among species and treatments
kruskal.test(av.a.diver_mirl_4.f.sub_Sa_Tco$mean ~ interaction(sample_data(av.a.diver_mirl_4.f.sub_Sa_Tco)$Light))
#  Kruskal-Wallis rank sum test - Kruskal-Wallis chi-squared = 2.6686, df = 1, p-value = 0.1023



#### Testing whether alpha diversity cycles in 24 or 12-hour periods
library(MetaCycle)
# Tco
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Tco.csv", outdir = "MetaCycle_alphadiv_Tco02", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Tco", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Tco.csv", outdir = "MetaCycle_alphadiv_Tco0212", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Tco12", parallelize = FALSE, nCores = 1, inDF = NULL)
# Tta
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Tta.csv", outdir = "MetaCycle_alphadiv_Tta", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Tta", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Tta.csv", outdir = "MetaCycle_alphadiv_Tta12", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Tta12", parallelize = FALSE, nCores = 1, inDF = NULL)
# Mde
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Mde.csv", outdir = "MetaCycle_alphadiv_Mde", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Mde", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Mde.csv", outdir = "MetaCycle_alphadiv_Mde12", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Mde12", parallelize = FALSE, nCores = 1, inDF = NULL)
# Mhi
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Mhi.csv", outdir = "MetaCycle_alphadiv_Mhi", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Mhi", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("SuppMatt_R_script/Period_search_methodologies/mirl_alphadiv_forMetacycle_Mhi.csv", outdir = "MetaCycle_alphadiv_Mhi12", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK", "LS"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Alphadiv_Mhi12", parallelize = FALSE, nCores = 1, inDF = NULL)
#none were significant



##### Rarefying with Gunifrac for composition, beta diversity and identification of diel rhythms in ASVs #####
# Disclaimer: Rarefaction may lead to differences in abundances, causing different results. Use gunifrac_micro_coral_4.f.sub_otu.csv to get the same result.

library(GUniFrac)
#Transposing ASV table to use GUniFrac
mc4.f.sub_otu_transp <- transpose(as.data.frame(otu_table(micro_coral_4.f.sub)))
colnames(mc4.f.sub_otu_transp) <- rownames(as.data.frame(otu_table(micro_coral_4.f.sub)))
rownames(mc4.f.sub_otu_transp) <- colnames(as.data.frame(otu_table(micro_coral_4.f.sub)))

#Rarefying
guni_micro_coral_4.f.sub <- Rarefy(mc4.f.sub_otu_transp, depth = 20371)

#Transposing ASV table back
guni_micro_coral_4.f.sub_otu <- transpose(as.data.frame(guni_micro_coral_4.f.sub$otu.tab.rff)) 
colnames(guni_micro_coral_4.f.sub_otu) <- rownames(guni_micro_coral_4.f.sub$otu.tab.rff)
rownames(guni_micro_coral_4.f.sub_otu) <- colnames(guni_micro_coral_4.f.sub$otu.tab.rff)

#Exporting rarefied dataset 
write.csv(guni_micro_coral_4.f.sub_otu, "./gunifrac_micro_coral_4.f.sub_otu.csv", row.names=TRUE)

# Creating phyloseq object with rarefied dataset
guni_micro_coral_4.f.sub_otu <- as.matrix(guni_micro_coral_4.f.sub_otu)
OTU2 = otu_table(guni_micro_coral_4.f.sub_otu, taxa_are_rows = TRUE)

guni_4fsub <- phyloseq(OTU2, TAX, samples)
sample_variables(guni_4fsub)
print(guni_4fsub)

# Removing ASVs with zero abundance in all samples
sum(taxa_sums(guni_4fsub) == 0, na.rm=TRUE) # counting how many they are
# 1263
guni_4fsub_1 <- prune_taxa(taxa_sums(guni_4fsub) > 0, guni_4fsub)
any(taxa_sums(micro_coral_1) == 0) #this line must return FALSE
sum(taxa_sums(guni_4fsub_1) == 0, na.rm=TRUE) # this line must return zero
ntaxa(guni_4fsub_1) # calculating number of ASVs with abundance
# 24570

##### Beta diversity #####

##DCA
ordu.dca = ordinate(guni_4fsub, method = "DCA")
plot_dca <- plot_ordination(guni_4fsub, ordu.dca, 
                            color = "Species", 
                            shape = "Light")

plot_dca <- plot_dca + 
  scale_fill_manual(values = c("#CBD588", "#5F7FC7", 
                               "#DA5724", "#508578")) + 
  ggtitle("DCA - absolute frequency") + 
  geom_point(size = 3) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),      
    axis.title = element_text(size = 14),                   
    axis.text = element_text(size = 12),                    
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12, face = "italic")   # Legend text in italics
  )

print(plot_dca)
ggsave("./beta_diversity_guni_4fsub.pdf", height = 6, width = 8)


##DCA separate for zooxanthelate an azooxanthelate species
guni_4fsub_azoo <-subset_samples(guni_4fsub, 
                                 guni_4fsub@sam_data$Symbiont == "Azooxanthellate")
ordu.dca_azoo = ordinate(guni_4fsub_azoo, method = "DCA")
plot_dca_azoo <- plot_ordination(guni_4fsub_azoo, ordu.dca_azoo, 
                                 color = "Species", 
                                 shape = "Light")
plot_dca_azoo <- plot_dca_azoo + 
  scale_fill_manual(values = c("#5F7FC7", "#CBD588")) + 
  ggtitle("DCA - absolute frequency") + 
  geom_point(size = 3) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),      # Title
    axis.title = element_text(size = 14),                   # Axis labels
    axis.text = element_text(size = 12),                    # Axis text
    legend.title = element_text(size = 14, face = "italic"), # Legend title in italics
    legend.text = element_text(size = 12, face = "italic")   # Legend text in italics
  )

print(plot_dca_azoo)
ggsave("./dca_azoo.pdf", height = 6, width = 8)



guni_4fsub_zoo <-subset_samples(guni_4fsub, 
                                 guni_4fsub@sam_data$Symbiont == "Zooxanthellate")
ordu.dca_zoo = ordinate(guni_4fsub_zoo, method = "DCA")
plot_dca_zoo <- plot_ordination(guni_4fsub_zoo, ordu.dca_zoo, 
                                 color = "Species", 
                                 shape = "Light")
plot_dca_zoo <- plot_dca_zoo + 
  scale_fill_manual(values = c("#CBD588", "#DA5724")) + 
  ggtitle("DCA - absolute frequency") + 
  geom_point(size = 3) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),      
    axis.title = element_text(size = 14),                 
    axis.text = element_text(size = 12),                    
    legend.title = element_text(size = 14, face = "italic"), 
    legend.text = element_text(size = 12, face = "italic")   
  )
print(plot_dca_zoo)
ggsave("./dca_zoo.pdf", height = 6, width = 8)




## Calculating beta diversity

library(vegan)

## Comparing species
ps_tr <- microbiome::transform(guni_4fsub, "clr") #necessary for adonis.pairwise to work
# Extract the distance matrix from the phyloseq object
dist_matrix <- phyloseq::distance(guni_4fsub, method = "bray")


## PERMDISP
#https://rpubs.com/maddieSC/R_SOP_UCR_Jan_2018
#https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/
library(vegan)
library(phyloseq)

## Grouping by species

dist_matrix.disper.sp <- betadisper(dist_matrix, sample_data(guni_4fsub)$Species)
dist_matrix.disper.sp$distances #contains the distances of each sample to the centroid of its group
pdf("dist_matrix.disper.sp_dist_to_centroid.pdf")
boxplot(dist_matrix.disper.sp) #makes box plot
dev.off()
permutest(dist_matrix.disper.sp) #checking whether the distances differ statistically
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#           Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     3 1.0123 0.33744 14.227    999  0.001 ***
#Residuals 65 1.5417 0.02372                         

##ANOSIM
#Because permdisp was significant for Species as a factor, running anosim
anosim(dist_matrix, sample_data(guni_4fsub)$Species, permutations = 1000)
#ANOSIM statistic R: 0.7184 
#Significance: 0.000999 
#Permutation: free
#Number of permutations: 1000
#Species is significant

##Post-hoc test - Pairwise adonis2
library("pairwiseAdonis")
pairwise.adonis(dist_matrix, phyloseq::sample_data(ps_tr)$Species)
#  #                      pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
#  1    M_decactis vs M_hispida  1 1.9855837  5.279231 0.1455166   0.001      0.006   *
#  2   M_decactis vs T_coccinea  1 3.0717363  9.808410 0.2346038   0.001      0.006   *
#  3 M_decactis vs T_tagusensis  1 3.8724424 14.529417 0.3122630   0.001      0.006   *
#  4    M_hispida vs T_coccinea  1 4.0006872 14.726742 0.3085637   0.001      0.006   *
#  5  M_hispida vs T_tagusensis  1 4.7849028 21.132341 0.3903829   0.001      0.006   *
#  6 T_coccinea vs T_tagusensis  1 0.6790329  3.957257 0.1042556   0.005      0.030   .


       
## Comparing Zooxanthellate and Azooxanthellate
## PERMDISP
dist_matrix.disper.symb <- betadisper(dist_matrix, sample_data(guni_4fsub)$Symbiont)
dist_matrix.disper.symb$distances
plot(dist_matrix.disper.symb) #makes plot (PCoA)
permutest(dist_matrix.disper.symb) #checking whether the distances differ statistically
#          Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups     1 1.0749 1.07485 57.739    999  0.001 ***
#Residuals 67 1.2472 0.01862 

## ANOSIM
anosim(dist_matrix, sample_data(guni_4fsub)$Symbiont, permutations = 1000)
#ANOSIM statistic R: 0.7878 
#Significance: 0.000999 
#Permutation: free
#Number of permutations: 1000


##Testing each species separately

## Ma decactis by Time and day vs night
guni_4fsub_Mde <-subset_samples(guni_4fsub, #subsetting phyloseq object to have one species 
                                    guni_4fsub@sam_data$Species == "M_decactis")
# Extract the distance matrix from the phyloseq object
dist_matrix_Mde <- phyloseq::distance(guni_4fsub_Mde, method = "bray")

## By Time
## PERMDISP
dist_matrix.disper.MdeTime <- betadisper(dist_matrix_Mde, sample_data(guni_4fsub_Mde)$Time)
plot(dist_matrix.disper.MdeTime) #makes plot (PCoA)
permutest(dist_matrix.disper.MdeTime) #checking whether the distances differ statistically
#Permutation test for homogeneity of multivariate dispersions
#Response: Distances
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.013582 0.0027164 1.0175    999  0.443
#Residuals 10 0.026698 0.0026698 

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Mde <- data.frame(Time = sample_data(guni_4fsub_Mde)$Time)
# Run the PERMANOVA test
permanova_result_Mde1 <- adonis2(dist_matrix_Mde ~ Time, data = group_data_Mde, permutations = 999)
print(permanova_result_Mde1)
#         Df SumOfSqs      R2      F Pr(>F)
#Time      1   0.4042 0.06358 0.9505  0.623
#Residual 14   5.9539 0.93642              
#Total    15   6.3581 1.00000  


## By Light (day vs night)
## PERMDISP
dist_matrix.disper.MdeLight <- betadisper(dist_matrix_Mde, sample_data(guni_4fsub_Mde)$Light)
plot(dist_matrix.disper.MdeLight) #makes plot (PCoA)
boxplot(dist_matrix.disper.MdeLight) #makes box plot
permutest(dist_matrix.disper.MdeLight) #checking whether the distances differ statistically
#          Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.0004048 0.00040478 0.2102    999  0.654
#Residuals 14 0.0269594 0.00192567

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Mde2 <- data.frame(Light = sample_data(guni_4fsub_Mde)$Light)
# Run the PERMANOVA test
permanova_result_Mde2 <- adonis2(dist_matrix_Mde ~ Light, data = group_data_Mde2, permutations = 999)
print(permanova_result_Mde2)
adonis2(formula = dist_matrix_Mde ~ Light, data = group_data_Mde2, permutations = 999)
#         Df SumOfSqs      R2      F Pr(>F)
#Light     1   0.3839 0.06038 0.8996  0.839
#Residual 14   5.9742 0.93962              
#Total    15   6.3581 1.00000  



## Mu hispida by time and day vs night
guni_4fsub_Mhi <-subset_samples(guni_4fsub, #subsetting phyloseq object to have one species
                                guni_4fsub@sam_data$Species == "M_hispida")
# Extract the distance matrix from the phyloseq object
dist_matrix_Mhi <- phyloseq::distance(guni_4fsub_Mhi, method = "bray")
## By Time
## PERMDISP
dist_matrix.disper.MhiTime <- betadisper(dist_matrix_Mhi, sample_data(guni_4fsub_Mhi)$Time)
plot(dist_matrix.disper.MhiTime) #makes plot (PCoA)
permutest(dist_matrix.disper.MhiTime) #checking whether the distances differ statistically
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.00910 0.001820 0.0384    999      1
#Residuals 11 0.52084 0.047349    

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Mhi <- data.frame(Time = sample_data(guni_4fsub_Mhi)$Time)
# Run the PERMANOVA test
permanova_result_Mhi1 <- adonis2(dist_matrix_Mhi ~ Time, data = group_data_Mhi, permutations = 999)
print(permanova_result_Mhi1)
#         Df SumOfSqs      R2      F Pr(>F)
#Time      1   0.2478 0.04675 0.7357  0.776
#Residual 15   5.0535 0.95325              
#Total    16   5.3014 1.00000    


## By Light (day vs night)
## PERMDISP
dist_matrix.disper.MhiLight <- betadisper(dist_matrix_Mhi, sample_data(guni_4fsub_Mhi)$Light)
boxplot(dist_matrix.disper.MhiLight) #makes plot (PCoA)
permutest(dist_matrix.disper.MhiLight) #checking whether the distances differ statistically
#          Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
#Groups     1 0.01069 0.010686 0.345    999  0.544
#Residuals 15 0.46461 0.030974 

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Mhi2 <- data.frame(Light = sample_data(guni_4fsub_Mhi)$Light)
# Run the PERMANOVA test
permanova_result_Mhi2 <- adonis2(dist_matrix_Mhi ~ Light, data = group_data_Mhi2, permutations = 999)
print(permanova_result_Mhi2)
#         Df SumOfSqs      R2      F Pr(>F)
#Light     1   0.2488 0.04693 0.7387  0.754
#Residual 15   5.0526 0.95307              
#Total    16   5.3014 1.00000  


## T coccinea by time and day vs night
guni_4fsub_Tco <-subset_samples(guni_4fsub, #subsetting phyloseq object to have one species
                                guni_4fsub@sam_data$Species == "T_coccinea")

# Extract the distance matrix from the phyloseq object
dist_matrix_Tco <- phyloseq::distance(guni_4fsub_Tco, method = "bray")

## By Time
## PERMDISP
dist_matrix.disper.TcoTime <- betadisper(dist_matrix_Tco, sample_data(guni_4fsub_Tco)$Time)
plot(dist_matrix.disper.TcoTime) #makes plot (PCoA)
permutest(dist_matrix.disper.TcoTime) #checking whether the distances differ statistically
#          Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
#Groups     5 0.15571 0.031141 1.001    999  0.468
#Residuals 12 0.37332 0.031110  

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Tco <- data.frame(Time = sample_data(guni_4fsub_Tco)$Time)
# Run the PERMANOVA test
permanova_result_Tco1 <- adonis2(dist_matrix_Tco ~ Time, data = group_data_Tco, permutations = 999)
print(permanova_result_Tco1)
#         Df SumOfSqs      R2      F Pr(>F)
#Time      1   0.2304 0.06289 1.0738  0.358
#Residual 16   3.4331 0.93711              
#Total    17   3.6635 1.00000    

## By Light (day vs night)
## PERMDISP
dist_matrix.disper.TcoLight <- betadisper(dist_matrix_Tco, sample_data(guni_4fsub_Tco)$Light)
boxplot(dist_matrix.disper.TcoLight) #makes plot (PCoA)
permutest(dist_matrix.disper.TcoLight) #checking whether the distances differ statistically
#          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)  
#Groups     1 0.10915 0.10915 3.8231    999  0.058 .
#Residuals 16 0.45680 0.02855

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Tco2 <- data.frame(Light = sample_data(guni_4fsub_Tco)$Light)
# Run the PERMANOVA test
permanova_result_Tco2 <- adonis2(dist_matrix_Tco ~ Light, data = group_data_Tco2, permutations = 999)
print(permanova_result_Tco2)
#        Df SumOfSqs      R2     F Pr(>F)
#Light     1   0.3050 0.08325 1.453  0.106
#Residual 16   3.3585 0.91675             
#Total    17   3.6635 1.00000   


## T tagusensis by time and day vs night
guni_4fsub_Tta <-subset_samples(guni_4fsub, #subsetting phyloseq object to have one species
                                guni_4fsub@sam_data$Species == "T_tagusensis")

# Extract the distance matrix from the phyloseq object
dist_matrix_Tta <- phyloseq::distance(guni_4fsub_Tta, method = "bray")

## By time
## PERMDISP
dist_matrix.disper.TtaTime <- betadisper(dist_matrix_Tta, sample_data(guni_4fsub_Tta)$Time)
#Warning message: In betadisper(dist_matrix_Tta, sample_data(guni_4fsub_Tta)$Time) : some squared distances are negative and changed to zero
plot(dist_matrix.disper.TtaTime) #makes plot (PCoA)
permutest(dist_matrix.disper.TtaTime) #checking whether the distances differ statistically
#          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     5 0.02676 0.0053521 0.2403    999  0.933
#Residuals 12 0.26729 0.0222745  

##PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Tta <- data.frame(Time = sample_data(guni_4fsub_Tta)$Time)
# Run the PERMANOVA test
permanova_result_Tta1 <- adonis2(dist_matrix_Tta ~ Time, data = group_data_Tta, permutations = 999)
print(permanova_result_Tta1)
#         Df SumOfSqs      R2     F Pr(>F)
#Time      1  0.11219 0.05169 0.872   0.47
#Residual 16  2.05848 0.94831             
#Total    17  2.17067 1.00000   


## By Light (day vs night)
## PERMDISP
dist_matrix.disper.TtaLight <- betadisper(dist_matrix_Tta, sample_data(guni_4fsub_Tta)$Light)
boxplot(dist_matrix.disper.TtaLight) #faz plot (uma PCoA)
permutest(dist_matrix.disper.TtaLight) #checking whether the distances differ statistically
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.02696 0.026962 1.0638    999  0.324
#Residuals 16 0.40550 0.025344     

## PERMANOVA
# Create a data frame with the relevant grouping variable (e.g., Species)
group_data_Tta2 <- data.frame(Light = sample_data(guni_4fsub_Tta)$Light)
# Run the PERMANOVA test
permanova_result_Tta2 <- adonis2(dist_matrix_Tta ~ Light, data = group_data_Tta2, permutations = 999)
print(permanova_result_Tta2)
#         Df SumOfSqs      R2      F Pr(>F)
#Light     1  0.10196 0.04697 0.7886  0.575
#Residual 16  2.06871 0.95303              
#Total    17  2.17067 1.00000  


##### Rarefying each species separately #####
library(GUniFrac)
mc4.f.sub_otu_transp_Mde <- mc4.f.sub_otu_transp[1:16,] 
mc4.f.sub_otu_transp_Mhi <- mc4.f.sub_otu_transp[17:33,]
mc4.f.sub_otu_transp_Tco <- mc4.f.sub_otu_transp[34:51,]
mc4.f.sub_otu_transp_Tta <- mc4.f.sub_otu_transp[52:69,]

#Mde 20371, Mhi 39403, Tco 65675, Tta 30977

guni_micro_coral_4.f.sub_Mde <- Rarefy(mc4.f.sub_otu_transp_Mde, depth = 20371)
guni_micro_coral_4.f.sub_otu_Mde <- transpose(as.data.frame(guni_micro_coral_4.f.sub_Mde$otu.tab.rff)) 
colnames(guni_micro_coral_4.f.sub_otu_Mde) <- rownames(guni_micro_coral_4.f.sub_Mde$otu.tab.rff)
rownames(guni_micro_coral_4.f.sub_otu_Mde) <- colnames(guni_micro_coral_4.f.sub_Mde$otu.tab.rff)
write.csv(guni_micro_coral_4.f.sub_otu_Mde,"Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mde.csv", row.names=TRUE)

guni_micro_coral_4.f.sub_Mhi <- Rarefy(mc4.f.sub_otu_transp_Mhi, depth = 39403)
guni_micro_coral_4.f.sub_otu_Mhi <- transpose(as.data.frame(guni_micro_coral_4.f.sub_Mhi$otu.tab.rff)) 
colnames(guni_micro_coral_4.f.sub_otu_Mhi) <- rownames(guni_micro_coral_4.f.sub_Mhi$otu.tab.rff)
rownames(guni_micro_coral_4.f.sub_otu_Mhi) <- colnames(guni_micro_coral_4.f.sub_Mhi$otu.tab.rff)
write.csv(guni_micro_coral_4.f.sub_otu_Mhi,"Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mhi.csv", row.names=TRUE)

guni_micro_coral_4.f.sub_Tco <- Rarefy(mc4.f.sub_otu_transp_Tco, depth = 65675)
guni_micro_coral_4.f.sub_otu_Tco <- transpose(as.data.frame(guni_micro_coral_4.f.sub_Tco$otu.tab.rff)) 
colnames(guni_micro_coral_4.f.sub_otu_Tco) <- rownames(guni_micro_coral_4.f.sub_Tco$otu.tab.rff)
rownames(guni_micro_coral_4.f.sub_otu_Tco) <- colnames(guni_micro_coral_4.f.sub_Tco$otu.tab.rff)
write.csv(guni_micro_coral_4.f.sub_otu_Tco,"Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tco.csv", row.names=TRUE)

guni_micro_coral_4.f.sub_Tta <- Rarefy(mc4.f.sub_otu_transp_Tta, depth = 30977)
guni_micro_coral_4.f.sub_otu_Tta <- transpose(as.data.frame(guni_micro_coral_4.f.sub_Tta$otu.tab.rff)) 
colnames(guni_micro_coral_4.f.sub_otu_Tta) <- rownames(guni_micro_coral_4.f.sub_Tta$otu.tab.rff)
rownames(guni_micro_coral_4.f.sub_otu_Tta) <- colnames(guni_micro_coral_4.f.sub_Tta$otu.tab.rff)
write.csv(guni_micro_coral_4.f.sub_otu_Tta,"Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tta.csv", row.names=TRUE)

##### Differential abundance analysis with DeSeq2 using rarefied data for each species #####


library(DESeq2)
# Using rarefied datasets (gunifrac) as count data

## M_decactis
# Converting the 2 data frames to deseq2 object
ds_guni_por_spp_Mde <- DESeqDataSetFromMatrix(countData = guni_micro_coral_4.f.sub_otu_Mde,
                              colData = rbind(samples_df[1:8,],samples_df[10:11,],samples_df[13:18,]),
                              design = ~ Light) #using rbind to exclude from the metadata table samples that are not in the ASV table (because DeSeq2 does not handle NAs)

# Running DESeq2
deseq_guni_por_spp_Mde <- DESeq(ds_guni_por_spp_Mde)
results_guni_spp_Mde_Light <- results(deseq_guni_por_spp_Mde, contrast = c("Light", "Day", "Night"))
# Assuming you want to extract genes with an adjusted p-value < 0.05
DEGs_guni_spp_Mde_Light <- subset(results_guni_spp_Mde_Light, padj < 0.05)
write.csv(DEGs_guni_spp_Mde_Light, "./guni_4fsub_DESVs_guni_spp_Mde_Light.csv", row.names=TRUE)

# Joining DESeq2 results with information of each species (whether they were collected during Day or Night) 
ASV_abunds_0 <- guni_micro_coral_4.f.sub_otu_Mde[rownames(DEGs_guni_spp_Mde_Light),]
Light_data <- transpose(sample_data(guni_4fsub)[1:16,8])
colnames(Light_data) <- rownames(sample_data(guni_4fsub)[1:16,8])
ASV_abunds_to_plot <- rbind(Light_data, ASV_abunds_0)
write.csv(ASV_abunds_to_plot, "ASV_abunds_to_plot_guni_spp_Mde.csv")
# In excel, day was separated from night and only ASVs with more than one non consecutive time point with abundance >0 was considered significant



## M_hispida
# Converting the 2 data frames to deseq2 object
ds_guni_por_spp_Mhi <- DESeqDataSetFromMatrix(countData = guni_micro_coral_4.f.sub_otu_Mhi,
                                              colData = samples_df[19:35,],
                                              design = ~ Light) 

# Rodando DESeq2
deseq_guni_por_spp_Mhi <- DESeq(ds_guni_por_spp_Mhi)
results_guni_spp_Mhi_Light <- results(deseq_guni_por_spp_Mhi, contrast = c("Light", "Day", "Night"))
# Assuming you want to extract genes with an adjusted p-value < 0.05
DEGs_guni_spp_Mhi_Light <- subset(results_guni_spp_Mhi_Light, padj < 0.05)
write.csv(DEGs_guni_spp_Mhi_Light, "./guni_4fsub_DESVs_guni_spp_Mhi_Light.csv", row.names=TRUE)

# Joining DESeq2 results with information of each species (whether they were collected during Day or Night) 
ASV_abunds_0_Mh <- guni_micro_coral_4.f.sub_otu_Mhi[rownames(DEGs_guni_spp_Mhi_Light),]
Light_data_Mh <- transpose(sample_data(guni_4fsub)[17:33,8]) 
colnames(Light_data_Mh) <- rownames(sample_data(guni_4fsub)[17:33,8])
ASV_abunds_to_plot_Mh <- rbind(Light_data_Mh, ASV_abunds_0_Mh)
write.csv(ASV_abunds_to_plot_Mh, "ASV_abunds_to_plot_guni_spp_Mhi.csv")
# In excel, day was separated from night and only ASVs with more than one non consecutive time point with abundance >0 was considered significant



## T_coccinea
# Converting the 2 data frames to deseq2 object
ds_guni_por_spp_Tco <- DESeqDataSetFromMatrix(countData = guni_micro_coral_4.f.sub_otu_Tco,
                                              colData = samples_df[36:53,],
                                              design = ~ Light)
# Rodando DESeq2
deseq_guni_por_spp_Tco <- DESeq(ds_guni_por_spp_Tco)
results_guni_spp_Tco_Light <- results(deseq_guni_por_spp_Tco, contrast = c("Light", "Day", "Night"))
# Assuming you want to extract genes with an adjusted p-value < 0.05
DEGs_guni_spp_Tco_Light <- subset(results_guni_spp_Tco_Light, padj < 0.05)
write.csv(DEGs_guni_spp_Tco_Light, "./guni_4fsub_DESVs_guni_spp_Tco_Light.csv", row.names=TRUE)

#  Joining DESeq2 results with information of each species (whether they were collected during Day or Night) 
ASV_abunds_0_Tc <- guni_micro_coral_4.f.sub_otu_Tco[rownames(DEGs_guni_spp_Tco_Light),]
Light_data_Tc <- transpose(sample_data(guni_4fsub)[34:51,8])
colnames(Light_data_Tc) <- rownames(sample_data(guni_4fsub)[34:51,8])
ASV_abunds_to_plot_Tc <- rbind(Light_data_Tc, ASV_abunds_0_Tc)
write.csv(ASV_abunds_to_plot_Tc, "ASV_abunds_to_plot_guni_spp_Tco.csv")
# In excel, day was separated from night and only ASVs with more than one non consecutive time point with abundance >0 was considered significant


## T_tagusensis
# Converting the 2 data frames to deseq2 object
ds_guni_por_spp_Tta <- DESeqDataSetFromMatrix(countData = guni_micro_coral_4.f.sub_otu_Tta,
                                              colData = samples_df[54:71,],
                                              design = ~ Light)
# Rodando DESeq2
deseq_guni_por_spp_Tta <- DESeq(ds_guni_por_spp_Tta)
results_guni_spp_Tta_Light <- results(deseq_guni_por_spp_Tta, contrast = c("Light", "Day", "Night"))
# Assuming you want to extract genes with an adjusted p-value < 0.05
DEGs_guni_spp_Tta_Light <- subset(results_guni_spp_Tta_Light, padj < 0.05)
write.csv(DEGs_guni_spp_Tta_Light, "./guni_4fsub_DESVs_guni_spp_Tta_Light.csv", row.names=TRUE)

# Joining DESeq2 results with information of each species (whether they were collected during Day or Night) 
ASV_abunds_0_Tt <- guni_micro_coral_4.f.sub_otu_Tta[rownames(DEGs_guni_spp_Tta_Light),]
Light_data_Tt <- transpose(sample_data(guni_4fsub)[52:69,8])
colnames(Light_data_Tt) <- rownames(sample_data(guni_4fsub)[52:69,8])
ASV_abunds_to_plot_Tt <- rbind(Light_data_Tt, ASV_abunds_0_Tt)
write.csv(ASV_abunds_to_plot_Tt, "ASV_abunds_to_plot_guni_spp_Tta.csv")
# In excel, day was separated from night and only ASVs with more than one non consecutive time point with abundance >0 was considered significant


##### Running JTKcycle with the datasets from each species rarefied separately #####

library(MetaCycle)

## JTK with frequency data
# Data were transformed to proportion in excel (frequency data mentioned below), by dividing each value by the depth of each species (the depth used to rarefy each species)
# Mde
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mde_rel.csv", outdir = "MetaCycle_guni_Mde24_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mde24", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mde_rel.csv", outdir = "MetaCycle_guni_Mde12_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mde12", parallelize = FALSE, nCores = 1, inDF = NULL)

# Mhi
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mhi_rel.csv", outdir = "MetaCycle_guni_Mhi24_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mhi24", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mhi_rel.csv", outdir = "MetaCycle_guni_Mhi12_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mhi12", parallelize = FALSE, nCores = 1, inDF = NULL)

# Tco
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tco_rel.csv", outdir = "MetaCycle_guni_Tco24_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tco24", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tco_rel.csv", outdir = "MetaCycle_guni_Tco12_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tco12", parallelize = FALSE, nCores = 1, inDF = NULL)

# Tta
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tta_rel.csv", outdir = "MetaCycle_guni_Tta24_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tta24", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tta_rel.csv", outdir = "MetaCycle_guni_Tta12_guni_indiv", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tta12", parallelize = FALSE, nCores = 1, inDF = NULL)

## JTK with count data
# Mde
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mde_wNA.csv", outdir = "MetaCycle_guni_Mde24_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mde24_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mde_wNA.csv", outdir = "MetaCycle_guni_Mde12_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mde12_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)

# Mhi
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mhi_wNA.csv", outdir = "MetaCycle_guni_Mhi24_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mhi24_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Mhi_wNA.csv", outdir = "MetaCycle_guni_Mhi12_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Mhi12_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)

# Tco
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tco.csv", outdir = "MetaCycle_guni_Tco24_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tco24_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tco.csv", outdir = "MetaCycle_guni_Tco12_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tco12_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)

# Tta
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tta.csv", outdir = "MetaCycle_guni_Tta24_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 24, maxper = 24, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tta24_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)
meta2d("Period_search_methodologies/guni_micro_coral_4.f.sub_otu_Tta.csv", outdir = "MetaCycle_guni_Tta12_guni_indiv_vals_abs", filestyle="csv", timepoints=seq(1, 72, by=4), minper = 12, maxper = 12, cycMethod = c("JTK"), analysisStrategy = "auto", outputFile = TRUE, outIntegration = "both", adjustPhase = "predictedPer", combinePvalue = "bonferroni", weightedPerPha = FALSE, ARSmle = "auto", ARSdefaultPer = 24, outRawData = FALSE, releaseNote = TRUE, outSymbol = "Guni_indiv_Tta12_vals_abs", parallelize = FALSE, nCores = 1, inDF = NULL)



###Ploting daily variations on ASVs 
# The following input data was generated out of R. The file has abundance data of all ASVs that were considered significant from DESeq and JTK.
#The metadata file has four columns: sample (species+ASV with significant results), analise (type of analysis - JTK or DESEq), species, and ritmo (if variation was along 12h, 24h or day vs night). 

library(ggplot2)
library(dplyr)
library(tidyr)

# Read the abundance data
abundance_data <- read.csv("jtk-deseq.csv", row.names = 1)
sample_metadata <- read.csv("jtk-deseq-metadata.csv")
abundance_long <- abundance_data %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = -sample, names_to = "sample_id", values_to = "abundance")
abundance_long <- abundance_long %>%
  left_join(sample_metadata, by = c("sample" = "sample"))

abundance_long$ritmo <- factor(abundance_long$ritmo)
abundance_long$sample_id <- factor(abundance_long$sample_id, levels = unique(abundance_long$sample_id))
abundance_long$ritmo <- factor(abundance_long$ritmo)
abundance_long <- abundance_long %>%
  mutate(ritmo_fill = ifelse(abundance > 0, as.character(ritmo), NA))

# Plot for absolute abundance
absolute_plot <- ggplot(abundance_long, aes(x = sample_id, y = sample, size = abundance)) +
  geom_point(aes(color = species, fill = ifelse(abundance == 0, NA, species)), shape = 21, stroke = 1) +
  scale_size_continuous(range = c(1, 12), breaks = c(0, 5, 10, 20, 50, 100, 200, 300, 800, 1000, 2000, 2500), name = "Abundance") +
  scale_fill_manual(values = c("#619CFF", "#F8766D", "#00BA38", "#E69F00"), na.value = "white") +
  scale_color_manual(values = c("#619CFF", "#F8766D", "#00BA38", "#E69F00")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(x = "", y = "Sample") +
  ggtitle("Absolute Abundance") +
  scale_y_discrete(limits = rev(unique(abundance_long$sample))) +
  guides(fill = "none") 

absolute_plot




##### Composition using rarefied data #####

# Using the dataset that was rarefied including all species
library(phyloseq)
library(microbiome)
library(ggplot2)
library(pals)
library(Polychrome)

## Barplot composition

taxic <- as.data.frame(guni_4fsub@tax_table)  # this will help in setting large color options

otu.df <- as.data.frame(otu_table(guni_4fsub))  # make a dataframe for OTU information.
head(otu.df) # check the rows and columns

taxic$ASV <- row.names.data.frame(otu.df)  # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic)  # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic)  # convert it into a matrix.

# Identify the column indices for the columns you want to modify
columns_to_modify <- c("Phylum", "Class", "Order", "Family", "Genus")

# Replace NA values and "NA" string for specific columns in the taxmat matrix
for (column in columns_to_modify) {
  taxmat[is.na(taxmat[, column]) | taxmat[, column] == "NA", column] <- paste("Unclassified", column, sep = " ")
}
# Print the modified taxmat matrix
taxmat
new.tax <- tax_table(taxmat)  # convert into phyloseq compatible file
tax_table(guni_4fsub) <- new.tax  # incorporate into phyloseq Object


#merge at family level
guni_4fsub_fam <- aggregate_taxa(guni_4fsub, "Family")
fam_abundance <- taxa_sums(guni_4fsub_fam)
sorted_families <- sort(fam_abundance, decreasing = TRUE)
top_20_families <- names(sorted_families)[1:20]
guni_4fsub_fam <- prune_taxa(top_20_families, guni_4fsub_fam)


#merge at genus level
guni_4fsub_gen <- aggregate_taxa(guni_4fsub, "Genus")
gen_abundance <- taxa_sums(guni_4fsub_gen)
sorted_genus <- sort(gen_abundance, decreasing = TRUE)
top_20_gneus <- names(sorted_genus)[1:20]
guni_4fsub_gen <- prune_taxa(top_20_gneus, guni_4fsub_gen)

#merge at order level
guni_4fsub_ord <- aggregate_taxa(guni_4fsub, "Order")
ord_abundance <- taxa_sums(guni_4fsub_ord)
sorted_order <- sort(ord_abundance, decreasing = TRUE)
top_20_orders <- names(sorted_order)[1:20]
guni_4fsub_ord <- prune_taxa(top_20_orders, guni_4fsub_ord)

#merge at phylum level
guni_4fsub_phyl <- aggregate_taxa(guni_4fsub, "Phylum")
phyl_abundance <- taxa_sums(guni_4fsub_phyl)
sorted_phyl <- sort(phyl_abundance, decreasing = TRUE)
top_20_phyl <- names(sorted_phyl)[1:20]
guni_4fsub_phyl <- prune_taxa(top_20_phyl, guni_4fsub_phyl)

## Abundance per species of 20 most abundant families

guni_4fsub_fam_fraction <- merge_samples(guni_4fsub_fam, "Species")

plot_bar(guni_4fsub_fam_fraction, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_geral_guni_4fsub.pdf", height = 6, width = 8)

## Abundance per species of 20 most abundant genera

guni_4fsub_gen_fraction <- merge_samples(guni_4fsub_gen, "Species")

plot_bar(guni_4fsub_gen_fraction, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_geral_guni_4fsub.pdf", height = 6, width = 8)

## Abundance per species of 20 most abundant orders

guni_4fsub_ord_fraction <- merge_samples(guni_4fsub_ord, "Species")

plot_bar(guni_4fsub_ord_fraction, fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_geral_guni_4fsub.pdf", height = 6, width = 8)

## Abundance per species of 20 most abundant phyla

guni_4fsub_phyl_fraction <- merge_samples(guni_4fsub_phyl, "Species")

plot_bar(guni_4fsub_phyl_fraction, fill = "Phylum") + 
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Phylum_barplot_geral_guni_4fsub.pdf", height = 6, width = 8)


#Plot using frequencies rather than count data 

guni_4fsub.fam.rel <- microbiome::transform(guni_4fsub_fam_fraction, "compositional")
plot_bar(guni_4fsub.fam.rel, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_relativol_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.gen.rel <- microbiome::transform(guni_4fsub_gen_fraction, "compositional")
plot_bar(guni_4fsub.gen.rel, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_relativol_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.ord.rel <- microbiome::transform(guni_4fsub_ord_fraction, "compositional")
plot_bar(guni_4fsub.ord.rel, fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_relativol_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.phyl.rel <- microbiome::transform(guni_4fsub_phyl_fraction, "compositional")
plot_bar(guni_4fsub.phyl.rel, fill = "Phylum") + 
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Phylum_barplot_relativol_guni_4fsub.pdf", height = 6, width = 8)


## Abundance per species by Time - 20 most abundant families

plot_bar(guni_4fsub_fam, x="Time", fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_species1_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub_fam, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_species2_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.fam.rel_sep <- microbiome::transform(guni_4fsub_fam, "compositional")

plot_bar(guni_4fsub.fam.rel_sep, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_species_relativa_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub.fam.rel_sep, x="Time", fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_species1_relativa_guni_4fsub.pdf", height = 6, width = 8)

## Abundance per species by Time - 20 most abundant genera

plot_bar(guni_4fsub_gen, x="Time", fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_species1_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub_gen, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_species2_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.gen.rel_sep <- microbiome::transform(guni_4fsub_gen, "compositional")

plot_bar(guni_4fsub.gen.rel_sep, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_species_relativa_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub.gen.rel_sep, x="Time", fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_species1_relativa_guni_4fsub.pdf", height = 6, width = 8)


## Abundance per species by Time - 20 most abundant orders

plot_bar(guni_4fsub_ord, x="Time", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_species1_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub_ord, fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_species2_guni_4fsub.pdf", height = 6, width = 8)

guni_4fsub.ord.rel_sep <- microbiome::transform(guni_4fsub_ord, "compositional")

plot_bar(guni_4fsub.ord.rel_sep, fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_species_relativa_guni_4fsub.pdf", height = 6, width = 8)

plot_bar(guni_4fsub.ord.rel_sep, x="Time", fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  facet_wrap(~ Species, scales = "free_x") +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_species1_relativa_guni_4fsub.pdf", height = 6, width = 8)


## Abundance per species by Light (Day vs Night) - 20 most abundant ...

# ... Families
plot_bar(guni_4fsub.fam.rel_sep, fill = "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack") +
  facet_wrap(~ Species+Light, scales = "free_x", nrow = 4) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=rev(as.vector(polychrome(20))))
ggsave("./Family_barplot_luz_f_sub.pdf", height = 12, width = 8)

# ... Genera
plot_bar(guni_4fsub.gen.rel_sep, fill = "Genus") + 
  geom_bar(aes(fill=Genus), stat="identity", position="stack") +
  facet_wrap(~ Species+Light, scales = "free_x", nrow = 4) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Genus_barplot_luz_guni_4fsub.pdf", height = 12, width = 8)

# ... Orders
plot_bar(guni_4fsub.ord.rel_sep, fill = "Order") + 
  geom_bar(aes(fill=Order), stat="identity", position="stack") +
  facet_wrap(~ Species+Light, scales = "free_x", nrow = 4) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Order_barplot_luz_guni_4fsub.pdf", height = 12, width = 8)

# ... Phyla
plot_bar(guni_4fsub.phyl.rel_sep, fill = "Phylum") + 
  geom_bar(aes(fill=Phylum), stat="identity", position="stack") +
  facet_wrap(~ Species+Light, scales = "free_x", nrow = 4) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values=as.vector(cols25(20)))
ggsave("./Phylum_barplot_luz_guni_4fsub.pdf", height = 12, width = 8)




##### Microbial core each species #####

### ASV level
# The code below calculates the core microbiome for the selected subset of samples and identifies taxa present in 100% (1.0) of samples.
core_threshold <- 1.0  # Adjust this threshold as needed

# all species
# Calculate the core microbiome
abund_mat_all <- as.matrix(otu_table(guni_4fsub))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_all <- ifelse(rowSums(abund_mat_all > 0) / ncol(abund_mat_all) >= core_threshold, 1, 0)
write.csv(core_matrix_all, "./core_matrix_all.csv", row.names=TRUE) #export to csv and sum column of presence/absence of ASVs


# M decactis
# Calculate the core microbiome
abund_mat_Mde <- as.matrix(otu_table(guni_4fsub_Mde))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mde <- ifelse(rowSums(abund_mat_Mde > 0) / ncol(abund_mat_Mde) >= core_threshold, 1, 0)
write.csv(core_matrix_Mde, "./core_matrix_Mde.csv", row.names=TRUE) #export to csv and sum column of presence/absence of ASVs


# M hispida
# Calculate the core microbiome
abund_mat_Mhi <- as.matrix(otu_table(guni_4fsub_Mhi))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mhi <- ifelse(rowSums(abund_mat_Mhi > 0) / ncol(abund_mat_Mhi) >= core_threshold, 1, 0)
write.csv(core_matrix_Mhi, "./core_matrix_Mhi.csv", row.names=TRUE) #export to csv and sum column of presence/absence of ASVs


# T tagusensis
# Calculate the core microbiome
abund_mat_Tta <- as.matrix(otu_table(guni_4fsub_Tta))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tta <- ifelse(rowSums(abund_mat_Tta > 0) / ncol(abund_mat_Tta) >= core_threshold, 1, 0)
write.csv(core_matrix_Tta, "./core_matrix_Tta.csv", row.names=TRUE) #export to csv and sum column of presence/absence of ASVs


# T coccinea
# Calculate the core microbiome
abund_mat_Tco <- as.matrix(otu_table(guni_4fsub_Tco))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tco <- ifelse(rowSums(abund_mat_Tco > 0) / ncol(abund_mat_Tco) >= core_threshold, 1, 0)
write.csv(core_matrix_Tco, "./core_matrix_Tco.csv", row.names=TRUE) #export to csv and sum column of presence/absence of ASVs


# Few ASVs in the microbial core (Mde 1 - ASV_12, Mhi 2 - ASV_4 ASV_6, Tta 3 - ASV_1 ASV_2 ASV_4, Tco 4 - ASV_1 ASV_2 ASV_4 ASV_7)

## Microbial core in other taxonomic levels
## Genus level

# Subsetting data for each species
guni_4fsub_gen2 <- aggregate_taxa(guni_4fsub, "Genus")
guni_4fsub_gen2_Mde <-subset_samples(guni_4fsub_gen2, #phyloseq object leaving only 1 species
                                guni_4fsub_gen2@sam_data$Species == "M_decactis")
guni_4fsub_gen2_Mhi <-subset_samples(guni_4fsub_gen2, #phyloseq object leaving only 1 species
                                guni_4fsub_gen2@sam_data$Species == "M_hispida")
guni_4fsub_gen2_Tta <-subset_samples(guni_4fsub_gen2, #phyloseq object leaving only 1 species
                                guni_4fsub_gen2@sam_data$Species == "T_tagusensis")
guni_4fsub_gen2_Tco <-subset_samples(guni_4fsub_gen2, #phyloseq object leaving only 1 species
                                guni_4fsub_gen2@sam_data$Species == "T_coccinea")

core_threshold <- 1.0  # Adjust this threshold as needed

# M decactis
# Calculate the core microbiome
abund_mat_gen2_Mde <- as.matrix(otu_table(guni_4fsub_gen2_Mde))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_gen2_Mde <- ifelse(rowSums(abund_mat_gen2_Mde > 0) / ncol(abund_mat_gen2_Mde) >= core_threshold, 1, 0)
write.csv(core_matrix_gen2_Mde, "./core_matrix_genus_Mde.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# M hispida
# Calculate the core microbiome
abund_mat_gen2_Mhi <- as.matrix(otu_table(guni_4fsub_gen2_Mhi))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_gen2_Mhi <- ifelse(rowSums(abund_mat_gen2_Mhi > 0) / ncol(abund_mat_gen2_Mhi) >= core_threshold, 1, 0)
write.csv(core_matrix_gen2_Mhi, "./core_matrix_genus_Mhi.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# T tagusensis
# Calculate the core microbiome
abund_mat_gen2_Tta <- as.matrix(otu_table(guni_4fsub_gen2_Tta))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_gen2_Tta <- ifelse(rowSums(abund_mat_gen2_Tta > 0) / ncol(abund_mat_gen2_Tta) >= core_threshold, 1, 0)
write.csv(core_matrix_gen2_Tta, "./core_matrix_genus_Tta.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# T coccinea
# Calculate the core microbiome
abund_mat_gen2_Tco <- as.matrix(otu_table(guni_4fsub_gen2_Tco))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_gen2_Tco <- ifelse(rowSums(abund_mat_gen2_Tco > 0) / ncol(abund_mat_gen2_Tco) >= core_threshold, 1, 0)
write.csv(core_matrix_gen2_Tco, "./core_matrix_genus_Tco.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# All
# Calculate the core microbiome
abund_mat_gen2 <- as.matrix(otu_table(guni_4fsub_gen2))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_gen2 <- ifelse(rowSums(abund_mat_gen2 > 0) / ncol(abund_mat_gen2) >= core_threshold, 1, 0)
write.csv(core_matrix_gen2, "./core_matrix_genus_all.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

## Family level

# subsetting data for each species
guni_4fsub_fam2 <- aggregate_taxa(guni_4fsub, "Family")
guni_4fsub_fam2_Mde <-subset_samples(guni_4fsub_fam2, #phyloseq object so deixando 1 spp por vez
                                     guni_4fsub_fam2@sam_data$Species == "M_decactis")
guni_4fsub_fam2_Mhi <-subset_samples(guni_4fsub_fam2, #phyloseq object so deixando 1 spp por vez
                                     guni_4fsub_fam2@sam_data$Species == "M_hispida")
guni_4fsub_fam2_Tta <-subset_samples(guni_4fsub_fam2, #phyloseq object so deixando 1 spp por vez
                                     guni_4fsub_fam2@sam_data$Species == "T_tagusensis")
guni_4fsub_fam2_Tco <-subset_samples(guni_4fsub_fam2, #phyloseq object so deixando 1 spp por vez
                                     guni_4fsub_fam2@sam_data$Species == "T_coccinea")

core_threshold <- 1.0  # Adjust this threshold as needed

# M decactis
# Calculate the core microbiome
abund_mat_fam2_Mde <- as.matrix(otu_table(guni_4fsub_fam2_Mde))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_fam2_Mde <- ifelse(rowSums(abund_mat_fam2_Mde > 0) / ncol(abund_mat_fam2_Mde) >= core_threshold, 1, 0)
write.csv(core_matrix_fam2_Mde, "./core_matrix_family_Mde.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# M hispida
# Calculate the core microbiome
abund_mat_fam2_Mhi <- as.matrix(otu_table(guni_4fsub_fam2_Mhi))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_fam2_Mhi <- ifelse(rowSums(abund_mat_fam2_Mhi > 0) / ncol(abund_mat_fam2_Mhi) >= core_threshold, 1, 0)
write.csv(core_matrix_fam2_Mhi, "./core_matrix_family_Mhi.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# T tagusensis
# Calculate the core microbiome
abund_mat_fam2_Tta <- as.matrix(otu_table(guni_4fsub_fam2_Tta))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_fam2_Tta <- ifelse(rowSums(abund_mat_fam2_Tta > 0) / ncol(abund_mat_fam2_Tta) >= core_threshold, 1, 0)
write.csv(core_matrix_fam2_Tta, "./core_matrix_family_Tta.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# T coccinea
# Calculate the core microbiome
abund_mat_fam2_Tco <- as.matrix(otu_table(guni_4fsub_fam2_Tco))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_fam2_Tco <- ifelse(rowSums(abund_mat_fam2_Tco > 0) / ncol(abund_mat_fam2_Tco) >= core_threshold, 1, 0)
write.csv(core_matrix_fam2_Tco, "./core_matrix_family_Tco.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# All
# Calculate the core microbiome
abund_mat_fam2 <- as.matrix(otu_table(guni_4fsub_fam2))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_fam2 <- ifelse(rowSums(abund_mat_fam2 > 0) / ncol(abund_mat_fam2) >= core_threshold, 1, 0)
write.csv(core_matrix_fam2, "./core_matrix_family_all.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

#.Genus level (1349 genera in total)
#Mde 18
#Mhi 12
#Tta 8
#Tco 7
#all 1 (Endozoicomonas)

#.Family level (630 families in total)
#Mde 28
#Mhi 20
#Tta 11
#Tco 10
#all 5 (Flavobacteriaceae, Cyanobiaceae, Rhodobacteraceae, Vibrionaceae, Endozoicomonadaceae)

## Each species during day and night separately

#ASV level dia de cada spp
guni_4fsub_Mde_day <-subset_samples(guni_4fsub_Mde, #phyloseq object including day or night of one species
                                    guni_4fsub_Mde@sam_data$Light == "Day")
guni_4fsub_Mde_night <-subset_samples(guni_4fsub_Mde, #phyloseq object including day or night of one species
                                    guni_4fsub_Mde@sam_data$Light == "Night")

guni_4fsub_Mhi_day <-subset_samples(guni_4fsub_Mhi, #phyloseq object including day or night of one species
                                    guni_4fsub_Mhi@sam_data$Light == "Day")
guni_4fsub_Mhi_night <-subset_samples(guni_4fsub_Mhi, #phyloseq object including day or night of one species
                                    guni_4fsub_Mhi@sam_data$Light == "Night")

guni_4fsub_Tco_day <-subset_samples(guni_4fsub_Tco, #phyloseq object including day or night of one species
                                    guni_4fsub_Tco@sam_data$Light == "Day")
guni_4fsub_Tco_night <-subset_samples(guni_4fsub_Tco, #phyloseq object including day or night of one species
                                      guni_4fsub_Tco@sam_data$Light == "Night")

guni_4fsub_Tta_day <-subset_samples(guni_4fsub_Tta, #phyloseq object including day or night of one species
                                    guni_4fsub_Tta@sam_data$Light == "Day")
guni_4fsub_Tta_night <-subset_samples(guni_4fsub_Tta, #phyloseq object including day or night of one species
                                      guni_4fsub_Tta@sam_data$Light == "Night")

## M decactis only day and only night
# Calculate the core microbiome
abund_mat_Mde_day <- as.matrix(otu_table(guni_4fsub_Mde_day))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mde_day <- ifelse(rowSums(abund_mat_Mde_day > 0) / ncol(abund_mat_Mde_day) >= core_threshold, 1, 0)
write.csv(core_matrix_Mde_day, "./core_matrix_Mde_day.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# Calculate the core microbiome
abund_mat_Mde_night <- as.matrix(otu_table(guni_4fsub_Mde_night))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mde_night <- ifelse(rowSums(abund_mat_Mde_night > 0) / ncol(abund_mat_Mde_night) >= core_threshold, 1, 0)
write.csv(core_matrix_Mde_night, "./core_matrix_Mde_night.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

## M hispida only day and only night
# Calculate the core microbiome
abund_mat_Mhi_day <- as.matrix(otu_table(guni_4fsub_Mhi_day))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mhi_day <- ifelse(rowSums(abund_mat_Mhi_day > 0) / ncol(abund_mat_Mhi_day) >= core_threshold, 1, 0)
write.csv(core_matrix_Mhi_day, "./core_matrix_Mhi_day.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# Calculate the core microbiome
abund_mat_Mhi_night <- as.matrix(otu_table(guni_4fsub_Mhi_night))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Mhi_night <- ifelse(rowSums(abund_mat_Mhi_night > 0) / ncol(abund_mat_Mhi_night) >= core_threshold, 1, 0)
write.csv(core_matrix_Mhi_night, "./core_matrix_Mhi_night.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

## T coccinea only day and only night
# Calculate the core microbiome
abund_mat_Tco_day <- as.matrix(otu_table(guni_4fsub_Tco_day))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tco_day <- ifelse(rowSums(abund_mat_Tco_day > 0) / ncol(abund_mat_Tco_day) >= core_threshold, 1, 0)
write.csv(core_matrix_Tco_day, "./core_matrix_Tco_day.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# Calculate the core microbiome
abund_mat_Tco_night <- as.matrix(otu_table(guni_4fsub_Tco_night))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tco_night <- ifelse(rowSums(abund_mat_Tco_night > 0) / ncol(abund_mat_Tco_night) >= core_threshold, 1, 0)
write.csv(core_matrix_Tco_night, "./core_matrix_Tco_night.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

## T tagusensis only day and only night
# Calculate the core microbiome
abund_mat_Tta_day <- as.matrix(otu_table(guni_4fsub_Tta_day))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tta_day <- ifelse(rowSums(abund_mat_Tta_day > 0) / ncol(abund_mat_Tta_day) >= core_threshold, 1, 0)
write.csv(core_matrix_Tta_day, "./core_matrix_Tta_day.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

# Calculate the core microbiome
abund_mat_Tta_night <- as.matrix(otu_table(guni_4fsub_Tta_night))
# Create a binary matrix indicating presence/absence based on the threshold
core_matrix_Tta_night <- ifelse(rowSums(abund_mat_Tta_night > 0) / ncol(abund_mat_Tta_night) >= core_threshold, 1, 0)
write.csv(core_matrix_Tta_night, "./core_matrix_Tta_night.csv", row.names=TRUE) #export to csv and sum column of presence/absence of genera

#.ASVs level (25833 ASVs in total)
#Mde dia 2 (ASV_4 ASV_12)<
#Mde noite 2 (ASV_12 ASV_22)*
#Mhi dia 2 (ASV_4 ASV_6)<
#Mhi noite 7 (ASV_4 ASV_6 ASV_8 ASV_91 ASV_119 ASV_413 ASV_444)*
#Tta dia 8 (ASV_1 ASV_2 ASV_4 ASV_7 ASV_11 ASV_60 ASV_98 ASV_254)<
#Tta noite 4 (ASV_1 ASV_2 ASV_4 ASV_101)*
#Tco dia 8 (ASV_1 ASV_2 ASV_4 ASV_7 ASV_112 ASV_176 ASV_219 ASV_228)<
#Tco noite 4 (ASV_1 ASV_2 ASV_4 ASV_7)*



##### Redundancy Analysis (RDA) #####

library(vegan)
library(readxl)

## Preparing spreadsheets

## Microbial data preadsheets
# 1. Transposing the ASV table from gunifrac (each sample corresponds to a row now)
guni4fsub_otu <- transpose(as.data.frame(guni_micro_coral_4.f.sub_otu))
colnames(guni4fsub_otu) <- rownames(guni_micro_coral_4.f.sub_otu)
rownames(guni4fsub_otu) <- colnames(guni_micro_coral_4.f.sub_otu)
View(guni4fsub_otu)


# 2. Dividing spreadsheet by species
#Mde
otuMde <- guni4fsub_otu[grep("Mde",rownames(guni4fsub_otu)),] # this will be used for Rinko CTD dataset
otuMde2 <- otuMde[-c(1,2), ] # this will be used for RBR CTD dataset
View(otuMde2)
#Mhi
otuMhi <- guni4fsub_otu[grep("Mhi",rownames(guni4fsub_otu)),] # this will be used for Rinko CTD dataset
otuMhi2 <- otuMhi[-c(1,2), ] # this will be used for RBR CTD dataset
View(otuMhi2)
#Tco
otuTco <- guni4fsub_otu[grep("Tc",rownames(guni4fsub_otu)),] # this will be used for Rinko CTD dataset
otuTco2 <- otuTco[-c(1,2), ] # this will be used for RBR CTD dataset
View(otuTco2)
#Tta
otuTta <- guni4fsub_otu[grep("Tt",rownames(guni4fsub_otu)),] # this will be used for Rinko CTD dataset
otuTta2 <- otuTta[-c(1,2), ] # this will be used for RBR CTD dataset
View(otuTta2)


## Environmental data spreadsheets obtained with RBR CTD
#    Prepared on excel. For rbr dataset, times in UTC were changed to Time_BR, and three spreadsheets were created: 30 minutes, 1 hour and 2 hours before collection time. They were saved as RDA_samplesRBR2hBefore.xlsx, RDA_samplesRBR1hBefore.xlsx, RDA_samplesRBR30minBefore.xlsx, each with four tabs, each corresponding to a species. Each tab was exported as a csv file and imported below.

rbrMde <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR2hBefore_Mde.csv"))
rbrMhi <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR2hBefore_Mhi.csv"))
rbrTco <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR2hBefore_Tco.csv"))
rbrTta <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR2hBefore_Tta.csv"))

rbr1hMde <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR1hBefore_Mde.csv"))
rbr1hMhi <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR1hBefore_Mhi.csv"))
rbr1hTco <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR1hBefore_Tco.csv"))
rbr1hTta <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR1hBefore_Tta.csv"))

rbr30mMde <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR30minBefore_Mde.csv"))
rbr30mMhi <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR30minBefore_Mhi.csv"))
rbr30mTta <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR30minBefore_Tta.csv"))
rbr30mTco <- as.data.frame(read.csv("SuppMatt_R_script/RDA_samplesRBR30minBefore_Tco.csv"))

## Testing for environmental variable correlation
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# If pearson correlation between any pair of variables is > 0.8, one of the variables is removed

#Correlation matrix

#2h before
cormat <- round(cor(rbrMde[-c(1,2,3)]),2)
cormat #ChlA vs ficoerit
cormat1 <- round(cor(rbrMhi[-c(1,2,3)]),2) 
cormat1 #ChlA vs ficoerit
cormat2 <- round(cor(rbrTta[-c(1,2,3)]),2) 
cormat2 #ChlA vs ficoerit
cormat3 <- round(cor(rbrTco[-c(1,2,3)]),2) 
cormat3 #ChlA vs ficoerit
#chlA vs Ficoeritrina > 0.8 (Mde -0.87, Mhi -0.87, Tta -0.86, Tco -0.86). 
#2h - remove ficoerit and leave chlA


#1h before
cormat4 <- round(cor(rbr1hMde[-c(1,2,3)]),2)
cormat4
cormat5 <- round(cor(rbr1hMhi[-c(1,2,3)]),2) 
cormat5 #ficoerit vs chlA -0.83, tirar
cormat6 <- round(cor(rbr1hTta[-c(1,2,3)]),2) 
cormat6
cormat7 <- round(cor(rbr1hTco[-c(1,2,3)]),2) 
cormat7
#1h - ficoerit vs chlA > 0.8 in Mhi, remove ficoerit. In other species, both variables were kept because values such as -0.71 and -0.74.

#30min before
cormat8 <- round(cor(rbr30mMde[-c(1)]),2)
cormat8 #ficoerit vs chlA -0.85
cormat9 <- round(cor(rbr30mMhi[-c(1)]),2)
cormat9 #ficoerit vs chlA -0.91; chlA vs sal 0.82
cormat10 <- round(cor(rbr30mTco[-c(1)]),2)
cormat10 #ficoerit vs chlA -0.90; chlA vs sal 0.80
cormat11 <- round(cor(rbr30mTta[-c(1)]),2)
cormat11 #ficoerit vs chlA -0.90; chlA vs sal 0.80
#30min - remove chlA


# Removing variables
#2h - remove ficoerit and leave chlA
rbrMde <- rbrMde[-c(6)]
rbrMhi <- rbrMhi[-c(6)]
rbrTco <- rbrTco[-c(6)]
rbrTta <- rbrTta[-c(6)]

#1h - ficoerit vs chlA > 0.8 in Mhi, remove ficoerit
rbr1hMhi <- rbr1hMhi[-c(6)]

#30min - remove chlA
rbr30mMde <- rbr30mMde[-c(3)]
rbr30mMhi <- rbr30mMhi[-c(3)]
rbr30mTco <- rbr30mTco[-c(3)]
rbr30mTta <- rbr30mTta[-c(3)]

# Organizing spreadsheets so that they only have environmental data and for sample names to be column labels

rownames(rbrMde) <- rbrMde[,1]
rbrMde <- rbrMde[-c(1,2,3)]
rownames(rbrMhi) <- rbrMhi[,1]
rbrMhi <- rbrMhi[-c(1,2,3)]
rownames(rbrTco) <- rbrTco[,1]
rbrTco <- rbrTco[-c(1,2,3)]
rownames(rbrTta) <- rbrTta[,1]
rbrTta <- rbrTta[-c(1,2,3)]
rownames(rbr1hMde) <- rbr1hMde[,1]
rbr1hMde <- rbr1hMde[-c(1,2,3)]
rownames(rbr1hMhi) <- rbr1hMhi[,1]
rbr1hMhi <- rbr1hMhi[-c(1,2,3)]
rownames(rbr1hTco) <- rbr1hTco[,1]
rbr1hTco <- rbr1hTco[-c(1,2,3)]
rownames(rbr1hTta) <- rbr1hTta[,1]
rbr1hTta <- rbr1hTta[-c(1,2,3)]

#In these datasets below, because without Time and Date, only the first column was removed
rownames(rbr30mMde) <- rbr30mMde[,1]
rbr30mMde <- rbr30mMde[-c(1)]
rownames(rbr30mMhi) <- rbr30mMhi[,1]
rbr30mMhi <- rbr30mMhi[-c(1)]
rownames(rbr30mTco) <- rbr30mTco[,1]
rbr30mTco <- rbr30mTco[-c(1)]
rownames(rbr30mTta) <- rbr30mTta[,1]
rbr30mTta <- rbr30mTta[-c(1)]


# Normalizing distances and ordering Time as a factor
# Because Time must not be normalized, only the other columns of the dataset will be normalized, and Time column will be bound back afterwards
vst_rbrMde2h_other<-decostand(rbrMde[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrMde2h <- cbind(vst_rbrMde2h_other,rbrMde[c(6)])
vst_rbrMde2h$Time<-factor(vst_rbrMde2h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMde2h) #checking if Time actually became a factor

vst_rbrMhi2h_0<-decostand(rbrMhi[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrMhi2h <- cbind(vst_rbrMhi2h_0,rbrMhi[c(6)])
vst_rbrMhi2h$Time<-factor(vst_rbrMhi2h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMhi2h) #checking if Time actually became a factor

vst_rbrTco2h_0<-decostand(rbrTco[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrTco2h <- cbind(vst_rbrTco2h_0,rbrTco[c(6)])
vst_rbrTco2h$Time<-factor(vst_rbrTco2h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTco2h) #checking if Time actually became a factor

vst_rbrTta2h_0<-decostand(rbrTta[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrTta2h <- cbind(vst_rbrTta2h_0,rbrTta[c(6)])
vst_rbrTta2h$Time<-factor(vst_rbrTta2h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTta2h) #checking if Time actually became a factor

vst_rbrMde1h_0<-decostand(rbr1hMde[-c(7)],"standardize", na.rm = TRUE) 
vst_rbrMde1h <- cbind(vst_rbrMde1h_0,rbr1hMde[c(7)])
vst_rbrMde1h$Time<-factor(vst_rbrMde1h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMde1h) #checking if Time actually became a factor

vst_rbrMhi1h_0<-decostand(rbr1hMhi[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrMhi1h <- cbind(vst_rbrMhi1h_0,rbr1hMhi[c(6)])
vst_rbrMhi1h$Time<-factor(vst_rbrMhi1h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMhi1h) #checking if Time actually became a factor

vst_rbrTco1h_0<-decostand(rbr1hTco[-c(7)],"standardize", na.rm = TRUE) 
vst_rbrTco1h <- cbind(vst_rbrTco1h_0,rbr1hTco[c(7)])
vst_rbrTco1h$Time<-factor(vst_rbrTco1h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTco1h) #checking if Time actually became a factor

vst_rbrTta1h_0<-decostand(rbr1hTta[-c(7)],"standardize", na.rm = TRUE) 
vst_rbrTta1h <- cbind(vst_rbrTta1h_0,rbr1hTta[c(7)])
vst_rbrTta1h$Time<-factor(vst_rbrTta1h$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTta1h) #checking if Time actually became a factor

vst_rbrMde30m_0<-decostand(rbr30mMde[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrMde30m <- cbind(vst_rbrMde30m_0,rbr30mMde[c(6)])
vst_rbrMde30m$Time<-factor(vst_rbrMde30m$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMde30m) #checking if Time actually became a factor

vst_rbrMhi30m_0<-decostand(rbr30mMhi[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrMhi30m <- cbind(vst_rbrMhi30m_0,rbr30mMhi[c(6)])
vst_rbrMhi30m$Time<-factor(vst_rbrMhi30m$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrMhi30m) #checking if Time actually became a factor

vst_rbrTco30m_0<-decostand(rbr30mTco[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrTco30m <- cbind(vst_rbrTco30m_0,rbr30mTco[c(6)])
vst_rbrTco30m$Time<-factor(vst_rbrTco30m$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTco30m) #checking if Time actually became a factor

vst_rbrTta30m_0<-decostand(rbr30mTta[-c(6)],"standardize", na.rm = TRUE) 
vst_rbrTta30m <- cbind(vst_rbrTta30m_0,rbr30mTta[c(6)])
vst_rbrTta30m$Time<-factor(vst_rbrTta30m$Time, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rbrTta30m) #checking if Time actually became a factor


## Plotting RDAs comparing microbial abundance and environmental data spreadsheets obtained with RBR CTD

RDA_rbrMde2h <- rda(otuMde2 ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2, vst_rbrMde2h), na.action = na.omit)
anova(RDA_rbrMde2h)
#         Df Variance      F Pr(>F)
#Model     5  3869715 0.5052  0.961
#Residual  8 12256137  
pdf('RDA_rbrMde2h.pdf')
plot(RDA_rbrMde2h)
text(RDA_rbrMde2h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde2h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


RDA_rbrMhi2h <- rda(otuMhi2 ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2, vst_rbrMhi2h), na.action = na.omit)
anova(RDA_rbrMhi2h)
#         Df Variance      F Pr(>F)
#Model     5 15628810 0.9572  0.501
#Residual  9 29388809 
pdf('RDA_rbrMhi2h.pdf')
plot(RDA_rbrMhi2h)
text(RDA_rbrMhi2h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi2h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


RDA_rbrTco2h <- rda(otuTco2 ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2, vst_rbrTco2h), na.action = na.omit)
anova(RDA_rbrTco2h)
#         Df Variance      F Pr(>F)
#Model     5  8234373 0.5875  0.803
#Residual 10 28031993 
pdf('RDA_rbrTco2h.pdf')
plot(RDA_rbrTco2h)
text(RDA_rbrTco2h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco2h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTta2h <- rda(otuTta2 ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTta2, vst_rbrTta2h), na.action = na.omit)
anova(RDA_rbrTta2h)
#         Df Variance      F Pr(>F)
#Model     5  3473143 0.3259  0.975
#Residual 10 21312151             
pdf('RDA_rbrTta2h.pdf')
plot(RDA_rbrTta2h)
text(RDA_rbrTta2h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTta2h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrMde1h <- rda(otuMde2 ~ Temperature + ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2, vst_rbrMde1h), na.action = na.omit)
anova(RDA_rbrMde1h)
#         Df Variance     F Pr(>F)
#Model     6  5011028 0.526  0.921
#Residual  7 11114824 
pdf('RDA_rbrMde1h.pdf')
plot(RDA_rbrMde1h)
text(RDA_rbrMde1h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde1h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


RDA_rbrMhi1h <- rda(otuMhi2 ~ Temperature + ChlorophyllA_corr + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2, vst_rbrMhi1h), na.action = na.omit)
anova(RDA_rbrMhi1h)
#         Df Variance      F Pr(>F)
#Model     5 21103448 1.5884  0.101
#Residual  9 23914170 
pdf('RDA_rbrMhi1h.pdf')
plot(RDA_rbrMhi1h)
text(RDA_rbrMhi1h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi1h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTco1h <- rda(otuTco2 ~ Temperature + ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2, vst_rbrTco1h), na.action = na.omit)
anova(RDA_rbrTco1h)
#         Df Variance      F Pr(>F)
#Model     6 13202203 0.8586  0.606
#Residual  9 23064163  
pdf('RDA_rbrTco1h.pdf')
plot(RDA_rbrTco1h)
text(RDA_rbrTco1h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco1h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTta1h <- rda(otuTta2 ~ Temperature + ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTta2, vst_rbrTta1h), na.action = na.omit)
anova(RDA_rbrTta1h)
#         Df Variance      F Pr(>F)
#Model     6  8872048 0.8363   0.59
#Residual  9 15913246
pdf('RDA_rbrTta1h.pdf')
plot(RDA_rbrTta1h)
text(RDA_rbrTta1h, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTta1h, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

# ficoerit and par axes were placed near each other in the plot, removing one at a time
#   without ficoerit
RDA_rbrTta1h2 <- rda(otuTta2 ~ Temperature + ChlorophyllA_corr + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTta2, vst_rbrTta1h), na.action = na.omit)
anova(RDA_rbrTta1h2)
#         Df Variance      F Pr(>F)
#Model     5  2758050 0.2504  0.991
#Residual 10 22027244

#   without PAR
RDA_rbrTta1h3 <- rda(otuTta2 ~ Temperature + ChlorophyllA_corr + Phycoerythrin + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTta2, vst_rbrTta1h), na.action = na.omit)
anova(RDA_rbrTta1h3)
#         Df Variance      F Pr(>F)
#Model     5  4455584 0.4383  0.921
#Residual 10 20329710


RDA_rbrMde30m <- rda(otuMde2 ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2, vst_rbrMde30m), na.action = na.omit)
anova(RDA_rbrMde30m)
#         Df Variance      F Pr(>F)
#Model     5  4034848 0.5339  0.897
#Residual  8 12091004
pdf('RDA_rbrMde30m.pdf')
plot(RDA_rbrMde30m)
text(RDA_rbrMde30m, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde30m, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrMhi30m <- rda(otuMhi2 ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2, vst_rbrMhi30m), na.action = na.omit)
anova(RDA_rbrMhi30m)
#         Df Variance      F Pr(>F)
#Model     5 16845709 1.0763  0.405
#Residual  9 28171910   
pdf('RDA_rbrMhi30m.pdf')
plot(RDA_rbrMhi30m)
text(RDA_rbrMhi30m, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi30m, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


RDA_rbrTco30m <- rda(otuTco2 ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2, vst_rbrTco30m), na.action = na.omit)
anova(RDA_rbrTco30m)
#         Df Variance      F Pr(>F)
#Model     5 12210509 1.0152  0.459
#Residual 10 24055857
pdf('RDA_rbrTco30m.pdf')
plot(RDA_rbrTco30m)
text(RDA_rbrTco30m, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco30m, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTta30m <- rda(otuTta2 ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTta2, vst_rbrTta30m), na.action = na.omit)
anova(RDA_rbrTta30m)
#         Df Variance      F Pr(>F)
#Model     5  7434084 0.8569  0.584
#Residual 10 17351211 
pdf('RDA_rbrTta30m.pdf')
plot(RDA_rbrTta30m)
text(RDA_rbrTta30m, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTta30m, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTta30m_wTime <- rda(otuTta2 ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas + Time, data = cbind(otuTta2, vst_rbrTta30m), na.action = na.omit)
anova(RDA_rbrTta30m_wTime)
#         Df Variance      F Pr(>F)
#Model    10 15437148 0.8257  0.665
#Residual  5  9348147 


## Running RDA without outlyers

# Mde
otuMde2_NoOutl2h <- otuMde2[!(row.names(otuMde2) %in% c('Mde13','Mde37')),]
vst_rbrMde2h_NoOutl <- vst_rbrMde2h[!(row.names(vst_rbrMde2h) %in% c('Mde13','Mde37')),]
otuMde2_NoOutl1h <- otuMde2[!(row.names(otuMde2) %in% c('Mde13')),]
vst_rbrMde1h_NoOutl <- vst_rbrMde1h[!(row.names(vst_rbrMde1h) %in% c('Mde13')),]
otuMde2_NoOutl30m <- otuMde2[!(row.names(otuMde2) %in% c('Mde13','Mde37', 'Mde65')),]
vst_rbrMde30m_NoOutl <- vst_rbrMde30m[!(row.names(vst_rbrMde30m) %in% c('Mde13','Mde37', 'Mde65')),]

RDA_rbrMde2h_NoOutl <- rda(otuMde2_NoOutl2h ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2_NoOutl2h, vst_rbrMde2h_NoOutl), na.action = na.omit)
anova(RDA_rbrMde2h_NoOutl)
#         Df Variance     F Pr(>F)
#Model     5  3250391 0.951  0.606
#Residual  6  4101505 
pdf('RDA_rbrMde2h_NoOutl.pdf')
plot(RDA_rbrMde2h_NoOutl)
text(RDA_rbrMde2h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde2h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

# chlA e PAR axes were placed near each other in the plot, removing one at a time
#   without chlA
RDA_rbrMde2h_NoOutl2 <- rda(otuMde2_NoOutl2h ~ Temperature + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2_NoOutl2h, vst_rbrMde2h_NoOutl), na.action = na.omit)
anova(RDA_rbrMde2h_NoOutl2)
#         Df Variance      F Pr(>F)
#Model     4  2769627 1.0577  0.422
#Residual  7  4582269   

#   without PAR
RDA_rbrMde2h_NoOutl3 <- rda(otuMde2_NoOutl2h ~ Temperature + ChlorophyllA_corr + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2_NoOutl2h, vst_rbrMde2h_NoOutl), na.action = na.omit)
anova(RDA_rbrMde2h_NoOutl3)
#         Df Variance      F Pr(>F)
#Model     4  2493491 0.8982  0.612
#Residual  7  4858405


RDA_rbrMde1h_NoOutl <- rda(otuMde2_NoOutl1h ~ Temperature + ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2_NoOutl1h, vst_rbrMde1h_NoOutl), na.action = na.omit)
anova(RDA_rbrMde1h_NoOutl)
#         Df Variance      F Pr(>F)
#Model     6  4691286 0.7767  0.704
#Residual  6  6040245
pdf('RDA_rbrMde1h_NoOutl.pdf')
plot(RDA_rbrMde1h_NoOutl)
text(RDA_rbrMde1h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde1h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrMde30m_NoOutl <- rda(otuMde2_NoOutl30m ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMde2_NoOutl30m, vst_rbrMde30m_NoOutl), na.action = na.omit)
anova(RDA_rbrMde30m_NoOutl)
#         Df Variance      F Pr(>F)
#Model     5  3860427 1.2485  0.244
#Residual  5  3092035
pdf('RDA_rbrMde30m_NoOutl.pdf')
plot(RDA_rbrMde30m_NoOutl)
text(RDA_rbrMde30m_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMde30m_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


### Mhi 
otuMhi2_NoOutl <- otuMhi2[!(row.names(otuMhi2) %in% c('Mhi14','Mhi26')),]
vst_rbrMhi2h_NoOutl <- vst_rbrMhi2h[!(row.names(vst_rbrMhi2h) %in% c('Mhi14','Mhi26')),]
vst_rbrMhi1h_NoOutl <- vst_rbrMhi1h[!(row.names(vst_rbrMhi1h) %in% c('Mhi14','Mhi26')),]
vst_rbrMhi30m_NoOutl <- vst_rbrMhi30m[!(row.names(vst_rbrMhi30m) %in% c('Mhi14','Mhi26')),]

RDA_rbrMhi2h_NoOutl <- rda(otuMhi2_NoOutl ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2_NoOutl, vst_rbrMhi2h_NoOutl), na.action = na.omit)
anova(RDA_rbrMhi2h_NoOutl)
#         Df Variance      F Pr(>F)
#Model     5 13157244 0.9799  0.494
#Residual  7 18797493 
pdf('RDA_rbrMhi2h_NoOutl.pdf')
plot(RDA_rbrMhi2h_NoOutl)
text(RDA_rbrMhi2h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi2h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrMhi1h_NoOutl <- rda(otuMhi2_NoOutl ~ Temperature + ChlorophyllA_corr + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2_NoOutl, vst_rbrMhi1h_NoOutl), na.action = na.omit)
anova(RDA_rbrMhi1h_NoOutl)
#         Df Variance      F Pr(>F)
#Model     5 16657109 1.5244  0.154
#Residual  7 15297629   
pdf('RDA_rbrMhi1h_NoOutl.pdf')
plot(RDA_rbrMhi1h_NoOutl)
text(RDA_rbrMhi1h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi1h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrMhi30m_NoOutl <- rda(otuMhi2_NoOutl ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuMhi2_NoOutl, vst_rbrMhi30m_NoOutl), na.action = na.omit)
anova(RDA_rbrMhi30m_NoOutl)
#         Df Variance     F Pr(>F)
#Model     5 15458941 1.312  0.261
#Residual  7 16495796 
RDA_rbrMhi30m_NoOutl0 <- rda(otuMhi2_NoOutl ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas + Time, data = cbind(otuMhi2_NoOutl, vst_rbrMhi30m_NoOutl), na.action = na.omit)
anova(RDA_rbrMhi30m_NoOutl0)
#         Df Variance      F Pr(>F)
#Model    10 28903057 1.8942  0.189
#Residual  2  3051680      
pdf('RDA_rbrMhi30m_NoOutl.pdf')
plot(RDA_rbrMhi30m_NoOutl)
text(RDA_rbrMhi30m_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrMhi30m_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


### Tco
otuTco2_NoOutl <- otuTco2[!(row.names(otuTco2) %in% c('Tc19')),]
vst_rbrTco2h_NoOutl <- vst_rbrTco2h[!(row.names(vst_rbrTco2h) %in% c('Tc19')),]
vst_rbrTco1h_NoOutl <- vst_rbrTco1h[!(row.names(vst_rbrTco1h) %in% c('Tc19')),]
vst_rbrTco30m_NoOutl <- vst_rbrTco30m[!(row.names(vst_rbrTco30m) %in% c('Tc19')),]

RDA_rbrTco2h_NoOutl <- rda(otuTco2_NoOutl ~ Temperature + ChlorophyllA_corr + PAR_corr_certo + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco2h_NoOutl), na.action = na.omit)
anova(RDA_rbrTco2h_NoOutl)
#         Df Variance      F Pr(>F)
#Model     5  6329651 0.5892  0.759
#Residual  9 19338618       
pdf('RDA_rbrTco2h_NoOutl.pdf')
plot(RDA_rbrTco2h_NoOutl)
text(RDA_rbrTco2h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco2h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rbrTco1h_NoOutl <- rda(otuTco2_NoOutl ~ Temperature + ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco1h_NoOutl), na.action = na.omit)
anova(RDA_rbrTco1h_NoOutl)
#         Df Variance      F Pr(>F)
#Model     6 12329164 1.2324  0.325
#Residual  8 13339105
pdf("RDA_rbrTco1h_NoOutl.pdf")
plot(RDA_rbrTco1h_NoOutl)
text(RDA_rbrTco1h_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco1h_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

# Temp, PAR and depth axes were placed near each other in the plot, removing two at a time
#   without Temp e depth
RDA_rbrTco1h_NoOutl2 <- rda(otuTco2_NoOutl ~ ChlorophyllA_corr + Phycoerythrin + PAR_corrigido + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco1h_NoOutl), na.action = na.omit)
anova(RDA_rbrTco1h_NoOutl2)
#         Df Variance      F Pr(>F)
#Model     4  8865810 1.3191   0.29
#Residual 10 16802459

#   without Temp e PAR
RDA_rbrTco1h_NoOutl3 <- rda(otuTco2_NoOutl ~ ChlorophyllA_corr + Phycoerythrin + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco1h_NoOutl), na.action = na.omit)
anova(RDA_rbrTco1h_NoOutl3)
#         Df Variance      F Pr(>F)
#Model     4  6600706 0.8654  0.555
#Residual 10 19067563

#   without depth e PAR
RDA_rbrTco1h_NoOutl4 <- rda(otuTco2_NoOutl ~ Temperature + ChlorophyllA_corr + Phycoerythrin + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco1h_NoOutl), na.action = na.omit)
anova(RDA_rbrTco1h_NoOutl4)
#         Df Variance      F Pr(>F)
#Model     4  8799355 1.3041  0.331
#Residual 10 16868914


RDA_rbrTco30m_NoOutl <- rda(otuTco2_NoOutl ~ Temperature + Phycoerythrin + PAR_corrigido + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco30m_NoOutl), na.action = na.omit)
anova(RDA_rbrTco30m_NoOutl)
#         Df Variance      F Pr(>F)
#Model     5  9710175 1.0953    0.4
#Residual  9 15958094
pdf("RDA_rbrTco30m_NoOutl.pdf")
plot(RDA_rbrTco30m_NoOutl)
text(RDA_rbrTco30m_NoOutl, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rbrTco30m_NoOutl, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

# Depth e PAR axes were placed near each other in the plot, removing one at a time
# Without depth
RDA_rbrTco30m_NoOutl2 <- rda(otuTco2_NoOutl ~ Temperature + Phycoerythrin + PAR_corrigido + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco30m_NoOutl), na.action = na.omit)
anova(RDA_rbrTco30m_NoOutl2)
#         Df Variance     F Pr(>F)
#Model     4  8619884 1.264  0.288
#Residual 10 17048385

# Without PAR
RDA_rbrTco30m_NoOutl3 <- rda(otuTco2_NoOutl ~ Temperature + Phycoerythrin + Depth_CORRIGIDA + Salinity_CORRIGIDAmedianas, data = cbind(otuTco2_NoOutl, vst_rbrTco30m_NoOutl), na.action = na.omit)
anova(RDA_rbrTco30m_NoOutl3)
#         Df Variance      F Pr(>F)
#Model     4  7706558 1.0726  0.412
#Residual 10 17961711


## Environmental data spreadsheets obtained with Rinko CTD
#   Rinko was deployed at collection time (colltime) and 2 hours before (2hbefore) collection. Because Rinko CTD provided a profile of the water column, data from 11 m (rk_p11m) and 12 m (rk_p12m) deep were used in the analysis

## Environmental data spreadsheets

rk_p11m_colltime <- as.data.frame(read.csv("SuppMatt_R_script/tabela_profsRINKO_isa_modifIsa_p11m_colltime.csv"))
rk_p11m_2hbefore <- as.data.frame(read.csv("SuppMatt_R_script/tabela_profsRINKO_isa_modifIsa_p11m_2hbefore.csv"))
rk_p12m_colltime <- as.data.frame(read.csv("SuppMatt_R_script/tabela_profsRINKO_isa_modifIsa_p12m_colltime.csv"))
rk_p12m_2hbefore <- as.data.frame(read.csv("SuppMatt_R_script/tabela_profsRINKO_isa_modifIsa_p12m_2hbefore.csv"))

#rk_p11m_colltime
rk_p11m_colltime_Mde <- rk_p11m_colltime[c(-9,-12),c(1,5:12)]
row.names(rk_p11m_colltime_Mde) <- rk_p11m_colltime_Mde$Md
rk_p11m_colltime_Mhi <- rk_p11m_colltime[c(-14),c(2,5:12)]
row.names(rk_p11m_colltime_Mhi) <- rk_p11m_colltime_Mhi$Mh
rk_p11m_colltime_Tco <- rk_p11m_colltime[c(3,5:12)]
row.names(rk_p11m_colltime_Tco) <- rk_p11m_colltime_Tco$Tc
rk_p11m_colltime_Tta <- rk_p11m_colltime[c(4,5:12)]
row.names(rk_p11m_colltime_Tta) <- rk_p11m_colltime_Tta$Tt

#rk_p11m_2hbefore
rk_p11m_2hbefore_Mde <- rk_p11m_2hbefore[c(-7),c(1,5:12)]
row.names(rk_p11m_2hbefore_Mde) <- rk_p11m_2hbefore_Mde$Md
rk_p11m_2hbefore_Mhi <- rk_p11m_2hbefore[c(-11),c(2,5:12)]
row.names(rk_p11m_2hbefore_Mhi) <- rk_p11m_2hbefore_Mhi$Mh
rk_p11m_2hbefore_Tco <- rk_p11m_2hbefore[,c(3,5:12)]
row.names(rk_p11m_2hbefore_Tco) <- rk_p11m_2hbefore_Tco$Tc
rk_p11m_2hbefore_Tta <- rk_p11m_2hbefore[,c(4,5:12)]
row.names(rk_p11m_2hbefore_Tta) <- rk_p11m_2hbefore_Tta$Tt

#rk_p12m_colltime
rk_p12m_colltime_Mde <- rk_p12m_colltime[c(-9,-12),c(1,5:12)]
row.names(rk_p12m_colltime_Mde) <- rk_p12m_colltime_Mde$Md
rk_p12m_colltime_Mhi <- rk_p12m_colltime[c(-14),c(2,5:12)]
row.names(rk_p12m_colltime_Mhi) <- rk_p12m_colltime_Mhi$Mh
rk_p12m_colltime_Tco <- rk_p12m_colltime[,c(3,5:12)]
row.names(rk_p12m_colltime_Tco) <- rk_p12m_colltime_Tco$Tc
rk_p12m_colltime_Tta <- rk_p12m_colltime[,c(4,5:12)]
row.names(rk_p12m_colltime_Tta) <- rk_p12m_colltime_Tta$Tt

#rk_p12m_2hbefore
rk_p12m_2hbefore_Mde <- rk_p12m_2hbefore[c(-7),c(1,5:12)]
row.names(rk_p12m_2hbefore_Mde) <- rk_p12m_2hbefore_Mde$Md
rk_p12m_2hbefore_Mhi <- rk_p12m_2hbefore[c(-11),c(2,5:12)]
row.names(rk_p12m_2hbefore_Mhi) <- rk_p12m_2hbefore_Mhi$Mh
rk_p12m_2hbefore_Tco <- rk_p12m_2hbefore[,c(3,5:12)]
row.names(rk_p12m_2hbefore_Tco) <- rk_p12m_2hbefore_Tco$Tc
rk_p12m_2hbefore_Tta <- rk_p12m_2hbefore[,c(4,5:12)]
row.names(rk_p12m_2hbefore_Tta) <- rk_p12m_2hbefore_Tta$Tt


## Testing for environmental variable correlation
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# If pearson correlation between any pair of variables is > 0.8, one of the variables is removed

#11m collection time --- remove turb
cormatRk1 <- round(cor(rk_p11m_colltime_Mde[-c(1,2)]),2)
cormatRk1 #DO vs turb -0.81
cormatRk2 <- round(cor(rk_p11m_colltime_Mhi[-c(1,2)]),2)
cormatRk2 #DO vs turb -0.81
cormatRk3 <- round(cor(rk_p11m_colltime_Tco[-c(1,2)]),2)
cormatRk3 #DO vs turb -0.79
cormatRk4 <- round(cor(rk_p11m_colltime_Tta[-c(1,2)]),2)
cormatRk4 #DO vs turb -0.79

#11m 2h before --- remove DO
cormatRk5 <- round(cor(rk_p11m_2hbefore_Mde[-c(1,2)]),2)
cormatRk5 #DO vs temp 0.89; DO vs turb -0.78; turb vs  temp -0.76
cormatRk6 <- round(cor(rk_p11m_2hbefore_Mhi[-c(1,2)]),2)
cormatRk6 #DO vs temp 0.89; DO vs turb -0.78; turb vs temp -0.78; chl vs turb 0.77
cormatRk7 <- round(cor(rk_p11m_2hbefore_Tco[-c(1,2)]),2)
cormatRk7 #DO vs temp 0.89; DO vs turb -0.78; turb vs temp -0.76; chl vs turb 0.74
cormatRk8 <- round(cor(rk_p11m_2hbefore_Tta[-c(1,2)]),2)
cormatRk8 #DO vs temp 0.89; DO vs turb -0.78; turb vs temp -0.76; chl vs turb 0.74

#12m collection time --- remove sal
cormatRk9 <- round(cor(rk_p12m_colltime_Mde[-c(1,2)]),2)
cormatRk9 #sal vs temp 0.79
cormatRk10 <- round(cor(rk_p12m_colltime_Mhi[-c(1,2)]),2)
cormatRk10 #sal vs temp 0.79
cormatRk11 <- round(cor(rk_p12m_colltime_Tco[-c(1,2)]),2)
cormatRk11 #sal vs temp 0.80 <---
cormatRk12 <- round(cor(rk_p12m_colltime_Tta[-c(1,2)]),2)
cormatRk12 #sal vs temp 0.80 <---

#12m 2hbefore --- remove DO
cormatRk13 <- round(cor(rk_p12m_2hbefore_Mde[-c(1,2)]),2)
cormatRk13 #DO vs temp 0.79; DO vs turb -0.76
cormatRk14 <- round(cor(rk_p12m_2hbefore_Mhi[-c(1,2)]),2)
cormatRk14 #DO vs temp 0.79; DO vs turb -0.76
cormatRk15 <- round(cor(rk_p12m_2hbefore_Tco[-c(1,2)]),2)
cormatRk15 #DO vs temp 0.79; DO vs turb -0.76
cormatRk16 <- round(cor(rk_p12m_2hbefore_Tta[-c(1,2)]),2)
cormatRk16 #DO vs temp 0.79; DO vs turb -0.76

## Removing variables

#11m collection time --- remove turb
rk_p11m_colltime_Mde1 <- rk_p11m_colltime_Mde[,-c(6)]
rk_p11m_colltime_Mhi1 <- rk_p11m_colltime_Mhi[,-c(6)]
rk_p11m_colltime_Tco1 <- rk_p11m_colltime_Tco[,-c(6)]
rk_p11m_colltime_Tta1 <- rk_p11m_colltime_Tta[,-c(6)]

#11m 2h before --- remove DO
rk_p11m_2hbefore_Mde1 <- rk_p11m_2hbefore_Mde[,-c(8)]
rk_p11m_2hbefore_Mhi1 <- rk_p11m_2hbefore_Mhi[,-c(8)]
rk_p11m_2hbefore_Tco1 <- rk_p11m_2hbefore_Tco[,-c(8)]
rk_p11m_2hbefore_Tta1 <- rk_p11m_2hbefore_Tta[,-c(8)]

#12m collection time --- remove sal
rk_p12m_colltime_Mde1 <- rk_p12m_colltime_Mde[,-c(4)]
rk_p12m_colltime_Mhi1 <- rk_p12m_colltime_Mhi[,-c(4)]
rk_p12m_colltime_Tco1 <- rk_p12m_colltime_Tco[,-c(4)]
rk_p12m_colltime_Tta1 <- rk_p12m_colltime_Tta[,-c(4)]

#12m 2h before --- remove DO
rk_p12m_2hbefore_Mde1 <- rk_p12m_2hbefore_Mde[,-c(8)]
rk_p12m_2hbefore_Mhi1 <- rk_p12m_2hbefore_Mhi[,-c(8)]
rk_p12m_2hbefore_Tco1 <- rk_p12m_2hbefore_Tco[,-c(8)]
rk_p12m_2hbefore_Tta1 <- rk_p12m_2hbefore_Tta[,-c(8)]


# Normalizing distances and ordering Time as a factor
# Because Time must not be normalized, only the other columns of the dataset will be normalized, and Time column will be bound back afterwards

#11m collection time
vst_rk_p11m_colltime_Mde_other<-decostand(rk_p11m_colltime_Mde1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_colltime_Mde1 <- cbind(rk_p11m_colltime_Mde1[c(2)],vst_rk_p11m_colltime_Mde_other)
vst_rk_p11m_colltime_Mde1$horario<-factor(vst_rk_p11m_colltime_Mde1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p11m_colltime_Mde1) #checking if Time actually became a factor

vst_rk_p11m_colltime_Mhi_other<-decostand(rk_p11m_colltime_Mhi1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_colltime_Mhi1 <- cbind(rk_p11m_colltime_Mhi1[c(2)],vst_rk_p11m_colltime_Mhi_other)
vst_rk_p11m_colltime_Mhi1$horario<-factor(vst_rk_p11m_colltime_Mhi1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p11m_colltime_Mhi1) #checking if Time actually became a factor

vst_rk_p11m_colltime_Tco_other<-decostand(rk_p11m_colltime_Tco1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_colltime_Tco1 <- cbind(rk_p11m_colltime_Tco1[c(2)],vst_rk_p11m_colltime_Tco_other)
vst_rk_p11m_colltime_Tco1$horario<-factor(vst_rk_p11m_colltime_Tco1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p11m_colltime_Tco1) #checking if Time actually became a factor

vst_rk_p11m_colltime_Tta_other<-decostand(rk_p11m_colltime_Tta1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_colltime_Tta1 <- cbind(rk_p11m_colltime_Tta1[c(2)],vst_rk_p11m_colltime_Tta_other)
vst_rk_p11m_colltime_Tta1$horario<-factor(vst_rk_p11m_colltime_Tta1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p11m_colltime_Tta1) #checking if Time actually became a factor

#11m 2h before
vst_rk_p11m_2hbefore_Mde_other<-decostand(rk_p11m_2hbefore_Mde1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_2hbefore_Mde1 <- cbind(rk_p11m_2hbefore_Mde1[c(2)],vst_rk_p11m_2hbefore_Mde_other)
vst_rk_p11m_2hbefore_Mde1$horario<-factor(vst_rk_p11m_2hbefore_Mde1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p11m_2hbefore_Mde1) #checking if Time actually became a factor

vst_rk_p11m_2hbefore_Mhi_other<-decostand(rk_p11m_2hbefore_Mhi1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_2hbefore_Mhi1 <- cbind(rk_p11m_2hbefore_Mhi1[c(2)],vst_rk_p11m_2hbefore_Mhi_other)
vst_rk_p11m_2hbefore_Mhi1$horario<-factor(vst_rk_p11m_2hbefore_Mhi1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p11m_2hbefore_Mhi1) #checking if Time actually became a factor

vst_rk_p11m_2hbefore_Tco_other<-decostand(rk_p11m_2hbefore_Tco1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_2hbefore_Tco1 <- cbind(rk_p11m_2hbefore_Tco1[c(2)],vst_rk_p11m_2hbefore_Tco_other)
vst_rk_p11m_2hbefore_Tco1$horario<-factor(vst_rk_p11m_2hbefore_Tco1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p11m_2hbefore_Tco1) #checking if Time actually became a factor

vst_rk_p11m_2hbefore_Tta_other<-decostand(rk_p11m_2hbefore_Tta1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p11m_2hbefore_Tta1 <- cbind(rk_p11m_2hbefore_Tta1[c(2)],vst_rk_p11m_2hbefore_Tta_other)
vst_rk_p11m_2hbefore_Tta1$horario<-factor(vst_rk_p11m_2hbefore_Tta1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p11m_2hbefore_Tta1) #checking if Time actually became a factor

#12m collection time
vst_rk_p12m_colltime_Mde1_other<-decostand(rk_p12m_colltime_Mde1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_colltime_Mde1 <- cbind(rk_p12m_colltime_Mde1[c(2)],vst_rk_p12m_colltime_Mde1_other)
vst_rk_p12m_colltime_Mde1$horario<-factor(vst_rk_p12m_colltime_Mde1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p12m_colltime_Mde1) #checking if Time actually became a factor

vst_rk_p12m_colltime_Mhi1_other<-decostand(rk_p12m_colltime_Mhi1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_colltime_Mhi1 <- cbind(rk_p12m_colltime_Mhi1[c(2)],vst_rk_p12m_colltime_Mhi1_other)
vst_rk_p12m_colltime_Mhi1$horario<-factor(vst_rk_p12m_colltime_Mhi1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p12m_colltime_Mhi1) #checking if Time actually became a factor

vst_rk_p12m_colltime_Tco1_other<-decostand(rk_p12m_colltime_Tco1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_colltime_Tco1 <- cbind(rk_p12m_colltime_Tco1[c(2)],vst_rk_p12m_colltime_Tco1_other)
vst_rk_p12m_colltime_Tco1$horario<-factor(vst_rk_p12m_colltime_Tco1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p12m_colltime_Tco1) #checking if Time actually became a factor

vst_rk_p12m_colltime_Tta1_other<-decostand(rk_p12m_colltime_Tta1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_colltime_Tta1 <- cbind(rk_p12m_colltime_Tta1[c(2)],vst_rk_p12m_colltime_Tta1_other)
vst_rk_p12m_colltime_Tta1$horario<-factor(vst_rk_p12m_colltime_Tta1$horario, ordered = TRUE, levels = c('0','4','8','12','16','20'))
str(vst_rk_p12m_colltime_Tta1) #checking if Time actually became a factor

#12m 2h before
vst_rk_p12m_2hbefore_Mde_other<-decostand(rk_p12m_2hbefore_Mde1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_2hbefore_Mde1 <- cbind(rk_p12m_2hbefore_Mde1[c(2)],vst_rk_p12m_2hbefore_Mde_other)
vst_rk_p12m_2hbefore_Mde1$horario<-factor(vst_rk_p12m_2hbefore_Mde1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p12m_2hbefore_Mde1) #checking if Time actually became a factor

vst_rk_p12m_2hbefore_Mhi_other<-decostand(rk_p12m_2hbefore_Mhi1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_2hbefore_Mhi1 <- cbind(rk_p12m_2hbefore_Mhi1[c(2)],vst_rk_p12m_2hbefore_Mhi_other)
vst_rk_p12m_2hbefore_Mhi1$horario<-factor(vst_rk_p12m_2hbefore_Mhi1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p12m_2hbefore_Mhi1) #checking if Time actually became a factor

vst_rk_p12m_2hbefore_Tco_other<-decostand(rk_p12m_2hbefore_Tco1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_2hbefore_Tco1 <- cbind(rk_p12m_2hbefore_Tco1[c(2)],vst_rk_p12m_2hbefore_Tco_other)
vst_rk_p12m_2hbefore_Tco1$horario<-factor(vst_rk_p12m_2hbefore_Tco1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p12m_2hbefore_Tco1) #checking if Time actually became a factor

vst_rk_p12m_2hbefore_Tta_other<-decostand(rk_p12m_2hbefore_Tta1[-c(1,2)],"standardize", na.rm = TRUE) 
vst_rk_p12m_2hbefore_Tta1 <- cbind(rk_p12m_2hbefore_Tta1[c(2)],vst_rk_p12m_2hbefore_Tta_other)
vst_rk_p12m_2hbefore_Tta1$horario<-factor(vst_rk_p12m_2hbefore_Tta1$horario, ordered = TRUE, levels = c('2','6','10','14','18','22'))
str(vst_rk_p12m_2hbefore_Tta1) #checking if Time actually became a factor


## Plotting RDAs comparing microbial abundance and environmental data spreadsheets obtained with Rinko CTD

#11m collection time
RDA_rk_p11m_colltime_Mde1 <- rda(otuMde ~ temp + sal + chl + pH + DO + PAR, data = cbind(otuMde, vst_rk_p11m_colltime_Mde1), na.action = na.omit)
anova(RDA_rk_p11m_colltime_Mde1)
#         Df Variance      F Pr(>F)
#Model     6  7055871 0.7639  0.708
#Residual  9 13854619  
pdf("RDA_rk_p11m_colltime_Mde1.pdf")
plot(RDA_rk_p11m_colltime_Mde1)
text(RDA_rk_p11m_colltime_Mde1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_colltime_Mde1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_colltime_Mhi1 <- rda(otuMhi ~ temp + sal + chl + pH + DO + PAR, data = cbind(otuMhi, vst_rk_p11m_colltime_Mhi1), na.action = na.omit)
anova(RDA_rk_p11m_colltime_Mhi1)
#         Df Variance      F Pr(>F)
#Model     6 17963405 1.2244  0.296
#Residual 10 24451697 
pdf("RDA_rk_p11m_colltime_Mhi1.pdf")
plot(RDA_rk_p11m_colltime_Mhi1)
text(RDA_rk_p11m_colltime_Mhi1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_colltime_Mhi1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_colltime_Tco1 <- rda(otuTco ~ temp + sal + chl + pH + DO + PAR, data = cbind(otuTco, vst_rk_p11m_colltime_Tco1), na.action = na.omit)
anova(RDA_rk_p11m_colltime_Tco1)
#         Df Variance      F Pr(>F)
#Model     6  9159262 0.6269   0.77
#Residual 11 26786609   
pdf("RDA_rk_p11m_colltime_Tco1.pdf")
plot(RDA_rk_p11m_colltime_Tco1)
text(RDA_rk_p11m_colltime_Tco1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_colltime_Tco1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_colltime_Tta1 <- rda(otuTta ~ temp + sal + chl + pH + DO + PAR, data = cbind(otuTta, vst_rk_p11m_colltime_Tta1), na.action = na.omit)
anova(RDA_rk_p11m_colltime_Tta1)
#         Df Variance      F Pr(>F)
#Model     6  8930234 1.1555  0.374
#Residual 11 14168978
pdf("RDA_rk_p11m_colltime_Tta1.pdf")
plot(RDA_rk_p11m_colltime_Tta1)
text(RDA_rk_p11m_colltime_Tta1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_colltime_Tta1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

#11m 2h before
RDA_rk_p11m_2hbefore_Mde1 <- rda(otuMde[-c(1,6),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuMde[-c(1,6),], vst_rk_p11m_2hbefore_Mde1), na.action = na.omit)
#otuMde[-c(1,6),] was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p11m_2hbefore_Mde1)
#         Df Variance      F Pr(>F)
#Model     6  7067275 0.9177  0.578
#Residual  7  8984314 
pdf("RDA_rk_p11m_2hbefore_Mde1.pdf")
plot(RDA_rk_p11m_2hbefore_Mde1)
text(RDA_rk_p11m_2hbefore_Mde1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_2hbefore_Mde1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_2hbefore_Mhi1 <- rda(otuMhi[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuMhi[-c(1,6,12),], vst_rk_p11m_2hbefore_Mhi1), na.action = na.omit)
#otuMhi[-c(1,6,12),]  was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p11m_2hbefore_Mhi1)
#         Df Variance      F Pr(>F)
#Model     6 24435438 1.3446  0.224
#Residual  7 21202105  
pdf("RDA_rk_p11m_2hbefore_Mhi1.pdf")
plot(RDA_rk_p11m_2hbefore_Mhi1)
text(RDA_rk_p11m_2hbefore_Mhi1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_2hbefore_Mhi1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_2hbefore_Tco1 <- rda(otuTco[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuTco[-c(1,6,12),], vst_rk_p11m_2hbefore_Tco1), na.action = na.omit)
anova(RDA_rk_p11m_2hbefore_Tco1)
#         Df Variance      F Pr(>F)
#Model     6 12348421 0.7581  0.691
#Residual  8 21718190    
pdf("RDA_rk_p11m_2hbefore_Tco1.pdf")
plot(RDA_rk_p11m_2hbefore_Tco1)
text(RDA_rk_p11m_2hbefore_Tco1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_2hbefore_Tco1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p11m_2hbefore_Tta1 <- rda(otuTta[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuTta[-c(1,6,12),], vst_rk_p11m_2hbefore_Tta1), na.action = na.omit)
anova(RDA_rk_p11m_2hbefore_Tta1)
#         Df Variance      F Pr(>F)
#Model     6  9513203 0.7774  0.667
#Residual  8 16316104   
pdf("RDA_rk_p11m_2hbefore_Tta1.pdf")
plot(RDA_rk_p11m_2hbefore_Tta1)
text(RDA_rk_p11m_2hbefore_Tta1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p11m_2hbefore_Tta1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


#12m collection time
RDA_rk_p12m_colltime_Mde1 <- rda(otuMde ~ temp + chl + turb + pH + DO + PAR, data = cbind(otuMde, vst_rk_p12m_colltime_Mde1), na.action = na.omit)
anova(RDA_rk_p12m_colltime_Mde1)
#         Df Variance     F Pr(>F)
#Model     6  6963990 0.749  0.769
#Residual  9 13946500  
pdf("RDA_rk_p12m_colltime_Mde1.pdf")
plot(RDA_rk_p12m_colltime_Mde1)
text(RDA_rk_p12m_colltime_Mde1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_colltime_Mde1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_colltime_Mhi1 <- rda(otuMhi ~ temp + chl + turb + pH + DO + PAR, data = cbind(otuMhi, vst_rk_p12m_colltime_Mhi1), na.action = na.omit)
anova(RDA_rk_p12m_colltime_Mhi1)
#         Df Variance      F Pr(>F)
#Model     6 16199982 1.0299  0.444
#Residual 10 26215120  
pdf("RDA_rk_p12m_colltime_Mhi1.pdf")
plot(RDA_rk_p12m_colltime_Mhi1)
text(RDA_rk_p12m_colltime_Mhi1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_colltime_Mhi1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_colltime_Tco1 <- rda(otuTco ~ temp + chl + turb + pH + DO + PAR, data = cbind(otuTco, vst_rk_p12m_colltime_Tco1), na.action = na.omit)
anova(RDA_rk_p12m_colltime_Tco1)
#         Df Variance      F Pr(>F)
#Model     6  9715656 0.6791  0.739
#Residual 11 26230215 
pdf("RDA_rk_p12m_colltime_Tco1.pdf")
plot(RDA_rk_p12m_colltime_Tco1)
text(RDA_rk_p12m_colltime_Tco1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_colltime_Tco1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_colltime_Tta1 <- rda(otuTta ~ temp + chl + turb + pH + DO + PAR, data = cbind(otuTta, vst_rk_p12m_colltime_Tta1), na.action = na.omit)
anova(RDA_rk_p12m_colltime_Tta1)
#         Df Variance     F Pr(>F)
#Model     6  9605396 1.305  0.289
#Residual 11 13493815     
pdf("RDA_rk_p12m_colltime_Tta1.pdf")
plot(RDA_rk_p12m_colltime_Tta1)
text(RDA_rk_p12m_colltime_Tta1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_colltime_Tta1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

#12m 2h before
RDA_rk_p12m_2hbefore_Mde1 <- rda(otuMde[-c(1,6),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuMde[-c(1,6),], vst_rk_p12m_2hbefore_Mde1), na.action = na.omit)
#otuMde[-c(1,6),]  was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p12m_2hbefore_Mde1)
#         Df Variance      F Pr(>F)
#Model     6  8127068 1.1965  0.294
#Residual  7  7924520              
pdf("RDA_rk_p12m_2hbefore_Mde1.pdf")
plot(RDA_rk_p12m_2hbefore_Mde1)
text(RDA_rk_p12m_2hbefore_Mde1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_2hbefore_Mde1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_2hbefore_Mhi1 <- rda(otuMhi[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuMhi[-c(1,6,12),], vst_rk_p12m_2hbefore_Mhi1), na.action = na.omit)
#otuMhi[-c(1,6,12),] was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p12m_2hbefore_Mhi1)
#         Df Variance      F Pr(>F)
#Model     6 25684782 1.5018  0.179
#Residual  7 19952761
pdf("RDA_rk_p12m_2hbefore_Mhi1.pdf")
plot(RDA_rk_p12m_2hbefore_Mhi1)
text(RDA_rk_p12m_2hbefore_Mhi1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_2hbefore_Mhi1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_2hbefore_Tco1 <- rda(otuTco[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuTco[-c(1,6,12),], vst_rk_p12m_2hbefore_Tco1), na.action = na.omit)
#otuMhi[-c(1,6,12),] was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p12m_2hbefore_Tco1)
#         Df Variance      F Pr(>F)
#Model     6 10334582 0.5806  0.826
#Residual  8 23732028   
pdf("RDA_rk_p12m_2hbefore_Tco1.pdf")
plot(RDA_rk_p12m_2hbefore_Tco1)
text(RDA_rk_p12m_2hbefore_Tco1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_2hbefore_Tco1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()

RDA_rk_p12m_2hbefore_Tta1 <- rda(otuTta[-c(1,6,12),] ~ temp + sal + chl + turb + pH + PAR, data = cbind(otuTta[-c(1,6,12),], vst_rk_p12m_2hbefore_Tta1), na.action = na.omit)
#otuTta[-c(1,6,12),]  was used because some time points were not available in the Rinko CTD dataset
anova(RDA_rk_p12m_2hbefore_Tta1)
#         Df Variance      F Pr(>F)
#Model     6 11874152 1.1345  0.392
#Residual  8 13955154 
pdf("RDA_rk_p12m_2hbefore_Tta1.pdf")
plot(RDA_rk_p12m_2hbefore_Tta1)
text(RDA_rk_p12m_2hbefore_Tta1, display = "sites", cex = 0.8, col = "black", font = 2)
text(RDA_rk_p12m_2hbefore_Tta1, display = "species", cex = 0.8, col = "red", font = 2)
dev.off()


