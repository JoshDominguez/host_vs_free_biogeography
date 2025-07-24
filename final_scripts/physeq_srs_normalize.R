# onloading phyloseq object from dada2 pipeline and performing SRS normalization
# goal of script - investigate quality of run, normalize samples, filter out ASVs and samples that will not be used

# loading packages 
packages <- c("ggplot2", "gridExtra", "knitr","phyloseq", "tidyverse", "vegan", "MuMIn", "dplyr", "ecodist", "glmmTMB", "glue", "ggforce", "geosphere", "lmPerm", "otuSummary", "ANCOMBC", "DECIPHER", "phangorn", "ape", "msa","gdm", "GUniFrac", "lme4", 'ggpubr', 'SRS' ,'car', 'reshape2', 'ggsignif', 'phytools', 'picante')
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

#### load in phyloseq object #### 
physeq_raw <- readRDS('data/dada2/ps.tree.RDS')

# take out repeat samples 
physeq_raw <- subset_samples(physeq_raw, !sampleNames %in% c("TR_74_S130" , "TR_10_S131", "TR_41_S132", "TR_52_S133","TR_48_S135","TR_01_S82","TR_61_S137","TR_36_S138", "TR_57_S139","TR_49_S140","TR_14_S141","TR_65_S142", "01_16S_S1", "28_16S_S28", "33_16S_S33", "42_16S_S42", "45_16S_S45", "57_16S_S57", "63_16S_S63", "66_16S_S66", "74_16S_S74", "78_16S_S78", "YoTB_623_S306", "YoTB_624_S307", "YoTB_626_S309", "YoTB_627_S310", "YoTB_628_S311", "YoTB_629_S312", "YoTB_630_S313", "YoTB_631_S314", "YoTB_636_S319"))

#### Filtering Chloroplasts and unwanted samples ####
# Double check and filter out chloroplasts and mitochondria  
physeq_raw <-  subset_taxa(physeq_raw, Order!= "Chloroplast" | is.na(Order))

# Remove Blue sample 
physeq_raw <- subset_samples(physeq_raw, !site == "Blue")

# remove stickleback from ormund lake - no water samples 
physeq_raw <- subset_samples(physeq_raw, !site == "ormund")

# Remove any asvs that were unique to those pops 
physeq_raw <- prune_taxa(taxa_sums(physeq_raw) > 0, physeq_raw) # now we have 26391 taxa and 862 samples after filtering out repeats, blue, and ormund
# how many samples of each type pre rarefaction? 
table(sample_data(physeq_raw)$sample_type) # 184 trout, 541 stickle, 63 sierra, 74 water VI 
nrow(sample_data(physeq_raw)) # 862 samples

#### Read Depth Examination ####
sample_data(physeq_raw)$read_depth <- sample_sums(physeq_raw)
hist(sample_data(physeq_raw)$read_depth)
mean(sample_data(physeq_raw)$read_depth) # mean depth of 21531
median(sample_data(physeq_raw)$read_depth)
nrow(sample_data(physeq_raw)[sample_data(physeq_raw)$read_depth < 3e3,]) # 74 samplse below 3k 

#### SRS normalizing - 3k ####
physeq_otus <- t(otu_table(physeq_raw)) # asvs as rows, samples as columns 
normalized_asvs <- SRS(data.frame(physeq_otus), Cmin = 3e3, set_seed = TRUE, seed = 1999) 
# 74 sample(s) discarded due to low number of counts (number of counts < Cmin):  X16_16S_S16, X2020_MB_125, X2020_MB_171, X2020_MB_174, X2020_MB_197, X2020_MB_23, X2020_MB_337, X2020_MB_338, X2020_MB_431, X2020_MB_56, X2021_MB_115, X2021_MB_120, X2021_MB_131, X2021_MB_140, X2021_MB_15, X2021_MB_151, X2021_MB_155, X2021_MB_158, X2021_MB_159, X2021_MB_174, X2021_MB_189, X2021_MB_19, X2021_MB_191, X2021_MB_196, X2021_MB_200, X2021_MB_201, X2021_MB_205, X2021_MB_22, X2021_MB_23, X2021_MB_24, X2021_MB_31, X2021_MB_32, X2021_MB_33, X2021_MB_35, X2021_MB_37, X2021_MB_38, X2021_MB_40, X2021_MB_41, X2021_MB_42, X2021_MB_45, X2021_MB_46, X2021_MB_52, X2021_MB_54, X2021_MB_58, X2021_MB_60, X2021_MB_65, X2021_MB_71, X2021_MB_72, X2021_MB_83, X2021_MB_85, X2021_MB_87, X2021_MB_94, X2022_MB_256, X2022_MB_258, X2022_MB_282, X40_16S_S40, AMOR_EDNA_3, TR_03_S66, TR_08_S18, TR_104_S29, TR_12_S12, TR_123_S31, TR_13_S88, TR_15_S16, TR_16_S53, TR_20_S62, TR_26_S11, TR_27_S22, TR_44_S57, TR_51_S49, TR_53_S63, TR_88_S106, UPPCAMP_EDNA_1, UPPCAMP_EDNA_3
colnames(normalized_asvs) <- sub("^X", "", colnames(normalized_asvs)) # this line! need to make it to where it only gets rid of Xs at the start of the character element 
rownames(normalized_asvs) <- rownames(data.frame(physeq_otus))

# get rid of any asvs that are not present in any sample 
normalized_asvs <- normalized_asvs[rowSums(normalized_asvs) >0,] # left with 21516 ASVs

# extract ASVs and sampleNames to subset the original phyloseq object 
remaining_asvs <- rownames(normalized_asvs)

# extract remaining samples 
remaining_samples <- colnames(normalized_asvs) # 788 samples 

# subset physeq 
physeq_norm <- subset_taxa(physeq_raw, taxa_names(physeq_raw) %in% remaining_asvs) # take out asvs
physeq_norm <- subset_samples(physeq_norm, sample_names(physeq_norm) %in% remaining_samples)


physeq_subset <- phyloseq(otu_table(as.matrix(normalized_asvs), taxa_are_rows = TRUE),
                          sample_data(sample_data(physeq_norm)),
                          tax_table(tax_table(physeq_norm)),
                          refseq(physeq_norm),
                          phy_tree(physeq_norm))
saveRDS(physeq_subset, 'data/dada2/physeq_normal.RDS')

