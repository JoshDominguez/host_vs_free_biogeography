### Ordinations of Fish/Water Bacterial Metacommunities in Sierras and Vancouver Island ####

#### load packages #### 
packages <- c("ggplot2", "gridExtra", "knitr","phyloseq", "tidyverse", "vegan", "MuMIn", "dplyr", "ecodist", "glmmTMB", "glue", "ggforce", "geosphere", "lmPerm", "otuSummary", "ANCOMBC", "DECIPHER", "phangorn", "ape", "msa","gdm", "GUniFrac", "lme4", 'ggpubr', 'SRS' ,'car', 'reshape2', 'ggsignif', 'phytools', 'picante', 'RColorBrewer', 'VennDiagram', 'ggvenn')
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

#### Loading Normalized Phyloseq Object #### 
physeq_subset <- readRDS('data/dada2/physeq_normal.RDS')

# make habitat type column 
sample_data(physeq_subset)$type <- ifelse(sample_data(physeq_subset)$sample_type %in% c("water_sierra", "water_VI"), "Water", "Fish")

sample_data(physeq_subset)$system <- ifelse(sample_data(physeq_subset)$sample_type %in% c("stickle_VI", "water_VI"), "VI", "SN")

# split into separate phyloseq object 
physeq_stickle <- subset_samples(physeq_subset, sample_type == "stickle_VI")
physeq_stickle <- prune_taxa(taxa_sums(physeq_stickle) > 0, physeq_stickle)

physeq_trout <- subset_samples(physeq_subset, sample_type == "trout_sierra")
physeq_trout <- prune_taxa(taxa_sums(physeq_trout) > 0, physeq_trout) 

physeq_waterVI <- subset_samples(physeq_subset, sample_type == "water_VI")
physeq_waterVI <- prune_taxa(taxa_sums(physeq_waterVI) > 0, physeq_waterVI)

physeq_watersierra <- subset_samples(physeq_subset, sample_type == "water_sierra")
physeq_watersierra <- prune_taxa(taxa_sums(physeq_watersierra) > 0, physeq_watersierra) 

#### Unweighted Unifrac Ordination #### 

# unw unifrac pcoa 
#uni_pcoa_full <- ordinate(physeq_subset, method = "PCoA", distance = "unifrac")
#saveRDS(uni_pcoa_full, 'data/dada2/uwuni_pcoa_full.RDS')
uni_pcoa_full <- readRDS('data/dada2/uwuni_pcoa_full.RDS')

pcoa_scores <- as.data.frame(uni_pcoa_full$vectors[, c("Axis.1", "Axis.2")])
uni_pcoa_full_df <- cbind(pcoa_scores, sample_data(physeq_subset))
round(uni_pcoa_full$values[1]*100/sum(uni_pcoa_full$values[1]),2)[1:2,]
uni_pcoa_full_df$type <- ifelse(uni_pcoa_full_df$sample_type %in% c("water_sierra", "water_VI"), "Water", "Fish")
uni_pcoa_full_df$system <- ifelse(uni_pcoa_full_df$sample_type %in% c("stickle_VI", "water_VI"), "VI", "SN")


# plot unweighted uni pcoa 
ggplot(uni_pcoa_full_df, aes(x=Axis.1, y=Axis.2, fill = factor(sample_type), 
                             color = factor(system), shape = factor(type))) +
  geom_vline(xintercept = c(0), color = "grey70",linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70",linetype = 2)+
  geom_point(size = 2.2) + 
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#e16f00","dodgerblue","dodgerblue","#e16f00"), guide = "none") +
  scale_shape_manual(values = c(21, 24)) +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  stat_ellipse(show.legend = FALSE, geom = "polygon", alpha = 0.2) +
  xlab("PCoA1 (11.7%)") + 
  ylab("PCoA2 (7.8%)") +
  theme_classic() + 
  theme(legend.title = element_blank(), legend.position = "bottom") # alot more phylo variation in host-associated 
ggsave("figures/final_figures/unw_PCOA.tiff", width = 8, height = 6)

#### WEIGHTED UNIFRAC PCOA #### 
set.seed(1999)
#wuni_pcoa_full <- ordinate(physeq_subset, method = "PCoA", distance = "wunifrac")
#saveRDS(wuni_pcoa_full, 'data/dada2/wuni_pcoa_full.RDS')

wuni_pcoa_full <- readRDS('data/dada2/wuni_pcoa_full.RDS')
wuni_pcoa_scores <- as.data.frame(wuni_pcoa_full$vectors[, c("Axis.1", "Axis.2")])
wuni_pcoa_full_df <- cbind(wuni_pcoa_scores, sample_data(physeq_subset))
round(wuni_pcoa_full$values[1]*100/sum(wuni_pcoa_full$values[1]),2)[1:2,]

# plot weighted uni pcoa 
ggplot(wuni_pcoa_full_df, aes(x=Axis.1, y=Axis.2, fill = factor(sample_type), 
                              color = factor(system), shape = factor(type))) +
  geom_vline(xintercept = c(0), color = "grey70",linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70",linetype = 2)+
  geom_point(size = 2.5) + 
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#e16f00","dodgerblue","dodgerblue","#e16f00"), guide = "none") +
  scale_shape_manual(values = c(24, 21)) +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  stat_ellipse(show.legend = FALSE, geom = "polygon", alpha = 0.2) +
  xlab("PCoA1 (18.2%)") + 
  ylab("PCoA2 (12.3%)") +
  theme_classic() + 
  theme(legend.title = element_blank(), 
        axis.title = element_text(size = 20), 
        legend.text = element_text(size = 15), 
        axis.text = element_text(size = 14))
ggsave("figures/final_figures/wun_PCOA.tiff", width = 8, height = 6)

#### Bray Curtis PcOA by Sample Types ####
set.seed(1999)
#bray_dist_full <- vegdist(decostand(t(otu_table(physeq_subset)), method = "hellinger"), method = 'bray')
#saveRDS(bray_dist_full, "data/dada2/bc_dist_full_norm.RDS")
bray_dist_full <- readRDS("data/dada2/bc_dist_full_norm.RDS")

pcoa_full <- cmdscale(bray_dist_full, k = 3, eig = TRUE)
pcoa_coords <- pcoa_full$points
pcoa_coords <- merge(sample_data(physeq_subset), pcoa_coords, by = 'row.names', all = TRUE)
round(pcoa_full$eig*100/sum(pcoa_full$eig),2)[1:2]

# new columns for plot 
pcoa_coords$type <- ifelse(pcoa_coords$sample_type %in% c("water_sierra", "water_VI"), "Water", "Fish")
pcoa_coords$system <- ifelse(pcoa_coords$sample_type %in% c("stickle_VI", "water_VI"), "VI", "SN")

# plot 
ggplot(pcoa_coords, aes(x=V1, y=V2, fill = factor(sample_type), 
                        color = factor(system), shape = factor(type))) +
  geom_vline(xintercept = c(0), color = "grey70",linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70",linetype = 2)+
  geom_point(size = 2.5) + 
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#e16f00","dodgerblue","dodgerblue","#e16f00"), guide = "none") +
  scale_shape_manual(values = c(24,21)) +
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  stat_ellipse(show.legend = FALSE, geom = "polygon", alpha = 0.2) +
  xlab("PCoA1 (9.3%)") + 
  ylab("PCoA2 (7.4%)") +
  theme_classic() + 
  theme(legend.title = element_blank(), 
        axis.title = element_text(size = 20), 
        legend.text = element_text(size = 15), 
        axis.text = element_text(size = 14))
ggsave("figures/final_figures/braycurt_pcoa_sampletype.tiff", width = 8, height = 6)

#### Full dataset PERMANOVA for Bray Curtis #### 
PERM_sampletype <- adonis2(bray_dist_full ~ pcoa_coords$sample_type, permutations = 999)
PERM_sampletype # Sample Type: R2 of 0.17, and P = 0.001

PERM_year <- adonis2(bray_dist_full ~ pcoa_coords$year, permutations = 999)
PERM_year # R2 of 0.018, P = 0.001

PERM_system <- adonis2(bray_dist_full ~ pcoa_coords$system, permutations = 999)
PERM_system # significant, R2 = 0.054, P = 0.001

PERM_type_fish_water <- adonis2(bray_dist_full ~ pcoa_coords$type, permutations = 999)
PERM_type_fish_water # fish type: R2 = 0.088, p = 0.001

PERM_site <- adonis2(bray_dist_full ~ pcoa_coords$site, permutations = 999)
PERM_site # site - R2 of 0.17, P = 0.001

# MAIN PERMANOVA - Does habitat type (fish vs water) or system explain more variation?
PERM_type_sys <- adonis2(bray_dist_full ~ type*system, data = pcoa_coords, permutations = 999, by = "terms")
PERM_type_sys

# TYPE: F1,784 = 82.7, R2 = 0.088, SYSTEM: F1,784 = 51.6, R2 = 0.055 TYPE*SYSTEM: F1,784 = 8.73 R2 = 0.025 ALL P = 0.001

#### Full dataset PERMANOVA FOR WEIGHTED UNIFRAC #### 
# calculating distance: 
#wuni_dist <- phyloseq::distance(physeq_subset, method = "wunifrac")
#saveRDS(wuni_dist, "data/dada2/wuni_dist_mat.RDS")
wuni_dist <- readRDS("data/dada2/wuni_dist_mat.RDS")

sampletype_wuni <- adonis2(wuni_dist ~ wuni_pcoa_full_df$sample_type, permutations = 999)
sampletype_wuni # Sample Type (water VI, water sierra, stickle, trout), R2 = 0.22, P = 0.001

year_wuni <- adonis2(wuni_dist ~ wuni_pcoa_full_df$year, permutations = 999)
year_wuni # R2 of 0.023, P = 0.001

system_wuni <- adonis2(wuni_dist ~ wuni_pcoa_full_df$system, permutations = 999)
system_wuni # significant, R2 = 0.051, P = 0.001

wuni_type_fish_water <- adonis2(wuni_dist ~ wuni_pcoa_full_df$type, permutations = 999)
wuni_type_fish_water # fish type: R2 = 0.14, p = 0.001

wuni_site <- adonis2(wuni_dist ~ wuni_pcoa_full_df$site, permutations = 999)
wuni_site # site - R2 of 0.155, P = 0.001

wuni_stickle_trout_water <- adonis2(wuni_dist ~ wuni_pcoa_full_df$habitat_type)
wuni_stickle_trout_water # habitat type - R2 of 0.15, p = 0.001

# type * system
PERM_wuni_hab_system <- adonis2(wuni_dist ~ type*system, data = wuni_pcoa_full_df, permutations = 999, by = "terms")
PERM_wuni_hab_system

# TYPE: F1,784 = 144, R2 = 0.14, SYSTEM: F1,784 = 50.7, R2 = 0.05, TYPE*SYSTEM: F1,784 = 24.2 R2 = 0.024 ALL P = 0.001

#### Full data PCOA PERMANOVA FOR UNWEIGHTED UNIFRAC #### 
uni_dist <- phyloseq::distance(physeq_subset, method = "unifrac")
saveRDS(uni_dist, "data/dada2/uni_dist_mat.RDS")
uni_dist <- readRDS("data/dada2/uni_dist_mat.RDS")

# habitat and system
PERM_uwuni_hab_system <- adonis2(uni_dist ~ type*system, data = uni_pcoa_full_df, permutations = 999)
PERM_uwuni_hab_system
# TYPE R2 = 0.082, SYSTEM R2 = 0.035, TYPE*SYSTEM R2 = 0.018 ALL P = 0.001




#### System Specific PERMANOVAs - Bray Curtis #### 
# separate systems 
sierra_physeq <- subset_samples(physeq_subset, sample_type %in% c("water_sierra", "trout_sierra"))
sierra_physeq <- prune_taxa(taxa_sums(sierra_physeq) > 0, sierra_physeq)

vancouver_physeq <- subset_samples(physeq_subset, sample_type %in% c("water_VI", "stickle_VI"))
vancouver_physeq <- prune_taxa(taxa_sums(vancouver_physeq) > 0, vancouver_physeq)


# subset bray curtis dist 
sierra_bc <- as.matrix(bray_dist_full)[rownames(as.matrix(bray_dist_full)) %in% sample_names(sierra_physeq), colnames(as.matrix(bray_dist_full)) %in% sample_names(sierra_physeq)]

vancouver_bc <- as.matrix(bray_dist_full)[rownames(as.matrix(bray_dist_full)) %in% sample_names(vancouver_physeq), colnames(as.matrix(bray_dist_full)) %in% sample_names(vancouver_physeq)]

# take out metadata to make cleaner 
sierra_metadata <- data.frame(sample_data(sierra_physeq))
vancouver_metadata <- data.frame(sample_data(vancouver_physeq))

# fix year column 
sierra_metadata$year <- as.factor(sierra_metadata$year)
vancouver_metadata$year <- as.factor(vancouver_metadata$year)

set.seed(1999)
# Sierra PERMANOVA 
permsierra_type <- adonis2(sierra_bc ~ sierra_metadata$type, permutations = 999)
permsierra_type
# adonis2(formula = sierra_bc ~ sample_data(sierra_physeq)$type, permutations = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# sample_data(sierra_physeq)$sample_type   1   18.373 0.20878 60.164  0.001 ***
#   Residual                               228   69.626 0.79122                  
# Total                                  229   87.998 1.00000     


permsierra_lake <- adonis2(sierra_bc ~ sierra_metadata$site, permutations = 999)
permsierra_lake 
# adonis2(formula = sierra_bc ~ sample_data(sierra_physeq)$site, permutations = 999)
# Df SumOfSqs      R2     F Pr(>F)    
# sample_data(sierra_physeq)$site  14   15.049 0.17101 3.168  0.001 ***
#   Residual                        215   72.950 0.82899                 
# Total                           229   87.998 1.00000     

perm_sierra_full <- adonis2(sierra_bc ~ sierra_metadata$type*sierra_metadata$site*sierra_metadata$year, by = "terms", permutations = 999)
perm_sierra_full

# Effects: 
# Type: R2 = 0.21, f1,141 = 82.1, p = 0.001
# Site/lake: R2 = 0.16, f14,141 = 4.53, p = 0.001
# year: R2 = 0.033, f4,141 = 3.30, p = 0.001
# type*site: R2 = 0.10, f14,141 = 2.84, p = 0.001
# type*year: R2 = 0.007, f1,141 = 2.81, p = 0.006 
# all else was not significant 

# Vancouver PERMANOVA
perm_vancouver_year <- adonis2(vancouver_bc ~ vancouver_metadata$type * vancouver_metadata$year, by = "terms", permutations = 999)
perm_vancouver_year
# Effects: 
# Year: R2 = 0.03, F2,552 = 9.8, p = 0.001
# type*year R2 = 0.01, f2,552 = 3.85, p = 0.001

perm_vancouver_full <- adonis2(vancouver_bc ~ vancouver_metadata$type * vancouver_metadata$site, by = "terms", permutations = 999)
perm_vancouver_full # year is highly correlated with site, site nested in vancouver just to simplify not looking at year effect, independently looked at affect of year and reported  

# Effects: 
# Type: R2 = 0.08780, f1,520 = 60.9363, p = 0.001
# Site/lake: R2 = 0.10906, f18,520 = 4.2050, p = 0.001
# type*site: R2 = 0.05388, f18,520 = 2.0773, p = 0.001


#### System Specific PERMANOVA - Weighted Unifrac #### 

# subset weighted unifrac distance matrix
sierra_wuni <- as.matrix(wuni_dist)[rownames(as.matrix(wuni_dist)) %in% sample_names(sierra_physeq), colnames(as.matrix(wuni_dist)) %in% sample_names(sierra_physeq)]

vancouver_wuni <- as.matrix(wuni_dist)[rownames(as.matrix(wuni_dist)) %in% sample_names(vancouver_physeq), colnames(as.matrix(wuni_dist)) %in% sample_names(vancouver_physeq)]

set.seed(1999)
# Fit PERMANOVA for Sierra 
perm_sierra_wuni <- adonis2(sierra_wuni ~ sierra_metadata$type*sierra_metadata$site*sierra_metadata$year, by = "terms", permutations = 999)
perm_sierra_wuni

# Effects: 
# Type: R2 = 0.39, f1,141 = 184, p = 0.001
# Site/lake: R2 = 0.10278, f14,141 =  3.4528, p = 0.001
# year: R2 = 0.02701, f4,141 = 3.1759, p = 0.001
# type*site: R2 = 0.06616, f14,141 = 2.2229, p = 0.001
# type*year: R2 = 0.00951, f1,141 = 4.4733, p = 0.004 
# all else was not significant 

# Fit PERMANOVA for Vancouver Island

# type and year 
perm_vancouver_yearwuni <- adonis2(vancouver_wuni ~ vancouver_metadata$type * vancouver_metadata$year, permutations = 999, by = "terms")
perm_vancouver_yearwuni
# Effects:
# Year: R2 = 0.03, F2,552 = 9.8,
# Year*Type: R2 = 0.014, F = 4.5, p = 0.001

# type and site
perm_vancouver_wuni <- adonis2(vancouver_wuni ~ vancouver_metadata$type * vancouver_metadata$site, permutations = 999, by = "terms")
perm_vancouver_wuni

# Effects: 
# Type: R2 = 0.10523, f1,520 = 74.0807, p = 0.001
# Site/lake: R2 = 0.10622, f18,520 =  4.15, p = 0.001
# type*site: R2 = 0.04989, f18,520 =  1.95, p = 0.001
