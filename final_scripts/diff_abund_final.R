# Differential Abundance Analyses and Stacked Bar Plots of Fish vs Water

# Question - What taxa underly differences between host-associated and free-living? Are these taxa consistent across systems? 

#### load packages #### 
packages <- c("ggplot2", "gridExtra", "phyloseq", "tidyverse", "vegan", "dplyr", "ANCOMBC", "ggVennDiagram")
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

#### load in unnormalized physeq object #### 
physeq <- readRDS('data/dada2/physeq_filt_nonnormal.RDS')
physeq # full phyloseq object 


#### Stacked Bar Charts #### 
# stacked barcharts of fish vs water, stickle vs trout vs sierra vs VI 
# phylum and class levels 
#ASV_rel_abund <- transform_sample_counts(PS.fin, function(x) x/sum(x))
#pm <- psmelt(ASV_rel_abund)
# load normalized phyloseq object 
physeq_subset <- readRDS('data/dada2/physeq_normal.RDS')
#pm <- psmelt(physeq_subset)
#saveRDS(pm, 'data/dada2/pmsmelt.RDS')
pm <- readRDS('data/dada2/pmsmelt.RDS')
pm$habitat_type2 <- ifelse(pm$sample_type %in% c("trout_sierra", "stickle_VI"), "fish", "water")



# by phylum 
Abund <- aggregate(Abundance ~ sample_type + site + year + habitat_type + habitat_type2 + Phylum, data=pm, FUN=sum)
sam <- aggregate(Abundance ~ sample_type + site + year + habitat_type + habitat_type2, data=pm, FUN=sum) # this should include everything Abund has except for Class

colnames(sam)[ncol(sam)] <- "tot.abund" # renaming the last column 'tot.abund' so that we can differentiate that from rel.abund

Abund.final <- dplyr::full_join(Abund, sam)
Abund.final$rel_abund <- Abund.final$Abundance/Abund.final$tot.abund

# phyla that have rel.abund less than 5% will be named as other
Abund.final$Phylum[Abund.final$rel_abund <= 0.05] <- "Other (<5%)" 

Abund.final$Phylum<-as.factor(Abund.final$Phylum)
Abund.final$Phylum<-factor(Abund.final$Phylum, 
                           levels=c(levels(Abund.final$Phylum)))

# plot by sample type
ggplot(Abund.final, aes(x=sample_type, y=rel_abund, fill=Phylum, color=Phylum))+
  geom_bar(aes(fill=Phylum, color=Phylum), stat="identity", position="fill") +
  xlab("Sample Type") +
  ylab("Relative Abundance (%)") + 
  theme_classic() + 
  theme(axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=14), 
        axis.title.y=element_blank(),
        axis.title.x = element_blank(), 
        legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("red4", "orangered", "gold", "plum", "turquoise3", "seagreen4", "lightgreen", "gray", "royalblue3", "tan3", "lightsalmon", "purple3")) +
  scale_color_manual(values = c("red4", "orangered", "gold", "plum", "turquoise3", "seagreen4", "lightgreen", "gray", "royalblue3", "tan3", "lightsalmon", "purple3"))
ggsave("figures/final_figures/phyl_rel_abund.pdf", width = 10, height = 8)

# values of mean relative abundance for each phyla across sample types
rel_abund_values <- Abund.final %>%
  group_by(sample_type, Phylum) %>%       
  summarise(mean_rel_abund = mean(rel_abund, na.rm = TRUE)) %>% 
  arrange(sample_type, Phylum)
# make csv of values
write.csv(rel_abund_values, "data/dada2/rel_abund_phyla.csv")

# What are the most common phyla across each sample type? 
top_phyla <- rel_abund %>%
  group_by(sample_type) %>%
  arrange(desc(mean_rel_abund), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  select(sample_type, Phylum, mean_rel_abund)

# rel abund plot across time, lake, and sample type
Rel.abund.phy <- ggplot(Abund.final, aes(x = Date, y = rel_abund, fill= Class, color = Class)) +
  geom_bar(width = 0.8, aes(fill = Class, color = Class), stat ="identity", position = "fill") +
  facet_grid(cols = vars(lake), rows = vars(phy_group)) + 
  xlab("Sampling Time Point") +  
  ylab("Relative Abundance (%)") +
  theme_bw() +
  theme(axis.text.x=element_text(size=8, angle = 75, vjust = 0.5, hjust=1), axis.title.x=element_text(size=13)) +
  theme(axis.text.y=element_text(size=11), axis.title.y=element_text(size=13)) +
  theme(strip.text.x = element_text(size = 13), strip.text.y = element_text(size = 13)) +
  scale_fill_manual(values = c("red4", "orangered", "gold", "plum", "turquoise3", "seagreen4", "lightgreen", "lightskyblue", "royalblue3", "tan3", "lightsalmon")) +
  scale_color_manual(values = c("red4", "orangered", "gold", "plum", "turquoise3", "seagreen4", "lightgreen", "lightskyblue", "royalblue3", "tan3", "lightsalmon"))

### ANCOMBC ####
# add fish vs water column 
sample_data(physeq)$habitat_type2 <- as.factor(ifelse(sample_data(physeq)$sample_type %in% c("trout_sierra", "stickle_VI"), "fish", "water"))


# split into physeq by region 
physeq_sierra <- subset_samples(physeq, sample_type == "trout_sierra" | sample_type == "water_sierra")
physeq_sierra <- prune_taxa(taxa_sums(physeq_sierra) > 0, physeq_sierra)

physeq_vancouver <- subset_samples(physeq, sample_type == "stickle_VI" | sample_type == "water_VI")
physeq_vancouver <- prune_taxa(taxa_sums(physeq_vancouver) > 0, physeq_vancouver)

#### ANCOMBC Stickle vs Vancouver Water #### 
vancouver_domout = ancombc(data = physeq_vancouver, assay_name = "counts",tax_level = "Phylum", formula = c("sample_type"), struc_zero = T, group="sample_type",prv_cut = 0.1,lib_cut = 0)

# check output
str(vancouver_domout$res)

# save output for plotting
res_prim_vancouver <- vancouver_domout$res

res_prim_vanc_df <- res_prim_vancouver$lfc # take out log fold change 
res_prim_vanc_df$se <- as.vector(res_prim_vancouver$se[,3]) # add se of each taxa

# check for log fold changes to plot
res_prim_vanc_df %>% filter(sample_typewater_VI< -1 | sample_typewater_VI > 1)

# plot
res_prim_vanc_df <- res_prim_vanc_df %>% arrange(-sample_typewater_VI) %>% mutate(taxon = factor(taxon))

dif_abund_vancouver <- res_prim_vanc_df %>% filter(sample_typewater_VI < -1 | sample_typewater_VI > 1) %>%
  mutate(cat = ifelse(sample_typewater_VI > 1, "Higher in Vancouver Lakes", "Higher in Stickleback")) %>%
  ggplot(aes(y = reorder(taxon, -sample_typewater_VI), x = sample_typewater_VI, color=cat)) + geom_point() +
  geom_errorbar(aes(xmin=sample_typewater_VI- se , xmax=sample_typewater_VI+ se), width=.2,
                position=position_dodge(.9))  +
  #scale_color_manual(values=c("dodgerblue", "coral")) +
  xlab("Log2FoldChange") + ylab("Microbial Phyla") + labs(color="") +
  geom_vline(xintercept = 0, linetype="dashed",
             linewidth=0.5, color='gray') + xlim(-8, 8) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14))
ggsave('figures/final_figures/ANCOMBC_sticklevsvancouver.tiff', width = 8, height = 6)

#### ANCOMBC Trout vs Sierra Water #### 
sierra_domout = ancombc(data = physeq_sierra, assay_name = "counts",tax_level = "Phylum", formula = c("sample_type"), struc_zero = T, group="sample_type",prv_cut = 0.1,lib_cut = 0)

# check output
str(sierra_domout$res)

# save output for plotting
res_prim_sierra <- sierra_domout$res

res_prim_sierra_df <- res_prim_sierra$lfc # take out log fold change 
res_prim_sierra_df$se <- as.vector(res_prim_sierra$se[,3]) # add se of each taxa

# check for log fold changes to plot
res_prim_sierra_df %>% filter(sample_typewater_sierra < -1 | sample_typewater_sierra > 1)

# plot
res_prim_sierra_df <- res_prim_sierra_df %>% arrange(-sample_typewater_sierra) %>% mutate(taxon = factor(taxon))

dif_abund_sierra <- res_prim_sierra_df %>% filter(sample_typewater_sierra< -1 | sample_typewater_sierra > 1) %>%
  mutate(cat = ifelse(sample_typewater_sierra > 1, "Higher in sierra Lakes", "Higher in Stickleback")) %>%
  ggplot(aes(y = reorder(taxon, -sample_typewater_sierra), x = sample_typewater_sierra, color=cat)) + geom_point() +
  geom_errorbar(aes(xmin=sample_typewater_sierra- se , xmax=sample_typewater_sierra+ se), width=.2,
                position=position_dodge(.9))  +
  scale_color_manual(values=c("dodgerblue", "coral")) +
  xlab("Log2FoldChange") + ylab("Microbial Phyla") + labs(color="") +
  geom_vline(xintercept = 0, linetype="dashed",
             linewidth=0.5, color='gray') + xlim(-8, 8) + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14))
ggsave('figures/final_figures/ANCOMBC_troutvssierra.tiff', width = 8, height = 6)


#### combining sierra and vancouver island ancombc ####
str(res_prim_sierra_df)
str(res_prim_vanc_df)

res_prim_sierra_df$se_sierra <- res_prim_sierra_df$se
res_prim_vanc_df$se_VI <- res_prim_vanc_df$se 

dif_abund_comb <- merge(res_prim_sierra_df[,c(1,3,5)],res_prim_vanc_df[,c(1,3,5)], by = "taxon") # wide data frame - lets make it long 

# make long 
dif_abund_long <- pivot_longer(
  dif_abund_comb,
  cols = c(sample_typewater_sierra, sample_typewater_VI, se_sierra, se_VI),
  names_to = c(".value", "sample_type"),
  names_pattern = "(.+)_(.+)"
)

dif_abund_long %>% #filter(sample_typewater < -1 | sample_typewater > 1) %>%
  mutate(cat = ifelse(sample_type == "sierra" & sample_typewater > 0 , "Higher in Sierra Lakes", ifelse(
    sample_type == "VI" & sample_typewater > 0 , "Higher in Vancouver Lakes", ifelse(
      sample_type == "VI" & sample_typewater < 0, "Higher in Stickleback", "Higher in Brook Trout"
    )))) %>%
  ggplot(aes(y = reorder(taxon, -sample_typewater), x = sample_typewater, color=cat)) + geom_point(size = 2) +
  geom_errorbar(aes(xmin=sample_typewater- se , xmax=sample_typewater+ se), width=.3)  +
  scale_color_manual(values=c("dodgerblue", "dodgerblue", "#e16f00", "#e16f00")) +
  xlab("Log2FoldChange") + ylab("Microbial Phyla") + labs(color="") +
  geom_vline(xintercept = 0, linetype="dashed",
             linewidth=0.65, color='gray') + xlim(-8, 8) + 
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 10),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14))
ggsave('figures/final_figures/dif_abund_final.pdf', width = 7, height = 6)

# list of phyla: 
# [1] "Acidobacteriota"   "Actinobacteriota"  "Armatimonadota"    "Bacteroidota"     
# [5] "Bdellovibrionota"  "Chloroflexi"       "Crenarchaeota"     "Cyanobacteria"    
# [9] "Deinococcota"      "Dependentiae"      "Desulfobacterota"  "Elusimicrobiota"  
# [13] "Euryarchaeota"     "Firmicutes"        "Fusobacteriota"    "Gemmatimonadota"  
# [17] "Halobacterota"     "MBNT15"            "Myxococcota"       "Nanoarchaeota"    
# [21] "NB1-j"             "Patescibacteria"   "Planctomycetota"   "Proteobacteria"   
# [25] "Verrucomicrobiota" "WPS-2" 

#### Ven Diagram of ASVs Across Sample Types ####
physeq_normal <- readRDS('data/dada2/physeq_normal.RDS')

# extracting asvs present in each sample type 
ps_list <- split(sample_names(physeq_normal), sample_data(physeq_normal)$sample_type)
asv_sets <- lapply(ps_list, function(samples) {
  ps_sub <- prune_samples(samples, physeq_normal)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  taxa_names(ps_sub)
})

ggVennDiagram(asv_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "coral") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
ggsave('figures/final_figures/VENN_DIAG.pdf')

# pairwise shared ASVs between sample types 
all_asvs <- taxa_names(prune_taxa(taxa_sums(physeq_normal) > 0, physeq_normal))
total_asvs <- length(all_asvs)

# Pairwise percent shared ASVs matrix (relative to total ASVs)
sample_types <- names(asv_sets)
shared_total_matrix <- outer(sample_types, sample_types, Vectorize(function(i, j) {
  shared <- length(intersect(asv_sets[[i]], asv_sets[[j]]))
  round(100 * shared / total_asvs, 1)
}))

# Add row and column names
dimnames(shared_total_matrix) <- list(sample_types, sample_types)

# View the result
shared_total_matrix
