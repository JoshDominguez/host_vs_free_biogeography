#### Distance Decay of Similarity for Fish Gut and Water Microbiomes Across Systems ####

#### Load Packages #### 
packages <- c("ggplot2", "gridExtra", "knitr","phyloseq", "tidyverse", "vegan", "MuMIn", "dplyr", "ecodist", "glmmTMB", "glue", "ggforce", "geosphere", "lmPerm", "otuSummary", "ANCOMBC", "DECIPHER", "phangorn", "ape", "msa","gdm", "GUniFrac", "lme4", 'ggpubr', 'SRS' ,'car', 'reshape2', 'ggsignif', 'phytools', 'picante')
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

# ADD SAME YEAR COLUMN - SPATIAL VS TEMPORAL TURNOVER !!!!!
#### Load in Normalized Phyloseq Object and Split Up #### 
physeq_meta <- readRDS('data/dada2/physeq_normal.RDS')

# fix May site name 
sample_data(physeq_meta)$site <- ifelse(sample_data(physeq_meta)$site == "May ", "May", sample_data(physeq_meta)$site)

# input latitude and longitude coordinates in metadata 
site_data <- read.csv('data/site_coords.csv')

# fix site names in site_data 
site_data$site[9] <- 'May'
site_data$site[10] <- 'Spillway'
site_data$site[14] <- 'Wasco'

# matching up and inputting latitudes
sample_data(physeq_meta)$latitude <- NA 

# add latitudes in loop 
for (site_name in unique(sample_data(physeq_meta)$site)) {
  match_index <- match(site_name, site_data$site)
  if (!is.na(match_index)) {
    sample_data(physeq_meta)$latitude[sample_data(physeq_meta)$site == site_name] <- site_data$latitude[match_index]
  }
}

# inputting longitudes 
sample_data(physeq_meta)$longitude <- NA

# add longitudes in loop 
for (site_name in unique(sample_data(physeq_meta)$site)) {
  match_index <- match(site_name, site_data$site)
  if (!is.na(match_index)) {
    sample_data(physeq_meta)$longitude[sample_data(physeq_meta)$site == site_name] <- site_data$longitude[match_index]
  }
}

# splitting up into separate phyloseq objects 
physeq_trout <- subset_samples(physeq_meta, sample_type == "trout_sierra")
physeq_trout <- prune_taxa(taxa_sums(physeq_trout) > 0, physeq_trout) # 4237 taxa, 169 samples 

physeq_sierra <- subset_samples(physeq_meta, sample_type == "water_sierra")
physeq_sierra <- prune_taxa(taxa_sums(physeq_sierra) > 0, physeq_sierra) # 3407 taxa, 61 samples 

physeq_stickle <- subset_samples(physeq_meta, sample_type == "stickle_VI")
physeq_stickle <- prune_taxa(taxa_sums(physeq_stickle) > 0, physeq_stickle) # 16109 taxa, 511 samples 

physeq_VI <- subset_samples(physeq_meta, sample_type == "water_VI")
physeq_VI <- prune_taxa(taxa_sums(physeq_VI) > 0, physeq_VI) # 4182 taxa and 67 samples 

#### Calculating Geographic Distance Matrices #### 
# calculate geographic distances 

# trout geo dist 
trout_geo_mat <- distm(cbind(sample_data(physeq_trout)$longitude, sample_data(physeq_trout)$latitude), fun = distGeo)
# alter row and column names to include year, site, and sample name 
colnames(trout_geo_mat) <- paste(rownames(sample_data(physeq_trout)), sample_data(physeq_trout)$site, sample_data(physeq_trout)$year, sep = ";") 
rownames(trout_geo_mat) <- paste(rownames(sample_data(physeq_trout)), sample_data(physeq_trout)$site, sample_data(physeq_trout)$year, sep = ";") 

# sierra nevada water geo dist 
sierra_geo_mat <- distm(cbind(sample_data(physeq_sierra)$longitude, sample_data(physeq_sierra)$latitude), fun = distGeo)
colnames(sierra_geo_mat) <- paste(rownames(sample_data(physeq_sierra)), sample_data(physeq_sierra)$site, sample_data(physeq_sierra)$year, sep = ";") 
rownames(sierra_geo_mat) <- paste(rownames(sample_data(physeq_sierra)), sample_data(physeq_sierra)$site, sample_data(physeq_sierra)$year, sep = ";") 

# vancouver island geo dist 
vi_geo_mat<- distm(cbind(sample_data(physeq_VI)$longitude, sample_data(physeq_VI)$latitude), fun = distGeo)
colnames(vi_geo_mat) <- paste(rownames(sample_data(physeq_VI)), sample_data(physeq_VI)$site, sample_data(physeq_VI)$year, sep = ";") 
rownames(vi_geo_mat) <- paste(rownames(sample_data(physeq_VI)), sample_data(physeq_VI)$site, sample_data(physeq_VI)$year, sep = ";") 

# stickleback geo dist 
stickle_geo_mat<- distm(cbind(sample_data(physeq_stickle)$longitude, sample_data(physeq_stickle)$latitude), fun = distGeo)
rownames(stickle_geo_mat) <- paste(rownames(sample_data(physeq_stickle)), sample_data(physeq_stickle)$site, sample_data(physeq_stickle)$year, sep = ";") 
colnames(stickle_geo_mat) <- paste(rownames(sample_data(physeq_stickle)), sample_data(physeq_stickle)$site, sample_data(physeq_stickle)$year, sep = ";") 



# convert dataframes to km  
trout_geo_mat <- trout_geo_mat/1e3
sierra_geo_mat <- sierra_geo_mat/1e3
vi_geo_mat <- vi_geo_mat/1e3
stickle_geo_mat <- stickle_geo_mat/1e3

# lengthen dataframes 
trout_geo_dist <- matrixConvert(trout_geo_mat, colname = c("site1", "site2", "dist_km"))
sierra_geo_dist <- matrixConvert(sierra_geo_mat, colname = c("site1", "site2", "dist_km"))
vi_geo_dist <- matrixConvert(vi_geo_mat, colname = c("site1", "site2", "dist_km"))
stickle_geo_dist <- matrixConvert(stickle_geo_mat, colname = c("site1", "site2", "dist_km"))

#### Distance Decay - Bray Curtis Similarity  #### 
# calculating distance matrices - bray curtis dissimilarity 
set.seed(1999)
trout_bc <- vegdist(decostand(t(otu_table(physeq_trout)), method = "hellinger"), method = 'bray') 
trout_bc_sim <- apply(as.matrix(trout_bc), 2, function(x) 1 - x)

set.seed(1999)
sierra_bc <- vegdist(decostand(t(otu_table(physeq_sierra)), method = "hellinger"), method = 'bray') 
sierra_bc_sim <- apply(as.matrix(sierra_bc), 2, function(x) 1 - x)

set.seed(1999)
vanc_bc <- vegdist(decostand(t(otu_table(physeq_VI)), method = "hellinger"), method = 'bray')
vanc_bc_sim <- apply(as.matrix(vanc_bc), 2, function(x) 1 - x)

set.seed(1999)
stickle_bc <- vegdist(decostand(t(otu_table(physeq_stickle)), method = "hellinger"), method = 'bray')
stickle_bc_sim <- apply(as.matrix(stickle_bc), 2, function(x) 1 - x)

# lengthen dataframes and merge with geographic distance dataframes 
trout_long <- matrixConvert(trout_bc_sim, colname = c("site1", "site2", "bc_sim"))
sierra_long <- matrixConvert(sierra_bc_sim, colname = c("site1", "site2", "bc_sim"))
vanc_long <- matrixConvert(vanc_bc_sim, colname = c("site1", "site2", "bc_sim")) 
stickle_long <- matrixConvert(stickle_bc_sim, colname = c("site1", "site2", "bc_sim"))

# combine with geographic dist df 
trout_sim_df <- data.frame(cbind(trout_long, trout_geo_dist))
sierra_sim_df <- data.frame(cbind(sierra_long, sierra_geo_dist))
vanc_sim_df <- data.frame(cbind(vanc_long, vi_geo_dist))
stickle_sim_df <- data.frame(cbind(stickle_long, stickle_geo_dist))

# separate out sample names, lake names, and years 
# trout
trout_sim_df <- trout_sim_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)
trout_sim_df <- trout_sim_df[,c(1:3, 5:7,9:12)] # clean up df to include only cols we want 
trout_sim_df$samelake <- as.factor(ifelse(trout_sim_df$dist_km == 0, "Y", "N"))
trout_sim_df$sameyear <- as.factor(ifelse(trout_sim_df$year1 == trout_sim_df$year2, "Y", "N"))

# sierra lakes 
sierra_sim_df <- sierra_sim_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)
sierra_sim_df <- sierra_sim_df[,c(1:3, 5:7,9:12)] # clean up df to include only cols we want 
sierra_sim_df$samelake <- as.factor(ifelse(sierra_sim_df$dist_km == 0, "Y", "N"))
sierra_sim_df$sameyear <- as.factor(ifelse(sierra_sim_df$year1 == sierra_sim_df$year2, "Y", "N"))

# vancouver lakes 
vanc_sim_df <- vanc_sim_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)
vanc_sim_df <- vanc_sim_df[,c(1:3, 5:7,9:12)] # clean up df to include only cols we want 
vanc_sim_df$samelake <- as.factor(ifelse(vanc_sim_df$dist_km == 0, "Y", "N"))
vanc_sim_df$sameyear <- as.factor(ifelse(vanc_sim_df$year1 == vanc_sim_df$year2, "Y", "N"))

# stickleback 
stickle_sim_df <- stickle_sim_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)
stickle_sim_df <- stickle_sim_df[,c(1:3, 5:7,9:12)] # clean up df to include only cols we want 
stickle_sim_df$samelake <- as.factor(ifelse(stickle_sim_df$dist_km == 0, "Y", "N"))
stickle_sim_df$sameyear <- as.factor(ifelse(stickle_sim_df$year1 == stickle_sim_df$year2, "Y", "N"))

#### Eliminating Comparisons from Same Lake #### 
stickle_sim_nosame <- stickle_sim_df[stickle_sim_df$samelake == "N",]
trout_sim_nosame <- trout_sim_df[trout_sim_df$samelake == "N",]
sierra_sim_nosame <- sierra_sim_df[sierra_sim_df$samelake == "N",]
vanc_sim_nosame <- vanc_sim_df[vanc_sim_df$samelake == "N",]

#### custom mantel test to work on vectors and non-square matrices #### 
custom_mantel_test <- function(x, y, method = "spearman", nperm = 999) {
  obs_cor <- cor(x, y, method = method, use = "complete.obs") # mantel r - observed spearman correlation
  permuted_cor <- replicate(nperm, {
    cor(x, sample(y), method = method, use = "complete.obs")
  }) # permute correlations
  
  p_val <- (sum(abs(permuted_cor) >= abs(obs_cor)) + 1) / (nperm + 1)
  
  list(
    mantel_r = obs_cor,
    p_value = p_val,
    permuted_correlation_mean = mean(permuted_cor),
    permuted_correlation_sd = sd(permuted_cor),
    permuted_correlations = permuted_cor
  )
}

#### trout dist decay with no samelake ####
summary(lm(bc_sim ~ dist_km, data = trout_sim_nosame))
# LM Results 
# Intercept:0.2306553  
# Coef Est: -0.0024415 
# Std error: 
# Mult R2:  0.03172

# fit with sameyear as a random effect 
fit_trout_2 <- glmmTMB(bc_sim ~ (1|sameyear) + dist_km, data = trout_sim_nosame)
summary(fit_trout_2) # model results: Coef Est: -0.0024412  (0.0001169), int: 0.2303457 (0.0057185) R2m : 0.029, R2c: 0.015,  
r.squaredGLMM(fit_trout_2) # model is 0.032m, 0.036c, year does not matter much  

# use coefficient and intercept to get line from mixed effects model
trout_sim_nosame$predicted <- -0.0024412*trout_sim_nosame$dist_km + 0.2303457

# mantel test - tests for significance - excluding samelake comparisons 


trout_mantel <- custom_mantel_test(trout_sim_nosame$dist_km, trout_sim_nosame$bc_sim)
trout_mantel$mantel_r # -0.16
trout_mantel$p_value # p = 0.001

ggplot(trout_sim_nosame, aes(dist_km, bc_sim)) +
  geom_point(size = 2.4, color = "dodgerblue", fill = 'dodgerblue', shape = 24, alpha = 0.3) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  ylim(0, 0.85) + 
  labs(x = "Distance (km)", y= "Community Similarity (1 - Bray-Curtis)") +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14)) # max is 0.66
ggsave("figures/final_figures/fig3_troutdistdecay.pdf",width = 6, height = 6)

#### sierras no same lake #### 
# Fitting model 
summary(lm(bc_sim ~ dist_km, data = sierra_sim_nosame))
# LM Results 
# Intercept:  0.3400760 (0.0042453)
# Coef Est: -0.0038069 (0.0002203), higher than others 
# Mult R2: 0.1327

# fit with samelake and sameyear as random effects 
fit_sierra_2 <- glmmTMB(bc_sim ~ (1|sameyear) + dist_km, data = sierra_sim_nosame)
summary(fit_sierra_2)
# LMM Results: 
# coef: -0.0037730  0.0002309 std error!
# intercept: 0.3501845  0.0136435

r.squaredGLMM(fit_sierra_2) 
# 0.1290567 R2c: 0.1655312 

# mantel test - tests for significance 
sierra_mantel <- custom_mantel_test(sierra_sim_nosame$dist_km, sierra_sim_nosame$bc_sim)
sierra_mantel$mantel_r # -0.37
sierra_mantel$p_value # p = 0.001

# get predicted values from LMM: 
sierra_sim_nosame$predicted <- -0.0037730*sierra_sim_nosame$dist_km + 0.3501845

ggplot(sierra_sim_nosame, aes(dist_km, bc_sim)) +
  geom_point(size = 2.4, color = "dodgerblue", alpha = 0.3) + 
  geom_line(aes(y = predicted, x = dist_km), color = 'black') +
  scale_y_continuous(limits=c(0, 0.85)) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14)) # max is 0.611
ggsave("figures/final_figures/fig3_sierradistdecay.pdf", width = 6, height = 6)

#### vancouver island no same lake #### 
# Fitting model 
summary(lm(bc_sim ~ dist_km, data = vanc_sim_nosame))
# LM Results 
# Intercept: 0.4435304  (0.0036325)
# Coef Est: -1.412e-03  (0.0001386)
# lower rate of distance decay compared to Sierra water - less env hetero and greater within lake similarity  
# Mult R2:  0.12, better fit though 


# fit with samelake and year as random effects (would be same thing here) 
fit_vanc_2 <- glmmTMB(bc_sim ~ (1|sameyear) + dist_km, data = vanc_sim_nosame)
summary(fit_vanc_2)
# estimates: 
# Intercept: 4.069e-01   1.919e-02 
# coef: -1.335e-03 7.637e-05

r.squaredGLMM(fit_vanc_2)
# R2m of 0.10, R2c is 0.19, some driven by year 

# mantel test - tests for significance 
vanc_mantel <- custom_mantel_test(vanc_sim_nosame$dist_km, vanc_sim_nosame$bc_sim)
vanc_mantel$mantel_r # r -0.36
vanc_mantel$p_value # p = 0.001

# predicted column 
vanc_sim_nosame$predicted <- -1.335e-03*vanc_sim_nosame$dist_km + 4.069e-01

ggplot(vanc_sim_nosame, aes(dist_km, bc_sim)) +
  geom_point(size = 2.4, color = "#e16f00", alpha = 0.3) + 
  geom_line(aes(y = predicted, x = dist_km), color = 'black') +
  scale_y_continuous(limits=c(0, 0.85)) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14)) # max is 0.62
ggsave("figures/final_figures/fig3_vancouverdistdecay.pdf", width = 6, height = 6)

#### stickleback nosame lake distdecay #### 
# Fitting model 
summary(lm(bc_sim ~ dist_km, data = stickle_sim_nosame))
# LM Results 
# Intercept: 8.839e-02 
# Coef Est: -3.476e-04 
# Mult R2: 0.01156 

# fit with samelake and same year as a random effect 
fit_stickle_2 <- glmmTMB(bc_sim ~ (1|sameyear) + dist_km, data = stickle_sim_nosame)
summary(fit_stickle_2)
r.squaredGLMM(fit_stickle_2)
# dist estimate: -2.851e-04  9.449e-06, says significant but its likely not much 
# intercept: 1.052e-01 (1.575e-02 )
# R2m: 0.00695714
# Rc: 0.1126613 - year or lake effect? 

stickle_sim_nosame$predicted <- -2.851e-04*stickle_sim_nosame$dist_km + 9.154e-02 

# mantel test - tests for significance 
stickle_mantel <- custom_mantel_test(stickle_sim_nosame$dist_km, stickle_sim_nosame$bc_sim)
stickle_mantel$mantel_r # r -0.11
stickle_mantel$p_value # p = 0.001

ggplot(stickle_sim_nosame, aes(dist_km, bc_sim)) +
  geom_point(size = 2.4, color = "#e16f00",fill = "#e16f00", alpha = 0.3, shape = 24) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  labs(x = "Distance (km)", y= "Community Similarity (1 - Bray-Curtis)") +
  ylim(0,0.85) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14))
ggsave("figures/final_figures/fig3_stickledistdecay.pdf", width = 6, height = 6)
#### Weighted Unifrac Dist Decay #### 
# load weighted unifrac distance matrix 
wuni_dist_mat <- as.matrix(readRDS("data/dada2/wuni_dist_mat.RDS"))

# subset sample names of different sample types 
wuni_trout <- wuni_dist_mat[sample_names(physeq_trout), sample_names(physeq_trout)]
wuni_stickle <- wuni_dist_mat[sample_names(physeq_stickle), sample_names(physeq_stickle)]
wuni_VI <- wuni_dist_mat[sample_names(physeq_VI), sample_names(physeq_VI)]
wuni_sierra <- wuni_dist_mat[sample_names(physeq_sierra), sample_names(physeq_sierra)]

# convert to similarity 
wuni_trout_sim <- apply(as.matrix(wuni_trout), 2, function(x) 1 - x)
wuni_stickle_sim <- apply(as.matrix(wuni_stickle), 2, function(x) 1 - x)
wuni_VI_sim <- apply(as.matrix(wuni_VI), 2, function(x) 1 - x)
wuni_sierra_sim <- apply(as.matrix(wuni_sierra), 2, function(x) 1 - x)

# lengthen dataframes and merge with geographic distance dataframes 
wuni_long_trout <- matrixConvert(wuni_trout_sim, colname = c("site1", "site2", "uni_sim"))
wuni_long_sierra <- matrixConvert(wuni_sierra_sim, colname = c("site1", "site2", "uni_sim"))
wuni_long_VI <- matrixConvert(wuni_VI_sim, colname = c("site1", "site2", "uni_sim"))
wuni_long_stickle <- matrixConvert(wuni_stickle_sim, colname = c("site1", "site2", "uni_sim"))

# combine with geographic dist df 
trout_wuni_df <- data.frame(cbind(wuni_long_trout, trout_geo_dist))
sierra_wuni_df <- data.frame(cbind(wuni_long_sierra, sierra_geo_dist))
vanc_wuni_df <- data.frame(cbind(wuni_long_VI, vi_geo_dist))
stickle_wuni_df <- data.frame(cbind(wuni_long_stickle, stickle_geo_dist))

# cleaning up dataframes 
# trout 
trout_wuni_df <- trout_wuni_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)

# clean up df to include only cols we want 
trout_wuni_df <- trout_wuni_df[,c(1:3, 5:7,9:12)]

# sierra 
sierra_wuni_df <- sierra_wuni_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)

# clean up df to include only cols we want 
sierra_wuni_df <- sierra_wuni_df[,c(1:3, 5:7,9:12)]

# stickleback
stickle_wuni_df <- stickle_wuni_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)

# clean up df to include only cols we want 
stickle_wuni_df <- stickle_wuni_df[,c(1:3, 5:7,9:12)]

# vancouver island 
vanc_wuni_df <- vanc_wuni_df %>%
  mutate(across(c(site1.1, site2.1), as.character)) %>%
  separate(site1.1, into = c("sample1", "lake1", "year1"), sep = ";", remove = FALSE) %>%
  separate(site2.1, into = c("sample2", "lake2", "year2"), sep = ";", remove = FALSE)

# clean up df to include only cols we want 
vanc_wuni_df <- vanc_wuni_df[,c(1:3, 5:7,9:12)]

# removing observations from the same lake in all dataframes 
stickle_wuni_df <- stickle_wuni_df[!stickle_wuni_df$dist_km == 0,]
trout_wuni_df <- trout_wuni_df[!trout_wuni_df$dist_km == 0,]
sierra_wuni_df <- sierra_wuni_df[!sierra_wuni_df$dist_km == 0,]
vanc_wuni_df <- vanc_wuni_df[!vanc_wuni_df$dist_km == 0,]

# adding same year column 
stickle_wuni_df$sameyear <- as.factor(ifelse(stickle_wuni_df$year1 == stickle_wuni_df$year2, "Y", "N"))
trout_wuni_df$sameyear <- as.factor(ifelse(trout_wuni_df$year1 == trout_wuni_df$year2, "Y", "N"))
vanc_wuni_df$sameyear <- as.factor(ifelse(vanc_wuni_df$year1 == vanc_wuni_df$year2, "Y", "N"))
sierra_wuni_df$sameyear <- as.factor(ifelse(sierra_wuni_df$year1 == sierra_wuni_df$year2, "Y", "N"))

# ylims: .55 to 1 

#### Weighted Uni Dist Decay Trout #### 
summary(lm(uni_sim ~ dist_km, data = trout_wuni_df))
# LM Results 
# Intercept: 0.86  
# Coef Est: -3.95e-4
# Std error: 4.32e-5
# Mult R2:  0.0063

# fit with sameyear as a random effect 
fit_wuni_trout <- glmmTMB(uni_sim ~ (1|sameyear) + dist_km, data = trout_wuni_df)
summary(fit_wuni_trout) # model results: Coef Est:  -3.949e-04  (4.319e-05), int: 8.588e-01 (4.319e-05 )  
r.squaredGLMM(fit_wuni_trout) # model is 0.006289391m, 0.006289392c
# use coefficient and intercept to get line from mixed effects model
trout_wuni_df$predicted <- -3.949e-04*trout_wuni_df$dist_km + 8.588e-01


# mantel test - tests for significance 
trout_wuni_mantel <- custom_mantel_test(trout_wuni_df$dist_km, trout_wuni_df$uni_sim)
trout_wuni_mantel$mantel_r # r -0.06
trout_wuni_mantel$p_value # p = 0.001

ggplot(trout_wuni_df, aes(dist_km, uni_sim)) +
  geom_point(size = 2, color = "dodgerblue", fill = 'dodgerblue', shape = 24, alpha = 0.3) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  scale_y_continuous(limits=c(0.55, 1)) +
  labs(x = "Distance (km)", y= "Phylogenetic Similarity (1-WuniFrac)") +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14))
ggsave("figures/final_figures/fig4_wuni_trout.pdf",width = 6, height = 6)


### weighted uni distance decay sierra lakes #### 
# model with no random effects 
summary(lm(uni_sim ~ dist_km, data = sierra_wuni_df))
# Call:
#   lm(formula = uni_sim ~ dist_km, data = sierra_wuni_df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.094327 -0.024103  0.001402  0.025355  0.087206 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.8770364  0.0016783 522.591  < 2e-16 ***
#   dist_km     -0.0007213  0.0000926  -7.789 1.16e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.035 on 1726 degrees of freedom
# Multiple R-squared:  0.03396,	Adjusted R-squared:  0.0334 
# F-statistic: 60.67 on 1 and 1726 DF,  p-value: 1.155e-14

# model with year as random effect 
fit_wuni_sierra <- glmmTMB(uni_sim ~ (1|sameyear) + dist_km, data = sierra_wuni_df)
summary(fit_wuni_sierra)

# model output: 
# Family: gaussian  ( identity )
# Formula:          uni_sim ~ (1 | sameyear) + dist_km
# Data: sierra_wuni_df
# 
# AIC      BIC   logLik deviance df.resid 
# -6712.7  -6690.9   3360.4  -6720.7     1724 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance  Std.Dev.
# sameyear (Intercept) 4.645e-05 0.006815
# Residual             1.193e-03 0.034544
# Number of obs: 1728, groups:  sameyear, 2
# 
# Dispersion estimate for gaussian family (sigma^2): 0.00119 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  8.808e-01  5.131e-03  171.67  < 2e-16 ***
#   dist_km     -7.087e-04  9.143e-05   -7.75  9.1e-15 ***


r.squaredGLMM(fit_wuni_sierra)
# R2m     R2c
# [1,] 0.03241585 0.06867

sierra_wuni_df$predicted <- -7.087e-04*sierra_wuni_df$dist_km + 8.808e-01

# plot 
ggplot(sierra_wuni_df, aes(dist_km, uni_sim)) +
  geom_point(size = 2, color = "dodgerblue", fill = 'dodgerblue', alpha = 0.3) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  scale_y_continuous(limits=c(0.55, 1)) +
  labs(x = "Distance (km)", y= "Phylogenetic Similarity (1-WuniFrac)") +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))
ggsave("figures/final_figures/fig4_wuni_sierra.pdf",width = 6, height = 6)

# mantel test - tests for significance 
sierra_wuni_mantel <- custom_mantel_test(sierra_wuni_df$dist_km, sierra_wuni_df$uni_sim)
sierra_wuni_mantel$mantel_r # r -0.18
sierra_wuni_mantel$p_value # p = 0.001

#### Weighted Uni Stickleback ##### 
# mormal linear model 
summary(lm(uni_sim ~ dist_km, data = stickle_wuni_df))
# Call:
#   lm(formula = uni_sim ~ dist_km, data = stickle_wuni_df)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.192119 -0.026538  0.001224  0.028382  0.201963 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 7.943e-01  2.236e-04 3552.07   <2e-16 ***
#   dist_km     6.109e-05  7.119e-06    8.58   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04693 on 112188 degrees of freedom
# Multiple R-squared:  0.0006558,	Adjusted R-squared:  0.0006469 
# F-statistic: 73.62 on 1 and 112188 DF,  p-value: < 2.2e-16

# model with year as random effect 
fit_wuni_stickle <- glmmTMB(uni_sim ~ (1|sameyear) + dist_km, data = stickle_wuni_df)
summary(fit_wuni_stickle)

# Family: gaussian  ( identity )
# Formula:          uni_sim ~ (1 | sameyear) + dist_km
# Data: stickle_wuni_df
# 
# AIC       BIC    logLik  deviance  df.resid 
# -368714.4 -368675.9  184361.2 -368722.4    112186 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance  Std.Dev.
# sameyear (Intercept) 1.554e-05 0.003941
# Residual             2.188e-03 0.046781
# Number of obs: 112190, groups:  sameyear, 2
# 
# Dispersion estimate for gaussian family (sigma^2): 0.00219 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 7.952e-01  2.796e-03  284.38   <2e-16 ***
#   dist_km     7.936e-05  7.132e-06   11.13   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

r.squaredGLMM(fit_wuni_stickle)
# R2m        R2c
# [1,] 0.001105353 0.00814635

stickle_wuni_df$predicted <- 7.936e-05*stickle_wuni_df$dist_km + 7.952e-01

# plot 
ggplot(stickle_wuni_df, aes(dist_km, uni_sim)) +
  geom_point(size = 2.4, color = "#e16f00",fill = "#e16f00", alpha = 0.3, shape = 24) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  scale_y_continuous(limits=c(0.55, 1)) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))
ggsave("figures/final_figures/fig4_wuni_stickle.pdf",width = 6, height = 6)

# mantel test - tests for significance 
stickle_wuni_mantel <- custom_mantel_test(stickle_wuni_df$dist_km, stickle_wuni_df$uni_sim)
stickle_wuni_mantel$mantel_r # r 0.02
stickle_wuni_mantel$p_value # p = 0.001


#### Weighted Uni Vancouver Island Lakes #### 
# mormal linear model 
summary(lm(uni_sim ~ dist_km, data = vanc_wuni_df))

# Call:
#   lm(formula = uni_sim ~ dist_km, data = vanc_wuni_df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.10919 -0.02216  0.00375  0.02538  0.07455 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 8.721e-01  1.106e-03 788.821   <2e-16 ***
#   dist_km     1.336e-05  3.267e-05   0.409    0.683    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03391 on 2383 degrees of freedom
# Multiple R-squared:  7.014e-05,	Adjusted R-squared:  -0.0003495 
# F-statistic: 0.1672 on 1 and 2383 DF,  p-value: 0.6827

# model with year as random effect 
fit_wuni_vanc <- glmmTMB(uni_sim ~ (1|sameyear) + dist_km, data = vanc_wuni_df)
summary(fit_wuni_vanc)

# Family: gaussian  ( identity )
# Formula:          uni_sim ~ (1 | sameyear) + dist_km
# Data: vanc_wuni_df
# 
# AIC      BIC   logLik deviance df.resid 
# -9551.1  -9528.0   4779.6  -9559.1     2381 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance  Std.Dev.
# sameyear (Intercept) 0.0001028 0.01014 
# Residual             0.0010596 0.03255 
# Number of obs: 2385, groups:  sameyear, 2
# 
# Dispersion estimate for gaussian family (sigma^2): 0.00106 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) 8.749e-01  7.250e-03  120.68   <2e-16 ***
#   dist_km     4.217e-05  3.143e-05    1.34     0.18    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

r.squaredGLMM(fit_wuni_vanc)
# R2m        R2c
# [1,] 0.0006909391 0.08905077

vanc_wuni_df$predicted <- 4.217e-05*vanc_wuni_df$dist_km +  8.749e-01

# plot 
ggplot(vanc_wuni_df, aes(dist_km, uni_sim)) +
  geom_point(size = 2.4, color = "#e16f00",fill = "#e16f00", alpha = 0.3) +
  geom_line(aes(y = predicted, x = dist_km), color = "black") +
  scale_y_continuous(limits=c(0.55, 1)) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))
ggsave("figures/final_figures/fig4_wuni_vanc.pdf",width = 6, height = 6)

# mantel test - tests for significance 
vanc_wuni_mantel <- custom_mantel_test(vanc_wuni_df$dist_km, vanc_wuni_df$uni_sim)
vanc_wuni_mantel$mantel_r # r -0.02
vanc_wuni_mantel$p_value # p = 0.399

#### Elevation Partial Mantel Tests with Sierra Nevada Samples ####
# custom partial mantel test function for dealing with vectors instead of square matrices  
partial_mantel <- function(x, y, z, method = "spearman", nperm = 999, alternative = "two.sided") {
  # Residuals of x ~ z and y ~ z
  x_resid <- resid(lm(x ~ z))
  y_resid <- resid(lm(y ~ z))
  
  # Observed partial correlation
  obs_cor <- cor(x_resid, y_resid, method = method, use = "complete.obs")
  
  # Null distribution
  permuted_cor <- replicate(nperm, {
    cor(x_resid, sample(y_resid), method = method, use = "complete.obs")
  })
  
  # Compute p-value based on alternative hypothesis
  p_val <- switch(
    alternative,
    "two.sided" = (sum(abs(permuted_cor) >= abs(obs_cor)) + 1) / (nperm + 1),
    "greater"   = (sum(permuted_cor >= obs_cor) + 1) / (nperm + 1),
    "less"      = (sum(permuted_cor <= obs_cor) + 1) / (nperm + 1),
    stop("Invalid alternative. Choose 'two.sided', 'greater', or 'less'.")
  )
  
  list(
    partial_mantel_r = obs_cor,
    p_value = p_val,
    alternative = alternative,
    permuted_correlation_mean = mean(permuted_cor),
    permuted_correlation_sd = sd(permuted_cor),
    permuted_correlations = permuted_cor
  )
}

# elevation distance matrix 
# add elevation data to sites 
sierra_elev <- read.csv('site_data.csv')
sierra_elev$site <- trimws(sierra_elev$site, which = "right") 

# inputting elevation values
sample_data(physeq_trout)$elevation <- NA
sample_data(physeq_sierra)$elevation <- NA

# add elevation values to physeq metadata in loop 
for (site_name in unique(sample_data(physeq_trout)$site)) {
  match_index <- match(site_name, sierra_elev$site)
  if (!is.na(match_index)) {
    sample_data(physeq_trout)$elevation[sample_data(physeq_trout)$site == site_name] <- sierra_elev$elevation[match_index]
  }
}

for (site_name in unique(sample_data(physeq_sierra)$site)) {
  match_index <- match(site_name, sierra_elev$site)
  if (!is.na(match_index)) {
    sample_data(physeq_sierra)$elevation[sample_data(physeq_sierra)$site == site_name] <- sierra_elev$elevation[match_index]
  }
}

# create elevation distance matrices
trout_elev_mat <- as.matrix(dist(sample_data(physeq_trout)$elevation, method = 'euclidean'))
sierra_elev_mat <- as.matrix(dist(sample_data(physeq_sierra)$elevation, method = 'euclidean'))

# alter row and column names to include year, site, and sample name 
colnames(trout_elev_mat) <- paste(rownames(sample_data(physeq_trout)), sample_data(physeq_trout)$site, sample_data(physeq_trout)$year, sep = ";") 
rownames(trout_elev_mat) <- paste(rownames(sample_data(physeq_trout)), sample_data(physeq_trout)$site, sample_data(physeq_trout)$year, sep = ";") 

colnames(sierra_elev_mat) <- paste(rownames(sample_data(physeq_sierra)), sample_data(physeq_sierra)$site, sample_data(physeq_sierra)$year, sep = ";") 
rownames(sierra_elev_mat) <- paste(rownames(sample_data(physeq_sierra)), sample_data(physeq_sierra)$site, sample_data(physeq_sierra)$year, sep = ";") 

# lengthen elev df 
trout_elev_long <- matrixConvert(trout_elev_mat, colname = c("site1", "site2", "elev_dist"))
sierra_elev_long <- matrixConvert(sierra_elev_mat, colname = c("site1", "site2", "elev_dist"))

# add to bc similarity dataframes 
trout_sim_df$elev_dist <- trout_elev_long$elev_dist
sierra_sim_df$elev_dist <- sierra_elev_long$elev_dist

# remove same lake comparisons
trout_sim_nosame <- trout_sim_df[trout_sim_df$samelake == "N",]
sierra_sim_nosame <- sierra_sim_df[sierra_sim_df$samelake == "N",]

# mantel test - tests for significance 
trout_elev_mantel <- custom_mantel_test(trout_sim_nosame$elev_dist, trout_sim_nosame$bc_sim)
trout_elev_mantel$mantel_r # r -0.105
trout_elev_mantel$p_value # p = 0.001 

# partial mantel accounting for geographic distance 
partial_mantel(trout_sim_nosame$bc_sim, trout_sim_nosame$elev_dist, trout_sim_nosame$dist_km) # partial mantel r = 0.07, p = 0.001

geo_control_elev <- partial_mantel(trout_sim_nosame$bc_sim, trout_sim_nosame$dist_km, trout_sim_nosame$elev_dist,)
#$partial_mantel_r
#[1] -0.1496488

#$p_value
#[1] 0.001

# geography more important than elevation for trout! 

# mantel test - tests for significance 
sierra_elev_mantel <- custom_mantel_test(sierra_sim_nosame$elev_dist, sierra_sim_nosame$bc_sim)
sierra_elev_mantel$mantel_r # r -0.34
sierra_elev_mantel$p_value # p = 0.001 

# partial mantel accounting for geographic distance 
partial_mantel(sierra_sim_nosame$bc_sim, sierra_sim_nosame$elev_dist, sierra_sim_nosame$dist_km) # partial mantel r = 0.06, p = 0.02

geo_control_elev <- partial_mantel(sierra_sim_nosame$bc_sim, sierra_sim_nosame$dist_km, sierra_sim_nosame$elev_dist)
#$partial_mantel_r
#[1] -0.21

#$p_value
#[1] 0.001

# again geography more important than elevation here 