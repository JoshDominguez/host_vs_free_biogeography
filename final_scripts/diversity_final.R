# Script is for calculating alpha and beta diversity within and across systems and sample types 

### loading packages ####
packages <- c("ggplot2", "gridExtra", "knitr","phyloseq", "tidyverse", "vegan", "MuMIn", "dplyr", "ecodist", "glmmTMB", "glue", "ggforce", "geosphere", "lmPerm", "otuSummary", "ANCOMBC", "DECIPHER", "phangorn", "ape", "msa","gdm", "GUniFrac", "lme4", 'ggpubr', 'SRS' ,'car', 'reshape2', 'ggsignif', 'phytools', 'picante', 'dunn.test', 'emmeans', 'coin', 'rstatix', 'otuSummary', 'sjPlot') 
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

#### load libraries and physeq object #### 
physeq <- readRDS('data/dada2/physeq_normal.RDS')

# metadata 
metadata <- as.data.frame(sample_data(physeq))

#### Calculating Beta Div and Organizing Dataframes ####
# load bc dist that's previously been calculated 
set.seed(1999)
bray_dist_full <- as.matrix(readRDS("data/dada2/bc_dist_full_norm.RDS"))

# make rownames pasted with site, sample type, and year 
newnames <- paste(metadata$sampleNames, metadata$site, metadata$sample_type, metadata$year, sep = ";")

rownames(bray_dist_full) <- newnames
colnames(bray_dist_full) <- newnames

# lengthen df 
bc_dist_long <- matrixConvert(bray_dist_full, colname = c("site1", "site2", "bc_dist"))
bc_dist_long$site1 <- as.character(bc_dist_long$site1)
bc_dist_long$site2 <- as.character(bc_dist_long$site2)

# strsplit and rename 
# split site1 and site2 by semicolon
site1_parts <- strsplit(bc_dist_long$site1, ";")
site2_parts <- strsplit(bc_dist_long$site2, ";")

# extract parts and create new columns
bc_dist_long$site1_id <- sapply(site1_parts, "[", 1)
bc_dist_long$site1_lake <- sapply(site1_parts, "[", 2)
bc_dist_long$site1_sampletype <- sapply(site1_parts, "[", 3)
bc_dist_long$site1_year <- sapply(site1_parts, "[",4)

bc_dist_long$site2_id <- sapply(site2_parts, "[", 1)
bc_dist_long$site2_lake <- sapply(site2_parts, "[", 2)
bc_dist_long$site2_sampletype <- sapply(site2_parts, "[", 3)
bc_dist_long$site2_year <- sapply(site2_parts, "[",4)

# remove oiginal site columns 
bc_dist_long <- bc_dist_long[, !(names(bc_dist_long) %in% c("site1", "site2"))]

# same or diff lake: 
bc_dist_long$samelake <- ifelse(bc_dist_long$site1_lake == bc_dist_long$site2_lake, "samelake", "difflake")

# same or diff sample type?  
bc_dist_long$sametype <- ifelse(bc_dist_long$site1_sampletype == bc_dist_long$site2_sampletype, "sametype", "difftype")

# same or diff year? 
bc_dist_long$sameyear <- ifelse(bc_dist_long$site1_year == bc_dist_long$site2_year, "sameyear", "diffyear")

# save dataframe as csv file for future use 
write.csv(bc_dist_long, "data/dada2/bc_dist_long.csv")

# load bc_dist_long 
bc_dist_long <- read.csv("data/dada2/bc_dist_long.csv")

# are sample types more different between or within lakes?
# split up into different into different dataframes for boxplots
sierra_bc <- bc_dist_long[bc_dist_long$site1_sampletype == "water_sierra" & bc_dist_long$site2_sampletype == "water_sierra",]

vancouver_bc <- bc_dist_long[bc_dist_long$site1_sampletype == "water_VI" & bc_dist_long$site2_sampletype == "water_VI",]

trout_bc <- bc_dist_long[bc_dist_long$site1_sampletype == "trout_sierra" & bc_dist_long$site2_sampletype == "trout_sierra",]

stickle_bc <- bc_dist_long[bc_dist_long$site1_sampletype == "stickle_VI" & bc_dist_long$site2_sampletype == "stickle_VI",]

#### Alpha Diversity Differences #### 
richness_est <- estimate_richness(physeq)
rownames(richness_est) <- gsub("^X", "", rownames(richness_est))

# combined with metadata 
richness_est <- merge(richness_est, metadata, by = 'row.names')
richness_est$system <- ifelse(richness_est$sample_type %in% c("stickle_VI", "water_VI"), "Vancouver", "Sierra Nevada")

rich_plot <- ggplot(richness_est, aes(x = factor(sample_type, levels = c("trout_sierra", "water_sierra", "stickle_VI","water_VI")), y = Observed, fill = factor(system, levels = c("Vancouver", "Sierra Nevada")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e16f00", "dodgerblue")) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave('figures/final_figures/ASV_richness.pdf', width = 6, height = 4,dpi = 300)

# Mean and SD for each sample type 
rich_summary <- tapply(richness_est$Observed, richness_est$sample_type, function(x) {
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  c(mean, median, sd)
})

# $stickle_VI
# [1] 127.5421  75.0000 139.1112
# 
# $trout_sierra
# [1] 130.3432  60.0000 155.0270
# 
# $water_sierra
# [1] 192.1475 175.0000  58.0617
# 
# $water_VI
# [1] 324.3803 276.0000 170.8726

log_observed <- aov(log(Observed) ~ sample_type, data = richness_est)
summary(log_observed)
TukeyHSD(log_observed) # same result here 

shan_plot <- ggplot(richness_est, aes(x = factor(sample_type, levels = c("trout_sierra", "water_sierra", "stickle_VI","water_VI")), y = Shannon, fill = factor(system, levels = c("Vancouver", "Sierra Nevada")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#e16f00", "dodgerblue")) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave('figures/final_figures/shannon_boxplot.pdf', width = 6, height = 4,dpi = 300)

# Mean and SD Shannon diversity 
shan_summary <- tapply(richness_est$Shannon, richness_est$sample_type, function(x) {
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  c(mean, median, sd)
})

# $stickle_VI
# [1] 2.686622 2.643791 1.325837
# 
# $trout_sierra
# [1] 2.188160 2.055638 1.046312
# 
# $water_sierra
# [1] 3.6715842 3.6869855 0.3481709
# 
# $water_VI
# [1] 4.4296399 4.2898766 0.5213757


#### ANOVA and Post Hoc Tests of Differences in Observed and Shannon ####
### shannon 
shannon_aov <- lm(Shannon~ sample_type, data=richness_est)
Anova(shannon_aov, type=2) # sample type and site significant, no effect of year

#  posthoc looking at sample type 
shannon_test <-emmeans(shannon_aov, pairwise~sample_type)

#contrast                    estimate    SE  df t.ratio p.value
#stickle_VI - trout_sierra      0.498 0.104 784   4.791  <.0001
#stickle_VI - water_sierra     -0.985 0.158 784  -6.224  <.0001
#stickle_VI - water_VI         -1.743 0.148 784 -11.775  <.0001
#trout_sierra - water_sierra   -1.483 0.174 784  -8.523  <.0001
#trout_sierra - water_VI       -2.241 0.165 784 -13.601  <.0001
#water_sierra - water_VI       -0.758 0.203 784  -3.726  0.0012

### observed richness 
observed_aov <- lm(Observed ~ sample_type, data=richness_est)
Anova(observed_aov, type=2) 
plot(allEffects(observed_aov))

#  posthoc looking at sample type 
observed_test <-emmeans(observed_aov, pairwise~sample_type)

#contrast                    estimate   SE  df t.ratio p.value
#stickle_VI - trout_sierra       -2.8 12.6 784  -0.222  0.9961
#stickle_VI - water_sierra      -64.6 19.2 784  -3.363  0.0045
#stickle_VI - water_VI         -196.8 18.0 784 -10.953  <.0001
#trout_sierra - water_sierra    -61.8 21.1 784  -2.925  0.0186
#trout_sierra - water_VI       -194.0 20.0 784  -9.699  <.0001
#water_sierra - water_VI       -132.2 24.7 784  -5.354  <.0001

#### within lake beta diversity figure - plot of mean within lake bray curtis distance #### 

# mean dist within lakes for each sample type as a dataframe 
mean_dist_trout <- data.frame(bc_dist = diag(vegan::meandist(trout_bc, sample_data(physeq_trout)$site)), sample_type = "trout", lake = names(diag(vegan::meandist(trout_bc, sample_data(physeq_trout)$site))))

mean_dist_stickle <- data.frame(bc_dist = diag(vegan::meandist(stickle_bc, sample_data(physeq_stickle)$site)), sample_type = "stickle", lake = names(diag(vegan::meandist(stickle_bc, sample_data(physeq_stickle)$site))))

mean_dist_sierra <- data.frame(bc_dist = diag(vegan::meandist(sierra_bc, sample_data(physeq_sierra)$site)), sample_type = "sierra", lake = names(diag(vegan::meandist(sierra_bc, sample_data(physeq_sierra)$site))))

mean_dist_vancouver <- data.frame(bc_dist = diag(vegan::meandist(vanc_bc, sample_data(physeq_VI)$site)), sample_type = "vancouver", lake = names(diag(vegan::meandist(vanc_bc, sample_data(physeq_VI)$site))))

# combine dataframes 
within_lake_beta <- rbind(mean_dist_trout, mean_dist_stickle, mean_dist_vancouver, mean_dist_sierra)
within_lake_beta$system <- ifelse(within_lake_beta$sample_type %in% c("trout", "sierra"), "Sierra Nevada", "Vancouver Island")

# plot 
within_lake <- ggplot(within_lake_beta, aes(x = factor(sample_type, levels = c("trout", "sierra", "stickle", "vancouver")), y = bc_dist, fill = system)) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab(NULL) + 
  ylab(NULL) + 
  scale_fill_manual(values = c("dodgerblue", "#e16f00")) +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 12), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))
ggsave("figures/final_figures/withinlake_beta.pdf", width = 6, height = 4, dpi = 300)

# Mean and SD within lake beta diversity 
beta_summary <- tapply(within_lake_beta$bc_dist, within_lake_beta$sample_type, function(x) {
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  c(mean, median, sd)
})

# $sierra
# [1] 0.58786064 0.57399343 0.03932412
# 
# $stickle
# [1] 0.86723831 0.87582660 0.02925626
# 
# $trout
# [1] 0.68811620 0.67832771 0.03874487
# 
# $vancouver
# [1] 0.2705502 0.2320688 0.1153358

# AOV 
beta_aov <- lm(bc_dist ~ sample_type, data=within_lake_beta)
Anova(beta_aov, type=2) 

#  posthoc looking at sample type 
beta_test <-emmeans(beta_aov, pairwise~sample_type)

# $emmeans
# sample_type emmean     SE df lower.CL upper.CL
# sierra       0.588 0.0176 64    0.553    0.623
# stickle      0.867 0.0156 64    0.836    0.898
# trout        0.688 0.0176 64    0.653    0.723
# vancouver    0.271 0.0156 64    0.239    0.302
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast            estimate     SE df t.ratio p.value
# sierra - stickle      -0.279 0.0235 64 -11.863  <.0001
# sierra - trout        -0.100 0.0249 64  -4.027  0.0009
# sierra - vancouver     0.317 0.0235 64  13.474  <.0001
# stickle - trout        0.179 0.0235 64   7.606  <.0001
# stickle - vancouver    0.597 0.0221 64  26.974  <.0001
# trout - vancouver      0.418 0.0235 64  17.731  <.0001
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# combine plots 
plot_width <- 5
plot_height <- 4
grid.arrange(rich_plot, shan_plot, within_lake, ncol = 3, heights = unit(rep(plot_height, 3), "in"), widths = unit(plot_width, "in"))

# Save to PDF
pdf("figures/final_figures/combined_diversity_horizontal.pdf", width = 6, height = 4 * 3)  # Adjust height to fit all plots
grid.arrange(rich_plot, shan_plot, within_lake, ncol = 1, heights = unit(rep(4, 3), "in"), widths = unit(6, "in"))

# Close the PDF device
dev.off()
ggsave("figures/final_figures/combined_diversity.pdf", width = 6, height = 4)

#### unifrac within lake beta diversity ####
# weighted unifrac distance 
wunifrac_dist <- as.matrix(readRDS("data/dada2/wuni_dist_mat.RDS"))

# individual sample type weighted unifrac matrices 
trout_wuni <- as.dist(wunifrac_dist[rownames(wunifrac_dist) %in% sample_names(physeq_trout), colnames(wunifrac_dist) %in% sample_names(physeq_trout)])

stickle_wuni <- as.dist(wunifrac_dist[rownames(wunifrac_dist) %in% sample_names(physeq_stickle), colnames(wunifrac_dist) %in% sample_names(physeq_stickle)])

sierra_wuni <- as.dist(wunifrac_dist[rownames(wunifrac_dist) %in% sample_names(physeq_sierra), colnames(wunifrac_dist) %in% sample_names(physeq_sierra)])

vanc_wuni <- as.dist(wunifrac_dist[rownames(wunifrac_dist) %in% sample_names(physeq_VI), colnames(wunifrac_dist) %in% sample_names(physeq_VI)])

# calculate mean between group distances
mean_wuni_trout <- data.frame(wuni_dist = diag(vegan::meandist(trout_wuni, sample_data(physeq_trout)$site)), sample_type = "trout", lake = names(diag(vegan::meandist(trout_wuni, sample_data(physeq_trout)$site))))

mean_wuni_stickle <- data.frame(wuni_dist = diag(vegan::meandist(stickle_wuni, sample_data(physeq_stickle)$site)), sample_type = "stickle", lake = names(diag(vegan::meandist(stickle_wuni, sample_data(physeq_stickle)$site))))

mean_wuni_sierra <- data.frame(wuni_dist = diag(vegan::meandist(sierra_wuni, sample_data(physeq_sierra)$site)), sample_type = "sierra", lake = names(diag(vegan::meandist(sierra_wuni, sample_data(physeq_sierra)$site))))

mean_wuni_vancouver <- data.frame(wuni_dist = diag(vegan::meandist(vanc_wuni, sample_data(physeq_VI)$site)), sample_type = "vancouver", lake = names(diag(vegan::meandist(vanc_wuni, sample_data(physeq_VI)$site))))

# combine dataframes 
within_lake_wuni <- rbind(mean_wuni_trout, mean_wuni_stickle, mean_wuni_vancouver, mean_wuni_sierra)
within_lake_wuni$system <- ifelse(within_lake_beta$sample_type %in% c("trout", "sierra"), "Sierra Nevada", "Vancouver Island")

# plot weighted unifrac 
ggplot(within_lake_wuni, aes(x = factor(within_lake_wuni$sample_type, levels = c("trout", "sierra", "stickle", "vancouver")), y = wuni_dist, fill = system)) + 
  geom_boxplot() + 
  theme_classic() + 
  xlab(NULL) + 
  ylab(NULL) + 
  scale_fill_manual(values = c("dodgerblue", "#e16f00")) +
  theme(legend.title = element_blank(), 
        axis.text = element_text(size = 12))
ggsave("figures/final_figures/withinlake_wuni.pdf", width = 6, height = 4, dpi = 300)

# anova and post hoc tests
wuni_aov <- lm(wuni_dist~ sample_type, data=within_lake_wuni)
Anova(wuni_aov, type=2) # sample type significant

#  posthoc looking at sample type 
wuni_posttest <-emmeans(wuni_aov, pairwise~sample_type)
wuni_posttest

# $emmeans
# sample_type emmean      SE df lower.CL upper.CL
# sierra      0.1161 0.00396 64   0.1081   0.1240
# stickle     0.1912 0.00352 64   0.1842   0.1982
# trout       0.1288 0.00396 64   0.1209   0.1367
# vancouver   0.0377 0.00352 64   0.0306   0.0447
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast            estimate      SE df t.ratio p.value
# sierra - stickle     -0.0751 0.00529 64 -14.191  <.0001
# sierra - trout       -0.0127 0.00560 64  -2.276  0.1144
# sierra - vancouver    0.0784 0.00529 64  14.803  <.0001
# stickle - trout       0.0624 0.00529 64  11.785  <.0001
# stickle - vancouver   0.1535 0.00497 64  30.867  <.0001
# trout - vancouver     0.0911 0.00529 64  17.210  <.0001
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# sierra and brook trout not different! 

# Mean and standard deviation 
wuni_summary <- tapply(within_lake_wuni$wuni_dist, within_lake_wuni$sample_type, function(x) {
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  c(mean, median, sd)
})
wuni_summary

# $sierra
# [1] 0.11605424 0.11927201 0.02023746
# 
# $stickle
# [1] 0.19119539 0.19238382 0.01606008
# 
# $trout
# [1] 0.12879533 0.13263134 0.01010514
# 
# $vancouver
# [1] 0.03767279 0.03379614 0.01340476
