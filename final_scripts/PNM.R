#### Prokaryote Neutral Model on Bacterial Metacommunities  #### 

#### Packages ####
packages <- c("ggplot2", "gridExtra", "knitr","phyloseq", "tidyverse", "vegan", "MuMIn", "dplyr", "ecodist", "glmmTMB", "glue", "ggforce", "geosphere", "lmPerm", "otuSummary", "ANCOMBC", "DECIPHER", "phangorn", "ape", "msa","gdm", "GUniFrac", "lme4", 'ggpubr', 'SRS' ,'car', 'reshape2', 'ggsignif', 'phytools', 'picante', 'statix')
#sapply(packages, install.packages(packages), character.only = TRUE)
sapply(packages, require, character.only = TRUE)

#### PNM Model Function #### 
#model function
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}

#### Fitting Model using new physeq objects #### 
physeq_meta <- readRDS('data/dada2/physeq_normal.RDS')

# split into separate objects 
physeq_stickle <- subset_samples(physeq_meta, sample_type == "stickle_VI")
physeq_stickle <- prune_taxa(taxa_sums(physeq_stickle ) > 0, physeq_stickle)

physeq_trout <- subset_samples(physeq_meta, sample_type == "trout_sierra")
physeq_trout <- prune_taxa(taxa_sums(physeq_trout) > 0, physeq_trout) 

physeq_waterVI <- subset_samples(physeq_meta, sample_type == "water_VI")
physeq_waterVI <- prune_taxa(taxa_sums(physeq_waterVI) > 0, physeq_waterVI)

physeq_watersierra <- subset_samples(physeq_meta, sample_type == "water_sierra")
physeq_watersierra <- prune_taxa(taxa_sums(physeq_watersierra) > 0, physeq_watersierra)

### Fitting to Trout ####
sncm_trout <- sncm.fit(data.frame(otu_table(t(physeq_trout))))
sncm_trout # R2 = 0.29, m = 0.0075, ci is high - 0.0076 

# fit is better than binomial and poisson models 

# asv level of frequency and relative to fit 
sncm_trout_freq <- sncm.fit(spp=data.frame(otu_table(t(physeq_trout))), pool=NULL, stats=FALSE, taxon=NULL)

# above or below for each 
sncm_trout_freq$low<- ifelse(sncm_trout_freq$freq<sncm_trout_freq$pred.lwr, 0, 1)
sncm_trout_freq$high<-ifelse(sncm_trout_freq$freq>sncm_trout_freq$pred.upr, 1, 0)
sncm_trout_freq$vals<-sncm_trout_freq$high + sncm_trout_freq$low

# assigning ASVs with trout 
sncm_trout_freq$vals[sncm_trout_freq$vals == 0] <- "Below Prediction"   
sncm_trout_freq$vals[sncm_trout_freq$vals == 1] <- "Within Prediction"   
sncm_trout_freq$vals[sncm_trout_freq$vals == 2] <- "Above Prediction"  
sncm_trout_freq$vals<-factor(sncm_trout_freq$vals, levels =c("Above Prediction" ,  "Within Prediction" ,"Below Prediction"))

trout_PNM <- ggplot(sncm_trout_freq, aes(x = log(p), y = freq)) +
  geom_point(size = 2, aes(color = vals)) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  geom_line(aes(y = freq.pred, x = log(p)), color = "red") +
  geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), color = "black", alpha = 0, linetype = "dashed") + 
  scale_x_continuous(limits = c(-14.2, -1,5)) +
  labs(x = NULL, y = NULL) +
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.position = c(0.05, 1), 
        legend.justification = c(0.05, 1), 
        axis.text = element_text(size = 14)) 
ggsave("figures/final_figures/trout_PNM.pdf", width = 5, height = 5)

#### Fitting to Sierra Water #### 
sncm_sierra <- sncm.fit(data.frame(otu_table(t(physeq_watersierra))))
sncm_sierra # R2 = 0.63

# asv level of frequency and relative to fit 
sncm_watersierra_freq <- sncm.fit(spp=data.frame(otu_table(t(physeq_watersierra))), pool=NULL, stats=FALSE, taxon=NULL)

# above or below for each 
sncm_watersierra_freq$low<- ifelse(sncm_watersierra_freq$freq<sncm_watersierra_freq$pred.lwr, 0, 1)
sncm_watersierra_freq$high<-ifelse(sncm_watersierra_freq$freq>sncm_watersierra_freq$pred.upr, 1, 0)
sncm_watersierra_freq$vals<-sncm_watersierra_freq$high + sncm_watersierra_freq$low

# assigning ASVs with watersierra 
sncm_watersierra_freq$vals[sncm_watersierra_freq$vals == 0] <- "Below Prediction"   
sncm_watersierra_freq$vals[sncm_watersierra_freq$vals == 1] <- "Within Prediction"   
sncm_watersierra_freq$vals[sncm_watersierra_freq$vals == 2] <- "Above Prediction"  
sncm_watersierra_freq$vals<-factor(sncm_watersierra_freq$vals, levels =c("Above Prediction" ,  "Within Prediction" ,"Below Prediction"))

sierra_PNM <-ggplot(sncm_watersierra_freq, aes(x = log(p), y = freq)) +
  geom_point(size = 2, aes(color = vals)) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  geom_line(aes(y = freq.pred, x = log(p)), color = "red") +
  geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), color = "black", alpha = 0, linetype = "dashed") + 
  scale_x_continuous(limits = c(-14.2, -1,5)) +
  theme_classic() + 
  labs(x = NULL, y = NULL) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.position = c(0.05, 1), 
        legend.justification = c(0.05, 1), 
        axis.text = element_text(size = 14))
ggsave("figures/final_figures/PNM_sierra.pdf", width = 5, height = 5)

#### Fitting to Stickleback #### 
sncm_stickle <- sncm.fit(data.frame(otu_table(t(physeq_stickle))))
sncm_stickle # R2 = 0.34, m = 0.002
# asv level of frequency and relative to fit 
sncm_stickle_freq <- sncm.fit(spp=data.frame(otu_table(t(physeq_stickle))), pool=NULL, stats=FALSE, taxon=NULL)

# above or below for each 
sncm_stickle_freq$low<- ifelse(sncm_stickle_freq$freq<sncm_stickle_freq$pred.lwr, 0, 1)
sncm_stickle_freq$high<-ifelse(sncm_stickle_freq$freq>sncm_stickle_freq$pred.upr, 1, 0)
sncm_stickle_freq$vals<-sncm_stickle_freq$high + sncm_stickle_freq$low

# assigning ASVs with stickle 
sncm_stickle_freq$vals[sncm_stickle_freq$vals == 0] <- "Below Prediction"   
sncm_stickle_freq$vals[sncm_stickle_freq$vals == 1] <- "Within Prediction"   
sncm_stickle_freq$vals[sncm_stickle_freq$vals == 2] <- "Above Prediction"  
sncm_stickle_freq$vals<-factor(sncm_stickle_freq$vals, levels =c("Above Prediction" ,  "Within Prediction" ,"Below Prediction"))

stickle_PNM <- ggplot(sncm_stickle_freq, aes(x = log(p), y = freq)) +
  geom_point(size = 2, aes(color = vals)) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) + 
  geom_line(aes(y = freq.pred, x = log(p)), color = "red") +
  geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), color = "black", alpha = 0, linetype = "dashed") + 
  scale_x_continuous(limits = c(-14.2, -1,5)) +
  theme_classic() + 
  labs(x = NULL, y = NULL) +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.position = c(0.05, 1), 
        legend.justification = c(0.05, 1),
        axis.text = element_text(size = 14))
ggsave('figures/final_figures/PNM_stickle.pdf', width = 5, height = 5)


#### Fitting to VI Water #### 
sncm_waterVI <- sncm.fit(data.frame(otu_table(t(physeq_waterVI))))
sncm_waterVI # R2 = 0.74, m = 0.098, ci is 0.006 

# fit is better than binomial and poisson models 

# asv level of frequency and relative to fit 
sncm_waterVI_freq <- sncm.fit(spp=data.frame(otu_table(t(physeq_waterVI))), pool=NULL, stats=FALSE, taxon=NULL)

# above or below for each 
sncm_waterVI_freq$low<- ifelse(sncm_waterVI_freq$freq<sncm_waterVI_freq$pred.lwr, 0, 1)
sncm_waterVI_freq$high<-ifelse(sncm_waterVI_freq$freq>sncm_waterVI_freq$pred.upr, 1, 0)
sncm_waterVI_freq$vals<-sncm_waterVI_freq$high + sncm_waterVI_freq$low

# assigning ASVs with waterVI 
sncm_waterVI_freq$vals[sncm_waterVI_freq$vals == 0] <- "Below Prediction"   
sncm_waterVI_freq$vals[sncm_waterVI_freq$vals == 1] <- "Within Prediction"   
sncm_waterVI_freq$vals[sncm_waterVI_freq$vals == 2] <- "Above Prediction"  
sncm_waterVI_freq$vals<-factor(sncm_waterVI_freq$vals, levels =c("Above Prediction" ,  "Within Prediction" ,"Below Prediction"))

PNM_vancouver <- ggplot(sncm_waterVI_freq, aes(x = log(p), y = freq)) +
  geom_point(size = 2, aes(color = vals)) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  geom_line(aes(y = freq.pred, x = log(p)), color = "red") +
  geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), color = "black", alpha = 0, linetype = "dashed") + 
  scale_x_continuous(limits = c(-14.2, -1,5)) +
  labs(x = NULL, y = NULL) +
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.position = c(0.05, 1), 
        legend.justification = c(0.05, 1), 
        axis.text = element_text(size = 14))
ggsave('figures/final_figures/PNM_vancouver.pdf', width = 5, height = 5)

# x lims 
# max 
# stickle: -2.6, VI: -3.2, Sierra: -2.8, Trout: -1.8 
# min 
# stickle: -14.2, VI: -12.3, Sierra: -12.1, Trout: -13.1

#### Fig 5 - Combing PNM Plots #### 
pdf("figures/final_figures/Fig_5.pdf", width = 10, height = 10)  # Adjust width and height as needed
grid.arrange(trout_PNM, sierra_PNM, stickle_PNM, PNM_vancouver, ncol = 2, nrow = 2)
dev.off()

#### What are the over, under, and neutral taxa in Trout? #### 


#### Redundant Taxa - Fish - What is over and under represented? #### 
# First
# combine tax table with Within Above and Below 
sncm_trout_freq$ASV <- rownames(sncm_trout_freq)
trout_asv_tax <- merge(sncm_trout_freq, data.frame(tax_table(physeq_trout)), by = 'row.names')

sncm_stickle_freq$ASV <- rownames(sncm_stickle_freq)
stickle_asv_tax <- merge(sncm_stickle_freq, data.frame(tax_table(physeq_stickle)), by = 'row.names')


# calculate residuals/deviance from model for each asv 
trout_asv_tax$dev <- trout_asv_tax$freq - trout_asv_tax$freq.pred
stickle_asv_tax$dev <- stickle_asv_tax$freq - stickle_asv_tax$freq.pred

# make new column that differentiates neutral vs nonneutral 
trout_asv_tax$neutrality <- ifelse(trout_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")
stickle_asv_tax$neutrality <- ifelse(stickle_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")

# which phyla have the highest frequency ? 
mean_phyla_freq_trout <- tapply(trout_asv_tax$freq, factor(trout_asv_tax$Phylum), mean)
mean_phyla_freq_trout[order(-mean_phyla_freq_trout)] # Cyano, Fuso, Actino, Synergistota, Chloroflexi, Firmicutes, Planctomycetota

mean_phyla_freq_stickle <- tapply(stickle_asv_tax$freq, factor(stickle_asv_tax$Phylum), mean)
mean_phyla_freq_stickle[order(-mean_phyla_freq_stickle)] # Deinococcota, Cyano, Thermo, Spiro, Euryarchaeota, Halo 

# which have the greatest mean rel abundance ? 
mean_phyla_abundance_trout <- tapply(log(trout_asv_tax$p), factor(trout_asv_tax$Phylum), mean)
mean_phyla_abundance_trout[order(-mean_phyla_abundance_trout)] # Fuso, Cyano, Synergistota, Firmicutes, Actinobacteriota

mean_phyla_abundance_stickle <- tapply(log(stickle_asv_tax$p), factor(stickle_asv_tax$Phylum), mean)
mean_phyla_abundance_stickle[order(-mean_phyla_abundance_stickle)] # Fuso, Thermoplasmatota, Deinococcota, Euryarchaeota, Cyanobacteria, Firmicutes

# which have the greatest mean deviance from prediction? 
mean_phyla_dev_trout <- tapply(abs(trout_asv_tax$dev), factor(trout_asv_tax$Phylum), 
                         mean)
mean_phyla_dev_trout[order(-mean_phyla_dev_trout)] # Fuso, Cyano, Actino, Synergistota, Chloroflexi 

mean_phyla_dev_stickle <- tapply(abs(stickle_asv_tax$dev), factor(stickle_asv_tax$Phylum), 
                               mean)
mean_phyla_dev_stickle[order(-mean_phyla_dev_stickle)] # Deinococcota, Spirochaetota, Thermo

# max value (absolute)
max_phyla_dev_trout <- tapply(abs(trout_asv_tax$dev), factor(trout_asv_tax$Phylum), 
                        max)
max_phyla_dev_trout[order(-max_phyla_dev_trout)] # Fuso, Firmicutes, Actino, Proteo, Plancto

max_phyla_dev_stickle <- tapply(abs(stickle_asv_tax$dev), factor(stickle_asv_tax$Phylum), 
                              max)
max_phyla_dev_stickle[order(-max_phyla_dev_stickle)] # Firmicutes, Spiro, Proteo, Actino 

# which phyla have the most non-neutral asvs? 
table(trout_asv_tax$vals, trout_asv_tax$Phylum)
table(trout_asv_tax$neutrality, trout_asv_tax$Phylum)

table(stickle_asv_tax$vals, stickle_asv_tax$Phylum)
table(stickle_asv_tax$neutrality, stickle_asv_tax$Phylum)

# Shared taxa/asvs between stickle and trout? 
trout_asv_tax$sample_type <- rep('trout', length.out = nrow(trout_asv_tax))
stickle_asv_tax$sample_type <- rep('stickle', length.out = nrow(stickle_asv_tax))

common_fish_asvs <- intersect(trout_asv_tax$ASV, stickle_asv_tax$ASV) 
length(common_fish_asvs) # 2154 shared asvs between stickle and trout

trout_common <- trout_asv_tax[trout_asv_tax$ASV %in% common_fish_asvs, ]
stickle_common <- stickle_asv_tax[stickle_asv_tax$ASV %in% common_fish_asvs, ]

# rbind them 
common_fish_df <- rbind(trout_common, stickle_common)
table(common_fish_df$neutrality) # 2440 neutral, 2116 noneutral 
table(common_fish_df$neutrality, common_fish_df$sample_type) # lots of non-neutral stickle, less trout - just more samples? 

# max and min deviances 
tapply(common_fish_df$dev, common_fish_df$sample_type, function(x) {
  min <- min(x)
  max <- max(x)
  c(min,max)
} 
) # stickle has less high deviance compared to trout 

# I need to eat some food but I want to keep coding, I think i'll go eat food and code 
trout_common$freq_trout <- trout_common$freq
trout_common$predfreq_trout <-trout_common$freq.pred
trout_common$p_trout <- trout_common$p
trout_common$dev_trout <- trout_common$dev
trout_common$neutral_trout <- trout_common$neutrality
rownames(trout_common) <- trout_common$ASV

# new col names stickle 
stickle_common$freq_stickle <- stickle_common$freq
stickle_common$predfreq_stickle <- stickle_common$freq.pred
stickle_common$p_stickle <- stickle_common$p
stickle_common$dev_stickle <- stickle_common$dev
stickle_common$neutral_stickle <- stickle_common$neutrality
rownames(stickle_common) <- stickle_common$ASV

# merge the objects by ASVs - want stickle freq, trout freq, stickle p, trout p, and taxa stuff
common_wide <- merge(trout_common[,13:27], stickle_common[,24:28], by = 'row.names')

# is there a correlation in the abundance of asvs in trout and stickle? 
plot(log(common_wide$p_stickle), log(common_wide$p_trout))
PNM_lm <- lm(log(p_trout) ~ log(p_stickle), data = common_wide)
summary(PNM_lm)
cor(log(common_wide$p_trout), log(common_wide$p_stickle), method = "spearman")
# if its abundant in trout - its abundant in stickle! some sort of correlation - 0.35 cor coefficient with spearman 
# it seems like low abundance ones are not correlated but the higher abundance ones are 

# what about frequency across the metacommunity? are those correlated? 
plot(common_wide$freq_trout, common_wide$freq_stickle) # looks like a pretty strong correlation ? So many low frequency ones looks better on a log scale 
cor(common_wide$freq_trout, common_wide$freq_stickle, method = "spearman")
# stronger correlation of 0.39, but the data looks strange/weird ? 
plot(log(common_wide$freq_trout), log(common_wide$freq_stickle))

# Conclusion - more frequent and higher abundance asvs are common across both 

# nicer plot of abundance in trout vs abundance in stickle 
ggplot(common_wide, aes(log(p_stickle), log(p_trout))) +
  geom_point(size = 1.2, alpha = 1, shape = 1) + 
  geom_smooth(method = lm, level = 0.95, color = "black", alpha = 0.3, size = 0.8) +
  labs(x = "log Mean Relative Abundance in Stickleback", y= "log Mean Relative Abundance in Brook trout)") +
  theme_classic()
ggsave("figures/dada2/stickle_trout_abund.tiff", width = 7, height = 5)

# dev in one very abundance in the other? - no relationship 
plot(log(common_wide$p_stickle), common_wide$dev_trout)
Dev_lm <- lm(log(p_stickle) ~ dev_trout, data = common_wide)

# plot dev in one vs dev in the other? 
plot(common_wide$dev_stickle, common_wide$dev_trout) # looks like stuff that has high deviance in trout have high deviance in stickle? 

# plotting trout - overrepresented in stickle and trout 



#### Similar Taxa between Water ? ##### 
# combine tax table with Within Above and Below 
sncm_watersierra_freq$ASV <- rownames(sncm_watersierra_freq)
watersierra_asv_tax <- merge(sncm_watersierra_freq, data.frame(tax_table(physeq_watersierra)), by = 'row.names')

sncm_waterVI_freq$ASV <- rownames(sncm_waterVI_freq)
waterVI_asv_tax <- merge(sncm_waterVI_freq, data.frame(tax_table(physeq_waterVI)), by = 'row.names')


# calculate residuals/deviance from model for each asv 
watersierra_asv_tax$dev <- watersierra_asv_tax$freq - watersierra_asv_tax$freq.pred
waterVI_asv_tax$dev <- waterVI_asv_tax$freq - waterVI_asv_tax$freq.pred

# make new column that differentiates neutral vs nonneutral 
watersierra_asv_tax$neutrality <- ifelse(watersierra_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")
waterVI_asv_tax$neutrality <- ifelse(waterVI_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")

# which phyla have the highest frequency ? 
mean_phyla_freq_watersierra <- tapply(watersierra_asv_tax$freq, factor(watersierra_asv_tax$Phylum), mean)
mean_phyla_freq_watersierra[order(-mean_phyla_freq_watersierra)] # Actinobacteria, bacteroidota, Proteobacteria, deinococcota, armatimonadota, everything else less than 5%

mean_phyla_freq_waterVI <- tapply(waterVI_asv_tax$freq, factor(waterVI_asv_tax$Phylum), mean)
mean_phyla_freq_waterVI[order(-mean_phyla_freq_waterVI)] # Actino, Proteo, Bacteroidota, Armatimonadota, Verruco, Cyano, Campy, Fuso, RCP, Bdell, rest under 5%

# which have the greatest mean deviance from prediction? 
mean_phyla_dev_watersierra <- tapply(abs(watersierra_asv_tax$dev), factor(watersierra_asv_tax$Phylum), 
                               mean)
mean_phyla_dev_watersierra[order(-mean_phyla_dev_watersierra)] # Deino, Actino, Bacteroi, Proteo

mean_phyla_dev_waterVI <- tapply(abs(waterVI_asv_tax$dev), factor(waterVI_asv_tax$Phylum), 
                                 mean)
mean_phyla_dev_waterVI[order(-mean_phyla_dev_waterVI)] # Armat, Cyano, Actino, Campy, Bacteroidota

# max value (absolute)
max_phyla_dev_watersierra <- tapply(abs(watersierra_asv_tax$dev), factor(watersierra_asv_tax$Phylum), 
                              max)
max_phyla_dev_watersierra[order(-max_phyla_dev_watersierra)] # Bacteroidota, Proteo, Verruco

max_phyla_dev_waterVI <- tapply(abs(waterVI_asv_tax$dev), factor(waterVI_asv_tax$Phylum), 
                                max)
max_phyla_dev_waterVI[order(-max_phyla_dev_waterVI)] # Plancto, Bacteroidota, Verruco, Cyano

# which phyla have the most non-neutral asvs? 
table(watersierra_asv_tax$vals, watersierra_asv_tax$Phylum)
table(watersierra_asv_tax$neutrality, watersierra_asv_tax$Phylum)

table(waterVI_asv_tax$vals, waterVI_asv_tax$Phylum)
table(waterVI_asv_tax$neutrality, waterVI_asv_tax$Phylum)

# Shared taxa/asvs between waterVI and watersierra? 
watersierra_asv_tax$sample_type <- rep('watersierra', length.out = nrow(watersierra_asv_tax))
waterVI_asv_tax$sample_type <- rep('waterVI', length.out = nrow(waterVI_asv_tax))

common_water_asvs <- intersect(watersierra_asv_tax$ASV, waterVI_asv_tax$ASV) # 1044 shared asvs between waterVI and watersierra

watersierra_common <- watersierra_asv_tax[watersierra_asv_tax$ASV %in% common_water_asvs, ]
waterVI_common <- waterVI_asv_tax[waterVI_asv_tax$ASV %in% common_water_asvs, ]

# rbind them 
common_water_df <- rbind(watersierra_common, waterVI_common)
table(common_water_df$neutrality) # 824 neutral, 242 non neutral 
table(common_water_df$neutrality, common_water_df$sample_type) # more non-neutral waterVI ASVs

# max and min deviances 
tapply(common_water_df$dev, common_water_df$sample_type, function(x) {
  min <- min(x)
  max <- max(x)
  c(min,max)
} 
) # waterVI has less high deviance compared to watersierra, but lower deviance 

# I need to eat some food but I want to keep coding, I think i'll go eat food and code 
watersierra_common$freq_watersierra <- watersierra_common$freq
watersierra_common$predfreq_watersierra <-watersierra_common$freq.pred
watersierra_common$p_watersierra <- watersierra_common$p
watersierra_common$dev_watersierra <- watersierra_common$dev
watersierra_common$neutral_watersierra <- watersierra_common$neutrality
rownames(watersierra_common) <- watersierra_common$ASV

# new col names waterVI 
waterVI_common$freq_waterVI <- waterVI_common$freq
waterVI_common$predfreq_waterVI <- waterVI_common$freq.pred
waterVI_common$p_waterVI <- waterVI_common$p
waterVI_common$dev_waterVI <- waterVI_common$dev
waterVI_common$neutral_waterVI <- waterVI_common$neutrality
rownames(waterVI_common) <- waterVI_common$ASV

# merge the objects by ASVs - want waterVI freq, watersierra freq, waterVI p, watersierra p, and taxa stuff
common_wide_water <- merge(watersierra_common[,13:27], waterVI_common[,24:28], by = 'row.names')

# is there a correlation in the abundance of asvs in watersierra and waterVI? 
plot(log(common_wide_water$p_waterVI), log(common_wide_water$p_watersierra)) # yeah? is this real what the fuck 
PNM_lm_water <- lm(log(p_watersierra) ~ log(p_waterVI), data = common_wide_water)
summary(PNM_lm_water)
cor(log(common_wide_water$p_watersierra), log(common_wide_water$p_waterVI), method = "spearman")
# if its abundant in watersierra - its abundant in waterVI! decent correlation of 0.45
# it seems like low abundance ones are not correlated but the higher abundance ones are 

# what about frequency across the metacommunity? are those correlated? 
plot(common_wide_water$freq_watersierra, common_wide_water$freq_waterVI) # looks like a pretty strong correlation ? So many low frequency ones looks better on a log scale 
cor(common_wide_water$freq_watersierra, common_wide_water$freq_waterVI, method = "spearman")
# stronger correlation of 0.5, but the data looks strange/weird ? 
plot(log(common_wide_water$freq_watersierra), log(common_wide_water$freq_waterVI))

# Conclusion - more frequent and higher abundance asvs are common across both water 

# nicer plot of abundance in trout vs abundance in stickle 
ggplot(common_wide_water, aes(log(p_waterVI), log(p_watersierra))) +
  geom_point(size = 1.2, alpha = 1, shape = 1) + 
  geom_smooth(method = lm, level = 0.95, color = "black", alpha = 0.3, size = 0.8) +
  labs(x = "log Mean Relative Abundance in Vancouver Island Water", y= "log Mean Relative Abundance in Sierra Nevada Water") +
  theme_classic()
ggsave("figures/dada2/sierra_VI_abund.tiff", width = 7, height = 5)

# dev in one very abundance in the other?
plot(log(common_wide$p_watersierra), common_wide$dev_waterVI)
Dev_lm <- lm(log(p_stickle) ~ dev_trout, data = common_wide)

# plot dev in one vs dev in the other? 
plot(common_wide_water$dev_waterVI, common_wide$dev_watersierra) # semi similar result with wider distribution, stuff that is overrepresented is overrepresented in both? 

#### Plotting Taxa that are Overrepresented in Trout #### 
phyla_to_include <- c("Firmicutes", "Fusobacteriota", "Proteobacteria", "Cyanobacteria", "Planctomycetota")

trout_asv_tax %>% 
  mutate(Phylum = ifelse(Phylum %in% phyla_to_include, Phylum, "Other")) %>% 
  ggplot(aes(x = factor(Phylum, levels = c("Cyanobacteria", "Firmicutes", "Fusobacteriota", "Planctomycetota","Proteobacteria", "Other")), y = dev)) +
  geom_point(aes(color = vals), alpha = 0.5) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  labs(y = "Deviance from Model Prediction", x = "Phylum") +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 13))
ggsave('figures/final_figures/PNM_deviance_trout.tiff', width = 11, height = 7)

#### Presence/Absence Unifrac Dist between within, above, below for trout, stickle, sierra, and vancouver ####
# df will look like the entire ASV table with 12 samples - above, within, and below for each with ASVs as either present or absent 
# how I'm going to do this is uncertain - but I have an idea to try to it 
# take mean relative abundances of ASVs 

# end result - a dataframe with columns for each community (neutral_trout, above_trout, below_trout, etc.) and every asv across all of these things, we'll start with presence/absence then we can try 

# new otu table for pcoa
#PNM_asvs <- data.frame(ASV = colnames(data.frame(otu_table(t(physeq_meta)))))
#rownames(PNM_asvs) <- colnames(data.frame(otu_table(t(physeq_meta))))
#saveRDS(PNM_asvs, "data/dada2/PNM_asvs.RDS")
PNM_asvs <- readRDS("data/dada2/PNM_asvs.RDS") # loading ASV list  

# adding trout to new presence/absence otu table
neutral_trout_asvs <- sncm_trout_freq[sncm_trout_freq$vals == "Within Prediction",]$ASV
above_trout_asvs <- sncm_trout_freq[sncm_trout_freq$vals == "Above Prediction",]$ASV
below_trout_asvs <- sncm_trout_freq[sncm_trout_freq$vals == "Below Prediction",]$ASV

PNM_asvs$neutral_trout <- ifelse(PNM_asvs$ASV %in% neutral_trout_asvs, 1, 0) # 1 is present, 0 is absent 
PNM_asvs$above_trout <- ifelse(PNM_asvs$ASV %in% above_trout_asvs,1,0)
PNM_asvs$below_trout <- ifelse(PNM_asvs$ASV %in% below_trout_asvs,1,0) 

# adding stickleback to new presence/absence otu table 
neutral_stickle_asvs <- sncm_stickle_freq[sncm_stickle_freq$vals == "Within Prediction",]$ASV
above_stickle_asvs <- sncm_stickle_freq[sncm_stickle_freq$vals == "Above Prediction",]$ASV
below_stickle_asvs <- sncm_stickle_freq[sncm_stickle_freq$vals == "Below Prediction",]$ASV

PNM_asvs$neutral_stickle <- ifelse(PNM_asvs$ASV %in% neutral_stickle_asvs, 1, 0) # 1 is present, 0 is absent 
PNM_asvs$above_stickle <- ifelse(PNM_asvs$ASV %in% above_stickle_asvs,1,0)
PNM_asvs$below_stickle <- ifelse(PNM_asvs$ASV %in% below_stickle_asvs,1,0)

# adding sierra water 
neutral_watersierra_asvs <- sncm_watersierra_freq[sncm_watersierra_freq$vals == "Within Prediction",]$ASV
above_watersierra_asvs <- sncm_watersierra_freq[sncm_watersierra_freq$vals == "Above Prediction",]$ASV
below_watersierra_asvs <- sncm_watersierra_freq[sncm_watersierra_freq$vals == "Below Prediction",]$ASV

PNM_asvs$neutral_watersierra <- ifelse(PNM_asvs$ASV %in% neutral_watersierra_asvs, 1, 0) # 1 is present, 0 is absent 
PNM_asvs$above_watersierra <- ifelse(PNM_asvs$ASV %in% above_watersierra_asvs,1,0)
PNM_asvs$below_watersierra <- ifelse(PNM_asvs$ASV %in% below_watersierra_asvs,1,0)

# adding vancouver water 
neutral_waterVI_asvs <- sncm_waterVI_freq[sncm_waterVI_freq$vals == "Within Prediction",]$ASV
above_waterVI_asvs <- sncm_waterVI_freq[sncm_waterVI_freq$vals == "Above Prediction",]$ASV
below_waterVI_asvs <- sncm_waterVI_freq[sncm_waterVI_freq$vals == "Below Prediction",]$ASV

PNM_asvs$neutral_waterVI <- ifelse(PNM_asvs$ASV %in% neutral_waterVI_asvs, 1, 0) # 1 is present, 0 is absent 
PNM_asvs$above_waterVI <- ifelse(PNM_asvs$ASV %in% above_waterVI_asvs,1,0)
PNM_asvs$below_waterVI <- ifelse(PNM_asvs$ASV %in% below_waterVI_asvs,1,0)

# remove ASV column for calculating Uni distance matrix 
PNM_asvs <- PNM_asvs[,-1]

# throw into a new phyloseq object for calculating unifrac distance matrix
asv_table <- otu_table(as.matrix(t(PNM_asvs)), taxa_are_rows = FALSE)
physeq_uni <- phyloseq(asv_table, phy_tree(physeq_meta))

uni_dist <- phyloseq::distance(physeq_uni, method = "unifrac")
uni_PNM <- ordinate(physeq_uni, method = "PCoA", distance = "unifrac")
pcoa_scores <- as.data.frame(uni_PNM$vectors[, c("Axis.1", "Axis.2")])
round(uni_PNM$values[1]*100/sum(uni_PNM$values[1]),2)[1:2,]

# create columns for visualizing 
pcoa_scores$sampleNames <- rownames(pcoa_scores)
sampleNames_parts <- strsplit(pcoa_scores$sampleNames, "_")
pcoa_scores$selection <- sapply(sampleNames_parts, "[", 1)
pcoa_scores$sample_type <- sapply(sampleNames_parts, "[", 2)
pcoa_scores$habitat_type <- ifelse(pcoa_scores$sample_type %in% c("watersierra", "waterVI"), "water", "fish")
pcoa_scores$system <- ifelse(pcoa_scores$sample_type %in% c("watersierra", "trout"), "Sierra", "Vancouver")
pcoa_scores$habitat_select <- paste(pcoa_scores$habitat_type, pcoa_scores$selection, sep = "_")

# shape as fish vs water, outside as system, fill as selection 
PNM_pcoa <- ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2,
                        shape = factor(habitat_type), 
                        fill = factor(selection), 
                        color = factor(system))) +
  geom_vline(xintercept = c(0), color = "grey70",linetype = 2)+
  geom_hline(yintercept = c(0), color = "grey70",linetype = 2)+
  geom_point(size = 4, stroke = 1) + 
  guides(shape = guide_legend(override.aes = list(size = 4)), color = guide_legend(override.aes = list(size = 4))) +
  scale_shape_manual(values = c(24, 21)) +
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  xlab("PCoA1 (21.6%)") + 
  ylab("PCoA2 (14.5%)") +
  theme_classic() + 
  theme(legend.title = element_blank())
PNM_pcoa
ggsave('figures/final_figures/figure6_ordination.pdf', width = 6, height = 5)

# PERMANOVA 
set.seed(1999)
PERM_sampletype <- adonis2(uni_dist ~ pcoa_scores$sample_type, permutations = 999)
PERM_sampletype # sample type doesn't explain variation 

set.seed(1999)
PERM_selection <- adonis2(uni_dist ~ pcoa_scores$selection, permutations = 999)
PERM_selection # 30% of variation, P = 0.002

set.seed(1999)
PERM_system <- adonis2(uni_dist ~ pcoa_scores$system, permutations = 999)
PERM_system # system doesn't matter! 

set.seed(1999)
PERM_habitattype <- adonis2(uni_dist ~ pcoa_scores$habitat_type, permutations = 999)
PERM_habitattype # marginally significant - 0.14 R2 with p = 0.029

set.seed(1999)
PERM_select_habitattype <- adonis2(uni_dist ~ pcoa_scores$selection*pcoa_scores$habitat_type, by = "terms", permutations = 999)
PERM_select_habitattype 

# pcoa_scores$selection                           2   1.3500 0.30562 2.5875  0.001 ***
#   pcoa_scores$habitat_type                        1   0.6087 0.13780 2.3334  0.003 ** 
#   pcoa_scores$selection:pcoa_scores$habitat_type  2   0.8933 0.20223 1.7122  0.004 ** 
#   Residual                                        6   1.5652 0.35435                  
# Total                                          11   4.4171 1.00000  



#### bar plots - How many taxa are neutral, non-neutral etc.? #### 
# plot - simple whats above, below, within for each species - color code by the neutral model plots 

# modify PNM model outputs by adding sample type column to each 
sncm_trout_freq$sample_type <- rep("Trout", length = nrow(sncm_trout_freq))
sncm_stickle_freq$sample_type <- rep("Stickleback", length = nrow(sncm_stickle_freq))
sncm_watersierra_freq$sample_type <- rep("Sierra Nevada Lakes", length = nrow(sncm_watersierra_freq))
sncm_waterVI_freq$sample_type <- rep("Vancouver Island Lakes", length = nrow(sncm_waterVI_freq))

# merge into one dataframe and modify
PNM_full_df <- bind_rows(sncm_trout_freq, sncm_watersierra_freq, sncm_stickle_freq, sncm_waterVI_freq)
PNM_full_df$sample_type <- factor(PNM_full_df$sample_type, levels = c("Trout", "Sierra Nevada Lakes", "Stickleback", "Vancouver Island Lakes"))
PNM_full_df$vals <-  factor(as.character(PNM_full_df$vals), levels = c("Below Prediction", "Above Prediction", "Within Prediction"))
PNM_full_df$neutrality <- factor(ifelse(PNM_full_df$vals == "Within Prediction", "Neutral", "Non-Neutral"), levels = c("Non-Neutral", "Neutral"))


# plot - above, below, within 
ggplot(PNM_full_df, aes(x = sample_type, fill = vals)) +
  geom_bar(position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = c("#94FFB2","#FFBBE4", "black")) + 
  labs(y = "", x = "")+
  theme_classic() + 
  theme(legend.title = element_blank())
ggsave("figures/final_figures/PNM_barplot.pdf", width = 7, height = 7)

# plot - neutral vs non-neutral 
ggplot(PNM_full_df, aes(x = sample_type, fill = neutrality)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = c("#E94B3CFF", "#2D2926FF")) + 
  labs(y = "", x = "")+
  theme_classic() + 
  theme(legend.title = element_blank())
ggsave("figures/final_figures/PNM_neutral_barplot.pdf", width = 7, height = 5)

#### what proportion of ASVs are above/below/within and neutral/non-neutral for each sample type? #####

# above vs within vs below 
PNM_full_df %>%
  group_by(sample_type, vals) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(proportion = count / sum(count))
# Groups:   sample_type [4]
# sample_type            vals              count proportion
# <fct>                  <fct>             <int>      <dbl>
#   1 Trout                  Below Prediction     40    0.00943
# 2 Trout                  Above Prediction   1022    0.241  
# 3 Trout                  Within Prediction  3180    0.750  
# 4 Sierra Nevada Lakes    Below Prediction     66    0.0193 
# 5 Sierra Nevada Lakes    Above Prediction    323    0.0947 
# 6 Sierra Nevada Lakes    Within Prediction  3022    0.886  
# 7 Stickleback            Below Prediction     78    0.00513
# 8 Stickleback            Above Prediction   2354    0.155  
# 9 Stickleback            Within Prediction 12764    0.840  
# 10 Vancouver Island Lakes Below Prediction    124    0.0295 
# 11 Vancouver Island Lakes Above Prediction    484    0.115  
# 12 Vancouver Island Lakes Within Prediction  3589    0.855  

# neutral vs non-neutral 
PNM_full_df %>%
  group_by(sample_type, neutrality) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(proportion = count / sum(count))

# sample_type            neutrality  count proportion
# <fct>                  <fct>       <int>      <dbl>
#   1 Trout                  Non-Neutral  1062      0.250
# 2 Trout                  Neutral      3180      0.750
# 3 Sierra Nevada Lakes    Non-Neutral   389      0.114
# 4 Sierra Nevada Lakes    Neutral      3022      0.886
# 5 Stickleback            Non-Neutral  2432      0.160
# 6 Stickleback            Neutral     12764      0.840
# 7 Vancouver Island Lakes Non-Neutral   608      0.145
# 8 Vancouver Island Lakes Neutral      3589      0.855
# put in results as descriptive 

# how many taxa are below 50% frequency across the four sample types? 
PNM_full_df$freq_perc <- factor(ifelse(PNM_full_df$freq < 0.5, "Low", "High"))
PNM_full_df %>% 
  group_by(sample_type, freq_perc) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_type) %>%
  mutate(proportion = count / sum(count))

#### Code for 7B - What Taxa Have High Absolute Deviance Across Different Phyla? ##### 
# make a barplot of deviance on Y axis, different phyla on the x axis 
# what taxa represent the above/below/within groups across each community? 

# modify dataframes to combine with taxa 
sncm_trout_freq$ASV <- rownames(sncm_trout_freq)
trout_asv_tax <- merge(sncm_trout_freq, data.frame(tax_table(physeq_trout)), by = 'row.names')

sncm_stickle_freq$ASV <- rownames(sncm_stickle_freq)
stickle_asv_tax <- merge(sncm_stickle_freq, data.frame(tax_table(physeq_stickle)), by = 'row.names')

sncm_watersierra_freq$ASV <- rownames(sncm_watersierra_freq)
watersierra_asv_tax <- merge(sncm_watersierra_freq, data.frame(tax_table(physeq_watersierra)), by = 'row.names')

sncm_waterVI_freq$ASV <- rownames(sncm_waterVI_freq)
waterVI_asv_tax <- merge(sncm_waterVI_freq, data.frame(tax_table(physeq_waterVI)), by = 'row.names')

# calculate residuals/deviance from model for each asv 
trout_asv_tax$dev <- trout_asv_tax$freq - trout_asv_tax$freq.pred
stickle_asv_tax$dev <- stickle_asv_tax$freq - stickle_asv_tax$freq.pred
watersierra_asv_tax$dev <- watersierra_asv_tax$freq - watersierra_asv_tax$freq.pred
waterVI_asv_tax$dev <- waterVI_asv_tax$freq - waterVI_asv_tax$freq.pred

# make new column that differentiates neutral vs nonneutral 
trout_asv_tax$neutrality <- ifelse(trout_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")
stickle_asv_tax$neutrality <- ifelse(stickle_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")
watersierra_asv_tax$neutrality <- ifelse(watersierra_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")
waterVI_asv_tax$neutrality <- ifelse(waterVI_asv_tax$vals == "Within Prediction", "Neutral", "Non-Neutral")

# start with taxa dataframes for each sample type - add sample type columns 
trout_asv_tax$sample_type <- rep("Brook Trout", length = nrow(trout_asv_tax))
watersierra_asv_tax$sample_type <- rep("Sierra Nevada Lakes", length = nrow(watersierra_asv_tax))
stickle_asv_tax$sample_type <- rep("Stickleback", length = nrow(stickle_asv_tax))
waterVI_asv_tax$sample_type <- rep("Vancouver Island Lakes", length = nrow(waterVI_asv_tax))

# combine into same dataframe 
tax_full_df <- bind_rows(trout_asv_tax, watersierra_asv_tax, waterVI_asv_tax, stickle_asv_tax)

# absolute deviance as another way to visualize 
tax_full_df$abs_dev <- abs(tax_full_df$dev)

# filter so we only have the phyla that are in the stacked barchart figure, the rest are turned into "other" 
phyla <- c("Actinobacteriota", "Armatimonadota", "Bacteroidota", "Cyanobacteria", "Desulfobacterota", "Firmicutes", "Fusobacteriota", "Planctomycetota", "Proteobacteria", "Spirochaetota", "Verrucomicrobiota")
tax_full_df$phyla <- factor(ifelse(tax_full_df$Phylum %in% phyla, tax_full_df$Phylum, "Other"), levels = c("Actinobacteriota","Armatimonadota","Bacteroidota","Cyanobacteria","Desulfobacterota","Firmicutes","Fusobacteriota","Planctomycetota", "Proteobacteria","Spirochaetota","Verrucomicrobiota", "Other"))


#### Plotting Fig. 7B #### 
# reorder sample types 
tax_full_df$sample_type <- factor(tax_full_df$sample_type, levels = c("Brook Trout", "Stickleback", "Sierra Nevada Lakes", "Vancouver Island Lakes"))

# plot - phyla on y axis and facet wrap by sample type 
ggplot(tax_full_df,aes(x = factor(phyla), y = dev)) +
  geom_point(aes(color = vals), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  facet_wrap(~sample_type, ncol = 1) +
  scale_color_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  labs(y = "Deviance from Model Prediction", x = "Phylum") +
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.position = "bottom",
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.spacing = unit(0.8, "lines"), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Adds borders to facets
        strip.background = element_blank(),
        strip.text = element_text(size = 11)) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5))
ggsave("figures/final_figures/PNM_taxa.pdf", width = 6.5, height = 7.5)

#### Fig 7C and 7D - Are deviances correlated between fish spp and lake systems? #### 
# create separate dataframes
fish_tax_df <- tax_full_df[tax_full_df$sample_type %in% c("Brook Trout", "Stickleback"),]
water_tax_df <- tax_full_df[!tax_full_df$sample_type %in% c("Brook Trout", "Stickleback"),]

asvs_in_both <- fish_tax_df %>%
  filter(sample_type %in% c("Brook Trout", "Stickleback")) %>%
  group_by(ASV) %>%
  filter(n_distinct(sample_type) == 2) %>%  # ASVs that appear in both types
  pull(ASV) %>% 
  unique()
fish_tax_df <- fish_tax_df %>%
  filter(ASV %in% asvs_in_both)

df_list <- split(tax_full_df, tax_full_df$sample_type)

# find common asvs 
fish_asvs <- intersect(df_list$`Brook Trout`$ASV, df_list$Stickleback$ASV)
water_asvs <-intersect(df_list$`Sierra Nevada Lakes`$ASV, df_list$'Vancouver Island Lakes'$ASV)

# combine devs based on ASVs 
fish_df <- data.frame(trout_dev = df_list$`Brook Trout`[df_list$`Brook Trout`$ASV %in% fish_asvs,]$dev, stickle_dev = df_list$Stickleback[df_list$Stickleback$ASV %in% fish_asvs,]$dev, phyla = df_list$`Brook Trout`[df_list$`Brook Trout`$ASV %in% fish_asvs,]$phyla)

water_df <- data.frame(sierra_dev = df_list$`Sierra Nevada Lakes`[df_list$`Sierra Nevada Lakes`$ASV %in% water_asvs,]$dev, vi_dev = df_list$'Vancouver Island Lakes'[df_list$'Vancouver Island Lakes'$ASV %in% water_asvs,]$dev, phyla = df_list$'Vancouver Island Lakes'[df_list$'Vancouver Island Lakes'$ASV %in% water_asvs,]$phyla)

# plot - color by phyla ? 
ggplot(fish_df, aes(x = trout_dev, y = stickle_dev)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "lightgreen", alpha = 0.3) +   # Top-right
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "lightgreen", alpha = 0.3) +  # Bottom-left
  geom_point(size = 3, shape = 17, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  stat_cor(method = "spearman", label.x = 0, label.y = 0.5, color = "black", size = 4) +
  theme_classic() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 10))
ggsave("figures/final_figures/deviance_corr_fish.pdf", width = 5, height = 5)
  
ggplot(water_df, aes(x = sierra_dev, y = vi_dev)) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "lightgreen", alpha = 0.3) +   # Top-right
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "lightgreen", alpha = 0.3) +  # Bottom-left
  geom_point(size = 3, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
  stat_cor(method = "spearman", label.x = 0, label.y = 0.5, color = "black", size = 4) +
  theme_classic() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 10))
ggsave("figures/final_figures/deviance_corr_lakes.pdf", width = 5, height = 5)

# What proportion of ASVs are parallel within fish and water? 
lapply(list(water_df, fish_df), function(x) {
  same_sign <- (x[[1]] > 0 & x[[2]] > 0) | (x[[1]] < 0 & x[[2]] < 0)
  sum(same_sign)/nrow(x)
})

# 87% of ASVs are parallel in water, 98% in fish 

##### What Taxa Are Correlated with PCOA1? #### 
str(pcoa_scores)
str(PNM_asvs) # combine with sncm outputs! 
# replace 1's and 0's with the mean relative abundance of each ASV 
# add taxa info of each ASV 
# agglomerate at different taxa levels, combine with PCoA1 and do a correlation of relative abundance across PcoA1 

freq_lookup <- list(
  trout = sncm_trout_freq,
  stickle = sncm_stickle_freq,
  watersierra = sncm_watersierra_freq,
  waterVI = sncm_waterVI_freq
)

cols_to_replace <- grep("trout|stickle|watersierra|waterVI", names(PNM_asvs), value = TRUE)
PNM_abund <- data.frame(matrix(0, nrow = nrow(PNM_asvs), ncol = length(cols_to_replace)))
asv_names <- rownames(PNM_asvs)
rownames(PNM_abund) <- asv_names
colnames(PNM_abund) <- cols_to_replace

# add mean rel abund to each observation where the ASV is present
for (col in cols_to_replace) {
  group <- names(freq_lookup)[sapply(names(freq_lookup), function(g) grepl(g, col))]
  
  if (length(group) == 0) next
  
  freq_df <- freq_lookup[[group]]
  freq_vec <- setNames(freq_df$p, freq_df$ASV)
  
  matched_freq <- freq_vec[asv_names]  # match freq to ASV rownames
  
  PNM_abund[[col]] <- ifelse(PNM_asvs[[col]] == 1, matched_freq, 0)
}

# add taxa
PNM_abund <- merge(PNM_abund, data.frame(tax_table(physeq_meta)), by = "row.names")

# calculate the means of each phyla across the 12 partitions 
phylum_means <- PNM_abund %>%
  pivot_longer(
    cols = matches("trout|stickle|watersierra|waterVI"),
    names_to = "sampleNames",
    values_to = "abundance"
  ) %>%
  group_by(sampleNames, Phylum) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(sampleNames, desc(mean_abundance))
phylum_means

phylum_pcoa_merged <- phylum_means %>%
  left_join(pcoa_scores, by = "sampleNames")

# do correlations between PcOA1 and the mean relative abundance of different Phyla 
cor_results <- phylum_pcoa_merged %>%
  group_by(Phylum) %>%
  summarise(
    cor_pcoa1 = cor(mean_abundance, Axis.1, method = "spearman"),
    p_pcoa1 = cor.test(mean_abundance, Axis.1, method = "spearman")$p.value,
    cor_pcoa2 = cor(mean_abundance, Axis.2, method = "spearman"),
    p_pcoa2 = cor.test(mean_abundance, Axis.2, method = "spearman")$p.value
  )
sig_cor1 <- cor_results[cor_results$p_pcoa1 < 0.05,]
sig_cor2 <- cor_results[cor_results$p_pcoa2 < 0.05,]
# desulfobacterota 

ggplot(phylum_pcoa_merged, aes(x = Axis.1, y = mean_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Phylum, scales = "free_y") +
  theme_minimal()

# log transformed 
phylum_pcoa_log <- phylum_pcoa_merged %>%
  mutate(log_p = log10(mean_abundance + 1e-6))
log_cor_results <- phylum_pcoa_log %>%
  group_by(Phylum) %>%
  summarise(
    cor_pcoa1 = cor(log_p, Axis.1, method = "spearman"),
    p_pcoa1 = cor.test(log_p, Axis.1, method = "spearman")$p.value, 
    cor_pcoa2 = cor(log_p, Axis.2, method = "spearman"),
    p_pcoa2 = cor.test(log_p, Axis.2, method = "spearman")$p.value
  )

sig_phyla<- log_cor_results[log_cor_results$p_pcoa1 < 0.05,]$Phylum
ggplot(phylum_pcoa_log %>% filter(Phylum %in% sig_phyla), 
       aes(x = Axis.1, y = log_p, shape = factor(habitat_type), 
           fill = factor(selection), color = factor(system))) +
  geom_point(size = 2) +
  stat_cor(aes(label = ..r.label..), method = "spearman") +
  geom_smooth(aes(x = Axis.1, y = log_p), method = "lm", se = FALSE, inherit.aes = FALSE, color = "gray50", linewidth = 0.7) +
  scale_shape_manual(values = c(24, 21)) +
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  facet_wrap(~Phylum) +
  theme_classic()
ggsave('figures/final_figures/PNM_PcOA1_corr.pdf', height = 10, width = 10)

sig_phyla2 <-log_cor_results[log_cor_results$p_pcoa2 < 0.05,]$Phylum
ggplot(phylum_pcoa_log %>% filter(Phylum %in% sig_phyla2), 
       aes(x = Axis.2, y = log_p, shape = factor(habitat_type), 
           fill = factor(selection), color = factor(system))) +
  geom_point(size = 2) +
  geom_smooth(aes(x = Axis.2, y = log_p), method = "lm", se = FALSE, inherit.aes = FALSE, color = "gray50", linewidth = 0.7) +
  scale_shape_manual(values = c(24, 21)) +
  scale_color_manual(values = c("dodgerblue", "#e16f00")) +
  scale_fill_manual(values = c("#FFBBE4", "black", "#94FFB2")) +
  facet_wrap(~Phylum) +
  theme_classic()
ggsave('figures/final_figures/PNM_PcOA2_corr.pdf', height = 8, width = 8)


##### T-tests - how different is mean of Phyla from Null Expectation ? (0) ####
# what taxa are over and under represented across each sample type?  
# How we accomplish that - a directional t-test examining whether mean deviance from model prediction is significantly higher or lower than 0 
# supplemental table - taxa with mean deviance, SD, and p-values from 0 for each sample type? - could get large 

# Test 1 - is it significantly different above or below 0? 
str(tax_full_df)

tax_filtered <- tax_full_df %>%
  group_by(sample_type, phyla) %>%
  filter(n() >= 3) %>%
  ungroup()

# Now run the full analysis
dev_stats <- tax_filtered %>%
  group_by(sample_type, phyla) %>%
  group_modify(~ {
    t_out <- t_test(.x, dev ~ 1, mu = 0)
    w_out <- wilcox_test(.x, dev ~ 1, mu = 0)
    eff_t <- cohens_d(.x, dev ~ 1)
    eff_w <- wilcox_effsize(.x, dev ~ 1)
    
    tibble(
      n = nrow(.x),
      mean_dev = mean(.x$dev, na.rm = TRUE),
      t_p = t_out$p,
      wilcox_p = w_out$p,
      cohen_d = eff_t$effsize,
      wilcox_r = eff_w$effsize
    )
  }) %>%
  ungroup() %>%
  mutate(
    t_fdr = p.adjust(t_p, method = "fdr"),
    wilcox_fdr = p.adjust(wilcox_p, method = "fdr")
  )

significant_stats <- dev_stats[dev_stats$wilcox_p < 0.05,]
saveRDS(significant_stats, "data/dada2/PNM_dev_t_tests.rds")
