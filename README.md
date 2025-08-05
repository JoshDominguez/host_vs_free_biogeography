## Scale-dependent differences in the biogeography of host-associated and free-living microbiomes across systems

We have a limited understanding of how host-associated and free-living microbes vary in their distribution, diversity, and composition. Here, we test whether host-associated (fish gut) and free-living (lake water) bacteria and archaea have different metacommunity structure and biogeography across two distinct lake ecosystems. 

# File Directory
The file directory contains folders housing data (csvs and R objects) and scripts (R) for analysis. 

# Scripts 
- 'diff_abund_final.R' = script for performing differential abundance analysis and taxa barplots for comparisons of host-associated and free-living microbiomes.
- 'dist_decay_final.R' = script for performing distance decay (spatial autocorrelation) analysis.
- 'diversity_final.R' = script for calculating and comparing alpha and beta-diversity metrics between host-associated and free-living microbiomes.
- 'ordinations_final.R' = script for calculating ordinations (Principle Coordinates Analysis) and PERMANOVA analysis to examine factors explaining variation in community structure.
- 'physeq_srs_normalize.R' = script for preprocessing of phyloseq object output downstream of DADA2 pipeline. Includes Scaling with Ranked Subsampling (SRS) to normalize read counts to account for differences in read depth among samples.
- 'PNM.R' = script for fitting the Prokaryote Neutral Model to sample categories. 

# Data 
- 'metacomm_TSCC.csv' = csv file containing all samples corresponding to fastq sequence files used in analysis.
- 'ps.tree.RDS' = RDS file of phyloseq object before filtering out unwanted samples and normalization (SRS).
- 'physeq_normal.RDS' = RDS file of phyloseq object after filtering and normalization.
- 'site_coords.csv' = csv file of coordinates of lakes sampled.
- 'site_data.csv' = csv file of coordinates and elevation data of lakes sampled.
- 'physeq_filt_nonnormal.RDS' = RDS file of phyloseq object after filtering out unwanted/unanalyzed samples but before SRS normalization. 
