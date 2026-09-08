######################################################
# Paper: Conserved upper thermal limits and small   #
# safety margins in soil copiotrophic bacteria.     #

#######################################################
########## Phylogenetic correlations  #################
#######################################################
## Code: Maria Rebolleda-Gomez (mreboll1@uci.edu)
## Last updated: July 1st 2026
#######################################################

### Load libraries
library(brms)
library(data.table)
library(ape)
library(phytools)
library(tidyverse)

### Set parent directory 
setwd("PATH TO DIRECTORY/Paper_PhylogeneticThermalConstraints")

# Load files (tree and trait data)
tree <- read.tree("05_TCS_phylogenetic_models/03_Output/tpc_tree_kbase_relaxed.tre") 
thermal_traits <- fread("05_TCS_phylogenetic_models/03_Output/tpc_tree_metadata_auc.csv.gz")

# Filter data to match tree species
tree_metadata <- thermal_traits %>% filter(as.character(seqID) %in% tree$tip.label) %>% droplevels()

# Find node of Entero-Pseudo divergnece and date using data from Timetree
P_E <- findMRCA(tree, c("2", "264", "57"))
calib_df <- makeChronosCalib(tree, node = P_E, age.min = 1350, age.max = 1527)
chronogram <- chronos(tree, model = "clock", lambda = 1, calibration = calib_df)
plot(chronogram)

# Calculate correlation matrix phylogenetic
A_cor <- vcv(chronogram, corr = TRUE)

#Align species names and subset only the species we are using
tree_metadata$seqID <- as.character(tree_metadata$seqID)
all(tree_metadata$seqID %in% rownames(A_cor))   # should be TRUE
rownames(tree_metadata) <- tree_metadata$seqID
A_cor <- A_cor[rownames(tree_metadata), rownames(tree_metadata)]

# Model ()
b.2.bf <- brmsformula(mvbind(rmax, topt, ctmin, ctmax) ~ (1|u|gr(seqID, cov = C))) + set_rescor(TRUE)
b.2 <- brm(b.2.bf,
           data   = tree_metadata,
           family = gaussian(), 
           data2  = list(C = A_cor),
           cores  = 4, 
           chains = 4
           )

summary(b.2)

# calculate page lamdda for phylo conservation
l_rmax <- b.2 %>% as_tibble() %>% 
        dplyr::select(sigma_phy_rmax = sd_seqID__rmax_Intercept, sigma_res_rmax = sigma_rmax) %>% 
        mutate(lambda_rmax = sigma_phy_rmax^2/(sigma_phy_rmax^2 + sigma_res_rmax^2)) %>% 
        pull(lambda_rmax) %>% quantile(probs=c(0,0.25,0.5,0.75,1))

l_topt <- b.2 %>% as_tibble() %>% 
        dplyr::select(sigma_phy_topt = sd_seqID__topt_Intercept, sigma_res_topt = sigma_topt) %>% 
        mutate(lambda_topt = sigma_phy_topt^2/(sigma_phy_topt^2 + sigma_res_topt^2)) %>% 
        pull(lambda_topt) %>% quantile(probs=c(0,0.25,0.5,0.75,1))

l_ctmax <- b.2 %>% as_tibble() %>% 
        dplyr::select(sigma_phy_ctmax = sd_seqID__ctmax_Intercept, sigma_res_ctmax = sigma_ctmax) %>% 
        mutate(lambda_ctmax = sigma_phy_ctmax^2/(sigma_phy_ctmax^2 + sigma_res_ctmax^2)) %>% 
        pull(lambda_ctmax) %>% quantile(probs=c(0,0.25,0.5,0.75,1))

l_ctmin <- b.2 %>% as_tibble() %>% 
        dplyr::select(sigma_phy_ctmin = sd_seqID__ctmin_Intercept, sigma_res_ctmin = sigma_ctmin) %>% 
        mutate(lambda_ctmin = sigma_phy_ctmin^2/(sigma_phy_ctmin^2 + sigma_res_ctmin^2)) %>% 
        pull(lambda_ctmin) %>% quantile(probs=c(0,0.25,0.5,0.75,1))

plot(b.2) 

fit_summary <- summary(b.2)

# Residual correlations (rescor) — between-trait residual covariance
res <- fit_summary$rescor %>% as.data.frame   # quick look

# Phylogenetic correlations (from |u| random effects) for Supplementary Figure 11 (Fig. S11)
phil <- fit_summary$random %>% as.data.frame
phil <- phil[5:10,]

correlations <- data.frame(cor = c("rmax_topt", "rmax_ctmin", "topt_ctmin", "rmax_ctmax", "topt_ctmax", "ctmin_ctmax"), 
    est = c(res$Estimate, phil$seqID.Estimate), CI_l = c(res[,3], phil[,3]),
    CI_u = c(res[,4], phil[,4]), type = rep(c("residual", "phylogenetic"), each = 6))

pdf("05_TCS_phylogenetic_models/03_Output/Phylogenetic_correlations.pdf", 4, 5)
ggplot(data = correlations, aes(x = est, y = cor, color = type))+
    geom_point(shape = 21, size = 1)+
    geom_errorbarh(aes(xmin = CI_l, xmax = CI_u), height = 0.2)+
    ylab("Correlation traits")+
    xlab("Estimate")+
    geom_vline(xintercept = 0, linetype = 2)+
    scale_color_manual(values = c("#0391fe", "black") )+
    theme_bw()

dev.off()