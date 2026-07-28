# MCMCglmm: Genome size~ bioclimatic variables
# with phylogenetic correction
#Quick overview generated with Claude--anthropic
#load packages
library(MCMCglmm)   
library(ape)        
library(coda)
library(phytools)
library(dplyr)
#Read in data
tree<-read.nexus("data/drostrees.nex")
#plot(tree[1])
dat<-read.csv("data/climate_data_drosophila_Sept17.csv")
####read in climate data generated from 03.WorldClimdata-dros-code.R##
dros<-read.csv("data/output/drosophila_bioclim.csv")
#remove data which is missing values
dros<-na.omit(dros)
#remove integer column
dros$X<-NULL
#make species a factor
dros$Species<-as.factor(dros$Species)
#we want to pull all relevant data for species which have bioclimatic information
#tips to keep in trees
to.keep<-unique(dros$Species)

#generating table of data
dat.mat<-matrix(,nrow=length(unique(dros$Species)), ncol=21)
colnames(dat.mat)<-c("Species", "MbDNA", "BIO1", "BIO2", "BIO3",
                     "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9",
                     "BIO10", "BIO11", "BIO12", "BIO13", "BIO14","BIO15",
                     "BIO16", "BIO17", "BIO18", "BIO19")

i<-1
#loops through for species we have geographic data for
#will randomly select one location if there are multiple locations to get 
#bioclimatic data for
#makes a bit table
for(i in 1:length(unique(dros$Species))){
  dat.mat[i,1]<-as.character(unique(dros$Species)[[i]])
  dat.mat[i,2]<-dat$MbDNA_Female[dat$Species==as.character(unique(dros$Species)[[i]])]
  currentspecies<-levels(dros$Species)[i]
  foo<-subset(dros, dros$Species == currentspecies)
  if(length(foo$Species)>1){
    random_row_index <- sample(nrow(foo), 1)
    output.test <- as.matrix(foo[random_row_index, ])
    
  dat.mat[i,3]<-as.numeric(output.test[5]) #bio1
  dat.mat[i,4]<-as.numeric(output.test[6]) #bio2
  dat.mat[i,5]<-as.numeric(output.test[7]) #bio3
  dat.mat[i,6]<-as.numeric(output.test[8]) #bio4
  dat.mat[i,7]<-as.numeric(output.test[9]) #bio5
  dat.mat[i,8]<-as.numeric(output.test[10]) #bio6
  dat.mat[i,9]<-as.numeric(output.test[11]) #bio7
  dat.mat[i,10]<-as.numeric(output.test[12]) #bio8
  dat.mat[i,11]<-as.numeric(output.test[13]) #bio9
  dat.mat[i,12]<-as.numeric(output.test[14]) #bio10
  dat.mat[i,13]<-as.numeric(output.test[15]) #bio11
  dat.mat[i,14]<-as.numeric(output.test[16]) #bio12
  dat.mat[i,15]<-as.numeric(output.test[17]) #bio13
  dat.mat[i,16]<-as.numeric(output.test[18]) #bio14
  dat.mat[i,17]<-as.numeric(output.test[19]) #bio15
  dat.mat[i,18]<-as.numeric(output.test[20]) #bio16
  dat.mat[i,19]<-as.numeric(output.test[21]) #bio17
  dat.mat[i,20]<-as.numeric(output.test[22]) #bio18
  dat.mat[i,21]<-as.numeric(output.test[23]) #bio19
  }
  else{
    output.test<-as.matrix(foo)
    dat.mat[i,3]<-as.numeric(output.test[5]) #bio1
    dat.mat[i,4]<-as.numeric(output.test[6]) #bio2
    dat.mat[i,5]<-as.numeric(output.test[7]) #bio3
    dat.mat[i,6]<-as.numeric(output.test[8]) #bio4
    dat.mat[i,7]<-as.numeric(output.test[9]) #bio5
    dat.mat[i,8]<-as.numeric(output.test[10]) #bio6
    dat.mat[i,9]<-as.numeric(output.test[11]) #bio7
    dat.mat[i,10]<-as.numeric(output.test[12]) #bio8
    dat.mat[i,11]<-as.numeric(output.test[13]) #bio9
    dat.mat[i,12]<-as.numeric(output.test[14]) #bio10
    dat.mat[i,13]<-as.numeric(output.test[15]) #bio11
    dat.mat[i,14]<-as.numeric(output.test[16]) #bio12
    dat.mat[i,15]<-as.numeric(output.test[17]) #bio13
    dat.mat[i,16]<-as.numeric(output.test[18]) #bio14
    dat.mat[i,17]<-as.numeric(output.test[19]) #bio15
    dat.mat[i,18]<-as.numeric(output.test[20]) #bio16
    dat.mat[i,19]<-as.numeric(output.test[21]) #bio17
    dat.mat[i,20]<-as.numeric(output.test[22]) #bio18
    dat.mat[i,21]<-as.numeric(output.test[23]) #bio19
  }
}

dat.tab<-as.data.frame(dat.mat)
str(dat.tab)

dat.tab2<-dat.tab[!(dat.tab$Species %in% "Chymomyza_procnemis"),]
dat.tab2 <- dat.tab2 %>%
  dplyr::mutate(across(c(MbDNA, BIO1, BIO2, BIO3, BIO4, BIO5,
                  BIO6, BIO7, BIO8, BIO9, BIO10, BIO11,
                  BIO12, BIO13, BIO14, BIO15, BIO16,
                  BIO17, BIO18, BIO19), as.numeric))

str(dat.tab2)
#making an inverse covariance matrix
# MCMCglmm needs the INVERSE of the phylogenetic covariance (relatedness) matrix.
# IMPORTANT: species in your data frame must be a subset of tree tip labels.
# Any mismatch will cause an error or silently drop taxa.

tree.trim<-keep.tip(tree[[1]], as.character(to.keep))
#need to root the tree
root.tree.trim<-root(tree.trim, outgroup=c("Chymomyza_amoena", "Chymomyza_procnemis"), resolve.root = T)
#having both chymomyza messes up the analysis ability
rooted<-drop.tip(root.tree.trim, "Chymomyza_procnemis")
is.rooted(rooted)
plot(rooted)

Ainv <- inverseA(rooted, nodes = "TIPS", scale = TRUE)$Ainv

# Sanity check: confirm all species in data appear in the tree
stopifnot(all(dat.tab2$Species %in% rownames(Ainv)))

# Scale continuous predictors ------------------------------------------
# Scaling (mean = 0, SD = 1) is strongly recommended because:
#   a) It puts all predictors on a common scale for prior specification
#   b) It makes posterior estimates directly comparable (effect sizes)
#   c) It improves MCMC mixing (faster convergence)

dat.tab2$BIO1_s <-scale(dat.tab2$BIO1)[,1]
dat.tab2$BIO2_s <-scale(dat.tab2$BIO2)[,1]
dat.tab2$BIO3_s <-scale(dat.tab2$BIO3)[,1]
dat.tab2$BIO4_s <-scale(dat.tab2$BIO4)[,1]
dat.tab2$BIO5_s <-scale(dat.tab2$BIO5)[,1]
dat.tab2$BIO6_s <-scale(dat.tab2$BIO6)[,1]
dat.tab2$BIO7_s <-scale(dat.tab2$BIO7)[,1]
dat.tab2$BIO8_s <-scale(dat.tab2$BIO8)[,1]
dat.tab2$BIO9_s <-scale(dat.tab2$BIO9)[,1]
dat.tab2$BIO10_s <-scale(dat.tab2$BIO10)[,1]
dat.tab2$BIO11_s <-scale(dat.tab2$BIO11)[,1]
dat.tab2$BIO12_s <-scale(dat.tab2$BIO12)[,1]
dat.tab2$BIO13_s <-scale(dat.tab2$BIO13)[,1]
dat.tab2$BIO14_s <-scale(dat.tab2$BIO14)[,1]
dat.tab2$BIO15_s <-scale(dat.tab2$BIO15)[,1]
dat.tab2$BIO16_s <-scale(dat.tab2$BIO16)[,1]
dat.tab2$BIO17_s <-scale(dat.tab2$BIO17)[,1]
dat.tab2$BIO18_s <-scale(dat.tab2$BIO18)[,1]
dat.tab2$BIO19_s <-scale(dat.tab2$BIO19)[,1]

# Log-transform genome size: MbDNA is typically right-skewed; log
# transformation normalises it and is biologically meaningful.
dat.tab2$logMbDNA <- log(dat.tab2$MbDNA)

# --- 5. Specify priors -------------------------------------------------------
# MCMCglmm uses inverse-Wishart priors for variance components.
# For a univariate model you need:
#   G: prior for the PHYLOGENETIC (random) variance = variation due to phylogeny
#   R: prior for the RESIDUAL variance = unexplained variation
#
# The parameter-expanded prior below (V=1, nu=0.002) is weakly informative
# and is the standard "default" recommendation for animal models.
# It places very little prior information, letting the data speak.
#
# If you have strong prior knowledge about the expected variance (e.g. from
# published heritability estimates), you can incorporate that here.

prior_simple <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),  # phylogenetic variance
  R = list(R1 = list(V = 1, nu = 0.002))   # residual variance
)

# --- 6. MCMC settings --------------------------------------------------------
# nitt    = total number of MCMC iterations (chain length)
# burnin  = iterations discarded at the start (chain "warms up")
# thin    = thinning interval (keep every nth sample to reduce autocorrelation)
# Effective stored samples = (nitt - burnin) / thin
#
# Start with the settings below for exploration. For a manuscript:
#   - Run nitt = 1,000,000 or more
#   - Check effective sample size (ESS) > 200 for all parameters
#   - Check trace plots for stationarity (no trend, good mixing)

nitt_val   <- 50000000
burnin_val <- 50000
thin_val   <- 500
# Effective samples = (500000 - 50000) / 500 = 900

# --- 7. Fit the full model ---------------------------------------------------
# Formula explanation:
#   fixed:  logMbDNA ~ BIO1_s + BIO12_s + BIO4_s + BIO15_s
#           These are your bioclimatic predictors (fixed effects)
#   random: ~animal
#           "animal" is a RESERVED WORD in MCMCglmm that triggers the
#           phylogenetic covariance structure via ginverse = list(animal = Ainv)
#   ginverse: links the "animal" random effect to your Ainv matrix
#   family: "gaussian" for a continuous, approximately normal response

cat("Fitting full model (this may take a few minutes)...\n")

model_full <- MCMCglmm(
  fixed   = logMbDNA ~ BIO1_s + BIO2_s +BIO3_s+ BIO4_s+ BIO5_s+ BIO6_s +BIO7_s+ BIO8_s+ BIO9_s+ BIO10_s+ BIO11_s+ BIO12_s+ BIO13_s+ BIO14_s+ BIO15_s+ BIO16_s+ BIO17_s+ BIO18_s+ BIO19_s,
  #fixed   = logMbDNA ~ BIO5_s+BIO9_s+BIO19_s+BIO10_s,
  random  = ~Species,
  ginverse = list(Species = Ainv),
  family  = "gaussian",
  data    = dat.tab2,
  prior   = prior_simple,
  nitt    = nitt_val,
  burnin  = burnin_val,
  thin    = thin_val,
  verbose = FALSE   # set TRUE to watch chain progress
)
cat("Full model fitted.\n")

#fitting null model
cat("Fitting null model...\n")

model_null <- MCMCglmm(
  fixed    = logMbDNA ~ 1,
  random   = ~Species,
  ginverse = list(Species = Ainv),
  family   = "gaussian",
  data     = dat.tab2,
  prior    = prior_simple,
  nitt     = nitt_val,
  burnin   = burnin_val,
  thin     = thin_val,
  verbose  = FALSE
)

cat("Null model fitted.\n")
#
# --- 9a. Trace plots ----------------------------------------------------------
# Each trace plot shows the sampled value at each MCMC iteration.
# A GOOD trace looks like: horizontal, stationary "hairy caterpillar"
#   - No upward/downward trends
#   - Good mixing (rapid oscillation around a stable mean)
# A BAD trace looks like: trending, autocorrelated, or "sticky" (flat regions)
# If bad: increase nitt, thin, or reconsider your prior / model structure.

par(mfrow = c(3, 2))   # arrange plots in a grid

# Fixed effects (regression coefficients)
plot(model_full$Sol, main = "Fixed effects — trace & density")

# Variance components (phylogenetic and residual)
plot(model_full$VCV, main = "Variance components — trace & density")

par(mfrow = c(1, 1))


# --- 9b. Autocorrelation ------------------------------------------------------
# Low autocorrelation at lag 1 (< 0.1) means your thinning is adequate.
# High autocorrelation means increase 'thin' or run longer chains.

cat("\n--- Autocorrelation (fixed effects) ---\n")
print(autocorr.diag(model_full$Sol))

cat("\n--- Autocorrelation (variance components) ---\n")
print(autocorr.diag(model_full$VCV))


# --- 9c. Effective sample size (ESS) -----------------------------------------
# ESS measures how many *independent* samples you effectively have after
# accounting for autocorrelation. Rule of thumb: ESS > 200 for all parameters.
# If ESS is low, run a longer chain or thin more aggressively.

cat("\n--- Effective sample size (fixed effects) ---\n")
print(effectiveSize(model_full$Sol))

cat("\n--- Effective sample size (variance components) ---\n")
print(effectiveSize(model_full$VCV))


# --- 9d. Heidelberger-Welch stationarity test --------------------------------
# Tests whether the posterior distribution has reached stationarity.
# You want "passed" for all parameters under "Stationarity test" and
# "Halfwidth test". Failures suggest the chain hasn't converged.

cat("\n--- Heidelberger-Welch test (fixed effects) ---\n")
print(heidel.diag(model_full$Sol))

cat("\n--- Heidelberger-Welch test (variance components) ---\n")
print(heidel.diag(model_full$VCV))


# =============================================================================
# 10. RESULTS — FIXED EFFECTS
# =============================================================================
# summary() gives the posterior distribution of each fixed effect:
#   post.mean   = posterior mean estimate (your "effect size")
#   l-95% CI    = lower 95% credible interval
#   u-95% CI    = upper 95% credible interval
#   eff.samp    = effective sample size for this parameter
#   pMCMC       = Bayesian p-value
#
# HOW TO INTERPRET pMCMC:
#   pMCMC is the proportion of the posterior that has the OPPOSITE sign
#   to the mean estimate, multiplied by 2 (two-tailed).
#   A commonly used threshold is pMCMC < 0.05, analogous to a frequentist
#   p-value, but note this is a posterior probability — not a null
#   hypothesis significance test.
#
# PREFER interpreting the 95% credible interval:
#   If the 95% CI does NOT overlap zero → strong evidence for an effect
#   If the 95% CI overlaps zero         → weak or no evidence for an effect
#
# NOTE: because genome size is log-transformed, back-transform for
# biological interpretation: exp(estimate) gives a multiplicative factor.
# e.g. estimate = 0.15 for BIO1_s → a 1 SD increase in temperature is
# associated with exp(0.15) = 1.16× larger genome size.

cat("\n=== FIXED EFFECTS SUMMARY ===\n")
print(summary(model_full)$solutions)


# --- Manual posterior summaries (more flexible) -----
cat("\n--- Posterior means and 95% credible intervals ---\n")
sol <- model_full$Sol
for (param in colnames(sol)) {
  cat(sprintf("%-20s  mean = %6.3f  95%% CI [%6.3f, %6.3f]  pMCMC = %.4f\n",
              param,
              mean(sol[, param]),
              quantile(sol[, param], 0.025),
              quantile(sol[, param], 0.975),
              2 * min(mean(sol[, param] > 0), mean(sol[, param] < 0))
  ))
}


