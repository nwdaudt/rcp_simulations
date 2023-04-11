##
## RCP simulations
## 02_fit-RCPs
##
## Nicholas W Daudt
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##


## For using Skynet (Otago's Maths & Stats supercomputer), 
## check script 'skynet'

### Libraries & Functions ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ecomix)

## For Skynet using only
# setwd("/home/winni797/RCP_simulations/")

# Helper functions for plotting
source("./scripts/source.R")

### Load data ####
files <- dir("./matrices", full.names = TRUE)

for (file in files) {
  load(file)
}

rm("file", "files")

### Prep to fit RCPs ####

# Number of RCP groups to test
vec_rcps <- c(1:15)

# RCP formula (rcp_form) will be specified inside the loop

# Sampling artifact formula
spp_form <- as.formula(~1)

# Some controls, according to Skip Woolley (Stats Methods Seminar EFI):
# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/tree/main/woolley_ecomix
control <- list(penalty = 0.01, penalty.tau = 10, penalty.gamma = 10,
                penalty.disp = c(10, sqrt(10)), quiet = TRUE)

# Number of random starts to run each 'vec_rcps'
nstarts <- 100

### Fit NegBin RCPs ####

Lists <- list(gr2_nspp.matrices,
              gr3_nspp.matrices,
              gr4_nspp.matrices)

for (List in Lists) {
  
  for (i in 1:length(List)) {
    
    # Get name for saving results afterwards
    name <- names(List[i])
    print(name)
    
    # Get number of spp to specify 'rcp_form' below
    nspp <- stringr::str_sub(names(List[i]), start = 9)
    
    # Get the matrix data
    rcp_data <- List[[i]]
    
    # rcp_form = formula for fitting RCPs (i.e. "spp names ~ 1 + env vars")
    rcp_form <- as.formula(paste0(
      "cbind(", paste0("sp", 1:as.numeric(nspp), collapse = ","),")~1+temp_scaled+long_scaled"))
    
    ### Store results in a list
    nRCPs_samp <- list()
    
    ## Run 'regional_mix.multifit' ------------------------------------------- #
    for(ii in vec_rcps) {
      
      nRCPs_samp[[ii]] <- 
        ecomix::regional_mix.multifit(
          rcp_formula = rcp_form,
          species_formula = spp_form,
          data = rcp_data,
          nRCP = ii,
          family = "negative.binomial",
          inits = "random2",
          control = control,
          nstart = nstarts)
      
      rm("ii")
    }
    
    ## Get BIC values for each model ----------------------------------------- #
    RCPsamp_BICs <- sapply(nRCPs_samp, function(x) sapply(x, function(y) y$BIC))
    
    ## Get smallest BIC value for each RCP number ---------------------------- #
    RCPsamp_minBICs <- apply(RCPsamp_BICs, 2, min)
    
    ## Get True/False for <= N sites allocated for a specific RCP ------------ #
    L <- list()
    
    for (ii in 2:length(nRCPs_samp)) {
      
      x1 <- nRCPs_samp[[ii]]
      
      for (iii in 1:length(x1)) {
        x <- x1[[iii]]$postProbs
        L <- append(L, list(as.matrix(x)))
      }
      
      rm("ii", "iii", "x1", "x")
    }
    
    RCPsamp_BICs_less5 <- sapply(L, FUN = FUN_lessNsites, n = 5)
    RCPsamp_BICs_only1 <- sapply(L, FUN = FUN_lessNsites, n = 1)
    
    rm("L")
    
    ## Plot results ------------------------------------------------------------ #
    
    grps <- vec_rcps # needed for ggplot & df's below
    
    # data.frames for plotting
    df2a <- data.frame(grps = grps, bic = c(RCPsamp_minBICs))
    df2b <- data.frame(grps = rep(grps, each = 100), 
                       bic = c(as.numeric(unlist(RCPsamp_BICs))),
                       less_5 = as.factor(c(rep(NA, times = 100), RCPsamp_BICs_less5)), 
                       only_1 = as.factor(c(rep(NA, times = 100), RCPsamp_BICs_only1))) 
    
    df2b <- 
      df2b %>%
      tidyr::pivot_longer(cols = c("less_5", "only_1"),
                          names_to = "rcp_sites",
                          values_to = "value") %>%
      dplyr::mutate(rcp_sites = as.factor(rcp_sites))
    
    # Set colour palette
    pal <- c("#000000",                                      # Black
             colorRampPalette(brewer.pal(7, "YlOrRd"))(14))  # Yellow < Red
    
    # Set right labels for facets in the plot
    site_labs <- c("Five or less sites", "Only one site")
    names(site_labs) <- c("less_5", "only_1")
    
    ## Plot BIC vs no. RCPs
    gg_multifitBIC <- 
      # plot min BIC
      ggplot(df2a[df2a$grps > 1,], aes(x = grps, y = bic)) +
      geom_point() + geom_line() +
      # plot the hundred fitted models
      geom_point(data = df2b[df2b$grps > 1, ], 
                 aes(x = grps, y = bic, colour = value),
                 alpha = 0.6) +
      scale_x_continuous("Number of groups", 
                         labels = as.character(grps), breaks = grps) +
      scale_color_manual(values = pal,
                         name = "# RCPs") +
      # facet it
      facet_wrap(vars(rcp_sites), labeller = labeller(rcp_sites = site_labs)) +
      ylab("BIC") +
      theme_bw() + 
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.5, "cm"))
    
    ## Save plot
    ggsave(plot = gg_multifitBIC,
           filename = 
             paste0("./results/03_NegBin_", as.character(name), "_RCPmultifit-plot.png"),
           width = 26, height = 12, units = "cm", dpi = 300)
    
    ## Save raw 'multifit' results ------------------------------------------- #
    save("nRCPs_samp",
         file = paste0("./results/03_NegBin_", as.character(name), "_RCPmultifit.rda"))
    
    # Also save objects used to generate the plot, just in case we need some adjustment
    save("grps", "df2a", "df2b", "site_labs", "pal",
         file = paste0("./results/03_NegBin_", as.character(name), "_RCPmultifit-plot-objs.rda"))
    
    ## Clean environment ----------------------------------------------------- #
    rm("nRCPs_samp", "rcp_data", "rcp_form", 
       "RCPsamp_BICs", "RCPsamp_minBICs", "RCPsamp_BICs_less5", "RCPsamp_BICs_only1", 
       "grps", "df2a", "df2b", "pal", "gg_multifitBIC", "site_labs")
    
    gc()
  }
  
  rm("List", "i", "name", "nspp")
  gc()
}
