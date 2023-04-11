##
## This is a 'mirrored' version of script '03_',
## for Skynet - Otago's Math & Stats supercomputer
##

## Using Skynet, you will need to redirect R to the right directory
setwd("/home/winni797/RCP_simulations/")

## Skynet doesn't use the 'root' symbol (full stop) prior to backslash
## and that's why I created this new script, only to accommodate this small changes

tictoc::tic()
for (List in Lists) {
  
  for (i in 1:length(List)) {
    
    name <- names(List[i])
    print(name)
    
    nspp <- stringr::str_sub(names(List[i]), start = 9)
    
    rcp_data <- List[[i]]
    
    rcp_form <- as.formula(paste0(
      "cbind(", paste0("sp", 1:as.numeric(nspp), collapse = ","),")~1+temp_scaled+long_scaled"))
    
    nRCPs_samp <- list()
    
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
    
    RCPsamp_BICs <- sapply(nRCPs_samp, function(x) sapply(x, function(y) y$BIC))
    
    RCPsamp_minBICs <- apply(RCPsamp_BICs, 2, min)
    
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
    
    grps <- vec_rcps 
    
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
    
    pal <- c("#000000", colorRampPalette(brewer.pal(7, "YlOrRd"))(14))  
    
    site_labs <- c("Five or less sites", "Only one site")
    names(site_labs) <- c("less_5", "only_1")
    
    gg_multifitBIC <- 
      ggplot(df2a[df2a$grps > 1,], aes(x = grps, y = bic)) +
      geom_point() + geom_line() +
      geom_point(data = df2b[df2b$grps > 1, ], 
                 aes(x = grps, y = bic, colour = value),
                 alpha = 0.6) +
      scale_x_continuous("Number of groups", 
                         labels = as.character(grps), breaks = grps) +
      scale_color_manual(values = pal,
                         name = "# RCPs") +
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
    
    ggsave(plot = gg_multifitBIC,
           filename = 
             paste0("/home/winni797/RCP_simulations/results/03_NegBin_", as.character(name), "_RCPmultifit-plot.png"),
           width = 26, height = 12, units = "cm", dpi = 300)
    
    save("nRCPs_samp",
         file = paste0("/home/winni797/RCP_simulations/results/03_NegBin_", as.character(name), "_RCPmultifit.rda"))
    save("grps", "df2a", "df2b", "site_labs", "pal",
         file = paste0("/home/winni797/RCP_simulations/results/03_NegBin_", as.character(name), "_RCPmultifit-plot-objs.rda"))
    
    rm("nRCPs_samp", "rcp_data", "rcp_form", 
       "RCPsamp_BICs", "RCPsamp_minBICs", "RCPsamp_BICs_less5", "RCPsamp_BICs_only1", 
       "grps", "df2a", "df2b", "pal", "gg_multifitBIC", "site_labs")
  }
}
tictoc::toc()
