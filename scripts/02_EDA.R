##
## RCP simulations
## 02_EDA
##
## Nicholas W Daudt
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##


### Libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(RColorBrewer)

### Load data ####
files <- dir("./matrices", full.names = TRUE)

for (file in files) {
  load(file)
}

rm("file", "files")

Lists <- list(gr2_nspp.matrices,
              gr3_nspp.matrices,
              gr4_nspp.matrices)

### EDA ------------------------------------------------------------------- ####

### 1. Number of observations by group/nspp
df_n_obs_by_grp_nspp <- data.frame()
### 2. ### Species counts (total number) by group/nspp
df_spp_counts_by_grp_nspp <- data.frame()

for (List in Lists) {
  
  for (i in 1:length(List)) {
    
    # Get name for saving results afterwards
    name <- names(List[i])
    print(name)
    
    df <- List[[i]]
    
    ## First plot
    df1 <- 
      df %>% 
      dplyr::group_by(true_group) %>% 
      dplyr::summarise(n = n())
    
    df1$grp <- stringr::str_sub(names(List[i]), end = 3)
    df1$nspp <- stringr::str_sub(names(List[i]), start = 9)
    
    df_n_obs_by_grp_nspp <- rbind(df_n_obs_by_grp_nspp, df1)
    
    ## Second plot
    # Get number of spp & create a vector with their names
    nspp <- stringr::str_sub(names(List[i]), start = 9)
    spp_cols <- paste0("sp", 1:as.numeric(nspp))
    
    df2 <- 
      df %>% 
      tidyr::pivot_longer(cols = all_of(spp_cols),
                          names_to = "spp",
                          values_to = "value")
    
    df2$spp <- factor(df2$spp, levels = spp_cols)
    df2$grp <- stringr::str_sub(names(List[i]), end = 3)
    df2$nspp <- stringr::str_sub(names(List[i]), start = 9)
    
    df_spp_counts_by_grp_nspp <- rbind(df_spp_counts_by_grp_nspp, df2)
    
    rm("i", "name", "df", "df1", 
       "df2", "nspp", "spp_cols")
  }
  rm("List")
}

df_n_obs_by_grp_nspp$true_group <- factor(df_n_obs_by_grp_nspp$true_group, 
                                      levels = c("1", "2", "3", "4"))

df_n_obs_by_grp_nspp$nspp <- factor(df_n_obs_by_grp_nspp$nspp, 
                                      levels = c("10", "50", "100"))

grp_labs <- c("2 groups", "3 groups", "4 groups")
names(grp_labs) <- c("gr2", "gr3", "gr4")

### 1 ---
gg_n_obs_by_grp_nspp <- 
  ggplot (df_n_obs_by_grp_nspp, aes(x = nspp, y = n, fill = true_group)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_brewer(palette = "Dark2", name = "Sites\ntrue group") + 
  facet_wrap(~grp, labeller = labeller(grp = grp_labs)) + 
  ylab("") + xlab("Number of simulated species") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"))

ggsave(gg_n_obs_by_grp_nspp,
       filename = "./results/02_EDA_n-obs-by-grp-nspp.png",
       height = 5, width = 12, unit = "cm", dpi = 300)

rm("df_n_obs_by_grp_nspp", "gg_n_obs_by_grp_nspp", "grp_labs")

### 2 ---

df_spp_counts_by_grp_nspp$true_group <- 
  factor(df_spp_counts_by_grp_nspp$true_group, 
         levels = c("1", "2", "3", "4"))
df_spp_counts_by_grp_nspp$nspp <- 
  factor(df_spp_counts_by_grp_nspp$nspp, 
         levels = c("10", "50", "100"))
df_spp_counts_by_grp_nspp$grp <- 
  factor(df_spp_counts_by_grp_nspp$grp, 
         levels = c("gr2", "gr3", "gr4"), 
         labels = c("2 groups", "3 groups", "4 groups")) 

gg_spp_counts_by_grp_nspp <-
  ggplot(df_spp_counts_by_grp_nspp, 
         aes(x = spp, y = value, fill = true_group)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_brewer(palette = "Dark2", name = "True group") +
  ggh4x::facet_nested(~grp + nspp) +
  xlab("Species") + ylab("") +
  coord_flip() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_blank(),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm"))

ggsave(gg_spp_counts_by_grp_nspp,
       filename = "./results/02_EDA_spp-counts-by-grp-nspp.png",
       height = 10, width = 18, unit = "cm", dpi = 300)

rm("df_spp_counts_by_grp_nspp", "gg_spp_counts_by_grp_nspp")

### END ?