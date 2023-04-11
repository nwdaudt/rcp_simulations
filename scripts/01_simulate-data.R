##
## RCP simulations
## 01_simulate-data
##
## Nicholas W Daudt
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

### Libraries for plotting ####
library(ggplot2)
library(ggExtra)
library(RColorBrewer)

### Parameters for species matrix & 'environmental' data ####

## Values for Poisson lambda will be generated from a log-normal distribution 

# "Common species" with mean of -0.25 and sigma of 0.1 -- i.e.
# curve(dlnorm(x,-0.25, 0.1), from = 0.1, to = 1.2, n = 1001)

# "Rare species" with mean of -1.5 and sigma of 0.1 -- i.e.
# curve(dlnorm(x,-1.5, 0.2), from = 0.1, to = 1.2, n = 1001)

## Number of sites ('samples')
nsite <- 500

## For reproducibility
set.seed(9)

### The simulation design ####

df_simulation_design <- data.frame(
  nspp10 = c("3_4","3_1","2_2"),
  nspp50 = c("20_10","12_14","10_10"),
  nspp100 = c("35_30","25_25","20_20")
)

rownames(df_simulation_design) <- c("gr2", "gr3", "gr4")

# write.csv(df_simulation_design, "./results/00_simulation_design.csv")

rm("df_simulation_design")

### Simulate 'environmental variables' ####

temperature <- rnorm(n = nsite, mean = 20, sd = 4)
longitude <- rnorm(n = nsite, mean = 155, sd = 5)

## Standardize 'environmental variables'
temp_norm <- temperature - mean(temperature)
long_norm <- longitude - mean(longitude)

# ---------------------------------------------------------------------------- #
# env <- data.frame(
#   temperature = temperature,
#   longitude = longitude
# )
# 
# env_plot <-
#   ggplot(env, aes(x = longitude, y = temperature, color = temperature)) +
#   geom_point(size = 0.5) +
#   scale_color_viridis_c(guide = NULL) +
#   xlab("'Longitude'") + ylab("'Temperature'") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7))
# 
# env_plot <- ggMarginal(env_plot, 
#                        type = "densigram",
#                        xparams = list(fill = "lightgrey", size = 0.3),
#                        yparams = list(fill = "lightgrey", size = 0.3))
# 
# ggsave(plot = env_plot,
#        filename = "00_simulated-environmental-data.png",
#        path = "./results/",
#        width = 5, height = 5, units = "cm")
# 
# rm("env", "env_plot")
# ---------------------------------------------------------------------------- #

### Simulate 2 'true' groups ---------------------------------------------- ####

logit12 = 0 + 0.7*temp_norm + 0.6*long_norm   # ('high temperature', 'offshore')
denom12 = 1 + exp(logit12)

## Get probabilities for each group ("1", "2")
gr2_prob.grp1 = exp(logit12) / denom12
gr2_prob.grp2 = 1 / denom12

rm("logit12", "denom12")

# ------------------ Check how the probs/group looks like -------------------- #
# df <-
#   as.data.frame(
#   rbind(cbind(temperature, longitude, prob = gr2_prob.grp1, grp = 1),
#         cbind(temperature, longitude, prob = gr2_prob.grp2, grp = 2))) |>
#   tidyr::pivot_longer(cols = c(temperature, longitude),
#                       names_to = "vars",
#                       values_to = "values")
# 
# df$grp <- factor(df$grp, levels = c("1", "2"))
# 
# plot_env_gr2 <-
#   ggplot(df, aes(x = values, y = prob, color = grp)) +
#   geom_point(size = 0.35, shape = 1, alpha = 0.5) +
#   scale_color_brewer(palette = "Dark2", name = "Group") +
#   ylab("Probability") + xlab("Environmental gradient") +
#   facet_grid(~ vars, scales = "free_x") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(size = 8),
#         legend.position = "right",
#         legend.title = element_text(size = 7),
#         panel.grid = element_blank())
# 
# ggsave(plot = plot_env_gr2,
#        filename = "01_multinomial-probs_gr2.png",
#        path = "./results/",
#        width = 13, height = 6, units = "cm")
# 
# rm("df", "plot_env_gr2")
# ---------------------------------------------------------------------------- #

### Probability matrix
gr2_prob.matrix <- cbind(gr2_prob.grp1, gr2_prob.grp2)

### Lambdas for different groups to be used as Poisson parameters
# nspp 10 = 3 exclusive spp, 3 low prob shared spp
gr2_lambda_nspp10 = list(lambda1 = c(rep.int(x = 0, times = 3), 
                                     rlnorm(4, -1.5, 0.2), 
                                     rlnorm(3, -0.25, 0.1)),
                         lambda2 = c(rlnorm(3, -0.25, 0.1), 
                                     rlnorm(4, -1.5, 0.2),
                                     rep.int(x = 0, times = 3)))

# nspp 50 = 20 exclusive spp, 10 low prob shared spp
gr2_lambda_nspp50 = list(lambda1 = c(rep.int(x = 0, times = 20), 
                                     rlnorm(10, -1.5, 0.2), 
                                     rlnorm(20, -0.25, 0.1)),
                         lambda2 = c(rlnorm(20, -0.25, 0.1), 
                                     rlnorm(10, -1.5, 0.2), 
                                     rep.int(x = 0, times = 20)))

# nspp 100 = 35 exclusive spp, 15 low prob shared spp
gr2_lambda_nspp100 = list(lambda1 = c(rep.int(x = 0, times = 35),
                                      rlnorm(30, -1.5, 0.2), 
                                      rlnorm(35, -0.25, 0.1)),
                          lambda2 = c(rlnorm(35, -0.25, 0.1), 
                                      rlnorm(30, -1.5, 0.2),
                                      rep.int(x = 0, times = 35)))

gr2_lambda <- list(nspp10 = gr2_lambda_nspp10,
                   nspp50 = gr2_lambda_nspp50,
                   nspp100 = gr2_lambda_nspp100)

rm("gr2_lambda_nspp10", "gr2_lambda_nspp50", "gr2_lambda_nspp100")

### Create species matrix and merge with normalized 'environmental' data
### to get final 'data.frame'

gr2_nspp.matrices <- list()

for (nspp in 1:length(gr2_lambda)) {
  
  nspec <- length(gr2_lambda[[nspp]]$lambda1)
  
  lambda_list <- gr2_lambda[[nspp]]
  
  # Assign a 'true' group for each 'nsite' inside the loop
  group <- rep(NA, nsite)
  
  # Create species matrix based on group and lambdas
  gr_y <- matrix(NA, nrow = nsite, ncol = nspec)
  
  # Loop through 
  for(i in 1:nsite){ 
    
    # Assign a group based on its probability
    group[i] <- sample(x = 1:2, size = 1, prob = gr2_prob.matrix[i, ])
    
    # Given the assigned group, use the according lambda
    gr_y[i, ] <- rpois(n = nspec, lambda = lambda_list[[group[i]]])
  } 
  
  ## Final dataset: Species, Environment, True group
  gr_y <- as.data.frame(gr_y) # It needs to be "df" for RCPs
  
  gr_y <- cbind(gr_y, temp_norm, long_norm, group)
  colnames(gr_y) <- c(paste0("sp", rep(1:nspec)), 
                       "temp_scaled", "long_scaled", "true_group")
  
  ## Name and save it
  name <- paste0("gr2_", as.character(names(gr2_lambda[nspp])))
  
  gr2_nspp.matrices[[name]] <- gr_y
  
  rm("nspp", "nspec", "lambda_list", "group", "gr_y", "i", "name")
}

### Save matrices and clean environment
save("gr2_nspp.matrices", file = "./matrices/gr2_nspp.matrices.rda")

rm("gr2_prob.grp1", "gr2_prob.grp2", "gr2_prob.matrix", "gr2_lambda",
   "gr2_nspp.matrices")

### Simulate 3 'true' groups ---------------------------------------------- ####

logit13 = 0 + 0.7*temp_norm + 0.6*long_norm       # ('high temp', 'offshore')
logit23 = -0.5 - 0.7*temp_norm + 0.6*long_norm    # ('low temp', 'offshore')
denom123 = 1 + exp(logit13) + exp(logit23)

gr3_prob.grp1 = exp(logit13) / denom123
gr3_prob.grp2 = exp(logit23) / denom123
gr3_prob.grp3 = 1 / denom123

rm("logit13", "logit23", "denom123")

# ------------------ Check how the probs/group looks like -------------------- #
# df <-
#   as.data.frame(
#   rbind(cbind(temperature, longitude, prob = gr3_prob.grp1, grp = 1),
#         cbind(temperature, longitude, prob = gr3_prob.grp2, grp = 2),
#         cbind(temperature, longitude, prob = gr3_prob.grp3, grp = 3))) |>
#   tidyr::pivot_longer(cols = c(temperature, longitude),
#                       names_to = "vars",
#                       values_to = "values")
# 
# df$grp <- factor(df$grp, levels = c("1", "2", "3"))
# 
# plot_env_gr3 <-
#   ggplot(df, aes(x = values, y = prob, color = grp)) +
#   geom_point(size = 0.35, shape = 1, alpha = 0.5) +
#   scale_color_brewer(palette = "Dark2", name = "Group") +
#   ylab("Probability") + xlab("Environmental gradient") +
#   facet_grid(~ vars, scales = "free_x") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(size = 8),
#         legend.position = "right",
#         legend.title = element_text(size = 7),
#         panel.grid = element_blank())
# 
# ggsave(plot = plot_env_gr3,
#        filename = "01_multinomial-probs_gr3.png",
#        path = "./results/",
#        width = 13, height = 6, units = "cm")
# 
# rm("df", "plot_env_gr3")
# ---------------------------------------------------------------------------- #

### Probability matrix
gr3_prob.matrix <- cbind(gr3_prob.grp1, gr3_prob.grp2, gr3_prob.grp3)

### Lambdas for different groups to be used as Poisson parameters
# nspp 10 = 3 exclusive spp, 1 low prob shared spp
gr3_lambda_nspp10 = list(lambda1 = c(rlnorm(3, -0.25, 0.1), 
                                     0, 0, 0,
                                     0, 0, 0,
                                     rlnorm(1, -1.5, 0.2)),
                         lambda2 = c(0, 0, 0, 
                                     rlnorm(3, -0.25, 0.1),
                                     0, 0, 0,
                                     rlnorm(1, -1.5, 0.2)),
                         lambda3 = c(0, 0, 0, 
                                     0, 0, 0,
                                     rlnorm(3, -0.25, 0.1),
                                     rlnorm(1, -1.5, 0.2)))

# nspp 50 = 12 exclusive spp, 14 low prob shared spp
gr3_lambda_nspp50 = list(lambda1 = c(rlnorm(12, -0.25, 0.1),
                                     rep.int(x = 0, times = 24), 
                                     rlnorm(14, -1.5, 0.2)), 
                         lambda2 = c(rep.int(x = 0, times = 12), 
                                     rlnorm(12, -0.25, 0.1),
                                     rep.int(x = 0, times = 12),
                                     rlnorm(14, -1.5, 0.2)),
                         lambda3 = c(rep.int(x = 0, times = 24),
                                     rlnorm(12, -0.25, 0.1), 
                                     rlnorm(14, -1.5, 0.2)))

# nspp 100 = 25 exclusive spp, 25 low prob shared spp
gr3_lambda_nspp100 = list(lambda1 = c(rlnorm(25, -0.25, 0.1),
                                      rep.int(x = 0, times = 50), 
                                      rlnorm(25, -1.5, 0.2)), 
                          lambda2 = c(rep.int(x = 0, times = 25), 
                                      rlnorm(25, -0.25, 0.1),
                                      rep.int(x = 0, times = 25),
                                      rlnorm(25, -1.5, 0.2)),
                          lambda3 = c(rep.int(x = 0, times = 50),
                                      rlnorm(25, -0.25, 0.1), 
                                      rlnorm(25, -1.5, 0.2)))

gr3_lambda <- list(nspp10 = gr3_lambda_nspp10,
                   nspp50 = gr3_lambda_nspp50,
                   nspp100 = gr3_lambda_nspp100)

rm("gr3_lambda_nspp10", "gr3_lambda_nspp50", "gr3_lambda_nspp100")

### Create species matrix and merge with normalized 'environmental' data
### to get final 'data.frame'

gr3_nspp.matrices <- list()

for (nspp in 1:length(gr3_lambda)) {
  
  nspec <- length(gr3_lambda[[nspp]]$lambda1)
  
  lambda_list <- gr3_lambda[[nspp]]
  
  # Assign a 'true' group for each 'nsite' inside the loop
  group <- rep(NA, nsite)
  
  # Create species matrix based on group and lambdas
  gr_y <- matrix(NA, nrow = nsite, ncol = nspec)
  
  # Loop through 
  for(i in 1:nsite){ 
    
    # Assign a group based on its probability
    group[i] <- sample(x = 1:3, size = 1, prob = gr3_prob.matrix[i, ])
    
    # Given the assigned group, use the according lambda
    gr_y[i, ] <- rpois(n = nspec, lambda = lambda_list[[group[i]]])
  } 
  
  ## Final dataset: Species, Environment, True group
  gr_y <- as.data.frame(gr_y) # It needs to be "df" for RCPs
  
  gr_y <- cbind(gr_y, temp_norm, long_norm, group)
  colnames(gr_y) <- c(paste0("sp", rep(1:nspec)), 
                       "temp_scaled", "long_scaled", "true_group")
  
  ## Name and save it
  name <- paste0("gr3_", as.character(names(gr3_lambda[nspp])))
  
  gr3_nspp.matrices[[name]] <- gr_y
  
  rm("nspp", "nspec", "lambda_list", "group", "gr_y", "i", "name")
}

### Save matrices and clean environment
save("gr3_nspp.matrices", file = "./matrices/gr3_nspp.matrices.rda")

rm("gr3_prob.grp1", "gr3_prob.grp2", "gr3_prob.grp3", "gr3_prob.matrix", "gr3_lambda",
   "gr3_nspp.matrices")

### Simulate 4 'true' groups ---------------------------------------------- ####

logit14 = -1 + 0.7*temp_norm + 0.5*long_norm # ('high temp', 'offshore')
logit24 = -1 - 0.7*temp_norm + 0.5*long_norm # ('low temp', 'offshore')
logit34 = -1 + 0.2*temp_norm - 0.5*long_norm # ('medium temp', 'inshore')
denom1234 = 1 + exp(logit14) + exp(logit24) + exp(logit34)

gr4_prob.grp1 = exp(logit14) / denom1234
gr4_prob.grp2 = exp(logit24) / denom1234
gr4_prob.grp3 = exp(logit34) / denom1234
gr4_prob.grp4 = 1 / denom1234

rm("logit14", "logit24", "logit34", "denom1234")

# ------------------ Check how the probs/group looks like -------------------- #
# df <-
#   as.data.frame(
#   rbind(cbind(temperature, longitude, prob = gr4_prob.grp1, grp = 1),
#         cbind(temperature, longitude, prob = gr4_prob.grp2, grp = 2),
#         cbind(temperature, longitude, prob = gr4_prob.grp3, grp = 3),
#         cbind(temperature, longitude, prob = gr4_prob.grp4, grp = 4))) |>
#   tidyr::pivot_longer(cols = c(temperature, longitude),
#                       names_to = "vars",
#                       values_to = "values")
# 
# df$grp <- factor(df$grp, levels = c("1", "2", "3", "4"))
# 
# plot_env_gr4 <-
#   ggplot(df, aes(x = values, y = prob, color = grp)) +
#   geom_point(size = 0.35, shape = 1, alpha = 0.5) +
#   scale_color_brewer(palette = "Dark2", name = "Group") +
#   ylab("Probability") + xlab("Environmental gradient") +
#   facet_grid(~ vars, scales = "free_x") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 8),
#         axis.text = element_text(size = 7),
#         strip.text = element_text(size = 8),
#         legend.position = "right",
#         legend.title = element_text(size = 7),
#         panel.grid = element_blank())
# 
# ggsave(plot = plot_env_gr4,
#        filename = "01_multinomial-probs_gr4.png",
#        path = "./results/",
#        width = 13, height = 6, units = "cm")
# 
# rm("df", "plot_env_gr4")
# ---------------------------------------------------------------------------- #

### Probability matrix
gr4_prob.matrix <- cbind(gr4_prob.grp1, gr4_prob.grp2, gr4_prob.grp3, gr4_prob.grp4)

### Lambdas for different groups to be used as Poisson parameters
# nspp 10 = 2 exclusive spp, 2 low prob shared spp
gr4_lambda_nspp10 = list(lambda1 = c(rlnorm(2, -0.25, 0.1),
                                     0, 0, 
                                     0, 0,
                                     0, 0,
                                     rlnorm(2, -1.5, 0.2)),
                         lambda2 = c(0, 0, 
                                     rlnorm(2, -0.25, 0.1),
                                     0, 0,
                                     0, 0,
                                     rlnorm(2, -1.5, 0.2)),
                         lambda3 = c(0, 0, 
                                     0, 0,
                                     rlnorm(2, -0.25, 0.1),
                                     0, 0,
                                     rlnorm(2, -1.5, 0.2)),
                         lambda4 = c(0, 0, 
                                     0, 0,
                                     0, 0,
                                     rlnorm(2, -0.25, 0.1),
                                     rlnorm(2, -1.5, 0.2)))

# nspp 50 = 10 exclusive spp, 10 low prob shared spp
gr4_lambda_nspp50 = list(lambda1 = c(rlnorm(10, -0.25, 0.1),
                                     rep.int(x = 0, times = 30), 
                                     rlnorm(10, -1.5, 0.2)), 
                         lambda2 = c(rep.int(x = 0, times = 10), 
                                     rlnorm(10, -0.25, 0.1),
                                     rep.int(x = 0, times = 20),
                                     rlnorm(10, -1.5, 0.2)),
                         lambda3 = c(rep.int(x = 0, times = 20),
                                     rlnorm(10, -0.25, 0.1),
                                     rep.int(x = 0, times = 10),
                                     rlnorm(10, -1.5, 0.2)),
                         lambda4 = c(rep.int(x = 0, times = 30),
                                     rlnorm(10, -0.25, 0.1),
                                     rlnorm(10, -1.5, 0.2)))

# nspp 100 = 20 exclusive spp, 20 low prob shared spp
gr4_lambda_nspp100 = list(lambda1 = c(rlnorm(20, -0.25, 0.1),
                                      rep.int(x = 0, times = 60), 
                                      rlnorm(20, -1.5, 0.2)), 
                          lambda2 = c(rep.int(x = 0, times = 20), 
                                      rlnorm(20, -0.25, 0.1),
                                      rep.int(x = 0, times = 40),
                                      rlnorm(20, -1.5, 0.2)),
                          lambda3 = c(rep.int(x = 0, times = 40),
                                      rlnorm(20, -0.25, 0.1),
                                      rep.int(x = 0, times = 20),
                                      rlnorm(20, -1.5, 0.2)),
                          lambda4 = c(rep.int(x = 0, times = 60),
                                      rlnorm(20, -0.25, 0.1),
                                      rlnorm(20, -1.5, 0.2)))

gr4_lambda <- list(nspp10 = gr4_lambda_nspp10,
                   nspp50 = gr4_lambda_nspp50,
                   nspp100 = gr4_lambda_nspp100)

rm("gr4_lambda_nspp10", "gr4_lambda_nspp50", "gr4_lambda_nspp100")

### Create species matrix and merge with normalized 'environmental' data
### to get final 'data.frame'

gr4_nspp.matrices <- list()

for (nspp in 1:length(gr4_lambda)) {
  
  nspec <- length(gr4_lambda[[nspp]]$lambda1)
  
  lambda_list <- gr4_lambda[[nspp]]
  
  # Assign a 'true' group for each 'nsite' inside the loop
  group <- rep(NA, nsite)
  
  # Create species matrix based on group and lambdas
  gr_y <- matrix(NA, nrow = nsite, ncol = nspec)
  
  # Loop through 
  for(i in 1:nsite){ 
    
    # Assign a group based on its probability
    group[i] <- sample(x = 1:4, size = 1, prob = gr4_prob.matrix[i, ])
    
    # Given the assigned group, use the according lambda
    gr_y[i, ] <- rpois(n = nspec, lambda = lambda_list[[group[i]]])
  } 
  
  ## Final dataset: Species, Environment, True group
  gr_y <- as.data.frame(gr_y) # It needs to be "df" for RCPs
  
  gr_y <- cbind(gr_y, temp_norm, long_norm, group)
  colnames(gr_y) <- c(paste0("sp", rep(1:nspec)), 
                      "temp_scaled", "long_scaled", "true_group")
  
  ## Name and save it
  name <- paste0("gr4_", as.character(names(gr4_lambda[nspp])))
  
  gr4_nspp.matrices[[name]] <- gr_y
  
  rm("nspp", "nspec", "lambda_list", "group", "gr_y", "i", "name")
}

### Save matrices and clean environment
save("gr4_nspp.matrices", file = "./matrices/gr4_nspp.matrices.rda")

rm("gr4_prob.grp1", "gr4_prob.grp2", "gr4_prob.grp3", "gr4_prob.grp4",
   "gr4_prob.matrix", "gr4_lambda", "gr4_nspp.matrices")

rm(list = ls())

### END