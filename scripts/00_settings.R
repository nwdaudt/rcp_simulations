##
## RCP simulations
## 00_settings
##
## Nicholas W Daudt
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

### Set directories
dir.create("./scripts")
dir.create("./matrices")
dir.create("./results")

### Install needed libraries

# Install {pacman}, a wrapper for checking/installing packages
install.packages("pacman")

needed_libraries <- c("dplyr", "tidyr", "ggplot2", "ggExtra", 
                      "RColorBrewer", "devtools", "ggh4x")
pacman::p_install(needed_libraries)

# Will also need {ecomix} which currently is only available on GitHub (14 March 2023)
pacman::p_install_gh("skiptoniam/ecomix@dev")

# In case {ggh4x} don't install through CRAN
# pacman::p_install_gh("teunbrand/ggh4x")

### END