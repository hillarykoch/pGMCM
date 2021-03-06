# For checking with Valgrind
library(pGMCM)

# Load data from pGMCM
load("/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/scripts/sim_studies/10dim_20class/hidim_fits.rda")
load("/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/scripts/sim_studies/10dim_20class/hidim_sim.rda")
load("/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/scripts/sim_studies/10dim_20class/hidim_trueassoc.rda")


# Process with LEMON processors otherwise used to prepare the LGF file
d <- 8
n <- nrow(sim$data)
fits <- fits[c(1,2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,25,26,27,28,31,32,33,36,37,40)]
        #fits[c(1,2,3,4,5,6,10,11,12,13,14,18,19,20,21,25,26,27,31,32,36)]
filepath = "/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/scripts/sim_studies/10dim_20class/lgf.txt"
red_class <- get_reduced_classes(fits, d, filepath=filepath)

prior_prop <- get_prior_prop(red_class, fits, d, n, dist_tol = 0)

trueidx <- get_true_assoc_idx(red_class, true_assoc[,1:8])
pal <- get_pals(1)
plot(prior_prop, col = pal[(seq_along(prior_prop) %in% trueidx)+1], pch = 16)
