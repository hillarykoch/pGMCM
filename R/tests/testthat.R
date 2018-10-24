# For checking with Valgrind

library(pGMCM)

# Load data from pGMCM
data(fits)
data(sim)
data(true_association)

# Process with LEMON processors otherwise used to prepare the LGF file
d <- 3
n <- nrow(sim$data)
filepath = "/Users/hillarykoch/Box Sync/School/research - Qunhua/Project_2/scripts/lgf.txt"
#red_class <- get_reduced_classes(fits, d, filepath=filepath)
red_class <- get_reduced_classes(fits, n, d, filepath=filepath, dist_tol = 0)
prior_prop <- get_prior_prop(red_class, fits, d, n)

trueidx <- get_true_assoc_idx(red_class, true_assoc)
pal <- get_pals(1)
plot(prior_prop, col = pal[(seq_along(prior_prop) %in% trueidx)+1], pch = 16)
