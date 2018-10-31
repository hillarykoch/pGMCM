# Functions needed to process fconstr_pGMCM to then dump into an LGF file

# Get the 2-mers from the fits
get_h <- function(fits) {
  mus <- purrr::map(fits, "mu") %>%
    lapply(function(X)
      replace(X, X > 0, 1)) %>%
    lapply(function(X)
      replace(X, X < 0, -1))
  mus
}

# Find out which associations exist in each dimension
# (may be all, but if filterable, can make problem easier)
filter_h <- function(h, d) {
  spl <- strsplit(names(h), split = "_") %>%
    purrr::map(as.numeric)

  filt <- list()
  for (i in seq(d)) {
    lapply(seq_along(spl), function(X) {
      idx <- spl[[X]] == i
      if (any(idx)) {
        h[[X]][, idx]
      }
    }) %>%
      unlist() %>%
      unique() %>%
      sort() -> filt[[i]]
  }
  names(filt) <- seq(d)
  filt
}

# Extract which results correspond to consecutive dimensions (e.g. 1_2, 2_3, but not 1_3)
get_consecutive <- function(h, non_consec = FALSE) {
  nms <- names(h)
  consec <- strsplit(names(h), split = "_") %>%
    purrr::map(as.numeric) %>%
    map(diff) %>%
    `==`(1)
  if (!non_consec) {
    h[consec]
  } else {
    h[!consec]
  }
}

# Make LGF file to be read by LEMON DigraphReader
write_LGF <- function(h, d, path) {
  cat("Writing LGF file...")
  # Make node section of LGF file
  filt_h <- filter_h(h, d)
  dims <- map(filt_h, length)
  node <-
    data.frame(
      "label" = c(0, seq_along(unlist(filt_h)), 3 * d + 1),
      "dim" = c(0, lapply(seq(dims), function(X)
        rep(X, times = dims[[X]])) %>% unlist(), d + 1),
      "assoc" = c(0, unlist(filt_h), 0)
    )
  readr::write_tsv(data.frame("@nodes"),
    path = path,
    col_names = FALSE
  )
  readr::write_tsv(node,
    path = path,
    col_names = TRUE,
    append = TRUE
  )
  cat("\n", file = path, append = TRUE)

  # Make arc section of LGF file
  consec <- get_consecutive(h)

  src <- trg <- list()
  for (i in seq_along(consec)) {
    src[[i]] <-
      sapply(consec[[i]][, 1], function(X)
        node$label[node$assoc == X & node$dim == i])
    trg[[i]] <-
      sapply(consec[[i]][, 2], function(X)
        node$label[node$assoc == X & node$dim == i + 1])
  }

  sources <- dplyr::filter(node, dim == 1)$label
  targets <- dplyr::filter(node, dim == d)$label
  arcs <- data.frame(
    "source" = c(
      rep(0, length(sources)),
      abind(src),
      targets
    ),
    "target" = c(
      sources,
      abind(trg),
      rep(3 * d + 1, length(targets))
    )
  )
  readr::write_tsv(
    data.frame("@arcs"),
    path = path,
    col_names = FALSE,
    append = TRUE
  )
  cat("\t\t -\n", file = path, append = TRUE)
  readr::write_tsv(
    arcs,
    path = path,
    col_names = FALSE,
    append = TRUE,
    na = ""
  )

  # Make attributes section of LGF file to note the sources
  cat("\n", file = path, append = TRUE)
  readr::write_tsv(
    data.frame("@attributes"),
    path = path,
    col_names = FALSE,
    append = TRUE
  )
  attrib <- data.frame(
    "type" = c("source", "target"),
    "label" = c(0, 3 * d + 1)
  )
  readr::write_tsv(attrib,
    path = path,
    col_names = FALSE,
    append = TRUE
  )
  cat("done!\n")
}

# get paths using LEMON, then convert to latent association vectors
get_paths <- function(filepath, len_filt_h, nonconsec, mus, labels, n, dist_tol) {
  cat("Finding latent classes...")
  path_build <- cgetPaths(filepath = filepath,
                          len_filt_h,
                          nonconsec,
                          mus,
                          labels,
                          n,
                          dist_tol)
  cat("done!\n")
  path_build
}

# Convert node labels back to associations
associate <- function(paths, filepath, filt_h) {
  node <- read_table2(
    file = filepath,
    skip = 1,
    col_names = TRUE,
    n_max = length(unlist(filt_h)) + 2
  )

  assoc_mx <-
    matrix(data = node[paths + 1, ]$assoc, nrow = nrow(paths))
  assoc_mx <- assoc_mx[, 2:(ncol(assoc_mx) - 1)]
}

# # Prune paths that are discordant with "non-consecutive" pairwise estimates
# prune_paths <- function(h, assoc_mx) {
#     nonconsec <- get_consecutive(h, non_consec = TRUE)
#     labs <- names(nonconsec)
#     keepers <- matrix(0, nrow = nrow(assoc_mx), ncol = length(nonconsec))
#
#     for (i in seq_along(nonconsec)) {
#         pair <- strsplit(labs[i], split = "_") %>%
#             `[[` (1) %>%
#             as.numeric
#         keepers[,i] <- crowMatch(assoc_mx[, pair], nonconsec[[i]])[,1]
#     }
#
#     # If a row is ever not a keeper (contains at least one 0), remove it
#     prunes <- apply(keepers, MARGIN = 1, function(X) any(X == 0))
#     assoc_mx[!prunes,]
# }

# Put everything together in one function here, get_reduced_classes
get_reduced_classes <- function(fits, n, d, filepath = "lgf.txt", dist_tol = 0) {
  h <- get_h(fits)
  filt <- filter_h(h, d)
  write_LGF(h, d, filepath)
  paths <- get_paths(
    filepath,
    length(unlist(filt)),
    get_consecutive(h, non_consec = TRUE),
    map(fits, "mu"),
    map(fits, "cluster") %>% simplify2array(),
    n,
    dist_tol
  )
  assoc <- associate(paths, filepath, filt_h = filt) # Still need to associate one last time
  assoc
  # prune_paths(h, assoc)
}

# Among all reduced classes, get the indices which correspond to the truth
# For simulated data only
get_true_assoc_idx <- function(red_class, true_assoc) {
  true_assoc <- as.matrix(true_assoc)
  trueidx <- cget_true_assoc_idx(red_class, true_assoc)
  sort(trueidx)
}
