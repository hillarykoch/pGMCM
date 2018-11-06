# Make some useful functions to extract info from the trace
# THESE NEED TO UPDATE WITH MY BETTER GIBBS SAMPLER

mean_traceplot <- function(param.tr, class_num, dim, ...) {
    plot(
        purrr::map(param.tr, "NIW") %>%
            purrr::map(class_num) %>%
            purrr::map("mu") %>%
            purrr::map(dim) %>%
            unlist,
        type = 'l',
        xlab = "iteration",
        ylab = substitute(mu[.(class_num, dim)]),
        main = "",
        ...
    )
}

var_traceplot <- function(param.tr, class_num, dim, ...) {
    plot(
        purrr::map(param.tr, "NIW") %>%
            purrr::map(class_num) %>%
            purrr::map("Sigma") %>%
            purrr::map(function(X)
                X[dim, dim]) %>%
            unlist,
        type = 'l',
        xlab = "iteration",
        ylab = substitute(sigma[.(class_num, dim)]),
        main = "",
        ...
    )
}

rho_traceplot <- function(param.tr, class_num, dim1, dim2, ...) {
    plot(
        purrr::map(param.tr, "NIW") %>%
            purrr::map(class_num) %>%
            purrr::map("Sigma") %>%
            purrr::map(function(X)
                X[dim1, dim2]) %>%
            unlist,
        type = 'l',
        xlab = "iteration",
        ylab = substitute(sigma[.(class_num, dim1, dim2)]),
        main = "",
        ...
    )
}

mixprop_traceplot <- function(param.tr, class_num = "all", ...) {
    if (class_num == "all") {
        nclass <- length(param.tr[[1]]$mix_prop)
        plot(
            purrr::map(param.tr, "mix_prop") %>%
                purrr::map(1) %>%
                unlist,
            type = 'l',
            xlab = "iteration",
            ylab = expression(pi),
            ylim = c(0, 1),
            main = "",
            ...
        )
        for (i in 2:nclass) {
            lines(
                purrr::map(param.tr, "mix_prop") %>%
                    purrr::map(i) %>%
                    unlist,
                type = 'l',
                ylim = c(0, 1),
                col = i
            )
        }
    } else {
        plot(
            purrr::map(param.tr, "mix_prop") %>%
                purrr::map(class_num) %>%
                unlist,
            type = 'l',
            xlab = "iteration",
            ylab = substitute(pi[.(class_num)]),
            ylim = c(0, 1),
            main = "",
            ...
        )
    }
  } else {
    plot(purrr::map(param.tr, "mix_prop") %>%
      purrr::map(class_num) %>%
      unlist(),
    type = "l",
    xlab = "iteration",
    ylab = substitute(pi[.(class_num)]),
    ylim = c(0, 1),
    main = "",
    ...
    )
  }
}

mean_post_density_plot <- function(param.tr, class_num, dim, ...) {
    plot(
        density(
            purrr::map(param.tr, "NIW") %>%
                purrr::map(class_num) %>%
                purrr::map("mu") %>%
                purrr::map(dim) %>%
                unlist
        ),
        xlab = substitute(mu[.(class_num, dim)]),
        ylab = "density",
        main = "",
        ...
    )
}

var_post_density_plot <- function(param.tr, class_num, dim, ...) {
    plot(
        density(
            purrr::map(param.tr, "NIW") %>%
                purrr::map(class_num) %>%
                purrr::map("Sigma") %>%
                purrr::map(function(X)
                    X[dim, dim]) %>%
                unlist
        ),
        xlab = substitute(sigma[.(class_num, dim)]),
        ylab = "density",
        main = "",
        ...
    )
}

rho_post_density_plot <-
    function(param.tr, class_num, dim1, dim2, ...) {
        plot(
            density(
                purrr::map(param.tr, "NIW") %>%
                    purrr::map(class_num) %>%
                    purrr::map("Sigma") %>%
                    purrr::map(function(X)
                        X[dim1, dim2]) %>%
                    unlist
            ),
            xlab = substitute(sigma[.(class_num, dim1, dim2)]),
            ylab = "density",
            main = "",
            ...
        )
    }

mixprop_post_density_plot <- function(param.tr, class_num, ...) {
    plot(
        density(
            purrr::map(param.tr, "mix_prop") %>%
                purrr::map(class_num) %>%
                unlist
        ),
        xlab = substitute(pi[.(class_num)]),
        ylab = "density",
        main = "",
        ...
    )
}

label_traceplot <- function(param.tr, obs_num, ...) {
    ylim <- c(0, max(param.tr[[1]]$z))
    plot(
        purrr::map(param.tr, "z") %>%
            map(obs_num) %>%
            unlist,
        xlab = "iteration",
        ylab = "label",
        main = "",
        ylim = ylim,
        ...
    )
}

# Which cluster label is most frequent for each observation
get_consensus_labels <- function(param.tr) {
  zees <- purrr::map(param.tr, "z") %>%
    abind(along = 2) %>%
    apply(1, sort) %>%
    apply(2, rle)
  maxidx <- purrr::map(zees, "lengths") %>%
    purrr::map(which.max)
  vals <- purrr::map(zees, "values")
  purrr::map2(.x = vals, .y = maxidx, `[`) %>% unlist()
}

# Match labels so that direct comparisons of results to truth can be made
relabel <- function(consensus_labels, red_class, true_assoc) {
  reorders <- rep(NA, nrow(red_class))

    for (i in seq(nrow(red_class))) {
        for (j in seq(nrow(true_assoc))) {
            if (all(red_class[i, ] == true_assoc[j, ])) {
                reorders[i] <- j
            }
        }
    }
  }

    # Account for classes that are not actually part of the true association patterns
    naidx <- is.na(reorders)
    if (sum(naidx) > 0) {
        reorders[naidx] <- (j + 1):(j + sum(naidx))
    }

    relabels <- rep(NA, length(consensus_labels))
    for (i in seq(nrow(red_class))) {
        relabels[consensus_labels == i] <- reorders[i]
    }

  relabels
}
