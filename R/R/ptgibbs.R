# Functions to call anything from ptgibbs.jl

prepare_julia <- function() {
    # Find julia v1.0.2 binary
    julia <- JuliaCall::julia_setup()
    ver <- as.numeric(stringr::str_split(string = julia$VERSION, pattern = "\\.")[[1]][1])
    if(ver < 1) {
        stop("Julia version > 1.0 required for this package to run.")
    }

    # load ptgibbs module
    JuliaCall::julia_command("using ptgibbs")
    JuliaCall::julia_command("using StatsBase")
    JuliaCall::julia_command("using DataFrames")
    JuliaCall::julia_command("using LinearAlgebra")
    JuliaCall::julia_command("using Distributions")
}


run_mcmc <- function(dat, kappa0, mu0, Psi0, alpha, nstep, reduced_classes, tune_df = 100, nw = 1) {
    if(is.null(reduced_classes)) {
        stop("Object `reduced_classes` output from get_reduced_classes must be provided.")
    }
    if(length(alpha) != nrow(reduced_classes)) {
        stop("length(alpha) must be the same as nrow(reduced_classes). These both
             correspond to the cluster number.")
    }

    dat <- data.frame(dat)

    prepare_julia()

    # Compute some values from the inputs
    dm <- ncol(dat)
    n <- nrow(dat)
    nm <- length(alpha)

    # Assign inputs names in Julia
    JuliaCall::julia_assign("dat", dat)
    JuliaCall::julia_assign("kappa0", kappa0)
    JuliaCall::julia_assign("mu0", mu0)
    JuliaCall::julia_assign("Psi0", Psi0)
    JuliaCall::julia_command("hyp = (kappa0, mu0, Psi0);")
    JuliaCall::julia_assign("alpha", alpha)
    JuliaCall::julia_assign("reduced_classes", reduced_classes)

    # Conver floats to integers when appropriate
    JuliaCall::julia_assign("nw", nw)
    JuliaCall::julia_assign("nw", JuliaCall::julia_eval("Int64(nw)"))

    JuliaCall::julia_assign("dm", dm)
    JuliaCall::julia_assign("dm", JuliaCall::julia_eval("Int64(dm)"))

    JuliaCall::julia_assign("n", n)
    JuliaCall::julia_assign("n", JuliaCall::julia_eval("Int64(n)"))

    JuliaCall::julia_assign("nm", nm)
    JuliaCall::julia_assign("nm", JuliaCall::julia_eval("Int64(nm)"))

    JuliaCall::julia_assign("nstep", nstep)
    JuliaCall::julia_assign("nstep", JuliaCall::julia_eval("Int64(nstep)"))

    JuliaCall::julia_assign("tune_df", tune_df)
    JuliaCall::julia_assign("tune_df", JuliaCall::julia_eval("Int64(tune_df)"))

    # Generate starting values
    JuliaCall::julia_assign("param",
                            JuliaCall::julia_eval("Array{Tuple{Array{Dict{String,Array{Float64,N} where N},1},Array{Float64,1},Array{Int64,1}}}(undef, (nw, 1));"))
    for(i in 1:nw) {
        JuliaCall::julia_assign("dictionary", JuliaCall::julia_eval("Dict{String,Array{Float64,N} where N}[];"))
        for(m in 1:nm) {
            JuliaCall::julia_assign("m", m)
            JuliaCall::julia_assign("m", JuliaCall::julia_eval("Int64(m)"))

            JuliaCall::julia_assign("h", reduced_classes[m,])
            JuliaCall::julia_assign("Sigma", JuliaCall::julia_eval("rand_constrained_Wish(Psi0[:,:,m], kappa0[m], h);"))
            JuliaCall::julia_assign("mu", JuliaCall::julia_eval("rand_constrained_MVN(Sigma/kappa0[m], mu0[m,:], reduced_classes[m,:]);"))

            JuliaCall::julia_command("push!(dictionary, Dict(\"mu\" => mu, \"Sigma\" => Sigma));")
        }
        JuliaCall::julia_assign("initprop", JuliaCall::julia_eval("alpha;"))
        JuliaCall::julia_assign("z", JuliaCall::julia_eval("rand(1:nm, n);"))

        JuliaCall::julia_assign("i", i)
        JuliaCall::julia_command("param[i,1] = (dictionary, initprop, z);")
    }

    labels <- apply(reduced_classes+1, 1, function(X) paste0(X, collapse = ""))
    JuliaCall::julia_assign("labels", labels)

    JuliaCall::julia_eval("ptgibbs.run_mcmc(dat, param, hyp, alpha, nstep, labels; tune_df = tune_df);")
}


get_mu_chain <- function(chain, walker_num, cluster_num, tempered = FALSE) {
    JuliaCall::julia_assign("walker_num", walker_num)
    JuliaCall::julia_assign("walker_num", JuliaCall::julia_eval("Int64(walker_num)"))

    JuliaCall::julia_assign("cluster_num", cluster_num)
    JuliaCall::julia_assign("cluster_num", JuliaCall::julia_eval("Int64(cluster_num)"))

    JuliaCall::julia_assign("chain", chain[[1]])

    if(!tempered) {
        JuliaCall::julia_eval("ptgibbs.get_gibbs_mu_chain(chain, walker_num, cluster_num)")
    } else {
        JuliaCall::julia_eval("ptgibbs.get_mu_chain(chain, walker_num, cluster_num)")
    }
}

get_Sigma_chain <- function(chain, walker_num, cluster_num, tempered = FALSE) {
    JuliaCall::julia_assign("walker_num", walker_num)
    JuliaCall::julia_assign("walker_num", JuliaCall::julia_eval("Int64(walker_num)"))

    JuliaCall::julia_assign("cluster_num", cluster_num)
    JuliaCall::julia_assign("cluster_num", JuliaCall::julia_eval("Int64(cluster_num)"))

    JuliaCall::julia_assign("chain", chain[[1]])

    if(!tempered) {
        JuliaCall::julia_eval("ptgibbs.get_gibbs_Sigma_chain(chain, walker_num, cluster_num)")
    } else {
        JuliaCall::julia_eval("ptgibbs.get_Sigma_chain(chain, walker_num, cluster_num)")
    }
}

get_prop_chain <- function(chain, walker_num, tempered = FALSE) {
    JuliaCall::julia_assign("walker_num", walker_num)
    JuliaCall::julia_assign("walker_num", JuliaCall::julia_eval("Int64(walker_num)"))

    JuliaCall::julia_assign("chain", chain[[1]])

    if(!tempered) {
        JuliaCall::julia_eval("ptgibbs.get_gibbs_prop_chain(chain, walker_num)")
    } else {
        JuliaCall::julia_eval("ptgibbs.get_prop_chain(chain, walker_num)")
    }
}

get_z_chain <- function(chain, walker_num, tempered = FALSE) {
    JuliaCall::julia_assign("walker_num", walker_num)
    JuliaCall::julia_assign("walker_num", JuliaCall::julia_eval("Int64(walker_num)"))

    JuliaCall::julia_assign("chain", chain[[1]])

    if(!tempered) {
        JuliaCall::julia_eval("ptgibbs.get_gibbs_z_chain(chain, walker_num)")
    } else {
        JuliaCall::julia_eval("ptgibbs.get_z_chain(chain, walker_num)")
    }
}
