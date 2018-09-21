// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <lemon/list_graph.h>
#include <lemon/dfs.h>
#include <lemon/lgf_reader.h>
#include <lemon/adaptors.h>
#include "findPath.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace lemon;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// Collect output from my algorithm
// [[Rcpp::export]]
arma::mat cgetPaths(std::string filepath) {
    ListDigraph gr;
    ListDigraph::NodeMap<int> dim(gr);
    ListDigraph::NodeMap<int> label(gr);
    ListDigraph::NodeMap<int> assoc(gr);
    Dfs<ListDigraph> dfs(gr);
    ListDigraph::ArcMap<bool> filter(gr);
    ListDigraph::Node src;
    ListDigraph::Node trg;
    ListDigraph::Node curr_node;
    std::vector<int> all_paths; // instantiate a resizable vector to contain all paths
    std::fstream infile;
    infile.open(filepath);


    // Read in Directed Graph from lgf.txt
    // "attributes" source and target are declared to be nodes and named src and trg
    // dim gives which "layer" the given node lies in
    digraphReader(gr, infile)
                    .nodeMap("label", label)
                    .nodeMap("dim", dim)
                    .nodeMap("assoc", assoc)
                    .node("source", src)
                    .node("target", trg)
                    .run();

    ListDigraph::NodeMap<int> layer(gr);
    for(ListDigraph::NodeIt u(gr); u != INVALID; ++u) {
        layer[u] = dim[u];
    }

    int d = dim[trg]-1;

    // Enumerate all paths using DFS and backtracking
    int num_paths = 0; // when not 0, enter a different section of "findPath" function
    PathEnumeration enumeration(gr, src, trg);

    findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter, curr_node, layer);

	// need to write this function
	num_paths = (all_paths.size())/(d+2);
	arma::mat out(d+2, num_paths, arma::fill::none);
	for(int i = 0; i < all_paths.size(); i++) {
	    out(i) = all_paths[i];
	}
	return out.t();
}


// Prune paths using all other [(d choose 2) - d + 1] combinations
// [[Rcpp::export]]
arma::uvec crowMatch(arma::mat assoc, arma::mat nonconsec) {
    // do matching
    int n_assoc = assoc.n_rows;
    int n_valid = nonconsec.n_rows;

    // 1 if keep, 0 otherwise
    arma::uvec keepers(n_assoc, arma::fill::zeros);
    for(int i = 0; i < n_assoc; i++) {
        for(int j = 0; j < n_valid; j++) {
                if(assoc(i,0) == nonconsec(j,0) & assoc(i,1) == nonconsec(j,1)) {
                    keepers(i) = 1;
                    break;
                }
        }
    }

    return keepers;
}
