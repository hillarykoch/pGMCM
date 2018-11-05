// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <lemon/list_graph.h>
#include <lemon/dfs.h>
#include <lemon/lgf_reader.h>
#include <lemon/adaptors.h>
#include <iostream>
#include <fstream>
#include <string>
#include <array>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace lemon;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Function to pass to .transform()
// [[Rcpp::export]]
int trans_func(double& x) {
    if(x < 0) {
        return(-1);
    } else if(x > 0) {
        return(1);
    } else {
        return(0);
    }
}

// Get names of an Rcpp::List
// [[Rcpp::export]]
std::vector<std::string> get_list_names(Rcpp::List L) {
    return L.names();
}

// C version of paste0
// [[Rcpp::export]]
std::string cpaste0(std::vector<std::string> str1) {
    std::string pasted;
    for(std::vector<std::string>::iterator it = str1.begin(); it != str1.end(); ++it) {
        pasted += *it;
    }
    
    return pasted;
}

// perform action just like R's str_split
// [[Rcpp::export]]
arma::mat cstr_split(std::vector<std::string> strings, std::string split) {
    int num_strings = strings.size();
    arma::mat dims(num_strings, 2, arma::fill::none);
    
    for(auto i = 0; i < num_strings; i++) {
        int num_substr = strings[i].length();
        std::vector<std::string> tmp1;
        std::vector<std::string> tmp2;
        bool found_split = FALSE;
        
        for(auto j=0; j < num_substr; j++) {
            if(found_split == FALSE) {
                if(strings[i].substr(j,1) != split) {
                    tmp1.push_back(strings[i].substr(j,1));
                } else {
                    found_split = TRUE;
                }
            } else {
                tmp2.push_back(strings[i].substr(j,1));
            }
        }
        
        // Still need to paste tmp1, tmp2 each into one individual int
        dims(i,0) = std::stoi(cpaste0(tmp1));
        dims(i,1) = std::stoi(cpaste0(tmp2));
    }
    
    return dims;
}

// Find which row in a matrix equals some vector
// [[Rcpp::export]]
arma::vec caccept(arma::mat x, arma::colvec y){
    int b = x.n_rows;
    arma::vec out(b, arma::fill::none);

    for(int i = 0; i < b; i++){
        bool vecmatch = arma::approx_equal(x.row(i), y.t(), "absdiff", 0.001);
        if(vecmatch) {
            out(i) = 0;
        } else {
            out(i) = 1;
        }
    }

    return out;
}

// Find which row in a matrix equals some vector
// [[Rcpp::export]]
arma::vec caccept2(arma::mat x, double y1, double y2){
    int b = x.n_rows;
    arma::vec out(b, arma::fill::none);
    arma::colvec y = { y1, y2 };
    
    for(int i = 0; i < b; i++){
        bool vecmatch = arma::approx_equal(x.row(i), y.t(), "absdiff", 0.001);
        if(vecmatch) {
            out(i) = 0;
        } else {
            out(i) = 1;
        }
    }

    return out;
}

// if the path needs to go, return true
//[[Rcpp::export]]
bool cprune_path(Rcpp::List nonconsec, std::vector<int> assoc) {
    bool out = false;
    int sz = nonconsec.size();

    std::vector<std::string> labs = get_list_names(nonconsec);
    arma::colvec keepers = arma::zeros<arma::colvec>(sz);
    arma::mat spl = cstr_split(labs, "_");
    arma::uvec idx;

    double pos1;
    double pos2;
    
    for(int i = 0; i < sz; i++)
    {   
        pos1 = spl(i,0)-1;
        pos2 = spl(i,1)-1;
        
        arma::mat testmat = as<arma::mat>(nonconsec[i]);
        arma::colvec m = caccept2(testmat, assoc[pos1], assoc[pos2]);
        
        idx = find(m == 0);
        if(!(idx.size() > 0))
        {
            out = true;
            return out;
        }
    }
    return out;
}

// test reading in a table
// [[Rcpp::export]]
std::vector<int> cassociate(std::vector<int> path,
                            std::string filepath,
                            int len_filt_h) {
    // path is the current putative path
    // filepath is location of node file
    // filt_h is only needed for its total length

    ifstream stream (filepath);
    std::string line;
    int sz = path.size()-2;
    int x, y, z;
    std::vector<int> c(len_filt_h + 2);
    std::vector<int> assoc(sz);

    unsigned int count = 0;
    while(std::getline(stream, line))
    {
        ++count;
        if (count > (len_filt_h + 3)) { break; }
        if(count < 2) { continue; }

        stream>>x>>y>>z;
        c[count-2] = z;
    }

    // this range is because we dont need dummy notes "src" and "trg"
    for(auto i = 0; i < sz; i++)
    {
        assoc[i] = c[path[i+1]];
    }

    stream.close(); // ok this might work

    return assoc;
}


// Prune one path during findPath, rather than after, based on empirical prior prop
// [[Rcpp::export]]
double cprune_path2(std::vector<int> assoc,
                            Rcpp::List mus,
                            arma::mat labels,
                            int d,
                            int n,
                            int dist_tol) {
    int n_pairs = mus.size();
    arma::mat tags(n, n_pairs, arma::fill::zeros);
    arma::colvec m;
    arma::uvec idx;
    arma::uvec tagidx;
    arma::colvec rSums(n);
    arma::uvec outidx;
    arma::colvec onevec;

    // convert std::vector to arma::vec for use in this function
    arma::colvec clsvec = arma::conv_to<arma::colvec>::from(assoc);

    // Get Combos
    arma::field<arma::mat> comb(n_pairs);
    for(int i = 0; i < n_pairs; i++) {
        comb(i) = Rcpp::as<arma::mat>(mus[i]);
        comb(i).transform([](double val) { return(trans_func(val)); } );
    }

    // Which fits model which pairs of dims
    arma::mat dims = cstr_split(get_list_names(mus), "_");

    // For each pair of dims
    for(int j = 0; j < n_pairs; j++) {
        m.set_size(comb(j).n_rows);
        idx.set_size(comb(j).n_rows);

        m = caccept(comb(j), { clsvec(dims(j,0)-1), clsvec(dims(j,1)-1) } );
        idx = find(m == 0);
        tagidx = find(labels.col(j) == (idx(0)+1));

        onevec.ones(tagidx.size());
        tags.elem(tagidx+(n*j)) = onevec;
    }

    rSums = sum(tags, 1);
    outidx = find(rSums >= (n_pairs - dist_tol));

    double emp_prop = static_cast<double>(outidx.size())/(double)(1.0*n);

    return emp_prop;
}