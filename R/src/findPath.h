#ifndef _findPath_H
#define _findPath_H

#include <RcppArmadillo.h>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/dfs.h>
#include <lemon/adaptors.h>
#include "pathenumeration.h"
#include "prune_enhancer.h"

using namespace lemon;
using namespace std;
using namespace Rcpp;
using namespace RcppArmadillo;

void findPath(ListDigraph& gr,
              ListDigraph::Node& src,
              ListDigraph::Node& trg,
              PathEnumeration& enumeration,
              int d,
              int& num_paths,
              vector<int>& all_paths,
              ListDigraph::ArcMap<bool>& filter,
              ListDigraph::Node& curr_node,
              ListDigraph::NodeMap<int>& layer, // everything after here has been added for my new pruner
              std::string& filepath,
              int& len_filt_h,
              Rcpp::List& nonconsec, // everything after here has been added for my second pruner
              Rcpp::List& mus,
              arma::mat& labels,
              int n,
              int dist_tol,
              bool prune_bool = false,
              double emp_prop = 0
              )
{
    //std::vector<int> assoc;

    if(num_paths == 0)
    {
        { // define a scope
            // FIND INITIAL PATH
        	Dfs<ListDigraph> dfs(gr);

            // Run DFS on full graph
        	dfs.run(src, trg);
        	for(ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
        		if(dfs.reached(n)) {
        			enumeration.push_first(n);
        		}
        	}
            num_paths++;

            // Add path to all_paths
            // d+2 = enumeration.len() when the enumeration is full
            for(int i=0; i < d+2; i++) {
                all_paths.push_back(gr.id(enumeration[i]));
            }
            

            // THEN CHECK IF PRUNE CHECK IS VALID USING NONCONSEC
            // if keep, prune_bool = false; true, prune_bool = false;
            std::vector<int> assoc;
            assoc = cassociate(all_paths, filepath, len_filt_h); // get the association value from node number
            prune_bool = cprune_path(nonconsec, assoc); // ask if path should be pruned

            // If the path is valid, it may still have 0 prior prop,
            // so check if we should prune it a second time
            if(!prune_bool)
            {
                emp_prop = cprune_path2(assoc, mus, labels, d, n, dist_tol);

                // if emp_prop is 0, we prune the path
                if(emp_prop == 0) {
                    prune_bool = true;
                }
            }

            if(prune_bool)
            {
                for(int i=0; i < d+2; i++)
                {
                    all_paths.pop_back();
                }
            }

            // move_curr_node one back in the path
            curr_node = enumeration[d];

            // POP TARGET NODE
            enumeration.pop_last(); // automatically filters between curr_node and target and decrements curr_node's outArcs

            // Delete stuff I don't need anymore
            std::vector<int>().swap(assoc);
        }

        findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter,
                    curr_node, layer, filepath, len_filt_h, nonconsec, mus,
                    labels, n, dist_tol, prune_bool, emp_prop);
    } else // IF NUM_PATHS > 0
    {
        //vector<int> prune_check;
        //std::vector<int> assoc;

        // STOPPING RULE
        while(!(enumeration.len() > 0 && enumeration[0] == src && enumeration.outArcs(enumeration[0]) == 0))
        {
            // WHILE THE CURRENT NODE STILL HAS FEASIBLE OUTGOING PATHS
            while(enumeration.outArcs(curr_node) > 0)
            {
                {   // define a scope
                    // update filter (might want to put this in main and pass by reference)
                    FilterArcs<ListDigraph> subgraph(gr, filter);
                    Dfs<FilterArcs<ListDigraph> > sub_dfs(subgraph);
                    vector<int> temp;
                    int sz;

                    // find new path based on filter
                    sub_dfs.run(curr_node, trg);

                    // Add to all_paths
                    sz = enumeration.len();
                    for(int i = 0; i < sz; i++) {
                        all_paths.push_back(gr.id(enumeration[i]));
                        //prune_check.push_back(gr.id(enumeration[i]));
                    }
                    for(ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
                        if(sub_dfs.reached(n) && gr.id(n) != gr.id(curr_node)) {
                            temp.push_back(gr.id(n));
                        }
                    }
                    sort(temp.begin(), temp.end());

                    for(auto i = 0; i < temp.size(); i++)
                    {
                        all_paths.push_back(temp[i]);
                        //prune_check.push_back(temp[i]);

                        for(ListDigraph::NodeIt n(gr); n != INVALID; ++n)
                        {
                            if(gr.id(n) == temp[i])
                            {
                                enumeration.push_last(n);
                            }
                        }
                    }
                    
                    vector<int> prune_check;
                    std::vector<int> assoc;
                    prune_check.assign(all_paths.end()-d-2, all_paths.end());

                    // THEN CHECK IF PRUNE CHECK IS VALID USING NONCONSEC
                    // if keep, prune_bool = false; true, prune_bool = false;
                    assoc = cassociate(prune_check, filepath, len_filt_h); // get the association value from node number
                    prune_bool = cprune_path(nonconsec, assoc); // ask if path should be pruned

                    // If the path is valid, it may still have 0 prior prop,
                    // so check if we should prune it a second time
                    if(!prune_bool)
                    {
                        emp_prop = cprune_path2(assoc, mus, labels, d, n, dist_tol);

                        // if emp_prop is 0, we prune the path
                        if(emp_prop == 0)
                        {
                            prune_bool = true;
                        }
                    }

                    if(prune_bool)
                    {
                        for(int i=0; i < d+2; i++)
                        {
                            all_paths.pop_back();
                        }
                    }

                    enumeration.pop_last();
                    curr_node = enumeration[d];

                    for(ListDigraph::ArcIt a(gr); a != INVALID; ++a) {
                        filter[a] = enumeration.filter(a);
                    }

                    // delete stuff I don't need anymore
                    std::vector<int>().swap(prune_check);
                    std::vector<int>().swap(assoc);
                    std::vector<int>().swap(temp);
                }
                
                findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter,
                    curr_node, layer, filepath, len_filt_h, nonconsec, mus,
                    labels, n, dist_tol, prune_bool, emp_prop);

            }
            // ONLY MOVE BACK 1 IF WE ARE NOT ALREADY AT THE SOURCE NODE
            // OTHERWISE, JUST EXIT
            if(curr_node != src)
            {
                {   // define a scope
                    // move_curr_node one back in the path
                    curr_node = enumeration[layer[curr_node]-1];

                    // pop all nodes after curr_node
                    for(int i = layer[curr_node]; i < enumeration.len()-1; i++)
                    {
                        enumeration.pop_last();
                    }

                    // Unfilter all nodes, except the one that curr_node was connected to
                    for(ListDigraph::ArcIt a(gr); a != INVALID; ++a)
                    {
                        if(!(gr.source(a) == curr_node) && layer[gr.source(a)] >= layer[curr_node])
                        {
                            enumeration.true_filter(a);
                        }
                    }

                    // Reset outArcs for every node except curr_node
                    for(ListDigraph::NodeIt u(gr); u != INVALID; ++u)
                    {
                        if(u != curr_node && layer[u] >= layer[curr_node])
                        {
                            enumeration.reset_outArcs(u);
                        }
                    }

                    //   UPDATE FILTER IN SUBGRAPH
                    for(ListDigraph::ArcIt a(gr); a != INVALID; ++a) {
                        filter[a] = enumeration.filter(a);
                    }

                    // delete stuff I don't need anymore
                    //std::vector<int>().swap(prune_check);
                    //std::vector<int>().swap(assoc);
                }

                findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter,
                    curr_node, layer, filepath, len_filt_h, nonconsec, mus,
                    labels, n, dist_tol, prune_bool, emp_prop);
            }
        }
    }
}


#endif
