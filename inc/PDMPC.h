#pragma once
#include "common.h"
#include "BasicGraph.h"
#include <boost/optional.hpp>

class PDMPC 
{
public:

    PDMPC(const BasicGraph& G);

    typedef boost::dynamic_bitset<> bitset_t;
    bitset_t get_reachable_set(const State &initial_state, const int &window) const;
    vector<bitset_t> get_coupling_graph(const vector<State> &agent_states , const int &window) const;

    vector<bitset_t> prioritize(const vector<bitset_t> &coupling_graph) const;
    vector<vector<int>> weigh(const vector<bitset_t> &directed_coupling_graph, const vector<State> &agent_states) const;

    vector<bitset_t> group(const vector<vector<int>> &directed_weighted_coupling_graph, const int &max_num_cls) const;

    boost::optional<vector<int>> get_computation_levels(const vector<bitset_t> &directed_coupling_graph) const;

    //debugging
    void print_reachable_set(bitset_t, int start_location) const;
    void print_coupling_graph(vector<vector<int>> graph) const;
    void print_coupling_graph(vector<bitset_t> graph) const;

private:
    const BasicGraph& G;

    static vector<bitset_t> to_unweighted_graph(vector<vector<int>> graph);

};