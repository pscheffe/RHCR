#pragma once
#include "common.h"
#include "BasicGraph.h"

class PDMPC 
{
public:

    PDMPC(const BasicGraph& G);

    typedef boost::dynamic_bitset<> bitset_t;
    bitset_t get_reachable_set(const State &initial_state, const int &window) const;
    vector<bitset_t> get_coupling_graph(const vector<State> &agent_states , const int &window) const;

    vector<bitset_t> prioritize(const vector<bitset_t> &coupling_graph) const;

    //debugging
    void print_reachable_set(bitset_t, int start_location) const;
    void print_coupling_graph(vector<bitset_t>) const;

private:
    const BasicGraph& G;
};