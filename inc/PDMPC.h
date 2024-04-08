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

private:
    const BasicGraph& G;
};