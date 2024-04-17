#pragma once
#include "common.h"
#include "BasicGraph.h"
#include <boost/optional.hpp>

using bitset_t = boost::dynamic_bitset<>;

using Graph = vector<bitset_t>;

template<typename T>
using WeightedGraph = vector<vector<T>>;

class PDMPC
{
public:

    explicit PDMPC(const BasicGraph& G);

    static boost::optional<vector<int>> get_computation_levels(const Graph &directed_coupling_graph);

    tuple<Graph, Graph, vector<bitset_t>> get_sequentially_planning_agents(const vector<State> &agent_states, const int &window, const int &max_num_CLs) const;

    //debugging
    void print_reachable_set(const bitset_t &reachable_set, const int &start_location=-1, const string &headline="reachable set") const;

    template<typename T>
    void print_coupling_graph(const WeightedGraph<T> &graph, const string &headline="weighted coupling graph") const;

    void print_coupling_graph(const Graph &graph, const string &headline="coupling graph") const;

private:
    const BasicGraph& G;

    bitset_t get_reachable_set(const State &initial_state, const int &window) const;

    pair<Graph, vector<bitset_t >> get_coupling_graph(const vector<State> &agent_states, const int &window) const;

    Graph prioritize(const Graph &coupling_graph) const;

    WeightedGraph<double> weigh(const Graph &directed_coupling_graph, const vector<State> &agent_states, const int &window) const;

    template<typename T>
    Graph group(const WeightedGraph<T> &directed_weighted_coupling_graph, const int &max_num_CLs) const;

    void print_headline(const string &headline) const;

};