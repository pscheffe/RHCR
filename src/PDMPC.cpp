#include "PDMPC.h"

#include <random>

/**
 * @brief Constructs a new PDMPC object.
 * @param G A BasicGraph object.
 */
PDMPC::PDMPC(const BasicGraph &G)
    : G {G} { }

/**
* @brief Gets the sequentially planning agents.
* @param agent_states A vector of agent states.
* @param window An integer representing the window.
* @param max_num_CLs An integer representing the maximum number of computation levels.
* @return A tuple containing the directed coupling graph and the directed sequential
 * coupling graph and a vector of reachable set objects.
*/
tuple<Graph, Graph, vector<bitset_t>> PDMPC::get_sequentially_planning_agents(const vector<State> &agent_states, const int &window, const int &max_num_CLs) const {
    auto[coupling, reachable_sets] = this->get_coupling_graph(agent_states, window);

    Graph directed_coupling = this->prioritize(coupling);

    WeightedGraph<double> directed_weighted_coupling = this->weigh(directed_coupling, agent_states, window);

    Graph directed_sequential_coupling = this->group(directed_weighted_coupling, max_num_CLs);

    return make_tuple(directed_coupling, directed_sequential_coupling, reachable_sets);
}

/**
 *
 * @brief Gets the computation levels for a given directed coupling graph.
 * @description Kahn's algorithm is used to find the computation levels of the agents.
 * In each step we assign an increasing computation level to the to-be cut vertices
 * starting from level 1 until we have assigned a level to all vertices.
 * @param directed_coupling_graph A Graph object.
 * @return An optional vector of integers with the agents computation level at the agents index.
 */
boost::optional<vector<int>> PDMPC::get_computation_levels(const Graph &directed_coupling_graph) {
    auto num_agents = directed_coupling_graph.size();
    Graph A(directed_coupling_graph);

    //Kahn's algorithm
    vector<int> computation_levels(num_agents);
    list<int> no_incoming_edges;
    bitset_t assigned_agents_mask(num_agents);

    int current_cl = 1;
    while(!assigned_agents_mask.all()) {
        //find unassigned nodes/agents without incoming edges
        for (auto i = 0; i < num_agents; i++) {
            if(!assigned_agents_mask[i]) {
                bool has_incoming_edge = false;
                for (auto j = 0; j < num_agents; j++) {
                    if(A[j][i]) {
                        has_incoming_edge = true;
                    }
                }
                if (!has_incoming_edge) {
                    no_incoming_edges.emplace_back(i);
                }
            }
        }

        if (no_incoming_edges.empty()) {
            return boost::optional<vector<int>>{};
        }

        for(auto agent : no_incoming_edges) {
            computation_levels[agent] = current_cl;
            assigned_agents_mask[agent] = true;
            //clear rows
            for(int i = 0; i < num_agents; i++) {
                A[agent][i] = false;
            }
        }

        no_incoming_edges.clear();
        current_cl = current_cl + 1;
    }

    return computation_levels;
}

/**
 * @brief Gets the reachable set for a given initial state and window.
 * @description The reachable set is calculated by performing a breadth-first search.
 * The reachable set is constructed by flattening the map (encoded in a BasicGraph)
 * into the bitset and setting the bits at reachable tiles (represent by the flattened indices) to 1
 * @param initial_state A State object representing the initial state.
 * @param window An integer representing the window.
 * @return A bitset_t object which represents the reachable set.
 */
bitset_t PDMPC::get_reachable_set(const State &initial_state, const int &window) const
{
    bitset_t reachable_set(G.size());

    deque<pair<int, int>> q;
    q.emplace_back(initial_state.location, 0);

    while(!q.empty()) {
        int current_location = q.front().first;
        int current_step = q.front().second;

        if (reachable_set[current_location] == 0) {
            reachable_set[current_location].flip();
        }

        if (current_step < window) {
            for (const auto &next_location: G.get_neighbors(current_location)) {
                if (!reachable_set[next_location]) {
                    q.emplace_back(next_location, current_step + 1);
                }
            }
        }

        q.pop_front();
    }

    return reachable_set;
}

/**
 * @brief Gets the coupling graph for a given vector of agent states and a window.
 * @description The coupling graph is constructed by checking if the reachable sets of the agents intersect.
 * If the reachable sets intersect, the agents are coupled, which means there is an edge between the agents.
 * @param agent_states A vector of State objects representing the agent states.
 * @param window An integer representing the window.
 * @return A pair containing a the coupling and a vector of reachable set objects.
 */
pair<Graph, vector<bitset_t>> PDMPC::get_coupling_graph(const vector<State> &agent_states, const int &window) const {

    auto num_agents = agent_states.size();

    Graph coupling(num_agents, bitset_t(num_agents));

    vector<bitset_t> reachable_sets(num_agents);

    for(int i = 0; i < num_agents; i++) {
        reachable_sets[i] = get_reachable_set(agent_states[i], window);
    }

    for(auto i = 0; i < num_agents; i++) {
        for(auto j = i + 1; j < num_agents; j++) {
            if(reachable_sets[i].intersects(reachable_sets[j])) {
                coupling[i][j].flip();
                coupling[j][i].flip();
            }
        }
    }

    return make_pair(coupling, reachable_sets);
}

/**
 * @brief Prioritizes a given coupling graph by ordering the agents randomly.
 * @param coupling_graph a coupling graph.
 * @return a directed graph with edges from higher to lower priority agents.
 */
Graph PDMPC::prioritize(const Graph &coupling_graph) const {
    auto num_agents = coupling_graph.size();

    vector<int> agent_indices(num_agents);

    for(auto i = 0; i < num_agents; i++) {
        agent_indices[i] = i;
    }

    std::shuffle(agent_indices.begin(), agent_indices.end(), std::mt19937(std::random_device()()));

    Graph directed_coupling_graph(coupling_graph);

    for(int &agent_index: agent_indices) {
        for(int i = 0; i < directed_coupling_graph[agent_index].size(); i++) {
            if(directed_coupling_graph[agent_index][i] && directed_coupling_graph[i][agent_index]) {
                directed_coupling_graph[i][agent_index].flip();
            }
        }
    }

    return directed_coupling_graph;
}

/**
 * @brief Weighs a given directed coupling graph by distance.
 * @description The weights are calculated by taking the inverse of the distance between the agents.
 * Therefore the distances are divided by a maximum distance of two agents, which is determined by 2*window
 * @param directed_coupling_graph A directed coupling graph.
 * @param agent_states A vector of State objects representing the agent states.
 * @param window An integer representing the window.
 * @return A WeightedGraph object with close agents having higher weights.
 */
WeightedGraph<double> PDMPC::weigh(const Graph &directed_coupling_graph, const vector<State> &agent_states, const int &window) const {
    WeightedGraph<double> directed_weighted_coupling_graph(directed_coupling_graph.size(), vector<double>(directed_coupling_graph.size()));

    double max_distance = 2 * window;

    for(int i = 0; i < directed_coupling_graph.size(); i++) {
        for(int j = 0; j < directed_coupling_graph[i].size(); j++) {
            int loc_a_i = agent_states[i].location;
            int loc_a_j = agent_states[j].location;

            double inverse_distance = 1 - (G.heuristics.at(loc_a_i).at(loc_a_j) / max_distance);

            directed_weighted_coupling_graph[i][j] = directed_coupling_graph[i][j] * inverse_distance;
        }
    }

    return directed_weighted_coupling_graph;
}

/**
 * @brief Group agents greedily (in decreasing order) by trying to sequentialize edges of closer agents first.
 * @param directed_weighted_coupling_graph A WeightedGraph object.
 * @param max_num_CLs An integer representing the maximum number of computation levels.
 * @return A Graph object with edges between sequentially computing agents.
 */
template<typename T>
Graph PDMPC::group(const WeightedGraph<T> &directed_weighted_coupling_graph, const int &max_num_CLs) const {
    auto num_agents = directed_weighted_coupling_graph.size();
    Graph directed_sequential_coupling_graph(num_agents, bitset_t(num_agents));

    if(max_num_CLs == 1) {
        return directed_sequential_coupling_graph;
    }

    //get list of weighted edges
    list<tuple<int, int, T>> weighted_edges;
    for(int i = 0; i < num_agents; i++) {
        for(int j = 0; j < num_agents; j++) {
            if(directed_weighted_coupling_graph[i][j]) {
                weighted_edges.push_back(make_tuple(i, j, directed_weighted_coupling_graph[i][j]));
            }
        }
    }

    //sort the edges in decreasing order
    weighted_edges.sort([](tuple<int, int, T> const &t1, tuple<int, int, T> const &t2) -> int {
        return std::get<2>(t1) > std::get<2>(t2);
    } );

    boost::optional<vector<int>> cls_optional = get_computation_levels(directed_sequential_coupling_graph);

    if(!cls_optional.has_value()) {
        return directed_sequential_coupling_graph;
    }

    vector<int> cls_of_agents = cls_optional.get();

    //find sequentializable edges greedily
    for(tuple<int, int, T> &node: weighted_edges) {
        int vertex_starting = std::get<0>(node);
        int vertex_ending = std::get<1>(node);
        int level_starting = cls_of_agents[vertex_starting];
        int level_ending = cls_of_agents[vertex_ending];

        //we can sequentialize edges if the priority of the ending vertex's level is lower than that of the starting vertex
        if(level_starting < level_ending) {
            directed_sequential_coupling_graph[vertex_starting][vertex_ending] = true;
            continue;
        }


        //otherwise we check if we can move the ending vertex to a level with lower priority
        directed_sequential_coupling_graph[vertex_starting][vertex_ending] = true;
        boost::optional<vector<int>> new_cls_optional = get_computation_levels(directed_sequential_coupling_graph);

        if(!new_cls_optional.has_value()) {
            //undo
            directed_sequential_coupling_graph[vertex_starting][vertex_ending] = false;
            continue;
        }

        vector<int> new_cls = new_cls_optional.get();
        if(*std::max_element(new_cls.begin(), new_cls.end()) <= max_num_CLs) {
            cls_of_agents = new_cls;
        } else {
            directed_sequential_coupling_graph[vertex_starting][vertex_ending] = false;
        }

    }

    return directed_sequential_coupling_graph;
}
template Graph PDMPC::group(const vector<vector<double>> &directed_weighted_coupling_graph, const int &max_num_CLs) const;


//debugging

/**
 * @brief Prints the map with '1' in place of the reachable set.
 * @param reachable_set A bitset_t object representing the reachable set.
 * @param start_location An integer representing the start location. Default is -1.
 * @param headline A string representing the headline. Default is "reachable set".
 */
void PDMPC::print_reachable_set(const bitset_t &reachable_set, const int &start_location, const string &headline) const {
    print_headline(headline);
    for(int i = 0; i < G.get_rows(); i++) {
        for(int j = 0; j < G.get_cols(); j++) {
            int id = i * G.get_cols() + j;
            if(id != start_location) {
                std::cout << reachable_set[id] << " ";
            } else {
                std::cout << "x ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * @brief Prints the matrix representing a coupling graph.
 * @param graph A WeightedGraph object.
 * @param headline A string representing the headline. Default is "weighted coupling graph".
 */
template<typename T>
void PDMPC::print_coupling_graph(const WeightedGraph<T> &graph, const string &headline) const {
    print_headline(headline);
    for(auto &row : graph) {
        for(auto &e : row) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
template void PDMPC::print_coupling_graph(const WeightedGraph<int> &graph, const string &headline) const;
template void PDMPC::print_coupling_graph(const WeightedGraph<double> &graph, const string &headline) const;

/**
 * @brief Prints the matrix representing a coupling graph.
 * @param graph A Graph object.
 * @param headline A string representing the headline. Default is "coupling graph".
 */
void PDMPC::print_coupling_graph(const Graph &graph, const string &headline) const {
    vector<vector<int>> g_v = vector<vector<int>>(graph.size(), vector<int>(graph.size()));
    for(int i = 0; i < g_v.size(); i++) {
        for(int j = 0; j < g_v.size(); j++) {
            g_v[i][j] = graph[i][j];
        }
    }
    print_coupling_graph(g_v, headline);
}

/**
 * @brief Prints a headline.
 * @param headline A string representing the headline.
 */
void PDMPC::print_headline(const std::string &headline) const {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << headline << std::endl;
    std::cout << std::endl;
}
