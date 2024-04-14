#include "PDMPC.h"

#include <random>

PDMPC::PDMPC(const BasicGraph &G) : G(G) {}

PDMPC::bitset_t PDMPC::get_reachable_set(const State &initial_state, const int &window) const
{
    PDMPC::bitset_t reachable_set(G.size());

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

vector<PDMPC::bitset_t> PDMPC::get_coupling_graph(const vector<State> &agent_states, const int &window) const {

    auto num_agents = agent_states.size();

    vector<PDMPC::bitset_t> coupling_matrix(num_agents, PDMPC::bitset_t(num_agents));

    vector<PDMPC::bitset_t> reachable_sets(num_agents);

    for(int i = 0; i < num_agents; i++) {
        reachable_sets[i] = get_reachable_set(agent_states[i], window);
    }

    for(auto i = 0; i < num_agents; i++) {
        for(auto j = i + 1; j < num_agents; j++) {
            if(reachable_sets[i].intersects(reachable_sets[j])) {
                coupling_matrix[i][j].flip();
                coupling_matrix[j][i].flip();
            }
        }
    }

    return coupling_matrix;
}

vector<PDMPC::bitset_t> PDMPC::prioritize(const vector<PDMPC::bitset_t>& coupling_graph) const {
    auto num_agents = coupling_graph.size();

    std::vector<int> range(num_agents);

    for(auto i = 0; i < num_agents; i++) {
        range[i] = i;
    }

    std::shuffle(range.begin(), range.end(), std::mt19937(std::random_device()()));

    std::cout << "priority order: ";
    for(int e: range) {
        std::cout << e + 1 << ", ";
    }
    std::cout << std::endl;

    vector<PDMPC::bitset_t> directed_coupling_graph(coupling_graph);

    for(auto &e: range) {
        for(int i = 0; i < directed_coupling_graph[e].size(); i++) {
            if(directed_coupling_graph[e][i] && directed_coupling_graph[i][e]) {
                directed_coupling_graph[i][e].flip();
            }
        }
    }

    return directed_coupling_graph;
}

vector<vector<int>> PDMPC::weigh(const vector<PDMPC::bitset_t> &directed_coupling_graph, const vector<State> &agent_states) const {
    vector<vector<int>> directed_weighted_coupling_graph(directed_coupling_graph.size(), vector<int>(directed_coupling_graph.size()));

    for(int i = 0; i < directed_coupling_graph.size(); i++) {
        for(int j = 0; j < directed_coupling_graph[i].size(); j++) {
            int loc_a_i = agent_states[i].location;
            int loc_a_j = agent_states[j].location;
            directed_weighted_coupling_graph[i][j] = directed_coupling_graph[i][j] * G.get_Manhattan_distance(loc_a_i, loc_a_j);
        }
    }

    return directed_weighted_coupling_graph;
}

vector<PDMPC::bitset_t> PDMPC::group(const vector<vector<int>> &directed_weighted_coupling_graph, const int &max_num_cls) const {
    auto num_agents = directed_weighted_coupling_graph.size();
    vector<PDMPC::bitset_t> directed_sequential_coupling_graph(num_agents, PDMPC::bitset_t(num_agents));

    if(max_num_cls == 1) {
        return directed_sequential_coupling_graph;
    }

    boost::optional<vector<int>> cls_optional = get_computation_levels(directed_sequential_coupling_graph);

    if(cls_optional.has_value()) {
        vector<int> cls = cls_optional.get();
        for(auto e : cls) {
            std::cout << e << " ,";
        }
        std::cout << std::endl;
    } else {
        std::cout << "optional has no value" << std::endl;
    }

    return directed_sequential_coupling_graph;
}

boost::optional<vector<int>> PDMPC::get_computation_levels(const vector<PDMPC::bitset_t> &directed_coupling_graph) const {
    auto num_agents = directed_coupling_graph.size();
    vector<PDMPC::bitset_t> A(directed_coupling_graph);

    //Kahn's algorithm
    vector<int> computation_levels(num_agents);
    std::list<int> no_incoming_edges;
    PDMPC::bitset_t assigned_agents_mask(num_agents);

    int current_cl_level = 1;
    while(!assigned_agents_mask.all()) {
        //get agents with no incoming edges
        for (auto i = 0; i < num_agents; i++) {
            if(!assigned_agents_mask[i]) {
                int num_incoming_edges = 0;
                for (auto j = 0; j < num_agents; j++) {
                    num_incoming_edges = num_incoming_edges + A[j][i];
                }
                if (num_incoming_edges == 0) {
                    no_incoming_edges.emplace_back(i);
                }
            }
        }

        if (no_incoming_edges.empty()) {
            return boost::optional<vector<int>>{};
        }

        for(auto agent : no_incoming_edges) {
            computation_levels[agent] = current_cl_level;
            assigned_agents_mask[agent] = true;
            //clear rows
            for(int i = 0; i < num_agents; i++) {
                A[agent][i] = false;
            }
        }

        no_incoming_edges.clear();
        current_cl_level = current_cl_level + 1;
    }

    return computation_levels;

}

//debugging

void PDMPC::print_reachable_set(PDMPC::bitset_t bs, int start_location) const {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "reachable set:" << std::endl;
    std::cout << std::endl;
    for(int i = 0; i < G.get_rows(); i++) {
        for(int j = 0; j < G.get_cols(); j++) {
            int id = i * G.get_cols() + j;
            if(id != start_location) {
                std::cout << bs[id] << " ";
            } else {
                std::cout << "x ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void PDMPC::print_coupling_graph(vector<vector<int>> graph) const {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "coupling graph:" << std::endl;
    std::cout << std::endl;
    for(auto &row : graph) {
        for(auto &e : row) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void PDMPC::print_coupling_graph(vector<PDMPC::bitset_t> graph) const {
    vector<vector<int>> g_v = vector<vector<int>>(graph.size(), vector<int>(graph.size()));
    for(int i = 0; i < g_v.size(); i++) {
        for(int j = 0; j < g_v.size(); j++) {
            g_v[i][j] = graph[i][j];
        }
    }
    print_coupling_graph(g_v);
}

vector<PDMPC::bitset_t> PDMPC::to_unweighted_graph(vector<vector<int>> graph) {
    vector<PDMPC::bitset_t> unweighted_graph(graph.size(), PDMPC::bitset_t(graph.size()));

    for(int i = 0; i < unweighted_graph.size(); i++) {
        for(int j = 0; j < unweighted_graph.size(); j++) {
            if(graph[i][j]) {
                unweighted_graph[i][j].flip();
            }
        }
    }

    return unweighted_graph;
}
