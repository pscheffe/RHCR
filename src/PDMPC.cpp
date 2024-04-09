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

    for(auto i = 0; i < range.size(); i++) {
        std::cout << range[i] << std::endl;
    }

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

void PDMPC::print_coupling_graph(vector<PDMPC::bitset_t> g) const {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "coupling graph:" << std::endl;
    std::cout << std::endl;
    for(auto & bs : g) {
        for(int j = 0; j < bs.size(); j++) {
            std::cout << bs[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
