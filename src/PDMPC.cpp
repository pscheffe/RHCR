#include "PDMPC.h"

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
