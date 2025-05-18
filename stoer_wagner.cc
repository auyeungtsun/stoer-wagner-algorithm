#include <iostream>
#include <vector>
#include <numeric>
#include <tuple>
#include <limits>
#include <algorithm>
#include <cassert>

using namespace std;

const int INF = numeric_limits<int>::max();

/**
 * @brief Computes the minimum cut of an undirected, weighted graph using the Stoer-Wagner algorithm.
 * 
 * @param edges A vector of tuples, where each tuple represents an edge (u, v, w), 
 *              where u and v are the vertices and w is the weight of the edge.
 * @param n The number of vertices in the graph.
 * @return The weight of the minimum cut in the graph.
 * 
 * @pre The graph must be undirected.
 * @pre Edge weights must be non-negative.
 * 
 * @note Time Complexity: O(V^3) where V is the number of vertices.
 * @note Space Complexity: O(V^2) for the adjacency matrix.
 */
int stoer_wagner(const vector<tuple<int, int, int>>& edges, int n) {
    if (n <= 1) {
        return 0;
    }
    
    // Adjacency Matrix
    vector<vector<int>> adj(n, vector<int>(n, 0));

    for (const auto& edge : edges) {
        int u, v, w;
        tie(u, v, w) = edge;
        if(u >= 0 && u < n && v >= 0 && v < n && w >= 0){
            adj[u][v] += w;
            adj[v][u] += w;
        }
    }

    int min_cut_weight = INF;
    vector<int> vertices(n); // the set of active vertices (nodes) in the graph during each phase
    iota(vertices.begin(), vertices.end(), 0); // Fill with 0, 1, ..., n-1

    for (int phase = 0; phase < n - 1; ++phase) {
        // --- Maximum Adjacency Search Phase ---
        // w(v_n, A_{n-1}) calculated at the end of the MAS phase is equal to the minimum v_{n-1}-v_n cut value.
        vector<bool> in_a(n, false);    // Tracks vertices added to set A
        vector<int> weights(n, 0);      // weights[v] = sum of weights of edges from v to A
        vector<int> order;              // Order in which vertices are added
        int num_active_in_phase = vertices.size();  // Number of vertices considered in this phase

        int prev_vertex = -1;           // Second-to-last vertex added
        int last_vertex = -1;           // Last vertex added
        
        for (int i = 0; i < num_active_in_phase; ++i) {
            int best_v = -1;
            int max_weight = -1;

            for (int v_idx : vertices) {
                if (!in_a[v_idx] && weights[v_idx] > max_weight) {
                    max_weight = weights[v_idx];
                    best_v = v_idx;
                }
            }

            in_a[best_v] = true;
            order.push_back(best_v);
            prev_vertex = last_vertex;
            last_vertex = best_v;

            for (int neighbor_idx : vertices) {
                if (!in_a[neighbor_idx]) {
                    weights[neighbor_idx] += adj[best_v][neighbor_idx];
                }
            }
        }
        
        int cut_of_the_phase = weights[last_vertex];
        min_cut_weight = min(min_cut_weight, cut_of_the_phase);

        // --- Contraction Phase ---
        // Merge last_vertex (t) into prev_vertex (s)
        int s = prev_vertex;
        int t = last_vertex;

        for (int v_idx : vertices) {
            if (v_idx != s && v_idx != t) {
                adj[s][v_idx] += adj[t][v_idx];
                adj[v_idx][s] = adj[s][v_idx];
            }
        }

        auto it = find(vertices.begin(), vertices.end(), t);
        if (it != vertices.end()) {
            vertices.erase(it);
        }

    }

    return min_cut_weight;
}

void test_stoer_wagner() {
    cout << "Running Stoer-Wagner Tests..." << endl;

    // Test Case 1: Trivial cases
    assert(stoer_wagner({}, 0) == 0);
    assert(stoer_wagner({}, 1) == 0);
    cout << "Test Case 1 Passed (Trivial cases n=0, n=1)" << endl;

    // Test Case 2: Two vertices
    vector<tuple<int, int, int>> edges2 = {{0, 1, 5}};
    assert(stoer_wagner(edges2, 2) == 5);
    vector<tuple<int, int, int>> edges2_multi = {{0, 1, 5}, {0, 1, 3}};
    assert(stoer_wagner(edges2_multi, 2) == 8);
    vector<tuple<int, int, int>> edges2_none = {};
    assert(stoer_wagner(edges2_none, 2) == 0);
     cout << "Test Case 2 Passed (Two vertices)" << endl;

    // Test Case 3: Triangle graph
    vector<tuple<int, int, int>> edges3 = {{0, 1, 2}, {1, 2, 3}, {0, 2, 4}};
    assert(stoer_wagner(edges3, 3) == 5);
    cout << "Test Case 3 Passed (Triangle graph)" << endl;

    // Test Case 4: Square graph
    vector<tuple<int, int, int>> edges4 = {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 0, 5}};
    assert(stoer_wagner(edges4, 4) == 5);
     cout << "Test Case 4 Passed (Square graph)" << endl;

    // Test Case 5: Square with diagonal
    vector<tuple<int, int, int>> edges4_diag = {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 0, 5}, {0, 2, 6}};
    assert(stoer_wagner(edges4_diag, 4) == 5);
     cout << "Test Case 5 Passed (Square with diagonal)" << endl;

    // Test Case 6: Disconnected graph
    vector<tuple<int, int, int>> edges_disconnected = {{0, 1, 5}, {2, 3, 6}};
    assert(stoer_wagner(edges_disconnected, 4) == 0);
    cout << "Test Case 6 Passed (Disconnected graph)" << endl;

    // Test Case 7: More complex graph (example often used for Stoer-Wagner)
    vector<tuple<int, int, int>> edges_complex = {
        {0, 1, 2}, {0, 4, 3}, {1, 2, 3}, {1, 4, 2}, {1, 5, 2},
        {2, 3, 4}, {2, 6, 2}, {3, 6, 2}, {3, 7, 2}, {4, 5, 3},
        {5, 6, 1}, {6, 7, 3}
    };
    assert(stoer_wagner(edges_complex, 8) == 4);
    cout << "Test Case 7 Passed (Complex graph)" << endl;

    // Test Case 8: Graph where merging order might matter
    vector<tuple<int, int, int>> edges_merge = {
        {0, 1, 1}, {1, 2, 10}, {2, 3, 1}, {0, 4, 10}, {1, 4, 1}, {2, 5, 10}, {3, 5, 1}
    };
    assert(stoer_wagner(edges_merge, 6) == 2);
    cout << "Test Case 8 Passed (Specific merge order test)" << endl;

    // Test Case 9: Empty graph with n > 1
    vector<tuple<int, int, int>> edges_empty_n3 = {};
    assert(stoer_wagner(edges_empty_n3, 3) == 0);
    cout << "Test Case 9 Passed (Empty graph n=3)" << endl;

    cout << "All Stoer-Wagner tests passed!" << endl;
}

void run_stoer_wagner_sample() {
    vector<tuple<int, int, int>> edges = {
        {0, 1, 10}, {0, 2, 5}, {1, 2, 3}, {1, 3, 2}, {2, 3, 6}
    };
    int n = 4;
    int min_cut = stoer_wagner(edges, n);
    cout << "Stoer-Wagner Sample - Minimum Cut: " << min_cut << endl;
}

int main() {
    test_stoer_wagner();
    run_stoer_wagner_sample();
    return 0;
}