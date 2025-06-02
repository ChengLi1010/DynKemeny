#pragma once
#include "../graph.hpp"
#include "Base.hpp"
#include <stack>

class SpanTree : public KCSolver {
public:
    Graph &graph;
    vector<vector<int>> &G;
    int T;  // Index Size
    int n, landmark;
    ll m;
    double eps;
    // res
    double kemeny_constant;
    
    vector<int> bfs_next;    

    SpanTree(Graph &g, double eps, int landmark, int T = 0) : graph(g), G(g.adj) {
        n = g.n, m = g.m;
        this->eps = eps;
        this->landmark = landmark;
        build_bfs_tree();
        
        if (T != 0) this->T = T;
        else this->T = ceil(log(n) / eps / eps);

        kemeny_constant = 0.;
    }
    void build_bfs_tree(int flag = FLAG_TREE_BFS) {
        int root = landmark;
        bfs_next = vector<int>(n + 1, -1);
        // BFS
        if (flag == FLAG_TREE_WILSON) {
            vector<int> next(n + 1, 0);
            vector<bool> intree(n + 1, false);
            int cnt = 1, cur;
            intree[root] = true;
            for (int i = 1; i <= n; i++) {
                for (int cur = i; !intree[cur]; cur = next[cur]) {
                    next[cur] = G[cur][rand() % G[cur].size()];
                }
                for (int cur = i; !intree[cur]; cur = next[cur]) {
                    intree[cur] = true;
                    bfs_next[cur] = next[cur];
                    cnt++;
                }
                if (cnt == n) break;
            }
        } else if (flag == FLAG_TREE_BFS) {
            vector<bool> visit(n + 1,false);
            queue<int> q;

            q.push(root);
            visit[root] = true;

            while (!q.empty()) {
                int cur = q.front();
                q.pop();

                for (int ne : G[cur]) {
                    if (!visit[ne]) {
                        q.push(ne);
                        visit[ne] = true;

                        bfs_next[ne] = cur;
                    }
                }
            }
        } else if (flag == FLAG_TREE_DFS) {
            vector<bool> intree(n + 1, false);

            function<void(int)> tree_dfs = [&](int cur) {
                for (int ne : G[cur]) {
                    if (!intree[ne]) {
                        intree[ne] = true;
                        bfs_next[ne] = cur;
                        tree_dfs(ne);
                    }
                }
            };

            intree[root] = true;
            tree_dfs(root);
        }
        return;
    }

    double spantree_sample_ust() {
        int cur, cnt = 1;
        double x1 = 0., x2 = 0.;
        vector<int> next(n + 1, 0);
        vector<bool> intree(n + 1, false);
        vector<double> ps(n + 1, 0.);
        Tree tree(n);
        
        // 1. sample ust by wilson
        intree[landmark] = true;        
        for (int i = 1; i <= n; i++) {
            if (intree[i]) continue;
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true;
                tree.add_edge(next[cur], cur);
            }
        }
        
        for (int i = 1; i <= n; i++) {
            double inc = (double) G[i].size() / 2. / m;
            for (cur = i; cur != landmark; cur = next[cur]) {
                ps[cur] += inc;
            }
        }

        vector<int> vis(n + 1, 0), fin(n + 1, 0);
        int timestamp = 0;
        function<void(int)> dfs = [&](int x) {
            vis[x] = ++timestamp;
            for (int i = tree.head[x]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }
            fin[x] = ++timestamp;
        };
        dfs(landmark);

        for (int i = 1; i <= n; i++) {
            for (cur = i; cur != landmark; cur = bfs_next[cur]) {
                if (next[cur] == bfs_next[cur] && vis[cur] <= vis[i] && fin[i] <= fin[cur]) {
                    x1 += (double)G[i].size();
                    x2 += (double)G[i].size() * ps[cur];
                }else if (next[bfs_next[cur]] == cur && vis[bfs_next[cur]] <= vis[i] && fin[i] <= fin[bfs_next[cur]]) {
                    x1 -= (double) G[i].size();
                    x2 -= (double) G[i].size() * ps[bfs_next[cur]];
                }
            }
        }

        return x1 - x2;       
    }

// Kemeny Constant
    double compute_kemeny_constant() {
        double kc = 0.;
        for (int i = 1; i <= T; i++) {
            kc = kc + spantree_sample_ust();
        }
        kemeny_constant = kc / T;
        return kemeny_constant;
    }

};


