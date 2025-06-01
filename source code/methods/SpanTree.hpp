#include "../graph.hpp"
#include <stack>

#define FLAG_ADD    true
#define FLAG_DEL    false

class SpanTreeIdx {
public:
    int n, landmark;
    ll m;
    vector<int> next;
    double x1, x2;

    SpanTreeIdx(Graph &graph, int landmark, vector<int> &bfs_next) {
        n = graph.n, m = graph.m;
        this->landmark = landmark;
        x1 = 0., x2 = 0.;
        sample_ust(graph, bfs_next);
    }

    void sample_ust(Graph &graph, vector<int> &bfs_next) {
        int cur;
        vector<vector<int>> &G = graph.adj, tree_rev(n + 1);
        vector<bool> intree(n + 1, false);
        vector<double> ps(n + 1, 0.);
        next = vector<int>(n + 1, -1);
        
        // double t1, t2, t3;
        // auto st = chrono::steady_clock::now();
        /* auto ed = chrono::steady_clock::now();
        t1 = chrono::duration<double>(ed - st).count();

        st = chrono::steady_clock::now(); */
        // 1. sample ust by wilson
        intree[landmark] = true;        
        for (int i = 1; i <= n; i++) {
            if (intree[i]) continue;
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true;
                tree_rev[next[cur]].push_back(cur);
            }
        }
        
        /* auto ed = chrono::steady_clock::now();
        t1 = chrono::duration<double>(ed - st).count();

        st = chrono::steady_clock::now(); */

        for (int i = 1; i <= n; i++) {
            double inc = (double) G[i].size() / 2. / m;
            for (cur = i; cur != landmark; cur = next[cur]) {
                ps[cur] += inc;
            }
        }

        /* ed = chrono::steady_clock::now();
        t2 = chrono::duration<double>(ed - st).count();
        
        st = chrono::steady_clock::now(); */

        vector<int> vis(n + 1, 0), fin(n + 1, 0);
        int timestamp = 0;
        function<void(int)> dfs = [&](int x) {
            vis[x] = ++timestamp;
            for (int ne : tree_rev[x]) {
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
       
    }

};

class SpanTree {
public:
    Graph &graph;
    vector<vector<int>> &G;
    int T;  // Index Size
    int n, landmark, delta;
    ll m;
    double eps;
    // res
    double kemeny_constant;
    
    vector<int> bfs_next;    

    vector<SpanTreeIdx> idxs;

    SpanTree(Graph &g, double eps, int landmark, int T = 0) : graph(g), G(g.adj) {
        n = g.n, m = g.m;
        this->eps = eps;
        this->landmark = landmark;
        build_bfs_tree();
        
        if (T != 0) this->T = T;
        else this->T = ceil(log(n) / eps / eps);

        kemeny_constant = 0.;
    }
// Init
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

    void build_index() {
        kemeny_constant = 0.;


        build_bfs_tree();
        for (int i = 0; i < T; i++) {
            idxs.emplace_back(graph, landmark, bfs_next);
        }
        kemeny_constant = estimate_sum_dii_LV_inv_ii() - estimate_d_LV_d();
    }

    double estimate_d_LV_d() {
        // 1/2m * d^T * (L_Vl)^+ * d
        double res = 0.;
        for (int i = 0; i < T; i++) {
            res += idxs[i].x2;
        }
        res = res / T;
        return res;
    }
    double estimate_sum_dii_LV_inv_ii() {
        // 1/2m * d^T * (L_Vl)^+ * d
        double res = 0.;
        for (int i = 0; i < T; i++) {
            res += idxs[i].x1;
        }
        res = res / T;
        return res;
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

    double spantree_sample_ust(double &t1, double &t2) {
        int cur, cnt = 1;
        double x1 = 0., x2 = 0.;
        vector<int> next(n + 1, 0);
        vector<bool> intree(n + 1, false);
        vector<double> ps(n + 1, 0.);
        Tree tree(n);
        
        auto st = chrono::steady_clock::now();
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
        
        auto ed = chrono::steady_clock::now();
        t1 = chrono::duration<double>(ed - st).count();

        st = chrono::steady_clock::now();

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

        ed = chrono::steady_clock::now();
        t2 = chrono::duration<double>(ed - st).count();

        return x1 - x2;       
    }

// Kemeny Constant
    double compute_kemeny_constant() {
        double kc = 0.;
        for (int i = 1; i <= T; i++) {
            kc = kc + spantree_sample_ust();
        }
        kemeny_constant = kc / T;
        return kc / T;
    }

    double compute_kemeny_constant(ofstream &of_res) {
        // static algorithm, with no index
        double kc = 0.;
        eps = 0.5;
        T = ceil(log(n) / eps / eps);
        double tot_time_build = 0., tot_time_query = 0.;

        for (int i = 1; i <= T; i++) {
            double time_build, time_query;
            kc = kc + spantree_sample_ust(time_build, time_query);
            tot_time_build += time_build;
            tot_time_query += time_query;
        }
        of_res << "Index_Building_Time: " << tot_time_build << endl;
        of_res << "Query_Time: " << tot_time_query << endl;
        return kc / T;
    }

};


