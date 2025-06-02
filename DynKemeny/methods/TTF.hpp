#pragma once

#include "../graph.hpp"
#include "Base.hpp"
#include <stack>

#define FLAG_ADD    true
#define FLAG_DEL    false
#define FLAG_BASIC  0
#define FLAG_IMPRV  1
#define FLAG_TREE_BFS 0
#define FLAG_TREE_WILSON 1
#define FLAG_TREE_DFS 2

class FixedTree {
public:
    vector<int> bfs_next; // bfs tree
    vector<int> dfs_in, dfs_out;
    Tree tree;
    FenwickTree<int64_t> bfs_fwt;
    FenwickTree<int64_t> f_cnt_fwt;
    int max_depth = 0;

    FixedTree(){
    }

    void build_bfs_tree(Graph &graph, const int root, int flag = FLAG_TREE_BFS) {
        int n = graph.n;
        vector<vector<int>> &G = graph.adj;
        
        bfs_next = vector<int>(n + 1, -1);
        dfs_in = vector<int> (n + 1), dfs_out = vector<int> (n + 1);
        tree = Tree(n);    
        bfs_fwt = FenwickTree<int64_t>(n);
        f_cnt_fwt = FenwickTree<int64_t>(n);
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
                    tree.add_edge(next[cur], cur);
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
                        tree.add_edge(cur, ne);
                    }
                }
            }
        } else if (flag == FLAG_TREE_DFS) {
            vector<bool> intree(n + 1, false);
            
            stack<int> s;
            s.push(root);
            intree[root] = true;

            while (!s.empty()) {
                int cur = s.top();
                s.pop();

                assert(cur >= 1 && cur < G.size());
                for (int ne : G[cur]) {
                    assert(ne >= 0 && ne < intree.size());
                    if (!intree[ne]) {
                        intree[ne] = true;
                        bfs_next[ne] = cur;
                        tree.add_edge(cur, ne);
                        s.push(ne);
                    }
                }
            }
           
        }
        // DFS order
        int time_stamp = 0;
        function<void(int)> dfs = [&](int x) {
            dfs_in[x] = ++time_stamp;

            for (int i = tree.head[x]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }
            dfs_out[x] = time_stamp;
        };
        dfs(root);

        return;
    }
};

class TreeIdx {
public:
    int n, root;
    ll m;

    vector<int> next;
    Tree tree;

    vector<ll> vol, vol_sum;
    vector<int> path_cnt;
    ll path_cnt_tot;

    double kc, weight;


    TreeIdx(Graph &graph, FixedTree &ft, int root) {
        n = graph.n, m = graph.m;
        kc = 0., weight = 1.;
        this->root = root;
        basic_sample_ust(graph, ft);
    }

    void basic_sample_ust(Graph &graph, FixedTree &ft) {
        kc = 0.;
        int cur, intree_cnt = 1, br_tot = 0;
        vector<vector<int>> &G = graph.adj;
        vector<bool> intree(n + 1, false);
        vector<Edge> branch(n);
        tree = Tree(n);
        next.resize(n + 1, 0), path_cnt.resize(n + 1, 0), vol.resize(n + 1, 0), vol_sum.resize(n + 1, 0);
        // 1. sample ust by wilson
        intree[root] = true;        
        for (int i = 1; i <= n; i++) {
            if (intree[i]) continue;
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true;
                intree_cnt++;
                tree.add_edge(next[cur], cur);
            }
            branch[br_tot++] = make_pair(i, cur);
            if (intree_cnt == n) break;
        }
        
        // 2. estimate f(\tau)
        // 2.1 calculate vol of each subtree
        for (int i = br_tot - 1; i >= 0; i--) {
            int s = branch[i].first, t = branch[i].second;
            for (cur = s; cur != t; cur = next[cur]) {
                vol[cur] += G[cur].size();
                vol[next[cur]] += vol[cur];
            }
        }
        // 2.2 DFS
        ll vol_tot = 2 * m;
        int f_cnt = 0;
        path_cnt_tot = 0;
        function<void(int)> dfs = [&](int cur) {
            int fa = next[cur];
            if (ft.bfs_next[cur] == fa) { // same direction
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], 1);
            } else if (ft.bfs_next[fa] == cur) { // reverse dirction
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -1);
            }

            path_cnt[cur] = ft.f_cnt_fwt.preSum(ft.dfs_in[cur]);
            path_cnt_tot += G[cur].size() * path_cnt[cur];
            vol_sum[cur] = ft.bfs_fwt.preSum(ft.dfs_in[cur]);
            kc += (path_cnt[cur] - (double)vol_sum[cur] / vol_tot) *  G[cur].size();

            for (int i = tree.head[cur]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }

            if (ft.bfs_next[cur] == fa) { // same direction
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -1);
            } else if (ft.bfs_next[fa] == cur) { // reverse dirction
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], 1);
            }
        };

        for (int i = tree.head[root]; i; i = tree.edgeNxt[i]) {
            dfs(tree.edgeTo[i]);
        }
        return;
    }

    void basic_update_sample(Graph &graph, FixedTree &ft, int u, int v) {
        n = graph.n, m = graph.m;
        kc = 0., weight = 1.;
        int cur, intree_cnt = 3;
        vector<vector<int>> &G = graph.adj;
        vector<bool> intree(n + 1, false);
        vector<Edge> branch;
        tree = Tree(n);
        next = vector<int>(n + 1, 0), path_cnt = vector<int> (n + 1, 0),
        vol = vector<ll>(n + 1, 0), vol_sum = vector<ll> (n + 1, 0);

        // 1. sample spanning tree contains (u, v) by wilson  & redirect edge
        intree[u] = intree[v] = true;
        next[u] = v;
        // root -> u(v) & redirect edge
        for (cur = root; !intree[cur]; cur = next[cur]) {
            next[cur] = G[cur][rand() % G[cur].size()];
        }
        intree[root] = true;
        int last = root, ne;
        for (cur = next[root]; !intree[cur]; cur = ne) {
            intree[cur] = true;
            intree_cnt++;
            tree.add_edge(last, cur);
            
            ne = next[cur];
            next[cur] = last;

            last = cur;
        }
        // cur = u or v, last -> cur -- ne  ==> last <- cur <- ne
        ne = (cur == u ? v : u);
        next[cur] = last; tree.add_edge(last, cur);
        next[ne] = cur; tree.add_edge(cur, ne);

        branch.emplace_back(ne, root);
        // normal wilson
        for (int i = 1; i <= n; i++) {
            if (intree[i]) continue;
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true; intree_cnt++;
                tree.add_edge(next[cur], cur);
            }
            branch.emplace_back(i, cur);
            if (intree_cnt == n) break;
        }

        // 2. estimate f(\tau)
        // 2.1 calculate vol of each subtree
        for (int i = branch.size() - 1; i >= 0; i--) {
            int s = branch[i].first, t = branch[i].second;
            for (cur = s; cur != t; cur = next[cur]) {
                vol[cur] += G[cur].size();
                vol[next[cur]] += vol[cur];
            }
        }
        // 2.2. DFS
        ll vol_tot = 2 * m;
        function<void(int)> dfs = [&](int cur) {
            int fa = next[cur];
            if (ft.bfs_next[cur] == fa) { // same direction
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], 1);
            } else if (ft.bfs_next[fa] == cur) { // reverse dirction
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -1);
            }

            path_cnt[cur] = ft.f_cnt_fwt.preSum(ft.dfs_in[cur]);
            path_cnt_tot += G[cur].size() * path_cnt[cur];
            vol_sum[cur] = ft.bfs_fwt.preSum(ft.dfs_in[cur]);
            kc += (path_cnt[cur] - (double)vol_sum[cur] / vol_tot) * (double) G[cur].size();

            for (int i = tree.head[cur]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }

            if (ft.bfs_next[cur] == fa) { // same direction
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -1);
            } else if (ft.bfs_next[fa] == cur) { // reverse dirction
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], 1);
            }
        };

        for (int i = tree.head[root]; i; i = tree.edgeNxt[i]) {
            dfs(tree.edgeTo[i]);
        }

        return;
    }

    // find unique path in tree from u to v, return the index of lca(u, v) in vector"path"
    int find_tree_path(int u, int v, vector<int> &path) {
        vector<int> path_u, path_v;
        for (int cur = u; cur != root; cur = next[cur]) {
            path_u.push_back(cur);
        }
        for (int cur = v; cur != root; cur = next[cur]) {
            path_v.push_back(cur);
        }
        
        int p1 = path_u.size() - 1, p2 = path_v.size() - 1;
        while (p1 >= 0 && p2 >= 0 && path_u[p1] == path_v[p2]) {
            p1--, p2--;
        }

        path_u.push_back(root);
        path.resize(p1 + p2 + 3);

        int i = 0;
        for (; i <= p1 + 1; i++) {
            path[i] = path_u[i];
        }
        for (int j = p2; j >= 0; j--, i++) {
            path[i] = path_v[j];
        }
        return p1 + 1;
    }

    // link-cut for ADD(u, v)
    void improved_update_link_cut(Graph &graph, FixedTree &ft, int u, int v) {
        m = graph.m;
        vector<vector<int>> &G = graph.adj;
        
        // 1. find tree cycle of (u,v)
        vector<int> path;
        int lca = find_tree_path(u, v, path);   
        int C = path.size();
        int del_edge_idx = rand() % (C - 1); // x = path[del_edge_idx], y = path[del_edge_idx + 1]
        if (lca <= del_edge_idx) {
            swap(u, v);
            reverse(path.begin(), path.end());
            lca = C - 1 - lca;
            del_edge_idx = C - 1 - del_edge_idx - 1;
        }
        int x = path[del_edge_idx], y = path[del_edge_idx + 1];
        // u -> ... -> x -> y ->... (lca) ... -> v

        // change u-> ... -> x tree edge direction 
        for (int i = 1; i <= del_edge_idx; i++) {
            int cur = path[i], ne = path[i - 1];
            tree.change_edge(next[cur], cur, ne, cur);
            next[cur] = ne;
        }
        tree.change_edge(next[u], u, v, u);
        next[u] = v;
        
        int dfs_st_y = -1, dfs_st_v = u, dfs_st_lca = -1;
        // 2. update vol
        // 2.1 y-> .. -> lca
        for (int i = del_edge_idx + 1; i < lca; i++) {
            int cur = path[i];
            vol[cur] -= vol[x];
            if (ft.bfs_next[cur] == next[cur] || ft.bfs_next[next[cur]] == cur) {
                dfs_st_y = cur; 
            }
        }
        // 2.2 u-> .. ->x
        for (int i = del_edge_idx; i >= 0; i--) {
            int cur = path[i];
            if (i > 0) vol[cur] -= vol[path[i-1]];
            if (i < del_edge_idx) vol[cur] += vol[path[i+1]];
        }
        vol[u]++; // G.add(u, v), d(u)+1, d(v)+1
        // 2.3 v -> lca
        for (int cur = v; cur != path[lca]; cur = next[cur]) {
            vol[cur] += vol[u] + 1;
            if (ft.bfs_next[cur] == next[cur] || ft.bfs_next[next[cur]] == cur) {
                dfs_st_v = cur;   
            }
        }
        // 2.4 lca -> root
        for (int cur = path[lca]; cur != root; cur = next[cur]) {
            vol[cur] += 2;
            if (ft.bfs_next[cur] == next[cur] || ft.bfs_next[next[cur]] == cur) {
                dfs_st_lca = cur;
            }
        }


        // 3. Update f(\tau) by dfs
        path_cnt_tot += path_cnt[u] + path_cnt[v];
        kc = kc - path_cnt_tot;
        kc = kc * m / (m + 1);
        ll vol_tot = 2 * m; 
        kc -= (double) (vol_sum[u] + vol_sum[v]) / vol_tot;

        vector<bool> inSub_u(n + 1, false);
        bool isPass_u = false;
        function<void(int)> dfs = [&](int cur) {
            if (cur == u) isPass_u = true;
            if (isPass_u) inSub_u[cur] = true;
            int fa = next[cur];
            if (ft.bfs_next[cur] == fa) {
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], 1);
            } else if (ft.bfs_next[fa] == cur) {
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -1);
            }

            path_cnt_tot -= path_cnt[cur] * G[cur].size();
            path_cnt[cur] = ft.f_cnt_fwt.preSum(ft.dfs_in[cur]);
            path_cnt_tot += path_cnt[cur] * G[cur].size();

            kc += (double) vol_sum[cur] * G[cur].size() / vol_tot;
            vol_sum[cur] = ft.bfs_fwt.preSum(ft.dfs_in[cur]);
            kc -= (double) vol_sum[cur] * G[cur].size() / vol_tot;

            for (int i = tree.head[cur]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }

            if (ft.bfs_next[cur] == fa) {
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur],  -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -1);
            } else if (ft.bfs_next[next[cur]] == cur) {
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], 1);
            }
            if (cur == u) isPass_u = false;
            return;
        };
        if (dfs_st_lca != -1) {
            dfs(dfs_st_lca);
        } else {
            if (dfs_st_y != -1) dfs(dfs_st_y);
            dfs(dfs_st_v);
        }
        ft.bfs_fwt.clear();
        
        kc += path_cnt_tot;

        // 4. compute weight 
        bool u_is_less_side = (vol[u] > m);
        ll d_tau_e = 0;
        for (int i = 1; i <= n; i++) {
            if (u_is_less_side && inSub_u[i]) {
                for (int ne : G[i]) {
                    if (!inSub_u[ne]) d_tau_e++;
                }
            }else if(!u_is_less_side && !inSub_u[i]){
                for (int ne : G[i]) {
                    if (inSub_u[ne]) d_tau_e++;
                }
            }
        }
        weight = weight * (C - 1) / d_tau_e;    // d_tau = (C - 1)

        return;
    }

    // cut-link for DEL(u, v)
    void improved_update_cut_link(Graph &graph, FixedTree &ft, int u, int v) {
        m = graph.m;
        vector<vector<int>> &G = graph.adj;
        // 1.1 find edges to link
        if (next[v] == u) swap(u, v); // let next[u] = v
        queue<int> q;
        vector<Edge> add_edges;
        vector<int> nodes;
        vector<bool> inSameTree(n + 1, false);
        if (vol[u] > m) {
            q.push(root);
        } else {
            q.push(u);
        }
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            nodes.push_back(cur);
            inSameTree[cur] = true;
            for (int i = tree.head[cur]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                if (ne == u) continue;
                q.push(ne);
            }
        } 
        for (int i : nodes) {
            for (int ne : G[i]) {
                if (!inSameTree[ne]) {
                    if ((i == u && ne == v) || (i == v && ne == u)) continue;
                    add_edges.emplace_back(i, ne);
                }
            }
        }
        // 1.2 cut & link
        int d_tau_e = add_edges.size(), d_tau = 0;
        int add_edge_idx = rand() % add_edges.size();
        int x = add_edges[add_edge_idx].first, y = add_edges[add_edge_idx].second;
        if (vol[u] > m) swap(x, y); // root belongs to {nodes}, first node in the same tree with root
        
            // change edge direction
        for (int cur = x, la = y, ne; cur != v; cur = ne) {
            ne = next[cur];
            next[cur] = la;
            tree.change_edge(ne, cur, la, cur);
            la = cur;
        }
        
        int dfs_st_y = x, dfs_st_v = -1, dfs_st_lca = -1;
        int cur, la;
        // 2. update vol
        // 2.1 v -> root
        vector<int> length_to_v(n + 1, -1);
        for (cur = v, la = 0; cur != root; cur = next[cur]) {
            vol[cur] -= vol[u] + 1; // -vol(u) - 1(d(v))
            length_to_v[cur] = length_to_v[la] + 1;
            la = cur;

            if (ft.bfs_next[cur] == next[cur] || ft.bfs_next[next[cur]] == cur) {
                dfs_st_v = cur;
            }
        }
        length_to_v[cur] = length_to_v[la];
        // 2.2 u -> x
        int length_to_u = 0;
        vol[u]--;
        for (cur = u, la = -1; cur != y; cur = next[cur]) {
            if (cur != x) vol[cur] -= vol[next[cur]];
            if (cur != u) vol[cur] += vol[la];
            la = cur;
            length_to_u++;
        }
        // 2.3 y -> root
        int lca = -1;
        for (cur = y; cur != root; cur = next[cur]) {
            vol[cur] += vol[x];
            if (length_to_v[cur] >= 0 && lca == -1) {   // cur = lca(u, v)
                d_tau = length_to_u + length_to_v[cur]; // d_tau = length of cycle
                lca = cur;
            }
            length_to_u++;

            if (ft.bfs_next[cur] == next[cur] || ft.bfs_next[next[cur]] == cur) {
                if (lca == -1) dfs_st_y = cur;
                else dfs_st_lca = cur;
            }
        }
        if (lca == -1) {
            lca = root;
            d_tau = length_to_u + length_to_v[root];
        }

        // 3. update f(\tau) by dfs
        path_cnt_tot -= path_cnt[u] + path_cnt[v];
        kc = (kc - path_cnt_tot) * (double) m / (m - 1);
        m = graph.m;
        ll vol_tot = 2 * m;
        kc += (double) (vol_sum[u] + vol_sum[v]) / vol_tot;

        function<void(int)> dfs = [&](int cur) {
            int fa = next[cur];
            if (ft.bfs_next[cur] == fa) {
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], 1);
            } else if (ft.bfs_next[fa] == cur) {
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], -1);
            }

            path_cnt_tot -= path_cnt[cur] * G[cur].size();
            path_cnt[cur] = ft.f_cnt_fwt.preSum(ft.dfs_in[cur]);
            path_cnt_tot += path_cnt[cur] * G[cur].size();

            kc += (double) vol_sum[cur] * G[cur].size() / vol_tot;
            vol_sum[cur] = ft.bfs_fwt.preSum(ft.dfs_in[cur]);
            kc -= (double) vol_sum[cur] * G[cur].size() / vol_tot;

            for (int i = tree.head[cur]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];
                dfs(ne);
            }

            if (ft.bfs_next[cur] == fa) {
                ft.bfs_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur],  -vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[cur], ft.dfs_out[cur], -1);
            } else if (ft.bfs_next[fa] == cur) {
                ft.bfs_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], vol[cur]);
                ft.f_cnt_fwt.addInt(ft.dfs_in[fa], ft.dfs_out[fa], 1);
            }
            return;
        };
        
        if (dfs_st_lca != -1) {
            dfs(dfs_st_lca);
        }else {
            if (dfs_st_v != -1) dfs(dfs_st_v);
            dfs(dfs_st_y);
        }
        ft.bfs_fwt.clear();
    
        kc += path_cnt_tot;
        weight = weight * d_tau_e / d_tau;

        return;
    }

};

class TTF : public KCSolver {
public:
    Graph &graph;
    vector<vector<int>> &G;
    ll m;
    int n;
    
    int root, T;
    double eps;
    double kemeny_constant;
    // bipush for effective resistance
    double rmax;
    int T_bp;
    FixedTree fix_tree;    
    vector<TreeIdx> samples;    // used in dynamic graphs

    TTF(Graph &g, double eps,int root, int build_tree_flag = FLAG_TREE_BFS, int T = 0) : graph(g), G(g.adj) {
        n = g.n, m = g.m;
        this->eps = eps;
        this->root = root;
        
        fix_tree.build_bfs_tree(graph, root, build_tree_flag);

        if (T != 0) this->T = T;
        else this->T = ceil(log(n) / eps / eps);

        rmax = 1e-3, T_bp = 50;
        kemeny_constant = 0.;
    }

    // copy
    TTF(const TTF&a, Graph& g_up) : graph(g_up), G(g_up.adj) {
        n = g_up.n, m = g_up.m;
        root = a.root, T = a.T, eps = a.T;
        rmax = a.rmax, T_bp = a.T_bp;
        fix_tree = a.fix_tree, samples = a.samples;
    }
    
    // Bipush (Calculate ER)
    double bipush_v_absorbed_random_walk(int s, vector<double> &r_s, vector<double> &r_t) {
        double res = 0.;
        if (T_bp == 0) return res;
        vector<vector<int>> &G = graph.adj;
        for (int i = 0; i < T_bp; i++) {
            for (int cur = s; cur != root; cur = G[cur][rand() % G[cur].size()]) {
                res += (r_s[cur] - r_t[cur]) / (double) G[cur].size();
            }
        }
        
        res = res / (double) T_bp;
        
        return res;
    }

    void bipush_v_absorbed_push(int s, vector<double> &r_s, vector<double> &X_s) {
        int n = graph.n;
        // init
        vector<bool> inQueue (n + 1, false);
        queue<int> q;

        r_s[s] = 1.;
        q.push(s);
        inQueue[s] = true;
        // push
        while (!q.empty()) {
            int cur = q.front();
            q.pop();
            inQueue[cur] = false;

            X_s[cur] += r_s[cur];
            
            double inc = r_s[cur] / (double) G[cur].size();
            for (int ne : G[cur]) {
                if (ne == root) continue;

                r_s[ne] += inc;

                if (!inQueue[ne] && r_s[ne] > rmax * G[ne].size()) {
                    q.push(ne);
                    inQueue[ne] = true;
                }
            }

            r_s[cur] = 0.;
        }

        return;
    }

    double single_pair_resistance_query_bipush(int s, int t) {
        int n = graph.n;
        vector<double> X_s(n + 1, 0.), X_t(n + 1, 0.);
        vector<double> R_s(n + 1, 0.), R_t(n + 1, 0.);
        
        double res = 0.;
        if (s == t) return res;
        if (s == root) {
            // r(v, t) = (L_v^-1)tt
            bipush_v_absorbed_push(t, R_t, X_t);
            res += X_t[t] / (double) G[t].size();
            res += bipush_v_absorbed_random_walk(t, R_t, R_s);
        } else if (t == root) {
            // r(s, v) = (L_v^-1)ss
            bipush_v_absorbed_push(s, R_s, X_s);
            res += X_s[s] / (double) G[s].size();
            res += bipush_v_absorbed_random_walk(s, R_s, R_t);
        } else {
            // r(s, t)
            bipush_v_absorbed_push(s, R_s, X_s);
            bipush_v_absorbed_push(t, R_t, X_t);
            res += (X_s[s] - X_t[s]) / (double) G[s].size() + (X_t[t] - X_s[t]) / (double) G[t].size();

            res += bipush_v_absorbed_random_walk(s, R_s, R_t);
            res += bipush_v_absorbed_random_walk(t, R_t, R_s);
        }
        // cout << "r(s, t) = " << res << endl;
        return res;
    }

// Kemeny Constant
// Static Graphs
    double sample_ust() {
        int cur, intree_cnt = 1, br_tot=0;
        double kc = 0.;
        vector<int> next(n + 1);
        vector<bool> intree(n + 1, false);
        vector<Edge> branch(n);
        vector<ll> vol(n + 1, 0);
        Tree tree(n);
        
        // 1. sample ust by wilson
        intree[root] = true;        
        for (int i = 1; i <= n; i++) {
            if (intree[i]) continue;
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true;
                intree_cnt++;
                tree.add_edge(next[cur], cur);
            }
            branch[br_tot++] = make_pair(i, cur);
            if (intree_cnt == n) break;
        }
        // 2. estimate f(\tau)
        // 2.1 calculate degree_sum
        for (int i = br_tot-1; i >= 0; i--) {
            int s = branch[i].first, t = branch[i].second;
            for (cur = s; cur != t; cur = next[cur]) {
                vol[cur] += G[cur].size();
                vol[next[cur]] += vol[cur];
            }
        }
        // 2.2 compare with fixed tree to estimate 1/2m * d^T * (L_Vl)+ * d
        ll total_vol = 2 * graph.m;
        function<void(int)> dfs = [&](int x) {

            for (int i = tree.head[x]; i; i = tree.edgeNxt[i]) {
                int ne = tree.edgeTo[i];

                if (fix_tree.bfs_next[ne] == x) { // same direction
                    fix_tree.bfs_fwt.addInt(fix_tree.dfs_in[ne], fix_tree.dfs_out[ne], total_vol - vol[ne]);
                } else if (fix_tree.bfs_next[x] == ne) { // reverse dirction
                    fix_tree.bfs_fwt.addInt(fix_tree.dfs_in[x], fix_tree.dfs_out[x], vol[ne] - total_vol);
                }
                
                kc += (double) fix_tree.bfs_fwt.preSum(fix_tree.dfs_in[ne]) / total_vol * G[ne].size();
                
                dfs(ne);

                if (fix_tree.bfs_next[ne] == x) { // same direction
                    fix_tree.bfs_fwt.addInt(fix_tree.dfs_in[ne], fix_tree.dfs_out[ne], vol[ne] - total_vol);
                } else if (fix_tree.bfs_next[x] == ne) { // reverse dirction
                    fix_tree.bfs_fwt.addInt(fix_tree.dfs_in[x], fix_tree.dfs_out[x], total_vol - vol[ne]);
                }
            }
        };
        dfs(root);
        return kc;
    }

    double compute_kemeny_constant() {
        double kc = 0.;

        for (int i = 1; i <= T; i++) {
            kc = kc + sample_ust();
        }

        kemeny_constant = kc / T;
        return kemeny_constant;
    }

// Evolving Graphs
// Init
    void init() {
        kemeny_constant = 0.;
        for (int i = 0; i < T; i++) {
            samples.emplace_back(graph, fix_tree, root);

            samples[i].weight = 1. / T;
            kemeny_constant += samples[i].weight * samples[i].kc;
        }
        // cout << kemeny_constant << endl;
    }
// Update 
    void add_edge_basic_update(int u, int v) {
        double r_uv, weight_sum = 0.;
        int T_new;
        if (!graph.add_edge(u, v)) {
            return;
        }
        m = graph.m;
        
        // compute r(u, v)
        r_uv = single_pair_resistance_query_bipush(u, v);
        T_new =ceil(r_uv * T);
            // cout << "r(s, t) = " << r_uv << ", T_new size = " << T_new << endl;
        // randomly choose T_new element
        vector<int> order(T);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), default_random_engine(time(NULL)));

            // double tot_build_time = 0, tot_query_time = 0;
        // Regenerate w*N idx with edge (s,t)
        for (int i = 0; i < T; i++) {
            int &idx = order[i];
            if (i < T_new) {
                samples[idx].basic_update_sample(graph, fix_tree, u, v);   // regenerate samples
                    // tot_build_time += samples[idx].building_time;
                    // tot_query_time += samples[idx].query_time;
            } else {
                weight_sum += samples[idx].weight;
            } 
        }
            // cout << "Bsc_Add_Build_Time: " << tot_build_time << endl;
            // cout << "Bsc_Add_Query_Time: " << tot_query_time << endl;
        kemeny_constant = 0.;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[order[i]];
            if (i < T_new) { // new samples
                sample.weight = r_uv / T_new;
            } else { // origin samples
                sample.weight = sample.weight / weight_sum * (1. - r_uv);
            }
            kemeny_constant += sample.weight * sample.kc;
        }

        return;
    }
    
    void del_edge_basic_update(int u, int v) {
        // update edge
        if (!graph.del_edge(u, v))  return;

        m = graph.m;
        int update_cnt = 0;
        double weight_sum = 0.;
        vector<bool> isUpdate(T, false);

        if (fix_tree.bfs_next[u] == v || fix_tree.bfs_next[v] == u) {
            fix_tree.build_bfs_tree(graph, root);
            init();
            return;
        }

            // double tot_build_time = 0, tot_query_time = 0;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[i];
            if (sample.next[u] == v || sample.next[v] == u) {
                weight_sum += sample.weight;
                sample = TreeIdx(graph, fix_tree, root);

                isUpdate[i] = true;
                update_cnt++;
                    // tot_build_time += sample.building_time;
                    // tot_query_time += sample.query_time;
            }
        }
            // cout << "Bsc_Del_Build_Time: " << tot_build_time << endl;
            // cout << "Bsc_Del_Query_Time: " << tot_query_time << endl;
        kemeny_constant = 0.;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[i];
            if (isUpdate[i]) {
                sample.weight = weight_sum / update_cnt;
            }
            kemeny_constant += sample.weight * sample.kc;
        }
            // cout << "update tree num :" << update_cnt << endl; 
        return;
    }

    void add_edge_improved_update(int u, int v) {
        double r_uv, weight_sum_org = 0., weight_sum_new = 0.;
        int T_new;
        if (!graph.add_edge(u, v)) {
            return;
        }
        m = graph.m;
        // compute r(u, v)
        r_uv = single_pair_resistance_query_bipush(u, v);
        T_new = ceil(r_uv * (double) T);
            // cout << "r(s, t) = " << r_uv << ", T_new size = " << T_new << endl;
        // randomly choose T_new element
        vector<int> order(T);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), default_random_engine(time(NULL)));

            // double tot_build_time = 0, tot_query_time = 0;
        // Regenerate w*N idx with edge (s,t)
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[order[i]];
            if (i < T_new) {
                sample.improved_update_link_cut(graph, fix_tree, u, v);   // regenerate samples
                weight_sum_new += sample.weight;
                    // tot_build_time += sample.building_time;
                    // tot_query_time += sample.query_time;
            } 
            else {
                weight_sum_org += sample.weight;
            } 
        }
            // cout << "Imp_Add_Build_Time: " << tot_build_time << endl;
            // cout << "Imp_Add_Query_Time: " << tot_query_time << endl;
        kemeny_constant = 0.;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[order[i]];
            if (i < T_new) {  // new samples
                sample.weight = sample.weight / weight_sum_new * r_uv;
            } else {         // origin samples
                sample.weight = sample.weight / weight_sum_org * (1. - r_uv);
            }
            kemeny_constant += sample.weight * sample.kc;
        }

        return;
    }

    void del_edge_improved_update(int u, int v) {
        // update edge
        if (!graph.del_edge(u, v))  return;

        m = graph.m;
        int update_cnt = 0;
        double weight_sum_org = 0., weight_sum_new = 0.;
        vector<bool> isUpdate(T, false);

        if (fix_tree.bfs_next[u] == v || fix_tree.bfs_next[v] == u) {
            fix_tree.build_bfs_tree(graph, root);
            init();
            return;
        }
        // double tot_build_time = 0, tot_query_time = 0;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[i];
            if (sample.next[u] == v || sample.next[v] == u) {
                weight_sum_org += sample.weight;

                sample.improved_update_cut_link(graph, fix_tree, u, v);
                weight_sum_new += sample.weight;

                update_cnt++;
                isUpdate[i] = true;
                // tot_build_time += sample.building_time;
                // tot_query_time += sample.query_time;
            } 
        }
        // cout << "Imp_Del_Build_Time: " << tot_build_time << endl;
        // cout << "Imp_Del_Query_Time: " << tot_query_time << endl;

        kemeny_constant = 0.;
        for (int i = 0; i < T; i++) {
            TreeIdx &sample = samples[i];
            if (isUpdate[i]) {
                sample.weight = sample.weight / weight_sum_new * weight_sum_org;
            }
            kemeny_constant += sample.weight * sample.kc;
        }
        // cout << "update tree num :" << update_cnt << endl; 
        
        return;
    }

};


