#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <chrono>
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <cassert>
#include <random>
#include <algorithm>

using namespace std;
using Edge = pair<int, int>;

typedef int ll;

class Graph {
public:

    vector<vector<int>> adj;
    vector<int> degree;
    string graph_file;
    int n;
    ll m;
    
    // vector<Edge> edges;
    // unordered_map<int, unordered_map<int, ll>> edge_idx;

    Graph() {
    } 

    Graph(int n, int m) {
        this->n = n, this->m = m;
        adj.resize(n + 1);
        degree.resize(n + 1);
    }

    Graph(string filename) {
        graph_file = filename;
        read_graph_undirect(filename);
    }

    bool add_edge(int u, int v) {
        if (find(adj[u].begin(), adj[u].end(), v) != adj[u].end()) {
            printf("Add Edge (%d, %d) Error!\n", u, v);
            printf("Edge is already existed!\n");
            return false;
        }

        adj[u].push_back(v);
        adj[v].push_back(u);
        
        m++;

        return true;
    }

    bool del_edge(int u, int v)
    {
        auto it_u = find(adj[u].begin(), adj[u].end(), v);
        auto it_v = find(adj[v].begin(), adj[v].end(), u);

        if (it_u == adj[u].end() || it_v == adj[v].end())
        {
            printf("Del Edge (%d, %d) Error!\n", u, v);
            printf("Edge is not existed!\n");
            return false;
        }

        adj[u].erase(it_u);
        adj[v].erase(it_v);

        m--;

        return true;
    }

    bool del_vertex(int v) {
        if (adj[v].size() == 0) {
            cout << v << " has no edges" << endl;
            return false;
        }
        
        for (int ne : adj[v]) {
            auto it = find(adj[ne].begin(), adj[ne].end(), v);
            adj[ne].erase(it);
        }
        adj[v].clear();

        return true;
    }

    void read_graph_undirect(string filename) {
        FILE *fin = fopen(filename.c_str(), "r");
        //printf("Start Reading...\n");
        if (fin == NULL) {
            printf("Reading Graph Failed: can not find file: %s\n", filename.c_str());
            assert(false);
        }
        int t1, t2;
        if (fscanf(fin, "%d%d", &t1, &t2) == EOF) {
            printf("Read Error: not start with (n, m).\n");
            return;
        }
        n = t1;
        adj = vector<vector<int>>(n + 1, vector<int>());

        m = 0;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            adj[t1].push_back(t2);
            adj[t2].push_back(t1);
            m++;
        }

        fclose(fin);
        //printf("Read Graph %s Successed!\n", filename.c_str());
        // cal_degree();
        return;
    }

    void cal_degree() {
        degree.resize(n + 1);
        for (int i = 1; i <= n; i++) {
            assert(adj[i].size() > 0);
            degree[i] = adj[i].size();
        }
        return;
    }
};

template<class T = int>
class FenwickTree {
    int n;
    vector<T> c;
public:

    FenwickTree() : n(-1) {}
    FenwickTree(int n) : n(n), c(n + 2) {}

    inline int lowbit(int x) {
        return x & (-x);
    }

    void add(int i, T k) {
        while (i <= n) {
            c[i] = c[i] + k;
            i += lowbit(i);
        }
    }

    void addInt(int l, int r, T k) {
        // assert(l <= r);
        add(l, k);
        add(r + 1, -k);
    }

    T preSum(int i) {
        T res = 0;
        while (i > 0) {
            res = res + c[i];
            i -= lowbit(i);
        }
        return res;
    }

    void clear() {
        if (n > 0) c = vector<T>(n + 2, 0);
    }

    size_t getMemoryUsage() {
        size_t res = 0;
        res += sizeof(c) + c.size() * sizeof(T);
        return res;
    }
};

class Tree {
public:
    int n, tot;
    vector<int> head, edgeTo, edgeNxt;
    Tree() {
        n = tot = -1;
    }

    Tree(int n) {
        this->n = n;
        tot = 0;
        head.resize(n+1, 0), edgeTo.resize(n, 0), edgeNxt.resize(n+1, 0);
    }

    void add_edge(int x, int y) {
        tot++;
        edgeTo[tot] = y;
        edgeNxt[tot] = head[x];
        head[x] = tot;
    }

    // change (u, v) to (x, y)
    void change_edge(int u, int v, int x, int y) {
        int idx = -1, last_idx = -1;
        // find idx of edge (u, v)
        for (int i = head[u]; i; i = edgeNxt[i]) {
            if (edgeTo[i] == v) {
                idx = i;
                break;
            } 
            last_idx = i;
        }
        if (idx == -1) {
            cout << "Error: no edge in Tree" << endl;
            assert(false);
        }
        // del edge (u, v)
        if (last_idx == -1) {
            head[u] = edgeNxt[idx];
        }else {
            edgeNxt[last_idx] = edgeNxt[idx];
        }
        // add edge (x, y)
        edgeTo[idx] = y;
        edgeNxt[idx] = head[x];
        head[x] = idx;

        return;
    }

};

#endif