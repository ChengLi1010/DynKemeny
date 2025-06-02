#pragma once

#include "../graph.hpp"
#include "Base.hpp"
#include <cmath>

class LEWalk : public KCSolver{
public:
    Graph &graph;
    vector<vector<int>> &G;
    int landmark, n;          // the largest degree vertex s
    int T;
    double kc, kemeny_constant;
    double eps, sig;
    std::mt19937 gen;
    std::discrete_distribution<> dist;

    LEWalk(Graph &g, double eps, int landmark, int T = 0): graph(g), G(g.adj) {
        this->eps = eps;
        this->landmark = landmark;
        n = g.n;
        if (T == 0)
            T = ceil(log(n) * log(n) / eps / eps);
        this->T = T;

        g.cal_degree();
        std::random_device rd;
        gen = std::mt19937(rd());
        dist = std::discrete_distribution<> (graph.degree.begin() + 1, graph.degree.end());
    }

    ll loop_erased_random_walk() {
        int cur, cnt = 1;
        ll y1 = 0;
        vector<bool> intree(n + 1, false);
        vector<int> next(n + 1, 0);
        intree[landmark] = true;
        
        for (int i = 1; i <= n; i++) {
            for (cur = i; !intree[cur]; cur = next[cur]) {
                next[cur] = G[cur][rand() % G[cur].size()];
                y1++;
            }
            for (cur = i; !intree[cur]; cur = next[cur]) {
                intree[cur] = true;
                cnt++;
            }
            if (cnt == n) break;
        }
        return y1;
    }

    ll hitting_time() {
        ll y2 = 0;
        int s = dist(gen) + 1, cur;

        cur = s;
        while(cur != landmark) {
            y2++;
            cur = G[cur][rand() % G[cur].size()];
        }
        return y2;
    }

    double compute_kemeny_constant() {
        kc = 0.;
        int n = graph.n, cur;
        
        // >>>>> Lewalk <<<<<<
        for (int i = 1; i <= T; i++) {
            ll est_i = loop_erased_random_walk() - hitting_time();
            kc = kc  + (double) est_i;
        }
        kemeny_constant = kc / T;
        return kemeny_constant;
    }

};