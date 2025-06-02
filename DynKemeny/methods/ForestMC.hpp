#pragma once

#include "../graph.hpp"
#include "Base.hpp"
#include <cmath>

class ForestMC : public KCSolver{
public:
    Graph &graph;
    vector<vector<int>> &G;
    int root, n;          // the largest degree vertex s
    int T;
    double kc, kemeny_constant;
    double eps, sig;

    ForestMC(Graph &g, double eps, int root, double sig = 0.9995, int T = 0): graph(g), G(g.adj) {
        this->eps = eps;
        this->root = root;
        this->sig = sig;

        n = g.n;

        if (T != 0)
            this->T = T;
        else
            this->T = ceil(log(n) * log(n) / eps / eps);
    }

    uint64_t truncated_random_walk(int start, int len) {
        int cur = start;
        uint64_t hits = 0;
        while(len--) {
            cur = G[cur][rand() % G[cur].size()];
            if (cur == start) hits++;
        }
        return hits;
    }

    double terminate_function (uint64_t n, double var, double sup, double delta) {
        return sqrt(2. * var * log(3. / delta) / n) + 3. * sup * log(3./delta) / n;
    }

    double compute_kemeny_constant() {
        kc = 0.;
        int n = graph.n, cur;
        uint64_t t_st = 0, i, rw_itr = 0, st_itr = 0;
        uint64_t l = ceil(log(2. / (eps - eps * sig)) / log(1. / sig));
        
        // >>>>> Lewalk <<<<<<
        for (st_itr = 1; st_itr <= T; st_itr++) {
            vector<bool> intree(n + 1, false);
            vector<int> next(n + 1, 0);
            intree[root] = true;

            for (int j = 1; j <= n; j++) {
                for (cur = j; !intree[cur]; cur = next[cur]) {
                    t_st++;
                    next[cur] = G[cur][rand() % G[cur].size()];
                }
                for (cur = j; !intree[cur]; cur = next[cur]) {
                    intree[cur] = true;
                }
            }
        }
        kc += (double)t_st / T;

        // >>>>> Return Time <<<<<
        double sum = 0, sqr_sum = 0;
        for (i = 1; i <= T; i++) {
            int hits = truncated_random_walk(root, l);
            
            sum += hits, sqr_sum += hits * hits;
            double var = (sqr_sum - sum * sum / (double) i) / (double) i;

            if (terminate_function(i, var, (double) l / 2., (double) 1. / n) <= (double) n * eps / 3.) break;
        }
        rw_itr = i;
        double t_s_l = sum / rw_itr; 

        kc -= (double) ((t_s_l + 1) * 2 * graph.m / G[root].size() - l);
        
        return kc;
    }

};