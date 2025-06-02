#include "methods/methods.hpp"

#define STR_EQUAL(s1, s2) strcmp(s1, s2) == 0

using namespace std;

vector<string> alg_names = {"ForestMC", "TTF", "SpanTree", "LEwalk"};


void read_updates(string updates_file, vector<pair<int, int>> &ups) {
    ifstream input(updates_file);
    int  u, v;
    
    ups.clear();
    string line;
    while (getline(input, line)) {
        stringstream ss(line);

        ss >> u >> v;
        ups.push_back(make_pair(u, v));
    }

    input.close();

    std::cout << "Updates Reading Finish" << endl;
    return;
}

double read_sig(string eigen_file) {
    ifstream input(eigen_file);
    if (!input.is_open()) {
        throw std::runtime_error("Unable to open file: " + eigen_file);
    }
    double sig;
    input >> sig;
    if (input.fail()) {
        throw std::runtime_error("Failed to read eigen value from file: " + eigen_file);
    }

    return sig;
}

int find_root(Graph &g) {
    vector<vector<int>> &G = g.adj;
    int n = g.n, root, max_deg = 0;
    for (int i = 1; i <= n; i++) {
        if (G[i].size() > max_deg) {
            max_deg = G[i].size();
            root = i;
        }
    }
    return root;
}

int main(int argc, char **argv) {

    string datafolder = "./datasets/", dataset_name = "powergrid", log_name = "log.txt", alg_name = "TTF";
    string graph_file, updates_ins_file, updates_del_file, log_file, eigen_file;
    int alg_flag = ALG_TTF, T = 0;
    bool static_flag = true, dyn_imp_flag = false, dyn_bsc_flag = false;
    double epsilon = 0.5, kc = 0.;
    
    vector<pair<int, int>> ins_edges, del_edges;
    std::ofstream of_log;
    chrono::time_point<std::chrono::steady_clock> st, ed;

    srand(time(0));

    for (int i = 1; i < argc; i++) {
        if (STR_EQUAL(argv[i], "-d")) {
            dataset_name = argv[++i];   // Dataset folder name
        } else if (STR_EQUAL(argv[i], "-a")) {
            alg_name = argv[++i];
            if (STR_EQUAL(argv[i], "ttf")) alg_flag = ALG_TTF;
            else if (STR_EQUAL(argv[i], "spantree")) alg_flag = ALG_SPAN_TREE;
            else if (STR_EQUAL(argv[i], "lewalk")) alg_flag = ALG_LEWALK;
            else if (STR_EQUAL(argv[i], "forestmc")) alg_flag = ALG_FORESTMC;
            else assert(false);
        } else if (STR_EQUAL(argv[i], "-log")) {
            log_name = argv[++i];       // log filename
        } else if (STR_EQUAL(argv[i], "-eps")) {
            epsilon = atof(argv[++i]);  // Set Epsilon (Default 0.5)
        } else if (STR_EQUAL(argv[i], "-imp")) {
            dyn_imp_flag = true;     // using ImprovedSampleMaintenance (ISM) algorithm
            alg_name = "ImprovedSM";
        } else if (STR_EQUAL(argv[i], "-basic")) {
            dyn_bsc_flag = true;    // using BasicSampleMaintenance (BSM) algorithm
            alg_name = "BasicSM";
        } else if (STR_EQUAL(argv[i], "-static")) {
            static_flag = true;         // Static Graphs
        } else if (STR_EQUAL(argv[i], "-dynamic")) {
            static_flag = false;        // Dynamic Graphs
        } else assert(false);
    }

    std::cout << "Dataset: " << dataset_name << endl;

    datafolder = datafolder + dataset_name + "/";
    log_file = "./" + log_name, updates_ins_file = datafolder + "add_edges.txt", updates_del_file = datafolder + "del_edges.txt", eigen_file = datafolder + "eigen.txt";
    graph_file = datafolder + (static_flag ? "graph.txt" : "base_graph.txt");
    double sig = read_sig(eigen_file);
    of_log.open(log_file, ios::app);
    
    // Read Graph & choose node with maximum degree as root
    Graph g = Graph(graph_file);
    vector<vector<int>> &G = g.adj;
    int n = g.n, root, max_deg = 0;
    for (int i = 1; i <= n; i++) {
        if (G[i].size() > max_deg) {
            max_deg = G[i].size();
            root = i;
        }
    }

    if (static_flag) {
        st = chrono::steady_clock::now();
        
        auto solver = KCSolverFactory::create_solver(alg_flag, g, epsilon, root, sig);
        double kc = solver->compute_kemeny_constant();
        
        ed = chrono::steady_clock::now();
        of_log << dataset_name << ", " << alg_names[alg_flag] << ", " << epsilon << ", " 
        << setprecision(8) << kc << ", " << chrono::duration<double>(ed - st).count() << endl;
    } else {
        read_updates(updates_ins_file, ins_edges);
        read_updates(updates_del_file, del_edges);

        if (alg_flag == ALG_TTF && (dyn_imp_flag || dyn_bsc_flag)) {
            TTF solver = TTF(g, epsilon, root);
            solver.init();

            for (auto [u, v] : del_edges) {
                Graph updated_graph = g; // copy
                TTF solver_up = TTF(solver, updated_graph);    // copy
                // updated_graph.del_edge(u, v);

                auto st = chrono::steady_clock::now();

                if (dyn_bsc_flag) solver_up.del_edge_basic_update(u, v);
                else if (dyn_imp_flag) solver_up.del_edge_improved_update(u, v);
                double kc = solver_up.kemeny_constant;

                auto ed = chrono::steady_clock::now();
                of_log << dataset_name << ", " << alg_name << "(Del), " 
                << setprecision(8) << kc << ", " << chrono::duration<double>(ed - st).count() << endl;
            }

            for (auto [u, v] : ins_edges) {
                Graph updated_graph = Graph(g); // copy
                TTF solver_up = TTF(solver, updated_graph);    // copy
                // updated_graph.add_edge(u, v);

                auto st = chrono::steady_clock::now();

                if (dyn_bsc_flag) solver_up.add_edge_basic_update(u, v);
                else if (dyn_imp_flag) solver_up.add_edge_improved_update(u, v);
                double kc = solver_up.kemeny_constant;

                auto ed = chrono::steady_clock::now();
                of_log << dataset_name << ", " << alg_name << "(Ins), " 
                << setprecision(8) << kc << ", " << chrono::duration<double>(ed - st).count() << endl;
            }

        }else {

            for (auto [u, v] : del_edges) {
                Graph updated_graph = g;
                updated_graph.del_edge(u, v);

                auto st = chrono::steady_clock::now();

                auto solver = KCSolverFactory::create_solver(alg_flag, g, epsilon, root, sig);
                double kc = solver->compute_kemeny_constant();
                
                auto ed = chrono::steady_clock::now();
                of_log << dataset_name << ", " << alg_names[alg_flag] << ", " 
                << setprecision(8) << kc << ", " << chrono::duration<double>(ed - st).count() << endl;
            }

            for (auto [u, v] : ins_edges) {
                Graph updated_graph = g;
                updated_graph.add_edge(u, v);

                auto st = chrono::steady_clock::now();

                auto solver = KCSolverFactory::create_solver(alg_flag, g, epsilon, root, sig);
                double kc = solver->compute_kemeny_constant();
                
                auto ed = chrono::steady_clock::now();
                of_log << dataset_name << ", " << alg_names[alg_flag] << ", " 
                << setprecision(8) << kc << ", " << chrono::duration<double>(ed - st).count() << endl;
            }
            
        }
    }
    
    of_log.close();
    return 0;
}