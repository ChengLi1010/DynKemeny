#include "TTF.hpp"
#include "SpanTree.hpp"
#include "ForestMC.hpp"
#include "LEWalk.hpp"

#define ALG_FORESTMC 0
#define ALG_TTF  1
#define ALG_SPAN_TREE  2
#define ALG_LEWALK 3

class KCSolverFactory {
public:
    static std::unique_ptr<KCSolver> create_solver(int alg_flag, Graph& g, double epsilon, int root, 
                                                double sig = 0.0) {
        switch (alg_flag) {
            case ALG_TTF:
                return std::make_unique<TTF>(g, epsilon, root);
            case ALG_SPAN_TREE:
                return std::make_unique<SpanTree>(g, epsilon, root);
            case ALG_LEWALK:
                return std::make_unique<LEWalk>(g, epsilon, root);
            case ALG_FORESTMC:
                return std::make_unique<ForestMC>(g, epsilon, root, sig);
            default:
                throw std::invalid_argument("Unsupported algorithm flag");
        }
    }
};