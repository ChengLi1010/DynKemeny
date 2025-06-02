#ifndef BASE_HPP
#define BASE_HPP

class KCSolver {
public:
    virtual ~KCSolver() = default;
    virtual double compute_kemeny_constant() = 0;
};  

#endif
