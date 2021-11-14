//
// Created by valentin on 14.11.21.
//

#ifndef NETWORKIT_SPGEMM_H
#define NETWORKIT_SPGEMM_H

#include <map>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/Globals.hpp>



namespace NetworKit {

class SPA {

    SPA();

    ~SPA();

    void accumulate(double value, index pos);

    index output(double & value, std::vector<index> & col);

private:

    std::map<index, double> spa;
};
class SpGEMM {
    SpGEMM();

    ~SpGEMM();

    CSRMatrix compute(const CSRMatrix & a, const CSRMatrix & b);

};

}


#endif //NETWORKIT_SPGEMM_H
