//
// Created by valentin on 14.11.21.
//

#include "SpGEMM.h"

namespace Networkit{

void SPA::accumulate(double value, index pos){
    if(!spa.contains(pos)){
        map.insert(value);
    }
    else{
        map[pos] += value;
    }
}

index SPA::output(std::vector<double> & val, std::vector<index> & col){
    index nzi = 0;
    for(std::iterator cptr = spa.begin(); cptr != spa.end();cptr++; nzi++){
        col.append(cptr->first);
        val.append(cptr->second);
    }
    return nzi;
}

CSRMatrix SpGEMM::compute(const CSRMatrix & a, const CSRMatrix & b) {
    CSRMatrix c = CSRMatrix();
    c.append(0);
    for(index i = 0; n < c.numberOfRows(), ++n){
        SPA spa = SPA();
        for(index k = a.rowIdx[i]; k < a.rowIdx[i+1]; ++k){
            for(index j = b.rowIdx[k]; j < b.rowIdx[k+1]; ++j){
                double value = a.nonZeros[k]*b.nonZeros[j];
                spa.accumulate(value, b.columnIdx[j])
            }
        }
        index nzi = spa.output(c.columnIdx, c.nonZeros);
        c.rowIdx.append(c.rowIdx[i]+nzi);
    }

}
}
