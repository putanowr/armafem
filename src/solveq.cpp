#include "solveq.h"

using namespace arma;

namespace armafem {

void solveq(const mat &K, const colvec &f, const uvec& fixId, const colvec &fixVal,
            colvec &d, colvec &Q) {
     auto nd = K.n_rows;
     uvec mask = ones<uvec>(nd);
     mask.elem(fixId).zeros();
     auto activeId = find(mask > 0);

     d=zeros(nd,1);
     Q=zeros(nd,1);
     
     mat Kr = K.submat(activeId, activeId);
     mat Kd = K.submat(activeId, fixId);
     colvec fr = f.elem(activeId) - Kd*fixVal;

     colvec s = solve(Kr, fr);

     d.elem(fixId)=fixVal;
     d.elem(activeId)=s;
     Q=K*d-f;
}

} // namespace armafem
