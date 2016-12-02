#include "solveq.h"

using namespace arma;

namespace armafem {

bool solveq(const mat &K, const colvec &f, const uvec& pdof, const colvec &dp,
            colvec &d, colvec &Q) {
     auto nd = K.n_rows;
     uvec mask = ones<uvec>(nd);
     mask.elem(pdof).zeros();
     auto fdof = find(mask > 0);

     d=zeros(nd,1);
     Q=zeros(nd,1);
     
     mat Kr = K.submat(fdof, fdof);
     mat Kd = K.submat(fdof, pdof);
     colvec fr = f.elem(fdof) - Kd*dp;

     colvec s = solve(Kr, fr);

     d.elem(pdof)=dp;
     d.elem(fdof)=s;
     Q=K*d-f;
  
     return true;
}

} // namespace armafem
