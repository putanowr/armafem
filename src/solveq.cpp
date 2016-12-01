bool solveq(const mat & K, const colvec &f, const uvec& pdof, const colvec &dp,
            colvec &d, colvec &Q) {
     int nd = K.n_rows;
     uvecfdof = linspace<uvec>(0, n-1, 1);

     d=zeros(n,1);
     Q=zeros(n,1);

     for (colvec::iterator i=pdof.begin(); 
                           i!= pdof.end(); i++) {
       fdof.shed_row(*i); 
     } 

     mat Kr = K.submat(fdof, fdof);
     mat Kd = K.submat(fdof, pdof);
     colvec fr = f.elem(fdof) - Kd*dp;

     colvec = solve(Kr, fr);

     d.elem(pdof)=dp;
     d.elem(fdof)=s;
     Q=K*d-f;
  
     return true;
}
