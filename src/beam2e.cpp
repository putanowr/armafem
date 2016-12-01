#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

bool solveq(const mat & K, const colvec &f, const uvec& pdof, const colvec &dp,
            colvec &d, colvec &Q) {
     int nd = K.n_rows;
     uvec adof = ones<uvec>(nd);

     uvec zdof = linspace<uvec>(0, nd-1, nd);

     d=zeros(nd,1);
     Q=zeros(nd,1);

     for (uvec::const_iterator i=  pdof.begin(); 
                           i!= pdof.end(); ++i) {
       adof(*i) = 0; 
     } 

     uvec fdof = zdof.elem(find(adof));

     mat Kr = K.submat(fdof, fdof);
     mat Kd = K.submat(fdof, pdof);
     colvec fr = f.elem(fdof) - Kd*dp;

     colvec s = solve(Kr, fr);

     d.elem(pdof)=dp;
     d.elem(fdof)=s;
     Q=K*d-f;
  
     return true;
}

bool assem(umat edof, mat &K, const mat &Ke, colvec &f, const colvec &fe) {
  int nie = edof.n_rows;
  int n = edof.n_cols;
  umat t = edof.cols(1,n-1);
  for (int i=0; i<nie; i++) {
    K.submat(t.row(i), t.row(i)) += Ke;
    f.elem(t.row(i)) += fe;
  }
  return true;
}


bool beam2e(const rowvec &ex, const rowvec &ey, const rowvec &ep, const rowvec &eq, 
            mat &Ke, colvec &fe) {
  rowvec b(2);
  b << (ex(1)-ex(0)) << (ey(1) - ey(0)) << endr; 
  double L = norm(b,2);
  rowvec n = b/L;
  double E=ep(0);  double A=ep(1);  double I=ep(2);
 
  double qx=eq(0); 
  double qy=eq(1);

  mat Kle = zeros(6,6);

  double EIL3 = E*I/(L*L*L);
  double EIL2 = E*I/(L*L);
  double EIL = E*I/L;
  double EAL = E*A/L; 
   
  Kle << EAL <<  0           <<  0   << -EAL <<  0       << 0       << endr
      <<  0   <<  12*EIL3  <<  6*EIL2 << 0    << -12*EIL3 << 6*EIL2  << endr
      <<   0   <<  6*EIL2   <<  4*EIL  << 0    << -6*EIL2  << 2*EIL   << endr
      << -EAL  <<  0           <<  0   <<  EAL <<  0       << 0       << endr
      <<   0   <<  -12*EIL3 << -6*EIL2 << 0    << 12*EIL3  << -6*EIL2 << endr
      <<   0   <<  6*EIL2   <<  2*EIL  << 0    << -6*EIL2  <<  4*EIL << endr;
   
  double qxf = qx*L/2;
  double qyf = qy*L/2;
  double qym = qy*L*L/12;
  
  colvec fle(6);

  fle << qxf << qyf << qym << qxf << qyf << -qym << endr;
  
  mat G = zeros(6,6);
  G(0,0) = G(3,3) = n(0);
  G(0,1) = G(3,4) = n(1);
  G(1,0) = G(4,3) = -n(1);
  G(1,1) = G(4,4) = n(0);
  G(2,2) = G(5,5) = 1.0;

  mat GT = trans(G);
   
  Ke=GT*Kle*G;   fe=GT*fle; 

  return true;
}

#ifdef TEST_IT
int main(int argc, char *argv[]) {
  rowvec ex(2);
  rowvec ey(2);
  rowvec ep(2);
  rowvec eq(2);

  ex << 2 << 4 << endr;
  ey << 1 << 1 << endr;
  ep << 1 << 1 << 1 << endr;
  eq << 0 << 4 << endr; 

  mat Ke;
  colvec fe;
  beam2e(ex,ey,ep,eq, Ke, fe);

  umat Edof(6,7);
  std::cout << "Ile : " << Edof.n_rows << endl;
  for (int i=0; i<Edof.n_rows; i++) {
     Edof << i;
  std::cout << "Ile : " << Edof.n_rows << endl;
     Edof << (3*i); 
     Edof << (3*i+1);
     Edof << (3*i+2);
     Edof << (3*i+3);
     Edof << (3*i+4);
     Edof << (3*i+5);
     Edof << endr; 
  }

  int n=3*(Edof.n_rows+1);
  std::cout << "Size: " << n << "\n";
  mat K = zeros(n,n);
  colvec f = zeros(n,1);

  for (int i=0; i<Edof.n_rows; i++) {
     assem(Edof.row(i), K, Ke, f, fe);
  }

  colvec s, Q;
  uvec bdof(3); 
  bdof << 0 << 1 << n-2;
  colvec bval = zeros(3);
  
  solveq(K, f, bdof, bval, s, Q);

  std::cout << s << endl;  

  return 0;
}
#endif
