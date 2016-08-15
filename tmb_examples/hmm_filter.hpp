// Utilities to discretize a SDE and calculate the likelihood function.

/*
  Finite volume discritize advection diffusion equation. Assuming
  equidistant grid and using central difference scheme for advection.
*/
template<template <typename> class sde_t, class Type>
struct fvade_t{
  matrix<Type> A;
  vector<Type> grid;
  fvade_t(sde_t<Type> sde, vector<Type> grid){
    this->grid = grid;
    int n    = grid.size();
    int nvol = n - 1;
    Type h   = grid[1] - grid[0];
    /* Setup Laplacian term */
    vector<Type> D(n);
    for (int i = 0; i < n; i++)
      D[i] = Type(.5) * pow(sde.dispersion(grid[i]),2);
    vector<Type> Dinner = D.segment(1,n-2);
    matrix<Type> L(nvol,nvol);
    L.setZero();
    L.diagonal(-1) = Dinner;
    L.diagonal( 1) = Dinner;
    L.diagonal()   = -L.colwise().sum();
    L = L / (h*h);
    /* Setup advection term */
    vector<Type> xm = Type(.5) * (grid.head(nvol) + grid.tail(nvol));
    vector<Type> v(nvol);
    for (int i = 0; i < nvol; i++)
      v[i] = sde.advection(xm[i]);
    matrix<Type> G(nvol,nvol);
    G.setZero();
    G.diagonal(-1) = -Type(.5) * v.head(nvol - 1);
    G.diagonal( 1) =  Type(.5) * v.tail(nvol - 1);
    G.diagonal()   = -G.colwise().sum();
    G = G / h;
    /* Setup generator */
    A = (-G + L);
  }
};
template<template <typename> class sde_t, class Type>
fvade_t<sde_t, Type> fvade(sde_t<Type> sde, vector<Type> grid){
  return fvade_t<sde_t, Type>(sde,grid);
}

/*
  Evaluate negative log-likelihood of hidden Markov model, specified
  through the (time-constant) generator matrix of the
  process. Assuming equidistant measurements in time.
*/
template<class Type>
struct hmm_filter{
  matrix<Type> P;  // Transition prob (transposed)
  matrix<Type> P0; // Observation error (transposed)
  Type dt;
  int n;
  vector<Type> grid;
  /* Input: 
     A:  Generator
     dt: Timestep between measurements
  */
  void init(matrix<Type> A, vector<Type> grid, Type dt){
    P = atomic::expm(matrix<Type>(A*dt));
    int n = A.rows();
    P0 = matrix<Type>(n, n);
    P0.setIdentity(); // Default: no obs error
    this->dt = dt;
    this->n = n;
    this->grid = grid;
  }
  hmm_filter(matrix<Type> A, vector<Type> grid, Type dt){
    init(A,grid,dt);
  }
  template<class fvade_t>
  hmm_filter(fvade_t fvol, Type dt){
    init(fvol.A,fvol.grid,dt);
  }
  void setGaussianError(Type sd){
    vector<Type> xm = Type(.5) * (grid.head(n) + grid.tail(n));
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	P0(i,j) = dnorm(xm[j],xm[i],sd,false);
    // Normalize
    vector<Type> cs = P0.colwise().sum();
    for(int i=0;i<n;i++) P0.col(i) /= cs[i];
  }
  /* Update step */
  vector<Type> multiply(matrix<Type> x, vector<Type> y){
    return atomic::matmul(x, matrix<Type>(y.matrix())).vec();
  }
  void update(vector<Type> &px, Type &nll, int yobs){
    // Update joint distribution (J) of state and measurement, and
    // extract 'yobs'-slice:
    vector<Type> Ppx = multiply(P, px);
    vector<Type> Jslice = Ppx * vector<Type>(P0.row(yobs));
    // Integrated slice is yobs-likelihood:
    Type Ly = Jslice.sum();
    nll -= log(Ly);
    // Normalized slice is updated px:
    px = Jslice / Ly;
  }
  /* Evaluate negative log likelihood of observations */
  Type operator()(vector<int> obs){
    vector<Type> px(n);
    px.setZero();
    px += 1.0 / Type(px.size()); // Uniform initial distribution
    px /= grid[1]-grid[0];       // prob -> density
    Type nll = 0;
    for(int k=0; k<obs.size(); k++){
      update(px, nll, obs[k]);
    }
    return nll;
  }
};
