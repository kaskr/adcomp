CppAD::vector<double> D_incpl_gamma_shape(CppAD::vector<double> x);
CppAD::vector<AD<double > > D_incpl_gamma_shape(CppAD::vector<AD<double> > x);
CppAD::vector<AD<AD<double> > > D_incpl_gamma_shape(CppAD::vector<AD<AD<double> > > x);
CppAD::vector<AD<AD<AD<double> > > > D_incpl_gamma_shape(CppAD::vector<AD<AD<AD<double> > > > x);

template <class Type>
class atomic_D_incpl_gamma_shape : public CppAD::atomic_base<Type> {
public:
  atomic_D_incpl_gamma_shape(const char* name) : CppAD::atomic_base<Type>(name){
    std::cout << "Constructing atomic " << "D_incpl_gamma_shape" << "\n" ;
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }
private:
  virtual bool forward(size_t p,
		       size_t q,
		       const CppAD::vector<bool>& vx,
		       CppAD::vector<bool>& vy,
		       const CppAD::vector<Type>& tx,
		       CppAD::vector<Type>& ty
		       )
  {
    if( vx.size() > 0 ){
      vy[0] = vx[0];
    }
    //ATOMIC_FORWARD;
    ty[0] = D_incpl_gamma_shape(tx)[0];
    return true;
  }
  virtual bool reverse(size_t q,
		       const CppAD::vector<Type>& tx,
		       const CppAD::vector<Type>& ty,
		       CppAD::vector<Type>& px,
		       const CppAD::vector<Type>& py
		       )
  {
    //ATOMIC_REVERSE;
    px[0] = exp(-tx[0])*pow(tx[0],tx[1]-Type(1.0))*pow(log(tx[0]),tx[2]) * py[0];
    CppAD::vector<Type> tx_(tx);
    tx_[2] = tx_[2] + Type(1.0);  // Add one to get partial wrt. tx[1]
    px[1] = D_incpl_gamma_shape(tx_)[0] * py[0];
    px[2] = Type(0);
    return true;
  }
  virtual bool rev_sparse_jac(size_t q,
			      const CppAD::vector<bool>& rt,
			      CppAD::vector<bool>& st)
  {
    for(size_t i=0;i<st.size();i++)st[i]=true;
    return true;
  }
  virtual bool rev_sparse_jac(size_t q,
			      const CppAD::vector< std::set<size_t> >& rt, 
			      CppAD::vector< std::set<size_t> >& st)
  {
    error("Should not be called");
  }
};

CppAD::vector<double> D_incpl_gamma_shape(CppAD::vector<double> vx){
  vector<double> vy(1);
  vy[0]=Rmath::D_incpl_gamma_shape(vx[0],vx[1],vx[2]);
  return vy;
}

atomic_D_incpl_gamma_shape<double> afun1_D_incpl_gamma_shape("atomic_" "D_incpl_gamma_shape"); 
CppAD::vector<AD<double > > D_incpl_gamma_shape(CppAD::vector<AD<double> > vx){
  CppAD::vector<AD<double> > vy(1);
  afun1_D_incpl_gamma_shape(vx,vy);
  return vy;
}

atomic_D_incpl_gamma_shape<AD<double> > afun2_D_incpl_gamma_shape("atomic_" "D_incpl_gamma_shape"); 
CppAD::vector<AD<AD<double> > > D_incpl_gamma_shape(CppAD::vector<AD<AD<double> > > vx){
  CppAD::vector<AD<AD<double> > > vy(1);
  afun2_D_incpl_gamma_shape(vx,vy);
  return vy;
}

atomic_D_incpl_gamma_shape<AD<AD<double> > > afun3_D_incpl_gamma_shape("atomic_" "D_incpl_gamma_shape");
CppAD::vector<AD<AD<AD<double> > > > D_incpl_gamma_shape(CppAD::vector<AD<AD<AD<double> > > > vx){
  CppAD::vector<AD<AD<AD<double> > > > vy(1);
  afun3_D_incpl_gamma_shape(vx,vy);
  return vy;
}

