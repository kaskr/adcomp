// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

/**
   \brief Automatic differentiation (gradient, hessian, jacobian)

   - Gives the TMB user easy access to the CppAD library
*/
namespace autodiff {
  /* Gradient */
  template<class Functor, class Type>
  struct gradient_t
  {
    Functor f;
    gradient_t(Functor f_) : f(f_) {}
    AD<Type> userfun(vector<AD<Type> > x){
      return f(x);
    }
    vector<Type> operator()(vector<Type> x0){
      CppAD::vector<AD<Type> > x( x0 );
      CppAD::vector<AD<Type> > y( 1 );
      CppAD::Independent(x);
      y[0] = userfun(x);
      CppAD::ADFun<Type> F(x, y);
      CppAD::vector<Type> x_eval(x0);
      return F.Jacobian(x_eval);
    }
  };

  /** \brief Calculate gradient of vector function with scalar values

      \param F Functor with scalar values
      \param x Vector evaluation point
      \return Gradient vector

      \note The evaluation method of the functor must be templated

      Example:

      \code
      #include <TMB.hpp>

      struct func {
        template <class T>
        T operator()(vector<T> x){  // Evaluate function
          return exp( -x.sum() );
        }
      };

      template<class Type>
      Type objective_function<Type>::operator() () {
        PARAMETER_VECTOR(theta);
	func f;
	// Calculate gradient
	vector<Type> g = autodiff::gradient(f, theta);
	REPORT(g);
	// Exit
	return 0;
      }
      \endcode
  */
  template<class Functor, class Type>
  vector<Type> gradient(Functor F, vector<Type> x){
    gradient_t<Functor, Type> f(F);
    return f(x);
  }

  /* Hessian */
  template<class Functor, class Type>
  struct hessian_t
  {
    Functor f;
    gradient_t<Functor,AD<Type> > gr;
    hessian_t(Functor f_) : f(f_), gr(f_) {}
    vector<AD<Type> > userfun(vector<AD<Type> > x){
      return gr(x);
    }
    matrix<Type> operator()(vector<Type> x0){
      CppAD::vector<AD<Type> > x( x0 );
      CppAD::vector<AD<Type> > y( x0 );
      CppAD::Independent(x);
      y = userfun(x);
      CppAD::ADFun<Type> F(x, y);
      CppAD::vector<Type> x_eval(x0);
      vector<Type> ans = F.Jacobian(x_eval);
      return asMatrix(ans, x.size(), x.size());
    }
  };

  /** \brief Calculate hessian of vector function with scalar values

      \param f Functor with scalar values
      \param x Vector evaluation point
      \return Hessian matrix

      \note The evaluation method of the functor must be templated

      Example:

      \code
      #include <TMB.hpp>

      struct func {
        template <class T>
        T operator()(vector<T> x){  // Evaluate function
          return exp( -x.sum() );
        }
      };

      template<class Type>
      Type objective_function<Type>::operator() () {
        PARAMETER_VECTOR(theta);
	func f;
	// Calculate hessian
	matrix<Type> h = autodiff::hessian(f, theta);
	REPORT(h);
	// Exit
	return 0;
      }
      \endcode
  */
  template<class Functor, class Type>
  matrix<Type> hessian(Functor f, vector<Type> x){
    hessian_t<Functor, Type> H(f);
    return H(x);
  }

  /* Jacobian */
  template<class Functor, class Type>
  struct jacobian_t
  {
    Functor f;
    jacobian_t(Functor f_) : f(f_) {}
    vector<AD<Type> > userfun(vector<AD<Type> > x){
      return f(x);
    }
    matrix<Type> operator()(vector<Type> x0){
      CppAD::vector<AD<Type> > x( x0 );
      CppAD::Independent(x);
      CppAD::vector<AD<Type> > y = userfun(x);
      CppAD::ADFun<Type> F(x, y);
      CppAD::vector<Type> x_eval(x0);
      vector<Type> ans = F.Jacobian(x_eval); // By row !
      return asMatrix(ans, x.size(), y.size()).transpose();
    }
  };

  /** \brief Calculate jacobian of vector function with vector values

      \param f Functor with vector values
      \param x Vector evaluation point
      \return Jacobian matrix

      \note The evaluation method of the functor must be templated

      Example:

      \code
      #include <TMB.hpp>

      struct func {
        template <class T>
        vector<T> operator()(vector<T> x){  // Evaluate function
	  vector<T> y(2);
	  y(0) = x.sum();
	  y(1) = x.prod();
	  return y;
        }
      };

      template<class Type>
      Type objective_function<Type>::operator() () {
        PARAMETER_VECTOR(theta);
	func f;
	// Calculate jacobian
	matrix<Type> j = autodiff::jacobian(f, theta);
	REPORT(j);
	// Exit
	return 0;
      }
      \endcode
  */
  template<class Functor, class Type>
  matrix<Type> jacobian(Functor f, vector<Type> x){
    jacobian_t<Functor, Type> J(f);
    return J(x);
  }


  /* Sparse Jacobian */
  template<class Functor, class Type>
  struct sparse_jacobian_t
  {
    Functor f;
    sparse_jacobian_t(Functor f_) : f(f_) {}
    vector<AD<Type> > userfun(vector<AD<Type> > x){
      return f(x);
    }
    Eigen::SparseMatrix<Type> operator()(vector<Type> x0) {
      /* Tape user functor F(x, y) */
      CppAD::vector<AD<Type> > x( x0 );
      CppAD::Independent(x);
      CppAD::vector<AD<Type> > y = userfun(x);
      CppAD::ADFun<Type> F(x, y);
      CppAD::vector<Type> x_eval(x0);
      /* Evaluate Sparse Jacobian */
      typedef vector<int> SizeVector;
      int nx = x.size();
      int ny = y.size();
      vector<bool> keep_col(nx); // Jacobian cols to get
      for(int i = 0; i < nx; i++) {
        keep_col[i] = true;
      }
      vector<Type> u(nx);
      // Prepare output
      typedef Eigen::Triplet<Type> T;
      std::vector<T> tripletList;
      // Do sweeps
      F.Forward(0, x_eval);
      SizeVector row_pattern(0);
      // Mark domain dependencies and initialize
      F.subgraph_reverse(keep_col);
      for(int r = 0; r < ny; r++) {
        // Do reverse sub sweep
        F.subgraph_reverse(1,           // order
                           r,           // range component
                           row_pattern, // valid rows
                           u);          // domain vector
        // Append value and integer pairs
        for (int k=0; k < row_pattern.size(); k++) {
          int c = row_pattern[k];
          tripletList.push_back( T(r, c, u[c]) );
        }
      }
      // Output
      Eigen::SparseMatrix<Type> mat(ny, nx);
      mat.setFromTriplets(tripletList.begin(), tripletList.end());
      return mat;
    }
  };

  /** \brief Calculate sparse jacobian of vector function with vector values

      \param f Functor with vector values
      \param x Vector evaluation point
      \return Sparse Jacobian matrix

      \note The evaluation method of the functor must be templated

      Example:

      \code
      #include <TMB.hpp>

      struct func {
        template <class T>
        vector<T> operator()(vector<T> x){  // Evaluate function
          vector<T> y(2);
          y(0) = x.sum();
          y(1) = x.prod();
          return y;
        }
      };

      template<class Type>
      Type objective_function<Type>::operator() () {
        PARAMETER_VECTOR(theta);
        func f;
        // Calculate jacobian
        Eigen::SparseMatrix<Type> j = autodiff::sparse_jacobian(f, theta);
        REPORT(j);
        // Exit
        return 0;
      }
      \endcode
  */
  template<class Functor, class Type>
  Eigen::SparseMatrix<Type> sparse_jacobian(Functor f, vector<Type> x){
    sparse_jacobian_t<Functor, Type> J(f);
    return J(x);
  }
}
