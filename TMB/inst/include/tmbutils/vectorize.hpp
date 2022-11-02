#ifdef TMBAD_FRAMEWORK

/** \file
    \brief Atomic vectorization
*/

/** \brief Vectorized atomic functions

    This namespace provides a vector class (`vectorize::vector`) for
    which 'native' AD vectorization is available.  A *native*
    vectorized operator is represented by a single node in the
    computational graph. Note that this is different from *syntactic*
    vectorization as implemented by the `Vectorize.hpp` macros.

    For an overview of vectorized methods see `TMBad/vectorize.hpp`.

    \note PROS/CONS: Native vectorization is fast and memory
    efficient. However, some optimization tricks may not be as
    effective (e.g. sparsity).
*/

namespace vectorize {

/** \brief Vector of consequtive values on the AD tape */
template <class Type>
struct vector : TMBad::ad_segment {
  typedef TMBad::ad_segment Base;
  using Base::Base;
  vector() {} // Surprisingly, needed on some systems.
  vector(tmbutils::vector<Type> x) :
    Base(x.data(), x.size()) { }
  vector(const Base &x) :
    Base(x) { }
  operator std::vector<Type>() const {
    std::vector<Type> x(this->size());
    for (size_t i=0; i<this->size(); i++)
      x[i] = (*this)[i];
    return x;
  }
  operator tmbutils::vector<Type>() const {
    tmbutils::vector<Type> x(this->size());
    for (size_t i=0; i<this->size(); i++)
      x[i] = (*this)[i];
    return x;
  }
  tmbutils::vector<Type> array() const {
    return *this;
  }
  Type sum() const {
    return TMBad::sum(*this);
  }
};
template <>
struct vector<double> : tmbutils::vector<double> {
  typedef tmbutils::vector<double> Base;
  using Base::Base;
  vector() {} // Surprisingly, needed on some systems.
  vector(tmbutils::vector<double> x) : Base(x) { }
};

/** \brief Demonstrate native vectorization of 'dnorm' */
template<class Type, class T1, class T2>
vector<Type> dnorm(vector<Type> x, T1 mean, T2 sd, int give_log=0)
{
  vector<Type> logres;
  x = (x - mean) / sd;
  logres = -log(Type(sqrt(2*M_PI)) * sd) - Type(.5) * x * x;
  if(give_log) return logres;
  else return exp(logres);
}

} // End namespace

#endif
