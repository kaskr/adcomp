/** \file
    \brief Atomic interpolation operator.
    \note Available for TMBad framework only.
*/

#ifdef TMBAD_FRAMEWORK

#include <memory>

namespace tmbutils {

/** \brief Configuration for `Interpol2D` */
template<class Type>
struct interpol2D_config {
  interpol2D_config(Type R=2) :
    R(asDouble(R)), safe_check(true) { }
  /** \brief Interpolation radius relative to integer lattice
      \details
      - `R` measures degree of smoothness. Default value: `R=2`.
      - Interpolation is only matching exact measurements if `R<=1`.
      - NaN returned outside boundaries
      - NA data values are allowed and will be given zero weight.
  */
  double R;
  /** \brief Enable safety check that data really are constant? */
  bool safe_check;
};

template<class Type>
struct interpol2Dtab {
  matrix<double> data;
  double xmin, xmax, ymin, ymax;
  interpol2D_config<Type> cfg;
  /** \brief Base kernel
      \details Properties
      - f(0) = 1, f(.5) = .5, f(1) = 0
      - Derivatives to 1st order are zero.
  */
  template<class T>
  T f(T x) { return .5 * ( 1 + cos(x * M_PI) ) ; }
  /** \brief Smooth kernel by iterating.
      \details Properties
      - f(0) = 1, f(.5) = .5, f(1) = 0
      - Derivatives to 3rd order are zero.
  */
  template<class T>
  T kernel(T x) { return f(1 - f(x)) ; }
  // helper
  template<class T>
  T sq(T x) { return x * x; }
  /** \brief Do interpolation at (x_,y_)
      \details Works with `double` or `tiny_ad`.
  */
  template<class T>
  T eval(T x_, T y_) {
    double R = cfg.R;
    int nrow = data.rows();
    int ncol = data.cols();
    T hx = (xmax - xmin) / (nrow - 1);
    T hy = (ymax - ymin) / (ncol - 1);
    // Transform (x_, y_) to 'lattice coordinates':
    T i_ = (x_ - xmin) / hx;
    T j_ = (y_ - ymin) / hy;
    // Check
    bool ok = (0 <= i_) && (i_ <= nrow-1) && (0 <= j_) && (j_ <= ncol-1);
    if (!ok) return R_NaN;
    // Envelope window: W := (i_, j_) + [-R, R]^2
    // Get sub-lattice within envelope: (i_min:i_max) x (j_min:j_max)
    int i_min = std::max(asDouble(i_) - R, (double) 0);
    int j_min = std::max(asDouble(j_) - R, (double) 0);
    int i_max = std::min(asDouble(i_) + R, (double) (nrow-1) );
    int j_max = std::min(asDouble(j_) + R, (double) (ncol-1) );
    // Loop through sub-lattice
    T FWsum = 0, Wsum = 0;
    for (int i = i_min; i <= i_max; i++) {
      for (int j = j_min; j <= j_max; j++) {
        T dist2 = sq( (T) i - i_ ) + sq( (T) j - j_ ) ;
        double tiny = 1e-100;
        T dist = sqrt( dist2 + tiny );
        if (dist <= R) {
          double F = data(i, j);
          if (! ISNA(F) ) {
            T W = kernel(dist / R);          
            FWsum += F * W;
            Wsum += W;
          }
        }
      }
    }
    return FWsum / Wsum;
  }
  /** \brief Helper to get derivatives */
  template<int order>
  double D_eval(double x_, double y_, int ny) {
    typedef atomic::tiny_ad::variable<order, 2> T;
    int i = (1 << ny) - 1;
    T x(x_, 0);
    T y(y_, 1);
    return eval(x, y).getDeriv()[i];
  }
  /** \brief Interpolation and derivatives up to 3rd order */
  double operator()(double x_, double y_, int nx=0, int ny=0) {
    int order = nx + ny;
    if (order == 0) {
      return eval(x_, y_);
    } else if (order == 1) {
      return D_eval<1>(x_, y_, ny);
    } else if (order == 2) {
      return D_eval<2>(x_, y_, ny);
    } else if (order == 3) {
      return D_eval<3>(x_, y_, ny);
    } else {
      Rf_error("Order not implemented");
    }
    return 0;
  }
};

/** \brief Get a smooth representation of a **data** matrix.
    \details This 2D interpolation can be put on the AD tape.
    Interpolation data is shared among instances so you can
    do interpolation many times without duplicating the data.

    In addition the operator is thread safe, whether used with manual
    parallelization (`parallel_accumulator` etc) or automatic
    parallization (`TMBad::autopar`).

    The actual interpolation is based on a kernel smoothing approach:
    To obtain the interpolated value at coordinate (x,y), the
    coordinate is first transformed to continous matrix index
    coordinates using the supplied `x_range` and `y_range`.  The
    distance to each matrix entry `d(i,j)` is then calculated, and a
    corresponding 'weight' is given as a function (kernel) of the
    distance divided by the 'interpolation radius' (R),
    i.e. `W(i,j)=kernel(d(i,j)/R)`.  The final interpolated value
    becomes `sum(W*F)/sum(W)` where `F` is the data matrix.

    We note that the interpolation radius (interpol2D_config::R) is
    relative to matrix index coordinates (and is thus not related to
    the `x_range` and `y_range`).

    \note An error is triggered if some data are not constant
*/
template<class Type>
struct interpol2D : TMBad::global::DynamicOperator<2, 1> {
  std::shared_ptr<interpol2Dtab<Type> > dtab;
  int xorder, yorder;
  matrix<double> asDoubleCheck(matrix<Type> x, bool do_check=true) {
    matrix<double> y(x.rows(), x.cols());
    for (int i=0; i<x.rows(); i++) {
      for (int j=0; j<x.cols(); j++) {
        if (do_check && CppAD::Variable(x(i,j)))
          Rf_error("Matrix values must be constants");
        y(i,j) = asDouble(x(i,j));
      }
    }
    return y;
  }
  /** \brief Construct interpolation object
      \param data Matrix of data values.
      \param x_range Vector of length 2 corresponding to first and last data row.
      \param y_range Vector of length 2 corresponding to first and last data column.
      \param cfg Control parameters
  */
  interpol2D(matrix<Type> data,
             vector<Type> x_range,
             vector<Type> y_range,
             interpol2D_config<Type> cfg=interpol2D_config<Type>() ) :
    dtab(std::make_shared<interpol2Dtab<Type> >(interpol2Dtab<Type>({
            asDoubleCheck(data, cfg.safe_check),
              asDouble(x_range[0]),
              asDouble(x_range[1]),
              asDouble(y_range[0]),
              asDouble(y_range[1]),
              cfg
              }))),
    xorder(0),
    yorder(0)
  { }
  TMBad::Index input_size() const {
    return 2;
  }
  TMBad::Index output_size() const {
    return 1;
  }
  typedef TMBad::ad_aug ad;
  /** \brief Scalar evaluation with double types
      \param x x-coordinate
      \param y y-coordinate
      \param nx x-derivative order
      \param ny y-derivative order
  */
  double operator()(double x, double y, int nx=0, int ny=0) {
    return (*dtab)(x, y, nx, ny);
  }
  /** \brief Scalar evaluation with ad types
      \param x x-coordinate
      \param y y-coordinate
      \param nx x-derivative order
      \param ny y-derivative order
  */
  ad operator()(ad x, ad y, int nx=0, int ny=0) {
    std::vector<ad> xy(2); xy[0]=x; xy[1]=y;
    interpol2D cpy(*this);
    cpy.xorder = nx; cpy.yorder=ny;
    std::vector<ad> z = TMBad::global::Complete<interpol2D>(cpy)(xy);
    return z[0];
  }
#define Type T
  VECTORIZE2_tt(operator())
#undef Type
  // Forward pass
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) {
    args.y(0) = (*this)(args.x(0), args.x(1), xorder, yorder);
  }
  // Reverse pass
  template<class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    T dy = args.dy(0);
    args.dx(0) += (*this)(args.x(0), args.x(1), xorder + 1, yorder) * dy;
    args.dx(1) += (*this)(args.x(0), args.x(1), xorder, yorder + 1) * dy;
  }
  // sources code writers are not implemented for this class
  void forward(TMBad::ForwardArgs<TMBad::Writer> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "IP2D"; }
};

}

#endif
