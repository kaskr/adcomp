#ifndef HAVE_VECTORIZE_HPP
#define HAVE_VECTORIZE_HPP
// Autogenerated - do not edit by hand !

namespace TMBad {

typedef global::ad_segment ad_segment;

template <class Type, bool S0 = 0, bool S1 = 0>
struct Vectorized {
  Type x;

  static constexpr bool stride(bool j) { return j == 0 ? S0 : S1; }
  operator Type() { return x; }
  Vectorized(Type x) : x(x) {}
  Vectorized() {}
};

template <class Type, bool S0, bool S1>
struct ForwardArgs<Vectorized<Type, S0, S1> > : ForwardArgs<Type> {
  typedef Vectorized<Type, S0, S1> T;
  typedef ForwardArgs<Type> Base;
  size_t k;
  /** \brief j'th input variable of this operator */
  Type x(bool j) const {
    return Base::values[Base::input(j) + k * T::stride(j)];
  }
  /** \brief j'th output variable of this operator */
  Type &y(Index j) { return Base::values[Base::output(j) + k]; }
  ForwardArgs(const Base &x) : Base(x) {}
};

template <class Type, bool S0, bool S1>
struct ReverseArgs<Vectorized<Type, S0, S1> > : ReverseArgs<Type> {
  typedef Vectorized<Type, S0, S1> T;
  typedef ReverseArgs<Type> Base;
  size_t k;
  /** \brief j'th input variable of this operator */
  Type x(bool j) const {
    return Base::values[Base::input(j) + k * T::stride(j)];
  }
  /** \brief j'th output variable of this operator */
  Type y(Index j) const { return Base::values[Base::output(j) + k]; }
  /** \brief Partial derivative of end result wrt. j'th input variable of
      this operator */
  Type &dx(bool j) const {
    return Base::derivs[Base::input(j) + k * T::stride(j)];
  }
  /** \brief Partial derivative of end result wrt. j'th output variable of
      this operator */
  Type dy(Index j) const { return Base::derivs[Base::output(j) + k]; }
  ReverseArgs(const Base &x) : Base(x) {}
};

struct VSumOp : global::DynamicOperator<1, 1> {
  static const bool is_linear = true;
  size_t n;
  VSumOp(size_t n);
  template <class Type>
  void forward(ForwardArgs<Type> &args) {
    const Type *x = args.x_ptr(0);
    Type &y = args.y(0);
    y = 0;
    for (size_t i = 0; i < n; i++) y += x[i];
  }
  template <class Type>
  void reverse(ReverseArgs<Type> &args) {
    Type *dx = args.dx_ptr(0);
    const Type &dy = args.dy(0);
    for (size_t i = 0; i < n; i++) dx[i] += dy;
  }

  void dependencies(Args<> &args, Dependencies &dep) const;
  static const bool have_dependencies = true;
  /** \brief This operator **has** implicit dependencies */
  static const bool implicit_dependencies = true;
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;
  void forward(ForwardArgs<Writer> &args);
  void reverse(ReverseArgs<Writer> &args);
  const char *op_name();
};

ad_aug sum(ad_segment x);

template <class dummy = void>
ad_segment operator/(ad_segment x, ad_segment y);
template <class dummy = void>
ad_segment operator*(ad_segment x, ad_segment y);
template <class dummy = void>
ad_segment operator+(ad_segment x, ad_segment y);
template <class dummy = void>
ad_segment operator-(ad_segment x, ad_segment y);
template <class dummy = void>
ad_segment operator-(ad_segment x);
template <class dummy = void>
ad_segment &operator+=(ad_segment &x, ad_segment y) {
  if (x.size() < y.size()) y = ad_segment(sum(y), 1);
  if (x.identicalZero())
    x = y;
  else
    x = x + y;
  return x;
}
template <class dummy = void>
ad_segment &operator-=(ad_segment &x, ad_segment y) {
  if (x.size() < y.size()) y = ad_segment(sum(y), 1);
  if (x.identicalZero())
    x = -y;
  else
    x = x - y;
  return x;
}

template <class Operator, bool S0 = false, bool S1 = false>
struct Vectorize : global::DynamicOperator<Operator::ninput, -1> {
  size_t n;
  static const bool have_input_size_output_size = true;
  Index input_size() const { return Operator::ninput; }
  Index output_size() const { return this->n; }
  Vectorize(size_t n) : n(n) {}
  void forward(ForwardArgs<Scalar> &args) {
    ForwardArgs<Vectorized<Scalar, S0, S1> > vargs(args);
    typename global::CPL<Operator>::type Op;
    for (vargs.k = 0; vargs.k < n; vargs.k++) {
      Op.forward(vargs);
    }
  }
  void forward(ForwardArgs<Replay> &args) {
    ad_segment x0(args.x_ptr(0), (S0 ? n : 1));
    ad_segment x1;
    if (Operator::ninput > 1) {
      x1 = ad_segment(args.x_ptr(1), (S1 ? n : 1));
    }
    global::Complete<Vectorize> F(*this);
    ad_segment y = F(x0, x1);
    for (size_t i = 0; i < y.size(); i++) args.y(i) = y[i];
  }
  void reverse(ReverseArgs<Scalar> &args) {
    ReverseArgs<Vectorized<Scalar, S0, S1> > vargs(args);
    typename global::CPL<Operator>::type Op;
    for (vargs.k = 0; vargs.k < n; vargs.k++) {
      Op.reverse(vargs);
    }
  }
  void reverse(ReverseArgs<Replay> &args) {
    std::vector<ad_segment> v;
    std::vector<ad_segment> d;
    std::vector<Index> i;

    v.push_back(ad_segment(args.x_ptr(0), (S0 ? n : 1)));
    d.push_back(ad_segment(args.dx_ptr(0), (S0 ? n : 1), true));
    i.push_back(i.size());
    if (Operator::ninput > 1) {
      v.push_back(ad_segment(args.x_ptr(1), (S1 ? n : 1)));
      d.push_back(ad_segment(args.dx_ptr(1), (S1 ? n : 1), true));
      i.push_back(i.size());
    }

    v.push_back(ad_segment(args.y_ptr(0), n));
    d.push_back(ad_segment(args.dy_ptr(0), n));

    ReverseArgs<ad_segment> vargs(i, v, d);

    vargs.ptr.first = 0;
    vargs.ptr.second = Operator::ninput;
    typename global::CPL<Operator>::type Op;
    Op.reverse(vargs);

    for (size_t i = 0; i < vargs.dx(0).size(); i++)
      args.dx_ptr(0)[i] = vargs.dx(0)[i];
    if (Operator::ninput > 1) {
      for (size_t i = 0; i < vargs.dx(1).size(); i++)
        args.dx_ptr(1)[i] = vargs.dx(1)[i];
    }
  }

  void dependencies(Args<> &args, Dependencies &dep) const {
    dep.add_segment(args.input(0), (S0 ? n : 1));
    if (Operator::ninput == 2) {
      dep.add_segment(args.input(1), (S1 ? n : 1));
    }
  }
  static const bool have_dependencies = true;
  /** \brief This operator **has** implicit dependencies */
  static const bool implicit_dependencies = true;
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;
  void forward(ForwardArgs<Writer> &args) { TMBAD_ASSERT(false); }
  void reverse(ReverseArgs<Writer> &args) { TMBAD_ASSERT(false); }
  const char *op_name() {
    global::Complete<Operator> Op;
    static const std::string name = std::string("V") + Op.op_name();
    return name.c_str();
  }
  Vectorize(const ad_segment &x, const ad_segment &y)
      : n(std::max(x.size(), y.size())) {}
};
template <class dummy>
ad_segment operator/(ad_segment x, ad_segment y) {
  size_t n = std::max(x.size(), y.size());
  if (x.size() > 1 && y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::DivOp, 1, 1> > F(n);
    return F(x, y);
  } else if (x.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::DivOp, 1, 0> > F(n);
    return F(x, y);
  } else if (y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::DivOp, 0, 1> > F(n);
    return F(x, y);
  } else {
    global::Complete<Vectorize<global::ad_plain::DivOp, 0, 0> > F(n);
    return F(x, y);
  }
  TMBAD_ASSERT(false);
  return ad_segment();
}
template <class dummy>
ad_segment operator*(ad_segment x, ad_segment y) {
  size_t n = std::max(x.size(), y.size());
  if (x.size() > 1 && y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::MulOp, 1, 1> > F(n);
    return F(x, y);
  } else if (x.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::MulOp, 1, 0> > F(n);
    return F(x, y);
  } else if (y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::MulOp, 0, 1> > F(n);
    return F(x, y);
  } else {
    global::Complete<Vectorize<global::ad_plain::MulOp, 0, 0> > F(n);
    return F(x, y);
  }
  TMBAD_ASSERT(false);
  return ad_segment();
}
template <class dummy>
ad_segment operator+(ad_segment x, ad_segment y) {
  size_t n = std::max(x.size(), y.size());
  if (x.size() > 1 && y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::AddOp, 1, 1> > F(n);
    return F(x, y);
  } else if (x.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::AddOp, 1, 0> > F(n);
    return F(x, y);
  } else if (y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::AddOp, 0, 1> > F(n);
    return F(x, y);
  } else {
    global::Complete<Vectorize<global::ad_plain::AddOp, 0, 0> > F(n);
    return F(x, y);
  }
  TMBAD_ASSERT(false);
  return ad_segment();
}
template <class dummy>
ad_segment operator-(ad_segment x, ad_segment y) {
  size_t n = std::max(x.size(), y.size());
  if (x.size() > 1 && y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::SubOp, 1, 1> > F(n);
    return F(x, y);
  } else if (x.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::SubOp, 1, 0> > F(n);
    return F(x, y);
  } else if (y.size() > 1) {
    global::Complete<Vectorize<global::ad_plain::SubOp, 0, 1> > F(n);
    return F(x, y);
  } else {
    global::Complete<Vectorize<global::ad_plain::SubOp, 0, 0> > F(n);
    return F(x, y);
  }
  TMBAD_ASSERT(false);
  return ad_segment();
}
template <class dummy = void>
ad_segment pow(ad_segment x, ad_segment y);
template <class dummy>
ad_segment pow(ad_segment x, ad_segment y) {
  size_t n = std::max(x.size(), y.size());
  if (x.size() > 1 && y.size() > 1) {
    global::Complete<Vectorize<PowOp, 1, 1> > F(n);
    return F(x, y);
  } else if (x.size() > 1) {
    global::Complete<Vectorize<PowOp, 1, 0> > F(n);
    return F(x, y);
  } else if (y.size() > 1) {
    global::Complete<Vectorize<PowOp, 0, 1> > F(n);
    return F(x, y);
  } else {
    global::Complete<Vectorize<PowOp, 0, 0> > F(n);
    return F(x, y);
  }
  TMBAD_ASSERT(false);
  return ad_segment();
}
template <class dummy>
ad_segment operator-(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<global::ad_plain::NegOp, 1, 0> > F(n);
  return F(x);
}

template <class dummy = void>
ad_segment fabs(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<AbsOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment sin(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<SinOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment cos(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<CosOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment exp(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<ExpOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment log(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<LogOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment sqrt(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<SqrtOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment tan(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<TanOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment sinh(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<SinhOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment cosh(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<CoshOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment tanh(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<TanhOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment expm1(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<Expm1, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment log1p(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<Log1p, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment asin(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<AsinOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment acos(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<AcosOp, 1, 0> > F(n);
  return F(x);
}
template <class dummy = void>
ad_segment atan(ad_segment x) {
  size_t n = x.size();
  global::Complete<Vectorize<AtanOp, 1, 0> > F(n);
  return F(x);
}
template <class T>
struct ScalarPack {
  static const int size = (sizeof(T) - 1) / sizeof(Scalar) + 1;
};

/** \brief Representation of a *specific* contiguous set of values on a
 * *specific* tape */
struct SegmentRef {
  global *glob_ptr;
  Index offset;
  Index size;
  Scalar *value_ptr();
  Scalar *deriv_ptr();
  SegmentRef();
  SegmentRef(const Scalar *x);
  SegmentRef(global *g, Index o, Index s);
  SegmentRef(const ad_segment &x);
  bool isNull();
  void resize(ad_segment &pack, Index n);
};

ad_segment pack(const ad_segment &x, bool up = false);
ad_segment unpack(const ad_segment &x);
void unpack(const ad_segment &x, ad_segment &y);

/** \brief Pack (`PackOp`) or unpack (`UnpkOp`) `n` consecutive values
    on the tape.

    The purpose of these operators is to provide an analogue of 'pass by
    reference' between AD tapes.

    An 'outer' tape which repeatedly (`N` times) calls an 'inner' tape
    (`R^n->R`) would normally require `O(N*n)` memory units because
    each of the `n` inputs must be referenced each time the inner tape
    is called. By packing the input array, the memory usage can be
    reduced to `O(N+n)`.
*/
template <bool UP>
struct PackOp : global::DynamicOperator<1, ScalarPack<SegmentRef>::size> {
  static const bool forward_updating = false;
  static const bool reverse_updating = UP;
  /** \brief Packed size (~2) */
  static const Index K = ScalarPack<SegmentRef>::size;
  /** \brief Unpacked size */
  Index n;
  PackOp(const Index n) : n(n) {}
  /** \brief Pack values */
  void forward(ForwardArgs<Scalar> &args) {
    SegmentRef *y = (SegmentRef *)args.y_ptr(0);
    y[0] = SegmentRef(args.glob_ptr, args.input(0), n);
  }
  /** \brief Pack values (replay) */
  void forward(ForwardArgs<Replay> &args) {
    ad_segment x(args.x_ptr(0), n);
    args.y_segment(0, K) = pack(x, UP);
  }
  /** \brief Unpack derivatives */
  void reverse(ReverseArgs<Scalar> &args) {
    SegmentRef tmp(args.dy_ptr(0));
    if (tmp.glob_ptr != NULL) {
      Scalar *dx = SegmentRef(args.y_ptr(0)).deriv_ptr();
      Scalar *dy = SegmentRef(args.dy_ptr(0)).deriv_ptr();
      for (Index i = 0; i < n; i++) dx[i] += dy[i];
    }
  }
  /** \brief Unpack derivatives (replay) */
  void reverse(ReverseArgs<Replay> &args) {
    ad_segment dy_packed(args.dy_ptr(0), K);

    if (SegmentRef(dy_packed).isNull()) {
      SegmentRef().resize(dy_packed, n);
    }
    if (UP) {
      TMBAD_ASSERT2(
          TMBad::global::ad_segment().is_contiguous(args.dx_ptr(0), n),
          "Internal inconsistency: "
          "Updatable derivative workspace *must* be 'on tape' and 'contiguous'")

          ;
      ad_segment dx(args.dx_ptr(0), n);
      unpack(dy_packed, dx);
    } else {
      ad_segment dy = unpack(dy_packed);
      ad_segment dx(args.dx_ptr(0), n, true);
      dx += dy;
      Replay *pdx = args.dx_ptr(0);
      for (Index i = 0; i < n; i++) pdx[i] = dx[i];
    }
  }
  const char *op_name() { return "PackOp"; }
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;
  static const bool have_dependencies = true;
  static const bool implicit_dependencies = true;
  void dependencies(Args<> &args, Dependencies &dep) const {
    dep.add_segment(args.input(0), n);
  }

  template <class T>
  void forward(ForwardArgs<T> &args) {
    TMBAD_ASSERT2(false, "PackOp: Invalid method!");
  }
  template <class T>
  void reverse(ReverseArgs<T> &args) {
    TMBAD_ASSERT2(false, "PackOp: Invalid method!");
  }
};

/** \copydoc PackOp */
template <bool UP>
struct UnpkOp : global::DynamicOperator<1 + UP, -1> {
  static const bool forward_updating = UP;
  static const bool reverse_updating = false;
  /** \brief Packed size (~2) */
  static const Index K = ScalarPack<SegmentRef>::size;
  /** \brief Unpacked size */
  Index noutput;
  static const bool have_input_size_output_size = true;
  Index input_size() const { return 1 + UP; }
  Index output_size() const { return UP ? 0 : noutput; }
  UnpkOp(const Index n) : noutput(n) {}
  /** \brief Unpack values */
  void forward(ForwardArgs<Scalar> &args) {
    Scalar *y = (UP ? args.x_ptr(1) : args.y_ptr(0));
    SegmentRef srx(args.x_ptr(0));
    if (srx.isNull()) {
      if (!UP)
        for (Index i = 0; i < noutput; i++) y[i] = 0;
      return;
    }
    Scalar *x = srx.value_ptr();
    if (UP)
      for (Index i = 0; i < noutput; i++) y[i] += x[i];
    else
      for (Index i = 0; i < noutput; i++) y[i] = x[i];

    ((SegmentRef *)args.x_ptr(0))->glob_ptr = NULL;
  }
  static const bool add_forward_replay_copy = true;
  /** \brief Pack derivatives */
  void reverse(ReverseArgs<Scalar> &args) {
    SegmentRef *dx = (SegmentRef *)args.dx_ptr(0);
    dx[0] =
        SegmentRef(args.glob_ptr, UP ? args.input(1) : args.output(0), noutput);
  }
  /** \brief Pack derivatives (replay) */
  void reverse(ReverseArgs<Replay> &args) {
    ad_segment dy(UP ? args.dx_ptr(1) : args.dy_ptr(0), noutput);
    ad_segment dy_packed = pack(dy, UP);
    Replay *pdx = args.dx_ptr(0);
    for (Index i = 0; i < dy_packed.size(); i++) pdx[i] = dy_packed[i];
  }
  const char *op_name() { return "UnpkOp"; }

  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;
  static const bool have_dependencies = true;
  static const bool implicit_dependencies = true;
  void dependencies(Args<> &args, Dependencies &dep) const {
    dep.add_segment(args.input(0), K);
  }
  void dependencies_updating(Args<> &args, Dependencies &dep) const {
    if (UP) dep.add_segment(args.input(1), noutput);
  }

  template <class T>
  void forward(ForwardArgs<T> &args) {
    TMBAD_ASSERT2(false, "UnpkOp: Invalid method!");
  }
  template <class T>
  void reverse(ReverseArgs<T> &args) {
    TMBAD_ASSERT2(false, "UnpkOp: Invalid method!");
  }
};

/** \brief Pack consecutive values on the tape */
ad_segment pack(const ad_segment &x, bool up);

/** \brief Unpack consecutive values on the tape */
ad_segment unpack(const ad_segment &x);
/** \brief Fused unpack and destination increment */
void unpack(const ad_segment &x, ad_segment &y);

/** \brief Unpack consecutive values on the tape */
template <class T>
ad_segment unpack(const std::vector<T> &x, Index j) {
  Index K = ScalarPack<SegmentRef>::size;
  ad_segment x_(x[j * K], K);
  return unpack(x_);
}
Scalar *unpack(const std::vector<Scalar> &x, Index j);

template <class T>
std::vector<T> repack(const std::vector<T> &x) {
  Index K = ScalarPack<SegmentRef>::size;
  size_t n = x.size() / K;
  std::vector<T> y;
  for (size_t j = 0; j < n; j++) {
    ad_segment x_(x[j * K], K);
    SegmentRef sr(x_);
    ad_segment orig(sr.offset, sr.size);
    ad_segment yj = pack(orig, true);
    for (size_t i = 0; i < K; i++) y.push_back(yj[i]);
  }
  return y;
}

std::vector<ad_aug> concat(const std::vector<ad_segment> &x);

}  // namespace TMBad
#endif  // HAVE_VECTORIZE_HPP
