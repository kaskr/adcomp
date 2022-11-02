/** \file
    \brief Object concatenation.
*/

namespace tmbutils {

template<class T>
struct concat_traits {
  typedef T value_type;
  T* ptr(T &x) { return &x; }
  size_t size(T &x) { return 1; }
};
template<class VT>
struct concat_traits_vector {
  typedef typename VT::value_type value_type;
  typedef value_type T;
  T* ptr(VT &x) { return x.data(); }
  size_t size(VT &x) { return x.size(); }
};
template<class T>
struct concat_traits<vector<T> > : concat_traits_vector<vector<T> > {};
template<class T>
struct concat_traits<matrix<T> > : concat_traits_vector<matrix<T> > {};
template<class T>
struct concat_traits<array<T> > : concat_traits_vector<array<T> > {};

/** \brief Serialized representation of objects of different types */
template<class S=void, class...Ts>
struct Concat {
  S& cur;
  Concat<Ts...> next;
  typedef typename concat_traits<S>::value_type value_type;
  void assign_from(const value_type* x) {
    value_type* dest = concat_traits<S>().ptr(cur);
    size_t n = concat_traits<S>().size(cur);
    for (size_t i=0; i<n; i++) {dest[i] = *x; x++;}
    next.assign_from(x);
  }
  void assign_to(value_type* x) {
    value_type* orig = concat_traits<S>().ptr(cur);
    size_t n = concat_traits<S>().size(cur);
    for (size_t i=0; i<n; i++) {*x = orig[i]; x++;}
    next.assign_to(x);
  }
  size_t size() {return concat_traits<S>().size(cur) + next.size(); }
  operator vector<value_type>() {
    vector<value_type> ans(size());
    assign_to(ans.data());
    return ans;
  }
  Concat& operator=(const vector<value_type> &x) {
    assign_from(x.data());
    return *this;
  }
};
template<>
struct Concat<void> {
  template<class T>
  void assign_from(T x) { }
  template<class T>
  void assign_to(T x) { }
  size_t size() { return 0; }
};
/** \brief Serialized representation of objects of different types

    \details This function provides an interface between a parameter
    list and its serialized (vector) representation.
    In R, this would correspond to 'unlist()/relist()'.

    Example:
    \code
    // Custom parameter list
    template<class Type>
    struct parameter_list {
      vector<Type> v;
      matrix<Type> m;
      array<Type> a;
      Type s;
      // To enable unlist/relist functionality, we must specify which members to include in representation:
      auto unlist() { return concat(v, m, a, s); }
    };

    // Declare a parameter list
    parameter_list<Type> pl;

    // To 'unlist' it do
    vector<Type> parms = pl.unlist();

    // To 'relist' it do
    pl.unlist() = parms;
    \endcode
    \note Return type 'auto' requires C++14.
*/
template<class...Ts>
Concat<Ts...> concat(Ts&... args) {
  return Concat<Ts...>({args...});
}

} // End namespace
