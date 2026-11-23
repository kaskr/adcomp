#include <iostream>
#include <iomanip>
#include <memory>
#include <map>
#include <cassert>
#define TMBAD_ASSERT(x) assert(x)
#define TMBAD_ASSERT2(x, msg) assert(x)
#define Rcout std::cout
#define Rcerr std::cerr
#include "TMBad/TMBad.hpp"
#include "TMBad/eigen_numtraits.hpp"
#include "TMBad/TMBad.cpp"
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
using TapeAD = TMBad::ad_aug;
using Tape = TMBad::ADFun<TapeAD>;
// Deliberately collides with a supported built-in name, but computes a cube.
struct FakeExp : TMBad::global::UnaryOperator {
  static const bool add_forward_replay_copy = true;
  static const bool have_eval = true;
  template <class T>
  T eval(T x) {
    return x * x * x;
  }
  template <class T>
  void reverse(TMBad::ReverseArgs<T> &a) {
    a.dx(0) += a.dy(0) * T(3) * a.y(0) / a.x(0);
  }
  const char *op_name() { return "ExpOp"; }
};
TapeAD fake_exp(TapeAD x) {
  return TMBad::get_glob()->add_to_stack<FakeExp>(TMBad::ad_plain(x));
}
struct SuppressedActive {
  std::vector<TapeAD> operator()(std::vector<TapeAD> x) {
    return {TMBad::get_glob()
                ->add_to_stack<TMBad::global::ad_plain::MulOp_<true, false> >(
                    TMBad::ad_plain(x[0]), TMBad::ad_plain(x[1]))};
  }
};
struct InnerExp {
  std::vector<TapeAD> operator()(std::vector<TapeAD> x) { return {exp(x[0])}; }
};
struct Term {
  int kind;
  std::vector<double> a;
};
struct Model {
  std::vector<Term> terms;
  std::vector<TapeAD> operator()(std::vector<TapeAD> x) {
    TapeAD y =
        x.back();  // inactive Hessian input, despite active objective input
    for (auto t : terms) {
      TapeAD z = 0;
      for (size_t i = 0; i < x.size(); ++i)
        if (t.a[i]) z += t.a[i] * x[i];
      if (t.kind == 0)
        y += exp(z);
      else if (t.kind == 1)
        y += sin(z);
      else if (t.kind == 2)
        y += tanh(z);
      else if (t.kind == 3)
        y += fake_exp(z);
      else {
        Tape inner(InnerExp{}, std::vector<double>{.1});
        y += inner.atomic_raw()(std::vector<TapeAD>{z})[0];
      }
    }
    return {y};
  }
};
std::vector<double> derivatives(int kind, double z) {
  if (kind == 0 || kind == 4) {
    double e = std::exp(z);
    return {e, e, e};
  }
  if (kind == 1) return {-std::sin(z), -std::cos(z), std::sin(z)};
  if (kind == 2) {
    double t = std::tanh(z), u = 1 - t * t;
    return {-2 * t * u, (-2 + 6 * t * t) * u, (16 * t - 24 * t * t * t) * u};
  }
  return {6 * z, 6, 0};
}
void check(const std::vector<double> &a, const std::vector<double> &b,
           double &err) {
  if (a.size() != b.size()) throw std::runtime_error("dimension mismatch");
  for (size_t i = 0; i < a.size(); ++i) {
    if (!std::isfinite(a[i]) || !std::isfinite(b[i]))
      throw std::runtime_error("nonfinite");
    err = std::max(err, std::abs(a[i] - b[i]));
    if (std::abs(a[i] - b[i]) > 1e-9 * (1 + std::abs(b[i])))
      throw std::runtime_error("analytic derivative mismatch");
  }
}
double analytic_suite() {
  try {
    double error = 0;
    int cases = 0;
    Tape suppressed(SuppressedActive{}, std::vector<double>{.2, .3});
    Tape::compress_helper::batch_transform guard(suppressed, 1, 4);
    if (guard.supported || guard.reason != "suppressed active derivative")
      throw std::runtime_error("suppressed active derivative accepted");
    std::fprintf(stderr, "GUARD,suppressed-active,PASS\n");
    Tape f;  // Replace the source repeatedly, including different operator
             // topologies.

    for (int model = 0; model < 7; ++model)
      for (int mask = 0; mask < 3; ++mask)
        for (int width : {0, 2, 4, 8})
          for (int remap : {0, 1}) {
            const size_t n = 12;
            Model m;
            auto term = [&](int kind, int first, int last) {
              Term t{kind, std::vector<double>(n, 0)};
              for (int j = first; j < last; ++j) t.a[j] = .3 + .07 * j;
              m.terms.push_back(t);
            };
            if (model == 0) {
              term(0, 0, 3);
              term(0, 3, 8);
              term(1, 8, 11);
            }  // separate partial groups
            if (model == 1)
              for (int j = 0; j < 11; ++j) term(0, j, j + 1);  // diagonal
            if (model == 2)
              for (int j = 0; j < 10; ++j)
                term(1, j, j + 2);  // overlapping band
            if (model == 3) {
              term(0, 0, 3);
              term(2, 3, 8);
              term(0, 8, 11);
            }  // mixed supported and fallback groups
            if (model == 4) {
              term(3, 0, 3);
              term(0, 3, 8);
              term(3, 8, 11);
            }  // spoofed operator name
            if (model == 6) {
              term(0, 0, 3);
              term(4, 3, 8);
              term(0, 8, 11);
            }  // nested atom fallback
            // model 5: linear objective and empty Hessian
            std::vector<double> x(n);
            for (size_t j = 0; j < n; ++j) x[j] = .02 * (j + 1);
            std::vector<bool> kx(n, true), ky(n, true);
            if (mask == 1)
              for (size_t j = 0; j < n; ++j) {
                kx[j] = (j % 3 != 1);
                ky[j] = (j % 2 == 0);
              }
            if (mask == 2) {
              kx.assign(n, false);
              ky.assign(n, false);
              kx[4] = kx[6] = true;
              ky[3] = ky[6] = true;
            }
            std::fprintf(stderr, "CASE,%d,%d,%d,%d\n", model, mask, width,
                         remap);
            f = Tape(m, x);
            if (model == 4) {
              Tape::compress_helper::batch_transform identity(f, n, 4);
              if (identity.supported || identity.reason != "ExpOp")
                throw std::runtime_error("spoofed built-in name accepted");
            }
            Tape g = f.JacFun();
            TMBad::SpJacFun_config cfg;
            cfg.batch_width = width;
            cfg.compress = false;
            cfg.index_remap = remap;
            auto plain = g.SpJacFun(kx, ky, cfg);
            cfg.compress = true;
            auto h = g.SpJacFun(kx, ky, cfg);
            if (h.i != plain.i || h.j != plain.j)
              throw std::runtime_error("sparse structure mismatch");
            std::vector<size_t> xi, yi;
            for (size_t j = 0; j < n; ++j) {
              if (kx[j]) xi.push_back(j);
              if (ky[j]) yi.push_back(j);
            }
            if (h.n != plain.n || h.m != plain.m)
              throw std::runtime_error("sparse dimensions mismatch");
            // Recorded weighted reverse exercises nested differentiation
            // through atoms.
            Tape nested = h.WgtJacFun();
            for (int point = 0; point < 4; ++point) {
              for (size_t j = 0; j < n; ++j) x[j] = .02 * (j + 1) - .03 * point;
              std::vector<double> q(h.Range()), v(n), eh(h.Range(), 0),
                  eg(n, 0), en(n + h.Range(), 0);
              for (size_t j = 0; j < q.size(); ++j)
                q[j] = point == 3 ? 0 : std::sin(.2 + j + point);
              for (size_t j = 0; j < n; ++j) v[j] = std::cos(.7 + j);
              for (auto t : m.terms) {
                double z = 0, av = 0, aq = 0;
                for (size_t j = 0; j < n; ++j) {
                  z += t.a[j] * x[j];
                  av += t.a[j] * v[j];
                }
                auto d = derivatives(t.kind, z);
                for (size_t j = 0; j < h.Range(); ++j) {
                  double aij = t.a[remap ? yi[h.i[j]] : h.i[j]] *
                               t.a[remap ? xi[h.j[j]] : h.j[j]];
                  eh[j] += d[0] * aij;
                  aq += q[j] * aij;
                  en[n + j] += d[1] * av * aij;
                }
                for (size_t j = 0; j < n; ++j) {
                  eg[j] += d[1] * t.a[j] * aq;
                  en[j] += d[2] * t.a[j] * av * aq;
                }
              }
              check(h(x), eh, error);
              check(h(x), plain(x), error);
              check(h.Jacobian(x, q), eg, error);
              auto input = x;
              input.insert(input.end(), q.begin(), q.end());
              check(nested(input), eg, error);
              check(nested.Jacobian(input, v), en, error);
              ++cases;
            }
          }
    return error;

  } catch (...) {
    throw;
  }
}
