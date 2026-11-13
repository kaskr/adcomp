#include "analytic_suite.hpp"
#include <memory>
struct Snapshot {
  TMBad::Sparse<Tape> h;
  Tape nested;
  Model model;
  std::vector<double> x;
};
void retape(Tape &f, Model m, const std::vector<double> &point) {
  f.glob.clear();
  std::vector<TapeAD> x(point.begin(), point.end());
  auto before = TMBad::get_glob();
  f.glob.ad_start();
  TMBad::Independent(x);
  auto y = m(x);
  TMBad::Dependent(y);
  f.glob.ad_stop();
  if (TMBad::get_glob() != before)
    throw std::runtime_error("retape did not restore recording context");
}
void verify(Snapshot &s, double &error) {
  size_t n = s.x.size();
  std::vector<double> q(s.h.Range()), hv(q.size(), 0), gr(n, 0);
  for (size_t j = 0; j < q.size(); ++j) q[j] = std::sin(.4 + j);
  for (auto t : s.model.terms) {
    double z = 0, aq = 0;
    for (size_t j = 0; j < n; ++j) z += t.a[j] * s.x[j];
    auto d = derivatives(t.kind, z);
    for (size_t j = 0; j < q.size(); ++j) {
      double a = t.a[s.h.i[j]] * t.a[s.h.j[j]];
      hv[j] += d[0] * a;
      aq += q[j] * a;
    }
    for (size_t j = 0; j < n; ++j) gr[j] += d[1] * aq * t.a[j];
  }
  check(s.h(s.x), hv, error);
  check(s.h.Jacobian(s.x, q), gr, error);
  auto in = s.x;
  in.insert(in.end(), q.begin(), q.end());
  check(s.nested(in), gr, error);
  // Retain only a copy of the nested derivative, then destroy all source tapes.
  Tape survivor = s.nested;
  Snapshot dead = std::move(s);
  s = Snapshot{};
  dead = Snapshot{};
  check(survivor(in), gr, error);
}
int main() {
  try {
    double error = analytic_suite();
    size_t verified = 0;
    auto root = TMBad::get_glob();
    for (int width : {0, 2, 4, 8}) {
      std::vector<Snapshot> retained;
      {
        Tape source, gradient;
        for (int cycle = 0; cycle < 80; ++cycle) {
          size_t n = 5 + cycle % 9;
          Model m;
          for (size_t start = 0; start < n - 1; start += 3) {
            Term t{(cycle % 3 == 0 && start == 3) ? 2 : 0,
                   std::vector<double>(n, 0)};
            for (size_t j = start; j < std::min(n - 1, start + 3); ++j)
              t.a[j] = .2 + .03 * j;
            m.terms.push_back(t);
          }
          std::vector<double> x(n, .1 + .001 * cycle);
          retape(source, m, x);
          gradient = source.JacFun();
          TMBad::SpJacFun_config cfg;
          cfg.batch_width = width;
          cfg.compress = true;
          cfg.index_remap = false;
          Snapshot s;
          s.h = gradient.SpJacFun({}, {}, cfg);
          s.nested = s.h.WgtJacFun();
          s.model = m;
          s.x = x;
          // Copies share atomic derivative tables; exercise independent
          // evaluations and destruction order.
          Snapshot copy = s;
          for (size_t j = 0; j < copy.x.size(); ++j)
            copy.x[j] += .007 * (j + 1);
          retained.push_back(std::move(copy));
          retained.push_back(s);
          if (cycle % 4 == 3) {
            verify(retained.front(), error);
            ++verified;
            retained.erase(retained.begin());
          }
          if (TMBad::get_glob() != root)
            throw std::runtime_error("recording context leaked");
        }
      }  // Destroy repeatedly retaped source and gradient before evaluating
         // retained Hessians.
      while (!retained.empty()) {
        auto s = std::move(retained.back());
        retained.pop_back();
        verify(s, error);
        ++verified;
      }
      if (TMBad::get_glob() != root)
        throw std::runtime_error("destruction changed recording context");
    }
    if (verified != 640) throw std::runtime_error("wrong lifecycle case count");
    std::cout << "PASS: 672 analytic cases; 320 in-place source retapes; "
              << verified << " retained/copy/destruction checks; max error "
              << std::setprecision(17) << error << "\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FAIL: " << e.what() << "\n";
    return 1;
  }
}
