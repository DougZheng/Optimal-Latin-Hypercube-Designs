#include "design.h"

#include <cmath>
#include <cassert>

#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <random>

/*
An efficient algorithm for constructing optimal design of computer experiments
https://sci-hub.yncjkj.com/https://www.sciencedirect.com/science/article/pii/S0378375804001922
*/

namespace SearchingAlgorithm {
class ESE {
 public:
  ESE(int n, int k);
  inline void SetMCol(int m_col) { m_col_ = m_col; }
  inline void SetJColPair(int j_col_pair) { j_col_pair_ = j_col_pair; }
  inline void SetAlpha1(double alpha1) { alpha1_ = alpha1; }
  inline void SetAlpha2(double alpha2) { alpha2_ = alpha2; }
  inline void SetAlpha3(double alpha3) { alpha3_ = alpha3; }
  inline void SetW(double w) { w_ = w; }
  inline void SetIterateCnt(int iterate_cnt) { iterate_cnt_ = iterate_cnt; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetPrintFrequence(int print_frequence) { print_frequence_ = print_frequence; }
  Design Search();
  Design IncrementalSearch(Design design);
 private:
  void InitDefaultParam();
  void ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m);
 private:
  int n_;
  int k_;
  int m_col_;
  int j_col_pair_;
  double alpha1_;
  double alpha2_;
  double alpha3_;
  double w_;
  int iterate_cnt_;
  std::mt19937 rng_;
  int print_frequence_;
};

ESE::ESE(int n, int k) : n_(n), k_(k) {
  InitDefaultParam();
}

Design ESE::Search() {
  Design design(n_, k_, rng_());
  return IncrementalSearch(design);
}

Design ESE::IncrementalSearch(Design design) {
  assert(design.GetN() == n_ && design.GetK() == k_);
  if (w_ == 0) design.DisableCorr();
  if (w_ == 1) design.DisableDis();
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int cnt = 0;
  double opt_corr = design.GetMaxAbsCorr();
  double opt_val = design.GetCritVal(w_, opt_corr, design.GetPhiP());
  Design::VecInt2D opt_design = design.GetDesign();
  double T_h = 0.005 * opt_val;
  std::vector<std::pair<int, int>> pair_list;
  pair_list.reserve(n_ * (n_ - 1) / 2);
  for (int i = 0; i < n_; ++i) {
    for (int j = i + 1; j < n_; ++j) {
      pair_list.emplace_back(i, j);
    }
  }
  double cur_val = opt_val;
  double cur_corr = opt_corr;
  int max_no_imp_cnt = 0;
  int explore_dir = 0;
  auto StopSearch = [this, &cnt, &max_no_imp_cnt]() -> bool {
    if (iterate_cnt_ >= 0) {
      return cnt >= iterate_cnt_;
    }
    return max_no_imp_cnt >= 1000;
  };
  auto PrintLog = [this, &cnt, &cur_val, &cur_corr, &design]() -> void {
    std::cerr << "Iteration: " << cnt
      << ", Val: " << cur_val 
      << ", PhiP: " << design.GetPhiP()
      << ", RhoMax: " << cur_corr << std::endl;
  };
  while (true) {
    bool stop = StopSearch();
    if (stop || cnt % print_frequence_ == 0) {
      PrintLog();
      if (stop) break;
    }
    ++cnt;

    double old_val = opt_val;
    int n_acpt = 0;
    int n_imp = 0;
    for (int i = 0; i < m_col_; ++i) {
      int col = i % k_;
      double inner_opt_val = __DBL_MAX__;
      double inner_opt_corr = __DBL_MAX__;
      ShuffleM(pair_list, j_col_pair_);
      std::pair<int, int> opt_pair;
      double corr_except = design.GetMaxAbsCorrExcept(col);
      for (int j = 0; j < j_col_pair_; ++j) {
        auto pair = pair_list[j];
        auto tmp_ret = design.PreSwapInCol(col, pair.first, pair.second);
        double tmp_corr = std::max(corr_except, tmp_ret.first);
        double tmp_val = design.GetCritVal(w_, tmp_corr, tmp_ret.second);
        if (tmp_val < inner_opt_val) {
          inner_opt_val = tmp_val;
          inner_opt_corr = tmp_corr;
          opt_pair = pair;
        }
      }
      if (inner_opt_val - cur_val <= T_h * uniform_dis(rng_)) {
        design.SwapInCol(col, opt_pair.first, opt_pair.second);
        ++n_acpt;
        if (inner_opt_val < opt_val) {
          opt_val = inner_opt_val;
          opt_design = design.GetDesign();
          ++n_imp;
        }
        cur_val = inner_opt_val;
        cur_corr = inner_opt_corr;
      }
    }

    bool flag_imp = old_val - opt_val > 0;
    max_no_imp_cnt = flag_imp ? 0 : max_no_imp_cnt + 1;
    double acpt_ratio = 1.0 * n_acpt / m_col_;
    double imp_ratio = 1.0 * n_imp / m_col_;
    if (flag_imp) {
      if (acpt_ratio > 0.1 && n_acpt > n_imp) {
        T_h *= alpha1_;
      } else if (acpt_ratio > 0.1 && n_acpt == n_imp) {
        ;
      } else {
        T_h /= alpha1_;
      }
    } else {
      if (acpt_ratio < 0.1) {
        T_h /= alpha3_;
        explore_dir = 1;
      } else if (acpt_ratio > 0.8) {
        T_h *= alpha2_;
        explore_dir = -1;
      } else if (explore_dir == 1) {
        T_h /= alpha3_;
      } else if (explore_dir == -1) {
        T_h *= alpha2_;
      }
    }
  }
  return Design(opt_design);
}

void ESE::InitDefaultParam() {
  int n_e = n_ * (n_ - 1) / 2;
  j_col_pair_ = std::min(50, (n_e + 4) / 5);
  m_col_ = k_;
  alpha1_ = 0.8;
  alpha2_ = 0.9;
  alpha3_ = 0.7;
  rng_.seed(0);
  w_ = 0.5;
  iterate_cnt_ = 1000;
  print_frequence_ = 100;
}

void ESE::ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m) {
  int n = pair_list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng_() % (n - i);
    std::swap(pair_list[i], pair_list[j]);
  }
}
} // namespace SearchingAlgorithm