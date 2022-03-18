#pragma once

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
Exploratory designs for computational experiments
https://sci-hub.yncjkj.com/10.1016/0378-3758(94)00035-t
*/

namespace SearchingAlgorithm {
class SA {
 public:
  SA(int n, int k);
  inline void SetInitTemp(double init_temp) { init_temp_ = init_temp; }
  inline void SetRate(double rate) { rate_ = rate; }
  inline void SetIMax(int i_max) { i_max_ = i_max; }
  inline void SetW(double w) { w_ = w; }
  inline void SetRepeatedCnt(int repeated_cnt) { repeated_cnt_ = repeated_cnt; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetPrintFrequence(int print_frequence) { print_frequence_ = print_frequence; }
  Design Search();
  Design IncrementalSearch(Design design);
 private:
  void InitDefaultParam();
  Design SearchOnce(Design Design);
 private:
  int n_;
  int k_;
  double init_temp_;
  double rate_;
  int i_max_;
  double w_;
  int repeated_cnt_;
  std::mt19937 rng_;
  int print_frequence_;
};

SA::SA(int n, int k) : n_(n), k_(k) {
  InitDefaultParam();
}

Design SA::Search() {
  Design design(n_, k_, rng_());
  return IncrementalSearch(design);
}

Design SA::IncrementalSearch(Design design) {
  assert(design.GetN() == n_ && design.GetK() == k_);
  std::vector<std::pair<double, double>> ret_val(repeated_cnt_);
  Design::VecInt2D opt_design = design.GetDesign();
  double opt_val = design.GetCritVal(w_, design.GetRhoMax(), design.GetPhiP());
  for (int i = 0; i < repeated_cnt_; ++i) {
    std::cerr << "Repeation: " << i << std::endl;
    Design search_ret = SearchOnce(design);
    ret_val[i] = {search_ret.GetRhoMax(), search_ret.GetPhiP()};
    double val = design.GetCritVal(w_, ret_val[i].first, ret_val[i].second);
    if (val < opt_val) {
      opt_val = val;
      opt_design = search_ret.GetDesign();
    }
  }
  for (int i = 0; i < repeated_cnt_; ++i) {
    std::cerr << "Repeation: " << i
      << ", Val: " << design.GetCritVal(w_, ret_val[i].first, ret_val[i].second)
      << ", PhiP: " << ret_val[i].second
      << ", RhoMax: " << ret_val[i].first << std::endl;
  }
  return Design(opt_design);
}

Design SA::SearchOnce(Design design) {
  if (w_ == 0) design.DisableCorr();
  if (w_ == 1) design.DisableDis();
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int cnt = 0;
  double opt_val = design.GetCritVal(w_, design.GetRhoMax(), design.GetPhiP());
  Design::VecInt2D opt_design = design.GetDesign();
  double cur_val = opt_val;
  double temp = init_temp_;
  auto PrintLog = [this, &cnt, &cur_val, &design]() -> void {
    std::cerr << "Iteration: " << cnt
      << ", Val: " << cur_val 
      << ", PhiP: " << design.GetPhiP()
      << ", RhoMax: " << design.GetRhoMax() << std::endl;
  };
  while (true) {
    int i = 0;
    bool imp_flag = false;
    while (i < i_max_) {
      if (cnt % print_frequence_ == 0) {
        PrintLog();
      }
      ++cnt;

      int col = rng_() % k_;
      int r1 = rng_() % n_;
      int r2 = rng_() % n_;
      while (r1 == r2) r2 = rng_() % n_;
      design.SwapInCol(col, r1, r2);
      double tmp_val = design.GetCritVal(w_, design.GetRhoMax(), design.GetPhiP());

      if (tmp_val < cur_val) {
        imp_flag = true;
        cur_val = tmp_val;
      } else {
        double acpt_prob = std::exp((cur_val - tmp_val) / temp);
        bool acpt = uniform_dis(rng_) <= acpt_prob;
        if (acpt) {
          imp_flag = true;
          cur_val = tmp_val;
        } else {
          design.SwapInCol(col, r1, r2); // restore
        }
      }

      if (cur_val < opt_val) {
        opt_val = cur_val;
        opt_design = design.GetDesign();
        i = 0;
      } else {
        ++i;
      }
    }
    temp *= rate_;
    if (!imp_flag) break;
  }
  return Design(opt_design);
}

void SA::InitDefaultParam() {
  init_temp_ = 0.5;
  rate_ = 0.5;
  i_max_ = n_ * (n_ - 1) / 2 * k_;
  repeated_cnt_ = 10;
  rng_.seed(0);
  print_frequence_ = 1000;
}
} // namespace SearchingAlgorithm