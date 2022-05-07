#include "SA.h"

#include <cassert>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>

namespace LHD {
SA::SA(int n, int k) : SearchAlgorithm(n, k) {
  SearchAlgorithm::InitDefaultParam();
  InitDefaultParam();
}

void SA::InitDefaultParam() {
  init_temp_ = 0.1;
  rate_ = 0.9;
  i_max_ = n_ * (n_ - 1) / 2 * k_ * 10;
}

Design::VecInt2D SA::Search() {
  return IncrementalSearch(Design(n_, k_, rng_()));
}

Design::VecInt2D SA::IncrementalSearch(Design design) {
  assert(design.GetN() == n_ && design.GetK() == k_);
  int swap_cnt = 0;
  design.InitCriteria(criteria_);
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  int cnt = 0;
  double opt_val = design.GetCriterion();
  Design::VecInt2D opt_design = design.GetDesign();
  double cur_val = opt_val;
  double temp = init_temp_;
  while (cnt < iterate_cnt_) {
    if (swap_cnt >= limit_swap_cnt_) break;
    int i = 0;
    bool imp_flag = false;
    while (i < i_max_ && cnt < iterate_cnt_) {
      if (swap_cnt >= limit_swap_cnt_) break;
      Log(cnt, design);
      ++cnt;

      int col = rng_() % k_;
      int r1 = rng_() % n_;
      int r2 = rng_() % n_;
      while (r1 == r2) r2 = rng_() % n_;
      design.SwapInCol(col, r1, r2);
      ++swap_cnt;
      double tmp_val = design.GetCriterion();

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
  return opt_design;
}
} // namespace LHD