#include "ESE.h"

#include <cassert>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>

namespace LHD {
ESE::ESE(int n, int k) : SearchAlgorithm(n, k) {
  SearchAlgorithm::InitDefaultParam();
  InitDefaultParam();
}

void ESE::InitDefaultParam() {
  int n_e = n_ * (n_ - 1) / 2;
  j_col_pair_ = std::max(1, std::min(50, n_e / 5));
  m_col_ = std::min(100, 2 * n_e * k_ / j_col_pair_);
  alpha1_ = 0.8;
  alpha2_ = 0.9;
  alpha3_ = 0.7;
  iterate_cnt_ = 1000;
}

Design::VecInt2D ESE::Search() {
  return IncrementalSearch(Design(n_, k_, rng_()).GetDesign());
}

Design::VecInt2D ESE::IncrementalSearch(const Design::VecInt2D& ori_design) {
  int swap_cnt = 0;
  Design design(ori_design);
  assert(design.GetN() == n_ && design.GetK() == k_);
  design.InitCriteria(criteria_);
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  int cnt = 0;
  double opt_val = design.GetCriterion();
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
  // int max_no_imp_cnt = 0;
  int explore_dir = 0;
  auto StopSearch = [this, &cnt/*, &max_no_imp_cnt*/]() -> bool {
    return cnt >= iterate_cnt_;
    // if (iterate_cnt_ >= 0) {
    //   return cnt >= iterate_cnt_;
    // }
    // return max_no_imp_cnt >= 1000;
  };
  while (!StopSearch()) {
    if (swap_cnt >= limit_swap_cnt_) break;
    Log(cnt, design);
    ++cnt;

    double old_val = opt_val;
    int n_acpt = 0;
    int n_imp = 0;
    for (int i = 0; i < m_col_; ++i) {
      int col = i % k_;
      double inner_opt_val = __DBL_MAX__;
      LHD::Utils::ShuffleM(pair_list, j_col_pair_, rng_);
      std::pair<int, int> opt_pair;
      for (int j = 0; j < j_col_pair_; ++j) {
        auto pair = pair_list[j];
        auto tmp_val = design.PreSwapInCol(col, pair.first, pair.second);
        ++swap_cnt;
        if (tmp_val < inner_opt_val) {
          inner_opt_val = tmp_val;
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
      }
    }

    bool flag_imp = old_val - opt_val > 0;
    // max_no_imp_cnt = flag_imp ? 0 : max_no_imp_cnt + 1;
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
  return opt_design;
}
} // namespace LHD