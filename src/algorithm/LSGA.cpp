#include "LSGA.h"

#include <algorithm>
#include <string>
#include <vector>
#include <set>

namespace LHD {
LSGA::LSGA(int n, int k) : SearchAlgorithm(n, k) {
  SearchAlgorithm::InitDefaultParam();
  InitDefaultParam();
}

void LSGA::InitDefaultParam() {
  population_num_ = 10;
  crossover_ratio_ = 0.5;
  mutation_prob_ = 0.2;
  select_max_prob_ = 0.3;
  select_min_prob_ = 0.01;
}

Design::VecInt2D LSGA::Search() {
  std::vector<Design> designs;
  designs.reserve(population_num_);
  for (int i = 0; i < population_num_; ++i) {
    designs.push_back(Design(n_, k_, rng_()));
    designs.back().InitCriteria(criteria_);
  }
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  int cnt = 0;
  int best_design_idx = 0;
  for (int i = 1; i < population_num_; ++i) {
    if (designs[i].GetCriterion() < designs[best_design_idx].GetCriterion()) {
      best_design_idx = i;
    }
  }
  std::vector<int> pos_list(n_);
  std::iota(pos_list.begin(), pos_list.end(), 0);
  auto AdjustCol = [this](Design& design, int col, const std::vector<int>& aim_col) {
    const auto& nums = design.GetDesignRef();
    std::vector<int> num_idx(n_ + 1);
    for (int i = 0; i < n_; ++i) {
      num_idx[nums[i][col]] = i;
    }
    for (int i = 0; i < n_; ++i) {
      int x = nums[i][col];
      int y = aim_col[i];
      if (x != y) {
        design.SwapInCol(col, i, num_idx[y]);
        num_idx[x] = num_idx[y];
        num_idx[y] = i;
      }
    }
  };
  auto GetCol = [this](const Design& design, int col) -> std::vector<int> {
    const auto& nums = design.GetDesignRef();
    std::vector<int> col_num(n_);
    for (int i = 0; i < n_; ++i) {
      col_num[i] = nums[i][col];
    }
    return col_num;
  };
  auto OX = [this](const std::vector<int>& par1_col, const std::vector<int>& par2_col, int r1, int r2) -> std::vector<int> {
    std::vector<int> aim_col = par1_col;
    std::set<int> sel_set;
    for (int i = r1; i <= r2; ++i) {
      sel_set.insert(par1_col[i]);
    }
    int par2_row_idx = 0;
    for (int i = 0, j = 0; i < n_; ++i) {
      if (i == r1) {
        i = r2;
        continue;
      }
      while (sel_set.count(par2_col[j])) {
        ++j;
      }
      aim_col[i] = par2_col[j];
    }
    return aim_col;
  };
  auto GenCrossoverPoint = [this](int& r1, int& c1, int& r2, int& c2) {
    int x, y;
    do {
      int x = rng_() % (k_ * n_);
      int y = rng_() % (k_ * n_);
    } while (std::abs(y - x) > n_ * k_ * crossover_ratio_);
    if (x > y) std::swap(x, y);
    r1 = x % n_, c1 = x / n_;
    r2 = y % n_, c2 = y / n_;
  };
  auto MOX = [&](const Design& par1, const Design& par2) -> Design {
    int r1, c1, r2, c2;
    GenCrossoverPoint(r1, c1, r2, c2);
    Design child = par2;
    if (c1 == c2) {
      if (r1 > r2) std::swap(r1, r2);
      AdjustCol(child, c1, OX(GetCol(par1, c1), GetCol(par2, c2), r1, r2));
      return child;
    }
    if (c1 > c2) {
      std::swap(c1, c2);
      std::swap(r1, r2);
    }
    AdjustCol(child, c1, OX(GetCol(par1, c1), GetCol(par2, c1), r1, n_ - 1));
    AdjustCol(child, c2, OX(GetCol(par1, c2), GetCol(par2, c2), 0, r2));
    for (int i = c1 + 1; i < c2; ++i) {
      AdjustCol(child, i, GetCol(par1, i));
    }
    return child;
  };
  while (cnt < iterate_cnt_) {
    Log(cnt, designs[best_design_idx]);
    ++cnt;
    const auto& best_design = designs[best_design_idx];
    auto new_designs = designs;
    for (int i = 0; i < population_num_; ++i) {
      if (i == best_design_idx) continue;
      new_designs[i] = MOX(best_design, designs[i]);
    }

    for (int i = 0; i < population_num_; ++i) {
      if (i == best_design_idx) continue;
      if (uniform_dis(rng_) < mutation_prob_) {
        for (int j = 0; j < k_; ++j) {
          Utils::ShuffleM(pos_list, 2, rng_);
          new_designs[i].SwapInCol(j, pos_list[0], pos_list[1]);
        }
      } else {
        int col = rng_() % k_;
        Utils::ShuffleM(pos_list, 2, rng_);
        new_designs[i].SwapInCol(col, pos_list[0], pos_list[1]);
      }
    }

    for (int i = 0; i < population_num_; ++i) {
      if (i == best_design_idx) continue;
      double select_prob = select_max_prob_ - 
        (select_max_prob_ - select_min_prob_) * cnt / iterate_cnt_;
      if (new_designs[i].GetCriterion() < designs[i].GetCriterion() ||
          uniform_dis(rng_) < select_prob) {
        designs[i] = new_designs[i];
      }
    }

    best_design_idx = 0;
    for (int i = 1; i < population_num_; ++i) {
      if (designs[i].GetCriterion() < designs[best_design_idx].GetCriterion()) {
        best_design_idx = i;
      }
    }

    for (int i = 0; i < k_; ++i) {
      Utils::ShuffleM(pos_list, 2, rng_);
      double val = designs[best_design_idx].GetCriterion();
      designs[best_design_idx].SwapInCol(i, pos_list[0], pos_list[1]);
      double new_val = designs[best_design_idx].GetCriterion();
      if (new_val >= val) {
        designs[best_design_idx].SwapInCol(i, pos_list[0], pos_list[1]);
      }
    }
  }
  return designs[best_design_idx].GetDesign();
}
} // namespace LHD