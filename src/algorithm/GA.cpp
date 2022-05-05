#include "GA.h"

#include <algorithm>
#include <string>
#include <vector>
#include <random>

namespace LHD {
GA::GA(int n, int k) : SearchAlgorithm(n, k) {
  SearchAlgorithm::InitDefaultParam();
  InitDefaultParam();
}

void GA::InitDefaultParam() {
  population_num_ = 10;
  mutation_prob_ = 1.0 / k_;
}

Design::VecInt2D GA::Search() {
  std::vector<Design> designs;
  designs.reserve(population_num_);
  for (int i = 0; i < population_num_; ++i) {
    designs.push_back(Design(n_, k_, rng_()));
    designs.back().InitCriteria(criteria_);
  }
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  int cnt = 0;
  std::vector<int> rank(population_num_);
  std::iota(rank.begin(), rank.end(), 0);
  std::sort(rank.begin(), rank.end(), [&](int x, int y) {
    return designs[x].GetCriterion() < designs[y].GetCriterion();
  });
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
  while (cnt < iterate_cnt_) {
    Log(cnt, designs[rank[0]]);
    ++cnt;
    const auto& best_design = designs[rank[0]];
    designs[rank[population_num_ / 2]] = best_design;
    for (int i = 1; i < population_num_ / 2; ++i) {
      int col = rng_() % k_;
      auto& design = designs[rank[i + population_num_ / 2]];
      const auto& aim_col = GetCol(best_design, col);
      design = designs[rank[i]];
      AdjustCol(design, col, aim_col);
    }
    for (int i = 1; i < population_num_ / 2; ++i) {
      int col = rng_() % k_;
      auto& design = designs[rank[i]];
      const auto& aim_col = GetCol(design, col);
      design = best_design;
      AdjustCol(design, col, aim_col);
    }
    for (int i = 0; i < population_num_; ++i) {
      if (i == rank[0]) continue;
      for (int j = 0; j < k_; ++j) {
        if (uniform_dis(rng_) < mutation_prob_) {
          LHD::Utils::ShuffleM(pos_list, 2, rng_);
          designs[i].SwapInCol(j, pos_list[0], pos_list[1]);
        }
      }
    }
    std::sort(rank.begin(), rank.end(), [&](int x, int y) {
      return designs[x].GetCriterion() < designs[y].GetCriterion();
    });
  }
  return designs[rank[0]].GetDesign();
}
} // namespace LHD