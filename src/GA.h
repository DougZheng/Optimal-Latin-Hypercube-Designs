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
Optimizing Latin hypercube designs by particle swarm
https://sci-hub.yncjkj.com/10.1007/s11222-012-9363-3
*/

namespace SearchingAlgorithm {
class GA {
 public:
  GA(int n, int k);
  inline void SetPopulation_num(int population_num) { population_num_ = population_num; }
  inline void SetMutationProb(double mutation_prob) { mutation_prob_ = mutation_prob; }
  inline void SetW(double w) { w_ = w; }
  inline void SetIterateCnt(int iterate_cnt) { iterate_cnt_ = iterate_cnt; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetPrintFrequence(int print_frequence) { print_frequence_ = print_frequence; }
  Design Search();
 private:
  void InitDefaultParam();
  void ShuffleM(std::vector<int>& pos_list, int m);
 private:
  int n_;
  int k_;
  int population_num_;
  double mutation_prob_;
  double w_;
  int iterate_cnt_;
  std::mt19937 rng_;
  int print_frequence_;
};

GA::GA(int n, int k) : n_(n), k_(k) {
  InitDefaultParam();
}

Design GA::Search() {
  std::vector<Design> designs;
  designs.reserve(population_num_);
  for (int i = 0; i < population_num_; ++i) {
    designs.push_back(Design(n_, k_, rng_()));
    if (w_ == 0) designs.back().DisableCorr();
    if (w_ == 1) designs.back().DisableDis();
  }
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int cnt = 0;
  std::vector<int> rank(population_num_);
  std::iota(rank.begin(), rank.end(), 0);
  std::sort(rank.begin(), rank.end(), [&](int x, int y) {
    return designs[x].GetCritVal(w_) < designs[y].GetCritVal(w_);
  });
  std::vector<int> pos_list(n_);
  std::iota(pos_list.begin(), pos_list.end(), 0);
  auto PrintLog = [this, &cnt, &rank, &designs]() -> void {
    const auto& gbest = designs[rank[0]];
    std::cerr << "Iteration: " << cnt
      << ", Val: " << gbest.GetCritVal(w_) 
      << ", PhiP: " << gbest.GetPhiP()
      << ", RhoMax: " << gbest.GetRhoMax() << std::endl;
  };
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
    if (cnt % print_frequence_ == 0) {
      PrintLog();
    }
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
    for (int i = 1; i < population_num_; ++i) {
      for (int j = 0; j < k_; ++j) {
        if (uniform_dis(rng_) < mutation_prob_) {
          ShuffleM(pos_list, 2);
          designs[i].SwapInCol(j, pos_list[0], pos_list[1]);
        }
      }
    }
    std::sort(rank.begin(), rank.end(), [&](int x, int y) {
      return designs[x].GetCritVal(w_) < designs[y].GetCritVal(w_);
    });
  }
  return designs[rank[0]];
}

void GA::InitDefaultParam() {
  rng_.seed(0);
  w_ = 0.5;
  population_num_ = 10;
  mutation_prob_ = 1.0 / k_;
  iterate_cnt_ = 1000;
  print_frequence_ = 100;
}

void GA::ShuffleM(std::vector<int>& pos_list, int m) {
  int n = pos_list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng_() % (n - i);
    std::swap(pos_list[i], pos_list[j]);
  }
}
} // namespace SearchingAlgorithm