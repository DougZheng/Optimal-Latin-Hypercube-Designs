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
class LaPSO {
 public:
  LaPSO(int n, int k);
  inline void SetParticleNum(int particle_num) { particle_num_ = particle_num; }
  inline void SetSameNumP(int same_num_p) { same_num_p_ = same_num_p; }
  inline void SetSameNumG(int same_num_g) { same_num_g_ = same_num_g; }
  inline void SetRatio(double ratio) { ratio_ = ratio; }
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
  int particle_num_;
  int same_num_p_;
  int same_num_g_;
  double ratio_;
  double w_;
  int iterate_cnt_;
  std::mt19937 rng_;
  int print_frequence_;
};

LaPSO::LaPSO(int n, int k) : n_(n), k_(k) {
  InitDefaultParam();
}

Design LaPSO::Search() {
  std::vector<Design> designs;
  designs.reserve(particle_num_);
  for (int i = 0; i < particle_num_; ++i) {
    designs.push_back(Design(n_, k_, rng_()));
    if (w_ == 0) designs.back().DisableCorr();
    if (w_ == 1) designs.back().DisableDis();
  }
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int cnt = 0;
  auto pbest = designs;
  auto gbest = designs[0];
  for (int i = 1; i < particle_num_; ++i) {
    if (designs[i].GetCritVal(w_) < gbest.GetCritVal(w_)) {
      gbest = designs[i];
    }
  }
  std::vector<int> pos_list(n_);
  std::iota(pos_list.begin(), pos_list.end(), 0);
  auto PrintLog = [this, &cnt, &gbest]() -> void {
    std::cerr << "Iteration: " << cnt
      << ", Val: " << gbest.GetCritVal(w_) 
      << ", PhiP: " << gbest.GetPhiP()
      << ", RhoMax: " << gbest.GetRhoMax() << std::endl;
  };
  while (cnt < iterate_cnt_) {
    if (cnt % print_frequence_ == 0) {
      PrintLog();
    }
    ++cnt;
    for (int i = 0; i < particle_num_; ++i) {
      for (int j = 0; j < k_; ++j) {
        std::vector<int> num_idx(n_ + 1);
        const auto& nums = designs[i].GetDesignRef();
        for (int o = 0; o < n_; ++o) {
          num_idx[nums[o][j]] = o;
        }

        ShuffleM(pos_list, same_num_p_);
        for (int o = 0; o < same_num_p_; ++o) {
          int p = pos_list[o];
          int x = designs[i].GetDesignRef()[p][j];
          int y = pbest[i].GetDesignRef()[p][j];
          if (x == y) continue;
          designs[i].SwapInCol(j, p, num_idx[y]);
          num_idx[x] = num_idx[y];
          num_idx[y] = p;
        }
        
        ShuffleM(pos_list, same_num_g_);
        for (int o = 0; o < same_num_g_; ++o) {
          int p = pos_list[o];
          int x = designs[i].GetDesignRef()[p][j];
          int y = gbest.GetDesignRef()[p][j];
          if (x == y) continue;
          designs[i].SwapInCol(j, p, num_idx[y]);
          num_idx[x] = num_idx[y];
          num_idx[y] = p;
        }

        if (uniform_dis(rng_) < ratio_) {
          ShuffleM(pos_list, 2);
          designs[i].SwapInCol(j, pos_list[0], pos_list[1]);
        }
      }
    }
    for (int i = 0; i < particle_num_; ++i) {
      if (designs[i].GetCritVal(w_) < pbest[i].GetCritVal(w_)) {
        pbest[i] = designs[i];
      }
      if (designs[i].GetCritVal(w_) < gbest.GetCritVal(w_)) {
        gbest = designs[i];
      }
    }
  }
  return gbest;
}

void LaPSO::InitDefaultParam() {
  rng_.seed(0);
  w_ = 0.5;
  particle_num_ = 10;
  same_num_p_ = 0;
  same_num_g_ = (n_ + 3) / 4;
  ratio_ = 1.0 / k_;
  iterate_cnt_ = 1000;
  print_frequence_ = 100;
}

void LaPSO::ShuffleM(std::vector<int>& pos_list, int m) {
  int n = pos_list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng_() % (n - i);
    std::swap(pos_list[i], pos_list[j]);
  }
}
} // namespace SearchingAlgorithm