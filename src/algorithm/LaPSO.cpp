#include "LaPSO.h"

#include <algorithm>
#include <string>
#include <vector>
#include <random>

namespace LHD {
LaPSO::LaPSO(int n, int k) : SearchAlgorithm(n, k) {
  SearchAlgorithm::InitDefaultParam();
  InitDefaultParam();
}

void LaPSO::InitDefaultParam() {
  particle_num_ = 10;
  same_num_p_ = 0;
  same_num_g_ = (n_ + 3) / 4;
  ratio_ = 1.0 / k_;
}

Design::VecInt2D LaPSO::Search() {
  std::vector<Design> designs;
  designs.reserve(particle_num_);
  for (int i = 0; i < particle_num_; ++i) {
    designs.push_back(Design(n_, k_, rng_()));
    designs.back().InitCriteria(criteria_);
  }
  std::uniform_real_distribution<double> uniform_dis(std::nextafter(0.0, 1.0), 1.0);
  int cnt = 0;
  auto pbest = designs;
  auto gbest = designs[0];
  for (int i = 1; i < particle_num_; ++i) {
    if (designs[i].GetCriterion() < gbest.GetCriterion()) {
      gbest = designs[i];
    }
  }
  std::vector<int> pos_list(n_);
  std::iota(pos_list.begin(), pos_list.end(), 0);
  while (cnt < iterate_cnt_) {
    Log(cnt, gbest);
    ++cnt;
    for (int i = 0; i < particle_num_; ++i) {
      for (int j = 0; j < k_; ++j) {
        std::vector<int> num_idx(n_ + 1);
        const auto& nums = designs[i].GetDesignRef();
        for (int o = 0; o < n_; ++o) {
          num_idx[nums[o][j]] = o;
        }

        LHD::Utils::ShuffleM(pos_list, same_num_p_, rng_);
        for (int o = 0; o < same_num_p_; ++o) {
          int p = pos_list[o];
          int x = designs[i].GetDesignRef()[p][j];
          int y = pbest[i].GetDesignRef()[p][j];
          if (x == y) continue;
          designs[i].SwapInCol(j, p, num_idx[y]);
          num_idx[x] = num_idx[y];
          num_idx[y] = p;
        }
        
        LHD::Utils::ShuffleM(pos_list, same_num_g_, rng_);
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
          LHD::Utils::ShuffleM(pos_list, 2, rng_);
          designs[i].SwapInCol(j, pos_list[0], pos_list[1]);
        }
      }
    }
    for (int i = 0; i < particle_num_; ++i) {
      if (designs[i].GetCriterion() < pbest[i].GetCriterion()) {
        pbest[i] = designs[i];
      }
      if (designs[i].GetCriterion() < gbest.GetCriterion()) {
        gbest = designs[i];
      }
    }
  }
  return gbest.GetDesign();
}
} // namespace LHD