#pragma once

#include "search_algorithm.h"

/*
Optimizing Latin hypercube designs by particle swarm
https://sci-hub.yncjkj.com/10.1007/s11222-012-9363-3
*/

namespace LHD {
class LaPSO : public SearchAlgorithm {
 public:
  LaPSO(int n, int k);
  inline void SetParticleNum(int particle_num) { particle_num_ = particle_num; }
  inline void SetSameNumP(int same_num_p) { same_num_p_ = same_num_p; }
  inline void SetSameNumG(int same_num_g) { same_num_g_ = same_num_g; }
  inline void SetRatio(double ratio) { ratio_ = ratio; }
  Design::VecInt2D Search() override;
 protected:
  void InitDefaultParam();
 private:
  int particle_num_;
  int same_num_p_;
  int same_num_g_;
  double ratio_;
};

} // namespace LHD