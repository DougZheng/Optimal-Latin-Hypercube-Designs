#pragma once

#include "search_algorithm.h"

/*
Optimizing Latin hypercube designs by particle swarm
https://sci-hub.yncjkj.com/10.1007/s11222-012-9363-3
*/

namespace LHD {
class GA : public SearchAlgorithm {
 public:
  GA(int n, int k);
  inline void SetPopulationNum(int population_num) { population_num_ = population_num; }
  inline void SetMutationProb(double mutation_prob) { mutation_prob_ = mutation_prob; }
  Design::VecInt2D Search() override;
 protected:
  void InitDefaultParam();
 private:
  int population_num_;
  double mutation_prob_;
};
} // namespace LHD