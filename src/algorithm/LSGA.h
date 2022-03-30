#pragma once

#include "search_algorithm.h"

/*
An efficient local search-based genetic algorithm for constructing optimal Latin hypercube design
https://sci-hub.yncjkj.com/10.1080/0305215X.2019.1584618
*/

namespace LHD {
class LSGA : public SearchAlgorithm {
 public:
  LSGA(int n, int k);
  inline void SetPopulation_num(int population_num) { population_num_ = population_num; }
  inline void SetCrossoverRatio(double crossover_ratio) { crossover_ratio_ = crossover_ratio; }
  inline void SetMutationProb(double mutation_prob) { mutation_prob_ = mutation_prob; }
  inline void SetSelectMaxProb(double select_max_prob) { select_max_prob_ = select_max_prob; }
  inline void SetSelectMinProb(double select_min_prob) { select_min_prob_ = select_min_prob; }
  Design::VecInt2D Search() override;
 protected:
  void InitDefaultParam();
 private:
  int population_num_;
  double crossover_ratio_;
  double mutation_prob_;
  double select_max_prob_;
  double select_min_prob_;
};
} // namespace LHD