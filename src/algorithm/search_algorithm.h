#pragma once

#include <iostream>
#include <string>
#include <random>

#include "comm/design.h"
#include "comm/utils.h"

namespace LHD {
class SearchAlgorithm {
 public:
  SearchAlgorithm(int n, int k) : n_(n), k_(k) { InitDefaultParam(); }
  inline void SetCriteria(const std::string& criteria) { criteria_ = criteria; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetLogFrequence(int log_frequence) { log_frequence_ = log_frequence; }
  inline void SetLogOp(bool log_op) { log_op_ = log_op; }
  inline void SetIterateCnt(int iterate_cnt) { iterate_cnt_ = iterate_cnt; }
  virtual Design::VecInt2D Search() = 0;
  void DisplayResult(const Design::VecInt2D& d) const {
    Design design(d);
    design.InitCriteria(criteria_);
    design.Display();
  }
 protected:
  void InitDefaultParam() {
    criteria_ = LHD::Utils::SMultiCriteria(LHD::Utils::SPhiPL2(0.5), LHD::Utils::SMaxAbsCor(0.5));
    rng_.seed(0);
    log_frequence_ = 100;
    log_op_ = false;
    iterate_cnt_ = 1000;
  }
  void Log(int iterate_cnt, const Design& design) const {
    if (!log_op_ || iterate_cnt % log_frequence_ != 0) return;
    std::cerr << "Iteration: " << iterate_cnt << ", " <<
      design.GetCriteriaLogString() << std::endl;
  }
 protected:
  int n_;
  int k_;
  std::string criteria_;
  std::mt19937 rng_;
  int log_frequence_;
  bool log_op_;
  int iterate_cnt_;
};
} // namespace LHD