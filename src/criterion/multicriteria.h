#pragma once

#include "criterion.h"

namespace LHD {
class MultiCriteria {
 public:
  MultiCriteria() = default;
  MultiCriteria(const MultiCriteria& o) { *this = o; }
  MultiCriteria& operator=(const MultiCriteria& o);
  double GetCriterion() const;
  void SwapInCol(int col, int r1, int r2);
  double PreSwapInCol(int col, int r1, int r2) const;
  void SetDesign(const std::vector<std::vector<int>>* design);
  void AddCritirion(std::shared_ptr<Criterion> criterion, double weight);
  std::string GetOneLineLogString() const;
  std::string GetDescriptionString() const;
  void Clear();
 private:
  std::vector<std::shared_ptr<Criterion>> criteria_;
  std::vector<double> weights_;
};
} // namespace LHD