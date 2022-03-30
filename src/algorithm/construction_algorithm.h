#pragma once

#include <vector>

namespace LHD {
class ConstructionAlgorithm {
 public:
  ConstructionAlgorithm() = default;
  virtual std::vector<std::vector<int>> Solve(int n, int k) = 0;
  virtual void DisplayResult(const std::vector<std::vector<int>>& d) const = 0;
};
} // namespace LHD