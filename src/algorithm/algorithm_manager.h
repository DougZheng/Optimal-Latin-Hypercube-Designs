#pragma once

#include <set>
#include <string>
#include <memory>

#include "search_algorithm.h"
#include "ESE.h"
#include "GA.h"
#include "LaPSO.h"
#include "LSGA.h"
#include "SA.h"
#include "construction_algorithm.h"
#include "Wang2018.h"

namespace LHD {
class AlgorithmManager {
 public:
  enum class AlgorithmType {
    kUnknown = 0,
    kSearchAlgorithm = 1,
    kConstructionAlgorithm = 2,
  };
  static AlgorithmType GetAlgorithmType(const std::string& algorithm) {
    static const std::set<std::string> search_algorithm_set{
      "ESE", "SA", "GA", "LaPSO", "LSGA",
    };
    static const std::set<std::string> construction_algorithm_set{
      "Wang2018",
    };
    if (search_algorithm_set.count(algorithm)) {
      return AlgorithmType::kSearchAlgorithm;
    } else if (construction_algorithm_set.count(algorithm)) {
      return AlgorithmType::kConstructionAlgorithm;
    } else {
      return AlgorithmType::kUnknown;
    }
  }
  static std::shared_ptr<SearchAlgorithm> GetSearchAlgorithm(int n, int k, const std::string& algorihtm = "Default") {
    if (algorihtm == "ESE") {
      return std::make_shared<ESE>(n, k);
    } else if (algorihtm == "SA") {
      return std::make_shared<SA>(n, k);
    } else if (algorihtm == "GA") {
      return std::make_shared<GA>(n, k);
    } else if (algorihtm == "LaPSO") {
      return std::make_shared<LaPSO>(n, k);
    } else if (algorihtm == "LSGA") {
      return std::make_shared<LSGA>(n, k);
    } else {
      return std::make_shared<ESE>(n, k);
    }
  }
  static std::shared_ptr<ConstructionAlgorithm> GetConstructionAlgorithm(const std::string& algorithm = "Default") {
    if (algorithm == "Wang2018") {
      return std::make_shared<Wang2018>();
    } else {
      return std::make_shared<Wang2018>();
    }
  }
};
} // namespace LHD