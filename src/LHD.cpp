#include <algorithm>
#include <string>

#include "algorithm/algorithm_manager.h"

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: ${exe} ${algorithm} ${n} ${k} [${iterate_cnt} ${criteria} ${seed}]\n";
    std::cerr << "Search algorithms: SA, ESE, GA, LaPSO, LSGA\n";
    std::cerr << "Construction algorithms: Wang2018\n";
    std::cerr << "Example1: ./LHD Wang2018 20 8\n";
    std::cerr << "Example2: ./LHD ESE 20 8 1000 0.5PhiP15L2+0.5MaxAbsCor" << std::endl;
    return 1;
  }
  std::string algorithm = argv[1];
  int n = std::stoi(argv[2]);
  int k = std::stoi(argv[3]);
  auto algorithm_type = LHD::AlgorithmManager::GetAlgorithmType(algorithm);
  if (algorithm_type == LHD::AlgorithmManager::AlgorithmType::kSearchAlgorithm) {
    auto solver = LHD::AlgorithmManager::GetSearchAlgorithm(n, k, algorithm);
    if (argc > 4) {
      int iterate_cnt = std::stoi(argv[4]);
      solver->SetIterateCnt(iterate_cnt);
    }
    if (argc > 5) {
      const std::string& criteria = argv[5];
      solver->SetCriteria(criteria);
    }
    if (argc > 6) {
      int seed = std::stoi(argv[6]);
      solver->SetSeed(seed);
    }
    solver->SetLogOp(true);
    auto design = solver->Search();
    solver->DisplayResult(design);
  } else if (algorithm_type == LHD::AlgorithmManager::AlgorithmType::kConstructionAlgorithm) {
    auto solver = LHD::AlgorithmManager::GetConstructionAlgorithm(algorithm);
    auto design = solver->Solve(n, k);
    solver->DisplayResult(design);
  } else {
    std::cerr << "Unkown algorithm" << std::endl;
    return 1;
  }
  return 0;
}