#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>

#include "algorithm/algorithm_manager.h"
#include "comm/utils.h"

void Display(const std::vector<std::vector<double>>& design) {
  std::cout << std::fixed << std::setprecision(2);
  double max_num = 0;
  for (int i = 0; i < design.size(); ++i) {
    max_num = std::max(max_num, *std::max_element(design[i].begin(), design[i].end()));
  }
  int w = std::floor(std::log10(max_num)) + 1 + 4;
  for (int i = 0; i < design.size(); ++i) {
    for (int j = 0; j < design[i].size(); ++j) {
      std::cout << std::setw(w) << design[i][j] << " \n"[j + 1 == design[i].size()];
    }
  }
  std::cout.flush();
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: ${exe} ${n} ${k} [${iterate_cnt}]\n";
    return 1;
  }
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  auto solver = LHD::AlgorithmManager::GetSearchAlgorithm(n, k); // Default: ESE
  { // Optional params
    if (argc > 3) {
      int iterate_cnt = std::stoi(argv[3]);
      solver->SetIterateCnt(iterate_cnt); // Default: 1000
    }
    solver->SetLogOp(true); // Debug info
  }
  auto ori_design = solver->Search(); // n x k Latin hypercube
  std::vector<std::pair<double, double>> ranges(k, {0, 1});
  auto design = LHD::Utils::ConvertLevel(ori_design, ranges); // n x k Latin hypercube design
  Display(design);
  return 0;
}