#include <iostream>
#include "../../mmcfilters/attributes/maxdist/EdtDIFT.hpp"

using namespace mmcfilters::maxdist;

void printAdj(const AdaptiveAdj &adj, const Point2D p)
{
  std::cout << "Neighbors for " << p << ":\n";
  for (const auto &[q, a] : adj.neighbors(p)) {
    std::cout << "\t" << q << " and " << a << "\n";
  }
}

int main()
{
  AdaptiveAdjBank bank;
  Point2D p{ 5, 5 };

  for (int i = 0; i < bank.size(); i++) {
    std::cout << "Adaptive Adj #" << i << "\n";
    printAdj(bank[i], p);
  }
  
  return 0;
}