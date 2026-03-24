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

  std::cout << "Adaptive Adj #0\n";
  printAdj(bank[0], p);

  std::cout << "Adaptive Adj #1\n";
  printAdj(bank[1], p);

  std::cout << "Adaptive Adj #2\n";
  printAdj(bank[2], p);

  std::cout << "Adaptive Adj #3\n";
  printAdj(bank[3], p);

  std::cout << "Adaptive Adj #4\n";
  printAdj(bank[4], p);

  return 0;
}