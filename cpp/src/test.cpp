#include "shtools_cpp.h"

int
main(int argc, char** argv)
{

  std::string infile = "../examples/ExampleDataFiles/MarsTopo719.shape";

  int lmax = 15;

  shtools::ShCoeff mars = shtools::sh_read(infile, lmax);

  std::cout << mars << std::endl;

  double val = shtools::make_grid_point(mars.begin(), mars.end(), 10.0, 30.0);
  std::cout << std::setprecision(16) << val << std::endl;
  std::cout << "diff to python " << val - 3395259.548270001 << std::endl;

  return 0;
}
