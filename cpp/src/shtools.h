#pragma once

#include "sh_tools_wrapper.h"

namespace shtools {

constexpr int n2 = 2;

inline int
deg_2_n(int deg)
{
  return std::sqrt(1. / 4. + 2. * deg) - 3. / 2.;
}

inline double
make_grid_point(const std::vector<double>& cilm,
                double lat,
                double lon,
                int norm = 1,
                int csphase = 1,
                int dealloc = 0)
{

  int n = cilm.size();
  int cilmd = std::sqrt(n / 2);
  int lmax = cilmd - 1;

  return cMakeGridPoint(
    &cilm[0], &cilmd, &lmax, &lat, &lon, &norm, &csphase, &dealloc);
}

inline void
sh_cindex_to_cilm(std::vector<double> cindex,
                  std::vector<double> cilm,
                  int degmax)
{

  int exitstatus;

  int n = cindex.size() / 2;
  int m = deg_2_n(n);

  cSHCindexToCilm(&cindex[0], &n2, &n, &cilm[0], &m, &degmax, &exitstatus);
}

inline std::vector<double>
sh_read(const std::string& filename,
        int degree,
        int skip = 0,
        std::vector<double>* error = nullptr)
{

  int exitstatus;
  int cilm_dim = degree + 1;

  std::vector<double> cilm(2 * cilm_dim * cilm_dim);

  double* error_ptr = nullptr;
  int error_dim = 0;
  if (error) {
    error_ptr = &error->at(0);
    error_dim = cilm_dim;
    error->resize(cilm.size());
  }

  double* header = nullptr;
  int header_d = 0;

  int s = filename.size();

  cSHRead(filename.c_str(),
          &s,
          &cilm[0],
          &cilm_dim,
          &degree,
          &skip,
          header,
          &header_d,
          error_ptr,
          &n2,
          &error_dim,
          &error_dim,
          &exitstatus);

  return cilm;
}

typedef Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::ColMajor> Cindex;
typedef Eigen::Tensor<double, 3, Eigen::ColMajor> Cilm;

}
