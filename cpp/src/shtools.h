#pragma once

#include "sh_tools_wrapper.h"

namespace shtools {

constexpr int n2 = 2;

// n is number of cosine coefficents or sine coefficents
inline int
n_to_deg(int n)
{
  return std::sqrt(1. / 4. + 2. * n) - 3. / 2.;
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

inline std::vector<double>
sh_cindex_to_cilm(const std::vector<double>& cindex, int degmax=-1)
{

  int exitstatus;
  
  int n = cindex.size() / 2;
  int cilm_dim = n_to_deg(n)+1;
  
  if(degmax < 0){
      degmax = cilm_dim-1;
  }
    
  std::vector<double> cilm(2 * cilm_dim * cilm_dim);
  cSHCindexToCilm(&cindex[0], &n2, &n, &cilm[0], &cilm_dim, &degmax, &exitstatus);
  
  return cilm;
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
