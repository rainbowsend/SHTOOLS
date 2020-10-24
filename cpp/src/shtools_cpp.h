#pragma once

#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

#include "shtools.h"

namespace shtools {

class ShCoeff
{
public:
  ShCoeff(int degree)
  {
    this->degree = degree;
    cilm = std::vector<double>(cilm_size(degree), 0.0);
  }

  ShCoeff(const std::vector<double>& cilm)
  {
    this->degree = std::sqrt(cilm.size() / 2) - 1;
    this->cilm = cilm;
  }

  // number of elements in cilm vector given the degree
  static constexpr size_t cilm_size(int degree)
  {
    return 2 * std::pow(degree + 1, 2);
  };
  size_t cilm_size() { return cilm_size(degree); };

  // degree given the number C or S coefficents
  static constexpr size_t calc_degree(size_t n_coefficents)
  {
    return std::sqrt(1. / 4. + 2. * n_coefficents) - 3. / 2.;
  };

  // size of second ond third dimension used in fortran cilm matrix
  int cilm_dim() const { return degree + 1; };

  int get_degree() const { return degree; };

  double* cilm_ptr() { return &cilm[0]; };

  auto begin() { return cilm.begin(); };
  auto end() { return cilm.end(); };

  double get_C(int degree, int order) const
  {
    if (degree > this->degree) {
      throw std::out_of_range("maximal degree is " +
                              std::to_string(this->degree));
    }
    if (order > degree) {
      throw std::out_of_range("order has to be samller than degree");
    }
    size_t c = (degree + cilm_dim() * order) * 2;
    return cilm[c];
  };

  double get_S(int degree, int order) const
  {
    if (degree > this->degree) {
      throw std::out_of_range("maximal degree is " +
                              std::to_string(this->degree));
    }
    if (order > degree) {
      throw std::out_of_range("order has to be samller than degree");
    }
    size_t c = (degree + cilm_dim() * order) * 2 + 1;
    return cilm[c];
  };

  std::string to_string_triangle() const
  {
    std::stringstream ss;
    for (int degree = 0; degree <= get_degree(); ++degree) {
      for (int order = 0; order <= degree; ++order) {
        ss << std::scientific << std::setw(20) << std::setprecision(12)
           << get_C(degree, order);
      }
      ss << std::endl;
    }
    ss << std::endl;
    for (int degree = 0; degree <= get_degree(); ++degree) {
      for (int order = 0; order <= degree; ++order) {
        ss << std::scientific << std::setw(20) << std::setprecision(12)
           << get_S(degree, order);
      }
      ss << std::endl;
    }
    return ss.str();
  }

  std::string to_string_list() const
  {
    std::stringstream ss;
    for (int degree = 0; degree <= get_degree(); ++degree) {
      for (int order = 0; order <= degree; ++order) {
        ss << std::setw(4) << degree << " " << std::setw(4) << order << " "
           << std::scientific << std::setw(20) << std::setprecision(12)
           << get_C(degree, order) << " " << std::scientific << std::setw(20)
           << std::setprecision(12) << get_S(degree, order) << std::endl;
      }
    }
    return ss.str();
  }

private:
  int degree = 0;
  // cilm has shape (2,cilm_dim,cilm_dim)
  std::vector<double> cilm;
};

inline std::ostream&
operator<<(std::ostream& os, const ShCoeff& shc)
{
  if (shc.get_degree() < 15) {
    os << shc.to_string_triangle();
  } else {
    os << shc.to_string_list();
  }
  return os;
}

template<class InputIt>
double
make_grid_point(const InputIt cilm_first,
                const InputIt cilm_last,
                double lat,
                double lon,
                int norm = 1,
                int csphase = 1,
                int dealloc = 0)
{

  int n = std::distance(cilm_first, cilm_last);
  int cilmd = std::sqrt(n / 2);
  int lmax = cilmd - 1;

  return MakeGridPoint(
    &*cilm_first, cilmd, lmax, lat, lon, &norm, &csphase, &dealloc);
}

template<class InputIt, class OutputIt>
constexpr OutputIt
sh_cindex_to_cilm(const InputIt cindex_first,
                  const InputIt cindex_last,
                  OutputIt cilm_first,
                  int degmax = -1)
{

  int exitstatus;

  int n = std::distance(cindex_first, cindex_last) / 2;
  int cilm_dim = ShCoeff::calc_degree(n) + 1;

  if (degmax < 0) {
    degmax = cilm_dim - 1;
  }
  SHCindexToCilm(
    &*cindex_first, n, &*cilm_first, cilm_dim, degmax, &exitstatus);

  return cilm_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt
sh_cilm_to_cindex(const InputIt cilm_first,
                  const InputIt cilm_last,
                  OutputIt cindex_first,
                  int degmax = -1)
{

  int exitstatus;

  int n = std::distance(cilm_first, cilm_last);
  int cilm_dim = sqrt(n / 2) - 1;

  if (degmax < 0) {
    degmax = cilm_dim - 1;
  }
  SHCilmToCindex(
    &*cilm_first, cilm_dim, &*cindex_first, n, degmax, &exitstatus);

  return cindex_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt
sh_vector_to_cilm(const InputIt vector_first,
                  const InputIt vector_last,
                  OutputIt cilm_first,
                  int degmax = -1)
{

  int exitstatus;

  int n = std::distance(vector_first, vector_last);
  int cilm_dim = sqrt(n);

  if (degmax < 0) {
    degmax = cilm_dim - 1;
  }
  SHVectorToCilm(&*vector_first, &*cilm_first, cilm_dim, degmax, &exitstatus);

  return cilm_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt
sh_cilm_to_vector(const InputIt cilm_first,
                  const InputIt cilm_last,
                  OutputIt vector_first,
                  int degmax = -1)
{

  int exitstatus;

  int n = std::distance(cilm_first, cilm_last);
  int cilm_dim = sqrt(n / 2);

  if (degmax < 0) {
    degmax = cilm_dim - 1;
  }
  SHCilmToVector(&*cilm_first, cilm_dim, &*vector_first, degmax, &exitstatus);

  return vector_first;
}

inline ShCoeff
sh_read(const std::string& filename,
        int degree,
        int skip = 0,
        std::vector<double>* error = nullptr)
{

  int exitstatus;
  int cilm_dim = degree + 1;

  ShCoeff result(degree);

  double* error_ptr = nullptr;
  if (error) {
    error_ptr = &error->at(0);
    error->resize(result.cilm_size());
  }

  double* header = nullptr;
  int header_d = 0;

  SHRead(filename.c_str(),
         filename.size(),
         result.cilm_ptr(),
         cilm_dim,
         &degree,
         &skip,
         header,
         header_d,
         error_ptr,
         &exitstatus);

  return result;
}

}
