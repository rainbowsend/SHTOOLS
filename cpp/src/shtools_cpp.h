#pragma once

#include "shtools.h"

namespace shtools {
    
    

// n is number of cosine coefficents or sine coefficents
inline int
n_to_deg(int n)
{
  return std::sqrt(1. / 4. + 2. * n) - 3. / 2.;
}



class ShCoeff{
    public:
        
        int clim_dim() const {return degree+1;};
        
        double get_C(int degree, int order) const{
            size_t c = clim_dim()*degree + order;
            return cilm[c];
        };
        
        double get_S(int degree, int order) const{
            size_t c = std::pow(clim_dim(),2) + clim_dim()*degree + order;
            return cilm[c];
        };
        
        
    private:
        int degree=0;
        // cilm has shape (2,clim_dim,clim_dim)
        std::vector<double> cilm;
        
};

inline std::ostream& operator<<(std::ostream& os, const ShCoeff& shc)
{

    for (int j = 0; j < shc.clim_dim(); ++j) {
        for (int k = 0; k < shc.clim_dim(); ++k) {
            os << std::setw(12) << std::setprecision(3)
                << shc.get_C(j,k) << " ";
        }
        std::cout << std::endl;
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
constexpr OutputIt sh_cindex_to_cilm( const InputIt cindex_first, const InputIt cindex_last, OutputIt cilm_first, int degmax=-1 )
{

  int exitstatus;
  
  int n = std::distance(cindex_first, cindex_last) / 2;
  int cilm_dim = n_to_deg(n)+1;
  
  if(degmax < 0){
      degmax = cilm_dim-1;
  }
  SHCindexToCilm(&*cindex_first, n, &*cilm_first, cilm_dim, degmax, &exitstatus);
  
  return cilm_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt sh_cilm_to_cindex( const InputIt cilm_first, const InputIt cilm_last, OutputIt cindex_first, int degmax=-1 )
{

  int exitstatus;
  
  int n = std::distance(cilm_first, cilm_last);
  int cilm_dim = sqrt(n/2)-1;
  
  if(degmax < 0){
      degmax = cilm_dim-1;
  }
  SHCilmToCindex( &*cilm_first, cilm_dim, &*cindex_first, n, degmax, &exitstatus);
  
  return cindex_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt sh_vector_to_cilm( const InputIt vector_first, const InputIt vector_last, OutputIt cilm_first, int degmax=-1 )
{

  int exitstatus;
  
  int n = std::distance(vector_first, vector_last);
  int cilm_dim = sqrt(n);
  
  if(degmax < 0){
      degmax = cilm_dim-1;
  }
  SHVectorToCilm( &*vector_first, &*cilm_first, cilm_dim, degmax, &exitstatus);
  
  return cilm_first;
}

template<class InputIt, class OutputIt>
constexpr OutputIt sh_cilm_to_vector( const InputIt cilm_first, const InputIt cilm_last, OutputIt vector_first, int degmax=-1 )
{

  int exitstatus;
  
  int n = std::distance(cilm_first, cilm_last);
  int cilm_dim = sqrt(n/2);
  
  
  if(degmax < 0){
      degmax = cilm_dim-1;
  }
  SHCilmToVector( &*cilm_first, cilm_dim, &*vector_first, degmax, &exitstatus);
  
  
  return vector_first;
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
  if (error) {
    error_ptr = &error->at(0);
    error->resize(cilm.size());
  }

  double* header = nullptr;
  int header_d = 0;


  SHRead(filename.c_str(),
         filename.size(),
          &cilm[0],
          cilm_dim,
          &degree,
          &skip,
          header,
          header_d,
          error_ptr,
          &exitstatus);

  return cilm;
}

typedef Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::ColMajor> Cindex;
typedef Eigen::Tensor<double, 3, Eigen::ColMajor> Cilm;

}
