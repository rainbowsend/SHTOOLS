The cpp interface is as consistent as possible to the Fortran interface. That means the cpp function and argument names are equivalent to the Fortran names. Moreover, we use the same order of arguments. However, additional arguments are required for Fortran arrays with assumed shape and character arrays. These additional arguments are inserted behind the corresponding array argument. 

# Arrays
For the cpp interface any data structure that satisfies the following conditions can be used as 'array' argument:

* The elements must be stored in column-major order
* The data types must agree
* The size must be the same (the total number of elements)

For example the Fortran Array

` real, dimension(5,10) :: matrix`

could be represented in cpp by:

* ` std::array<double,50> matrix; `
* ` std::vector<double> matrix(50); `
* ` std::valarray<double> matrix(50); `
* ` double matrix[50]; `
* ` Eigen::Matrix<double,5,10,Eigen::ColMajor>; `
* ...

In the first four cases it is the users responsibility to ensure column-major order.

Fortran arrays are passed by reference, thus in cpp we need to pass the address of the first element of the data structure. Additionally, we need to provide the dimension of the corresponding array. This is done with additional arguments. For example, the subroutine/function `Curve2Mask` has four additional arguments (`dhgrid_d0`,`dhgrid_d1`,`profile_d0`,`profile_d1`) specifying the size of the arrays `dhgrid` and `profile`.

``` fortran
subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP, &
                              extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(out) :: dhgrid(:,:)
            real(dp), intent(in) :: profile(:,:)
            integer, intent(in) :: n, sampling, nprofile, np
            integer, intent(in), optional :: extend
            integer, intent(out), optional :: exitstatus
```

``` cpp
void Curve2Mask(int* dhgrid,
                const int dhgrid_d0,
                const int dhgrid_d1,
                const int n,
                const int sampling,
                const double* profile,
                const int profile_d0,
                const int profile_d1,
                const int nprofile,
                int* NP,
                const int* extend = nullptr,
                int* exitstatus = nullptr);
```

# Optional arguments
Arguments that are optional in the Fortran interface are assigned with a default value in the cpp interface. Passing a null pointer in cpp is equivalent to not specifying the argument in Fortran. Thus, optional arguments are pointers in cpp. An exception are additional arguments for the dimension of optional arrays. Those are passed by value and are set to zero by default.

# Intent keyword
Fortran `intent(in)` arguments are `const` in the cpp interface.
