    subroutine cPlmBar(p,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmBar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmBar(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmBar

    subroutine cPlmBar_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmBar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmBar_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmBar_d1

    subroutine cPlBar(p,lmax,z,exitstatus)  bind(c, name="PlBar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlBar(p,lmax,z,exitstatus=exitstatus)
    end subroutine cPlBar

    subroutine cPlBar_d1(p,dp1,lmax,z,exitstatus)  bind(c, name="PlBar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), dimension(lmax+1),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlBar_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cPlBar_d1

    subroutine cPlmON(p,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmON")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmON(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmON

    subroutine cPlmON_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmON_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmON_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmON_d1

    subroutine cPlON(p,lmax,z,exitstatus)  bind(c, name="PlON")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlON(p,lmax,z,exitstatus=exitstatus)
    end subroutine cPlON

    subroutine cPlON_d1(p,dp1,lmax,z,exitstatus)  bind(c, name="PlON_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), dimension(lmax+1),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlON_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cPlON_d1

    subroutine cPlmSchmidt(p,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmSchmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmSchmidt(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmSchmidt

    subroutine cPlmSchmidt_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="PlmSchmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmSchmidt_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cPlmSchmidt_d1

    subroutine cPlSchmidt(p,lmax,z,exitstatus)  bind(c, name="PlSchmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlSchmidt(p,lmax,z,exitstatus=exitstatus)
    end subroutine cPlSchmidt

    subroutine cPlSchmidt_d1(p,dp1,lmax,z,exitstatus)  bind(c, name="PlSchmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), dimension(lmax+1),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlSchmidt_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cPlSchmidt_d1

    subroutine cPLegendreA(p,lmax,z,csphase,exitstatus)  bind(c, name="PLegendreA")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendreA(p,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine cPLegendreA

    subroutine cPLegendreA_d1(p,dp1,lmax,z,csphase,exitstatus)  bind(c, name="PLegendreA_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: p
        real(kind=c_double), dimension((lmax+1)*(lmax+2)/2),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendreA_d1(p,dp1,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine cPLegendreA_d1

    subroutine cPLegendre(p,lmax,z,exitstatus)  bind(c, name="PLegendre")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendre(p,lmax,z,exitstatus=exitstatus)
    end subroutine cPLegendre

    subroutine cPLegendre_d1(p,dp1,lmax,z,exitstatus)  bind(c, name="PLegendre_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre_d1
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: p
        real(kind=c_double), dimension(lmax+1),intent(out) :: dp1
        real(kind=c_double), value,intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendre_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cPLegendre_d1

    function cPlmIndex(l,m)  bind(c, name="PlmIndex")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmIndex
        implicit none
        integer(kind=c_int) :: cPlmIndex
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: m
        cPlmIndex=PlmIndex(l,m)
    end function cPlmIndex

    subroutine cSHExpandDH(grid,grid_d0,grid_d1,n,cilm,cilm_dim,lmax,norm,sampling&
                               ,csphase,lmax_calc,exitstatus)  bind(c, name="SHExpandDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDH
        implicit none
        integer(kind=c_int), value,intent(in) :: grid_d0
        integer(kind=c_int), value,intent(in) :: grid_d1
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandDH(grid,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                            ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cSHExpandDH

    subroutine cMakeGridDH(griddh,griddh_d0,griddh_d1,n,cilm,cilm_dim,lmax,norm,sampling&
                                 ,csphase,lmax_calc,extend,exitstatus)  bind(c, name="MakeGridDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: griddh_d0
        integer(kind=c_int), value,intent(in) :: griddh_d1
        real(kind=c_double), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridDH(griddh,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                              ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGridDH

    subroutine cSHExpandDHC(grid,grid_d0,grid_d1,n,cilm,cilm_dim,lmax,norm,sampling&
                                ,csphase,lmax_calc,exitstatus)  bind(c, name="SHExpandDHC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDHC
        implicit none
        integer(kind=c_int), value,intent(in) :: grid_d0
        integer(kind=c_int), value,intent(in) :: grid_d1
        complex(kind=c_double_complex), dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandDHC(grid,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                             ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cSHExpandDHC

    subroutine cMakeGridDHC(griddh,griddh_d0,griddh_d1,n,cilm,cilm_dim,lmax,norm,sampling&
                                  ,csphase,lmax_calc,extend,exitstatus)  bind(c, name="MakeGridDHC")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDHC
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: griddh_d0
        integer(kind=c_int), value,intent(in) :: griddh_d1
        complex(kind=c_double_complex), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridDHC(griddh,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                               ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGridDHC

    subroutine cSHGLQ(lmax,zero,w,plx,norm,csphase,cnorm,exitstatus)  bind(c, name="SHGLQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHGLQ
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(out) :: zero
        real(kind=c_double), dimension(lmax+1),intent(out) :: w
        real(kind=c_double), optional,dimension(lmax+1,(lmax+1)*(lmax+2)/2),intent(out) :: plx
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHGLQ(lmax,zero,w,plx=plx,norm=norm,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cSHGLQ

    subroutine cSHExpandGLQ(cilm,cilm_dim,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc&
                                ,exitstatus)  bind(c, name="SHExpandGLQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQ
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(in) :: w
        real(kind=c_double), dimension(lmax+1,2*lmax+1),intent(in) :: gridglq
        real(kind=c_double), optional,dimension(lmax+1,(lmax+1)*(lmax+2)/2),intent(in) :: plx
        real(kind=c_double), optional,dimension(lmax+1),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandGLQ(cilm,lmax,gridglq,w,plx=plx,zero=zero,norm=norm,csphase=csphase&
                             ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cSHExpandGLQ

    subroutine cMakeGridGLQ(gridglq,gridglq_d0,gridglq_d1,cilm,cilm_dim,lmax,plx,zero&
                                   ,norm,csphase,lmax_calc,extend,exitstatus)  bind(c&
                                   , name="MakeGridGLQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQ
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), optional,dimension(lmax+1,(lmax+1)*(lmax+2)/2),intent(in) :: plx
        real(kind=c_double), optional,dimension(lmax+1),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: gridglq_d0
        integer(kind=c_int), value,intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridGLQ(gridglq,cilm,lmax,plx=plx,zero=zero,norm=norm,csphase=csphase&
                                ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGridGLQ

    subroutine cSHExpandGLQC(cilm,cilm_dim,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc&
                                 ,exitstatus)  bind(c, name="SHExpandGLQC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQC
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1),intent(in) :: w
        complex(kind=c_double_complex), dimension(lmax+1,2*lmax+1),intent(in) :: gridglq
        real(kind=c_double), optional,dimension(lmax+1,(lmax+1)*(lmax+2)/2),intent(in) :: plx
        real(kind=c_double), optional,dimension(lmax+1),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandGLQC(cilm,lmax,gridglq,w,plx=plx,zero=zero,norm=norm,csphase=csphase&
                              ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cSHExpandGLQC

    subroutine cMakeGridGLQC(gridglq,gridglq_d0,gridglq_d1,cilm,cilm_dim,lmax,plx&
                                    ,zero,norm,csphase,lmax_calc,extend,exitstatus)  bind(c&
                                    , name="MakeGridGLQC")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQC
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), optional,dimension(lmax+1,(lmax+1)*(lmax+2)/2),intent(in) :: plx
        real(kind=c_double), optional,dimension(lmax+1),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: gridglq_d0
        integer(kind=c_int), value,intent(in) :: gridglq_d1
        complex(kind=c_double_complex), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridGLQC(gridglq,cilm,lmax,plx=plx,zero=zero,norm=norm,csphase=csphase&
                                 ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGridGLQC

    subroutine cGLQGridCoord(latglq,latglq_d0,longlq,longlq_d0,lmax,nlat,nlong,extend&
                                   ,exitstatus)  bind(c, name="GLQGridCoord")
        use, intrinsic :: iso_c_binding
        use shtools, only: GLQGridCoord
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        integer(kind=c_int), value,intent(in) :: latglq_d0
        real(kind=c_double), dimension(latglq_d0),intent(out) :: latglq
        integer(kind=c_int), value,intent(in) :: longlq_d0
        real(kind=c_double), dimension(longlq_d0),intent(out) :: longlq
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong,extend=extend,exitstatus=exitstatus)
    end subroutine cGLQGridCoord

    subroutine cSHExpandLSQ(cilm,cilm_dim,d,lat,lon,nmax,lmax,norm,chi2,csphase,weights&
                                ,exitstatus)  bind(c, name="SHExpandLSQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandLSQ
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        real(kind=c_double), dimension(nmax),intent(in) :: d
        real(kind=c_double), dimension(nmax),intent(in) :: lat
        real(kind=c_double), dimension(nmax),intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        real(kind=c_double), optional,intent(out) :: chi2
        real(kind=c_double), optional,dimension(nmax),intent(in) :: weights
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm=norm,chi2=chi2,csphase=csphase&
                             ,weights=weights,exitstatus=exitstatus)
    end subroutine cSHExpandLSQ

    subroutine cMakeGrid2d(grid,grid_d0,grid_d1,cilm,cilm_dim,lmax,interval,nlat,nlong&
                               ,norm,csphase,f,a,north,south,east,west,dealloc,exitstatus)  bind(c&
                               , name="MakeGrid2d")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGrid2d
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: interval
        integer(kind=c_int), value,intent(in) :: grid_d0
        integer(kind=c_int), value,intent(in) :: grid_d1
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(out) :: grid
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        real(kind=c_double), optional,intent(in) :: f
        real(kind=c_double), optional,intent(in) :: a
        real(kind=c_double), optional,intent(in) :: north
        real(kind=c_double), optional,intent(in) :: south
        real(kind=c_double), optional,intent(in) :: east
        real(kind=c_double), optional,intent(in) :: west
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm,csphase=csphase&
                            ,f=f,a=a,north=north,south=south,east=east,west=west,dealloc=dealloc&
                            ,exitstatus=exitstatus)
    end subroutine cMakeGrid2d

    function cMakeGridPoint(cilm,cilm_dim,lmax,lat,lon,norm,csphase,dealloc)  bind(c&
                                , name="MakeGridPoint")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPoint
        implicit none
        real(kind=c_double) :: cMakeGridPoint
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: lat
        real(kind=c_double), value,intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        cMakeGridPoint=MakeGridPoint(cilm,lmax,lat,lon,norm=norm,csphase=csphase,dealloc=dealloc)
    end function cMakeGridPoint

    function cMakeGridPointC(cilm,cilm_dim,lmax,lat,lon,norm,csphase,dealloc)  bind(c&
                                 , name="MakeGridPointC")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPointC
        implicit none
        complex(kind=c_double_complex) :: cMakeGridPointC
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: lat
        real(kind=c_double), value,intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        cMakeGridPointC=MakeGridPointC(cilm,lmax,lat,lon,norm=norm,csphase=csphase&
                                           ,dealloc=dealloc)
    end function cMakeGridPointC

    subroutine cSHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp,norm,csphase,exitstatus)  bind(c&
                                , name="SHMultiply")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiply
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax2
        integer(kind=c_int), value,intent(in) :: lmax1
        real(kind=c_double), dimension(2,lmax1+lmax2+1,lmax1+lmax2+1),intent(out) :: shout
        real(kind=c_double), dimension(2,lmax1+1,lmax1+1),intent(in) :: sh1
        real(kind=c_double), dimension(2,lmax2+1,lmax2+1),intent(in) :: sh2
        integer(kind=c_int), optional,intent(in) :: precomp
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp=precomp,norm=norm,csphase=csphase&
                             ,exitstatus=exitstatus)
    end subroutine cSHMultiply

    subroutine cSHRead(filename,filename_d1,cilm,cilm_dim,lmax,skip,header,header_d0&
                               ,error,exitstatus)  bind(c, name="SHRead")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        integer(kind=c_int), value,intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: header_d0
        real(kind=c_double), optional,dimension(header_d0),intent(out) :: header
        real(kind=c_double), optional,dimension(2,cilm_dim,cilm_dim),intent(out) :: error
        integer(kind=c_int), optional,intent(in) :: skip
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRead(filename,cilm,lmax,skip=skip,header=header,error=error,exitstatus=exitstatus)
    end subroutine cSHRead

    subroutine cSHRead2(filename,filename_d1,cilm,cilm_dim,lmax,gm,r0_pot,error,dot&
                                ,dot_d0,dot_d1,dot_d2,doystart,doyend,epoch,exitstatus)  bind(c&
                                , name="SHRead2")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead2
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        integer(kind=c_int), value,intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), intent(out) :: gm
        real(kind=c_double), intent(out) :: r0_pot
        real(kind=c_double), optional,dimension(2,cilm_dim,cilm_dim),intent(out) :: error
        integer(kind=c_int), value,intent(in) :: dot_d0
        integer(kind=c_int), value,intent(in) :: dot_d1
        integer(kind=c_int), value,intent(in) :: dot_d2
        real(kind=c_double), optional,dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real(kind=c_double), optional,intent(out) :: doystart
        real(kind=c_double), optional,intent(out) :: doyend
        real(kind=c_double), optional,intent(out) :: epoch
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRead2(filename,cilm,lmax,gm,r0_pot,error=error,dot=dot,doystart=doystart&
                             ,doyend=doyend,epoch=epoch,exitstatus=exitstatus)
    end subroutine cSHRead2

    subroutine cSHReadJPL(filename,filename_d1,cilm,cilm_dim,lmax,error,gm,formatstring&
                                  ,exitstatus)  bind(c, name="SHReadJPL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReadJPL
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        integer(kind=c_int), value,intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), optional,dimension(2,cilm_dim,cilm_dim),intent(out) :: error
        real(kind=c_double), optional,dimension(2),intent(out) :: gm
        character(kind=c_char), optional,intent(in) :: formatstring
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReadJPL(filename,cilm,lmax,error=error,gm=gm,formatstring=formatstring&
                               ,exitstatus=exitstatus)
    end subroutine cSHReadJPL

    subroutine cSHCilmToVector(cilm,cilm_dim,vector,lmax,exitstatus)  bind(c, name="SHCilmToVector")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToVector
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), dimension((lmax+1)**2),intent(out) :: vector
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCilmToVector(cilm,vector,lmax,exitstatus=exitstatus)
    end subroutine cSHCilmToVector

    subroutine cSHVectorToCilm(vector,cilm,cilm_dim,lmax,exitstatus)  bind(c, name="SHVectorToCilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHVectorToCilm
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), dimension((lmax+1)**2),intent(in) :: vector
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHVectorToCilm(vector,cilm,lmax,exitstatus=exitstatus)
    end subroutine cSHVectorToCilm

    subroutine cSHCilmToCindex(cilm,cilm_dim,cindex,cindex_d0,cindex_d1,degmax,exitstatus)  bind(c&
                                   , name="SHCilmToCindex")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToCindex
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: cindex_d0
        integer(kind=c_int), value,intent(in) :: cindex_d1
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCilmToCindex(cilm,cindex,degmax=degmax,exitstatus=exitstatus)
    end subroutine cSHCilmToCindex

    subroutine cSHCindexToCilm(cindex,cindex_d0,cindex_d1,cilm,cilm_dim,degmax,exitstatus)  bind(c&
                                     , name="SHCindexToCilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCindexToCilm
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: cindex_d0
        integer(kind=c_int), value,intent(in) :: cindex_d1
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCindexToCilm(cindex,cilm,degmax=degmax,exitstatus=exitstatus)
    end subroutine cSHCindexToCilm

    subroutine cSHrtoc(rcilm,cilm_dim,ccilm,degmax,convention,switchcs,exitstatus)  bind(c&
                            , name="SHrtoc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHrtoc
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: rcilm
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: ccilm
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(in) :: convention
        integer(kind=c_int), optional,intent(in) :: switchcs
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHrtoc(rcilm,ccilm,degmax=degmax,convention=convention,switchcs=switchcs&
                         ,exitstatus=exitstatus)
    end subroutine cSHrtoc

    subroutine cSHctor(ccilm,cilm_dim,rcilm,degmax,convention,switchcs,exitstatus)  bind(c&
                            , name="SHctor")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHctor
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: ccilm
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: rcilm
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(in) :: convention
        integer(kind=c_int), optional,intent(in) :: switchcs
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHctor(ccilm,rcilm,degmax=degmax,convention=convention,switchcs=switchcs&
                         ,exitstatus=exitstatus)
    end subroutine cSHctor

    subroutine cdjpi2(dj,lmax,exitstatus)  bind(c, name="djpi2")
        use, intrinsic :: iso_c_binding
        use shtools, only: djpi2
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(lmax+1,lmax+1,lmax+1),intent(out) :: dj
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call djpi2(dj,lmax,exitstatus=exitstatus)
    end subroutine cdjpi2

    subroutine cSHRotateCoef(x,cof,cof_d0,cof_d1,rcof,dj,lmax,exitstatus)  bind(c&
                              , name="SHRotateCoef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateCoef
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cof_d0
        integer(kind=c_int), value,intent(in) :: cof_d1
        real(kind=c_double), dimension(cof_d0,cof_d1),intent(in) :: cof
        real(kind=c_double), dimension(lmax+1,lmax+1,lmax+1),intent(in) :: dj
        real(kind=c_double), dimension(3),intent(in) :: x
        real(kind=c_double), dimension(2,(lmax+1)*(lmax+2)/2),intent(out) :: rcof
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateCoef(x,cof,rcof,dj,lmax,exitstatus=exitstatus)
    end subroutine cSHRotateCoef

    subroutine cSHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus)  bind(c, name="SHRotateRealCoef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateRealCoef
        implicit none
        integer(kind=c_int), value,intent(in) :: lmax
        real(kind=c_double), dimension(2,lmax+1,lmax+1),intent(in) :: cilm
        real(kind=c_double), dimension(3),intent(in) :: x
        real(kind=c_double), dimension(lmax+1,lmax+1,lmax+1),intent(in) :: dj
        real(kind=c_double), dimension(2,lmax+1,lmax+1),intent(out) :: cilmrot
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus=exitstatus)
    end subroutine cSHRotateRealCoef

    function cSHPowerL(cilm,cilm_dim,lmax)  bind(c, name="SHPowerL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerL
        implicit none
        real(kind=c_double) :: cSHPowerL
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        cSHPowerL=SHPowerL(cilm,lmax)
    end function cSHPowerL

    function cSHPowerDensityL(cilm,cilm_dim,lmax)  bind(c, name="SHPowerDensityL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityL
        implicit none
        real(kind=c_double) :: cSHPowerDensityL
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        cSHPowerDensityL=SHPowerDensityL(cilm,lmax)
    end function cSHPowerDensityL

    function cSHCrossPowerL(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax)  bind(c, name="SHCrossPowerL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerL
        implicit none
        real(kind=c_double) :: cSHCrossPowerL
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        real(kind=c_double), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        real(kind=c_double), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        cSHCrossPowerL=SHCrossPowerL(cilm1,cilm2,lmax)
    end function cSHCrossPowerL

    function cSHCrossPowerDensityL(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax)  bind(c, name="SHCrossPowerDensityL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityL
        implicit none
        real(kind=c_double) :: cSHCrossPowerDensityL
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        real(kind=c_double), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        real(kind=c_double), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        cSHCrossPowerDensityL=SHCrossPowerDensityL(cilm1,cilm2,lmax)
    end function cSHCrossPowerDensityL

    subroutine cSHPowerSpectrum(cilm,cilm_dim,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                    , name="SHPowerSpectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrum
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrum(cilm,lmax,spectra,exitstatus=exitstatus)
    end subroutine cSHPowerSpectrum

    subroutine cSHPowerSpectrumDensity(cilm,cilm_dim,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                           , name="SHPowerSpectrumDensity")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensity
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumDensity(cilm,lmax,spectra,exitstatus=exitstatus)
    end subroutine cSHPowerSpectrumDensity

    subroutine cSHCrossPowerSpectrum(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax,cspectra&
                                          ,cspectra_d0,exitstatus)  bind(c, name="SHCrossPowerSpectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrum
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        real(kind=c_double), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        real(kind=c_double), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cspectra_d0
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrum(cilm1,cilm2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cSHCrossPowerSpectrum

    subroutine cSHCrossPowerSpectrumDensity(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax,cspectra&
                                                 ,cspectra_d0,exitstatus)  bind(c&
                                                 , name="SHCrossPowerSpectrumDensity")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensity
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        real(kind=c_double), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        real(kind=c_double), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cspectra_d0
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensity(cilm1,cilm2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cSHCrossPowerSpectrumDensity

    subroutine cSHAdmitCorr(G,G_d0,G_d1,G_d2,T,T_d0,T_d1,T_d2,lmax,admit,admit_d0&
                             ,corr,corr_d0,admit_error,admit_error_d0,exitstatus)  bind(c&
                             , name="SHAdmitCorr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHAdmitCorr
        implicit none
        integer(kind=c_int), value,intent(in) :: G_d0
        integer(kind=c_int), value,intent(in) :: G_d1
        integer(kind=c_int), value,intent(in) :: G_d2
        real(kind=c_double), dimension(G_d0,G_d1,G_d2),intent(in) :: G
        integer(kind=c_int), value,intent(in) :: T_d0
        integer(kind=c_int), value,intent(in) :: T_d1
        integer(kind=c_int), value,intent(in) :: T_d2
        real(kind=c_double), dimension(T_d0,T_d1,T_d2),intent(in) :: T
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), value,intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), value,intent(in) :: admit_error_d0
        real(kind=c_double), optional,dimension(admit_error_d0),intent(out) :: admit_error
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error=admit_error,exitstatus=exitstatus)
    end subroutine cSHAdmitCorr

    function cSHConfidence(l_conf,r)  bind(c, name="SHConfidence")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHConfidence
        implicit none
        real(kind=c_double) :: cSHConfidence
        real(kind=c_double), value,intent(in) :: r
        integer(kind=c_int), value,intent(in) :: l_conf
        cSHConfidence=SHConfidence(l_conf,r)
    end function cSHConfidence

    function cSHPowerLC(cilm,cilm_dim,lmax)  bind(c, name="SHPowerLC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerLC
        implicit none
        real(kind=c_double) :: cSHPowerLC
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        cSHPowerLC=SHPowerLC(cilm,lmax)
    end function cSHPowerLC

    function cSHPowerDensityLC(cilm,cilm_dim,lmax)  bind(c, name="SHPowerDensityLC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityLC
        implicit none
        real(kind=c_double) :: cSHPowerDensityLC
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        cSHPowerDensityLC=SHPowerDensityLC(cilm,lmax)
    end function cSHPowerDensityLC

    function cSHCrossPowerLC(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax)  bind(c, name="SHCrossPowerLC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerLC
        implicit none
        complex(kind=c_double_complex) :: cSHCrossPowerLC
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        complex(kind=c_double_complex), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        complex(kind=c_double_complex), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        cSHCrossPowerLC=SHCrossPowerLC(cilm1,cilm2,lmax)
    end function cSHCrossPowerLC

    function cSHCrossPowerDensityLC(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax)  bind(c&
                                         , name="SHCrossPowerDensityLC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityLC
        implicit none
        complex(kind=c_double_complex) :: cSHCrossPowerDensityLC
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        complex(kind=c_double_complex), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        complex(kind=c_double_complex), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        cSHCrossPowerDensityLC=SHCrossPowerDensityLC(cilm1,cilm2,lmax)
    end function cSHCrossPowerDensityLC

    subroutine cSHPowerSpectrumC(cilm,cilm_dim,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                     , name="SHPowerSpectrumC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumC
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumC(cilm,lmax,spectra,exitstatus=exitstatus)
    end subroutine cSHPowerSpectrumC

    subroutine cSHPowerSpectrumDensityC(cilm,cilm_dim,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                            , name="SHPowerSpectrumDensityC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        complex(kind=c_double_complex), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumDensityC(cilm,lmax,spectra,exitstatus=exitstatus)
    end subroutine cSHPowerSpectrumDensityC

    subroutine cSHCrossPowerSpectrumC(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax,cspectra&
                                           ,cspectra_d0,exitstatus)  bind(c, name="SHCrossPowerSpectrumC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumC
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        complex(kind=c_double_complex), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        complex(kind=c_double_complex), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cspectra_d0
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumC(cilm1,cilm2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cSHCrossPowerSpectrumC

    subroutine cSHCrossPowerSpectrumDensityC(cilm1,cilm1_dim,cilm2,cilm2_dim,lmax&
                                                  ,cspectra,cspectra_d0,exitstatus)  bind(c&
                                                  , name="SHCrossPowerSpectrumDensityC")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        complex(kind=c_double_complex), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        complex(kind=c_double_complex), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: cspectra_d0
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensityC(cilm1,cilm2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cSHCrossPowerSpectrumDensityC

    subroutine cSHMultiTaperSE(mtse,mtse_d0,sd,sd_d0,cilm,cilm_dim,lmax,tapers,tapers_d0&
                                   ,tapers_d1,taper_order,taper_order_d0,lmaxt,k,alpha&
                                   ,alpha_d0,lat,lon,taper_wt,taper_wt_d0,norm,csphase&
                                   ,exitstatus)  bind(c, name="SHMultiTaperSE")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperSE
        implicit none
        integer(kind=c_int), value,intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), value,intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: lmaxt
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: alpha_d0
        real(kind=c_double), optional,dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), optional,intent(in) :: lat
        real(kind=c_double), optional,intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperSE(mtse,sd,cilm,lmax,tapers,taper_order,lmaxt,k,alpha=alpha&
                                ,lat=lat,lon=lon,taper_wt=taper_wt,norm=norm,csphase=csphase&
                                ,exitstatus=exitstatus)
    end subroutine cSHMultiTaperSE

    subroutine cSHMultiTaperCSE(mtse,mtse_d0,sd,sd_d0,cilm1,cilm1_dim,lmax1,cilm2&
                                    ,cilm2_dim,lmax2,tapers,tapers_d0,tapers_d1,taper_order&
                                    ,taper_order_d0,lmaxt,k,alpha,alpha_d0,lat,lon&
                                    ,taper_wt,taper_wt_d0,norm,csphase,exitstatus)  bind(c&
                                    , name="SHMultiTaperCSE")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperCSE
        implicit none
        integer(kind=c_int), value,intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), value,intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), value,intent(in) :: cilm1_dim
        real(kind=c_double), dimension(2,cilm1_dim,cilm1_dim),intent(in) :: cilm1
        integer(kind=c_int), value,intent(in) :: cilm2_dim
        real(kind=c_double), dimension(2,cilm2_dim,cilm2_dim),intent(in) :: cilm2
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: lmax1
        integer(kind=c_int), value,intent(in) :: lmax2
        integer(kind=c_int), value,intent(in) :: lmaxt
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: alpha_d0
        real(kind=c_double), optional,dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), optional,intent(in) :: lat
        real(kind=c_double), optional,intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperCSE(mtse,sd,cilm1,lmax1,cilm2,lmax2,tapers,taper_order,lmaxt&
                                 ,k,alpha=alpha,lat=lat,lon=lon,taper_wt=taper_wt&
                                 ,norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine cSHMultiTaperCSE

    subroutine cSHLocalizedAdmitCorr(tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                                           ,lwin,lat,lon,g,g_d0,g_d1,g_d2,t,t_d0,t_d1&
                                           ,t_d2,lmax,admit,admit_d0,corr,corr_d0&
                                           ,k,admit_error,admit_error_d0,corr_error&
                                           ,corr_error_d0,taper_wt,taper_wt_d0,mtdef&
                                           ,k1linsig,exitstatus)  bind(c, name="SHLocalizedAdmitCorr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHLocalizedAdmitCorr
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), value,intent(in) :: lat
        real(kind=c_double), value,intent(in) :: lon
        integer(kind=c_int), value,intent(in) :: g_d0
        integer(kind=c_int), value,intent(in) :: g_d1
        integer(kind=c_int), value,intent(in) :: g_d2
        real(kind=c_double), dimension(g_d0,g_d1,g_d2),intent(in) :: g
        integer(kind=c_int), value,intent(in) :: t_d0
        integer(kind=c_int), value,intent(in) :: t_d1
        integer(kind=c_int), value,intent(in) :: t_d2
        real(kind=c_double), dimension(t_d0,t_d1,t_d2),intent(in) :: t
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), value,intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), value,intent(in) :: admit_error_d0
        real(kind=c_double), optional,dimension(admit_error_d0),intent(out) :: admit_error
        integer(kind=c_int), value,intent(in) :: corr_error_d0
        real(kind=c_double), optional,dimension(corr_error_d0),intent(out) :: corr_error
        integer(kind=c_int), optional,intent(in) :: mtdef
        integer(kind=c_int), optional,intent(in) :: k1linsig
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,lmax,admit,corr&
                                        ,k,admit_error=admit_error,corr_error=corr_error&
                                        ,taper_wt=taper_wt,mtdef=mtdef,k1linsig=k1linsig&
                                        ,exitstatus=exitstatus)
    end subroutine cSHLocalizedAdmitCorr

    subroutine cSHReturnTapers(theta0,lmax,tapers,tapers_d0,tapers_d1,eigenvalues&
                                     ,eigenvalues_d0,taper_order,taper_order_d0,degrees&
                                     ,degrees_d0,exitstatus)  bind(c, name="SHReturnTapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapers
        implicit none
        real(kind=c_double), value,intent(in) :: theta0
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), value,intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(out) :: taper_order
        integer(kind=c_int), value,intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,degrees=degrees&
                                  ,exitstatus=exitstatus)
    end subroutine cSHReturnTapers

    subroutine cSHReturnTapersM(theta0,lmax,m,tapers,tapers_d0,tapers_d1,eigenvalues&
                                      ,eigenvalues_d0,shannon,degrees,degrees_d0,ntapers&
                                      ,exitstatus)  bind(c, name="SHReturnTapersM")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersM
        implicit none
        real(kind=c_double), value,intent(in) :: theta0
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: m
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), value,intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        real(kind=c_double), optional,intent(out) :: shannon
        integer(kind=c_int), value,intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: ntapers
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,shannon=shannon,degrees=degrees&
                                   ,ntapers=ntapers,exitstatus=exitstatus)
    end subroutine cSHReturnTapersM

    subroutine cComputeDm(dllm,dllm_d0,dllm_d1,lmax,m,theta0,degrees,degrees_d0,exitstatus)  bind(c&
                              , name="ComputeDm")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDm
        implicit none
        integer(kind=c_int), value,intent(in) :: dllm_d0
        integer(kind=c_int), value,intent(in) :: dllm_d1
        real(kind=c_double), dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        real(kind=c_double), value,intent(in) :: theta0
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: m
        integer(kind=c_int), value,intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDm(dllm,lmax,m,theta0,degrees=degrees,exitstatus=exitstatus)
    end subroutine cComputeDm

    subroutine cComputeDG82(dG82,dG82_d0,dG82_d1,lmax,m,theta0,exitstatus)  bind(c&
                                , name="ComputeDG82")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDG82
        implicit none
        integer(kind=c_int), value,intent(in) :: dG82_d0
        integer(kind=c_int), value,intent(in) :: dG82_d1
        real(kind=c_double), dimension(dG82_d0,dG82_d1),intent(out) :: dG82
        real(kind=c_double), value,intent(in) :: theta0
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: m
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDG82(dG82,lmax,m,theta0,exitstatus=exitstatus)
    end subroutine cComputeDG82

    function cSHFindLWin(theta0,m,alpha,taper_number)  bind(c, name="SHFindLWin")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHFindLWin
        implicit none
        integer(kind=c_int) :: cSHFindLWin
        real(kind=c_double), value,intent(in) :: theta0
        real(kind=c_double), value,intent(in) :: alpha
        integer(kind=c_int), value,intent(in) :: m
        integer(kind=c_int), optional,intent(in) :: taper_number
        cSHFindLWin=SHFindLWin(theta0,m,alpha,taper_number=taper_number)
    end function cSHFindLWin

    subroutine cSHBiasK(tapers,tapers_d0,tapers_d1,lwin,k,incspectra,incspectra_d0&
                              ,ldata,outcspectra,outcspectra_d0,taper_wt,taper_wt_d0&
                              ,save_cg,exitstatus)  bind(c, name="SHBiasK")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasK
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), value,intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: ldata
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt=taper_wt&
                           ,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cSHBiasK

    subroutine cSHMTCouplingMatrix(Mmt,Mmt_d0,Mmt_d1,lmax,tapers_power,tapers_power_d0&
                                      ,tapers_power_d1,lwin,k,taper_wt,taper_wt_d0&
                                      ,exitstatus)  bind(c, name="SHMTCouplingMatrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTCouplingMatrix
        implicit none
        integer(kind=c_int), value,intent(in) :: Mmt_d0
        integer(kind=c_int), value,intent(in) :: Mmt_d1
        real(kind=c_double), dimension(Mmt_d0,Mmt_d1),intent(out) :: Mmt
        integer(kind=c_int), value,intent(in) :: tapers_power_d0
        integer(kind=c_int), value,intent(in) :: tapers_power_d1
        real(kind=c_double), dimension(tapers_power_d0,tapers_power_d1),intent(in) :: tapers_power
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,taper_wt=taper_wt,exitstatus=exitstatus)
    end subroutine cSHMTCouplingMatrix

    subroutine cSHBiasAdmitCorr(sgt,sgt_d0,sgg,sgg_d0,stt,stt_d0,lmax,tapers,tapers_d0&
                                   ,tapers_d1,lwin,k,admit,admit_d0,corr,corr_d0,mtdef&
                                   ,taper_wt,taper_wt_d0,exitstatus)  bind(c, name="SHBiasAdmitCorr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasAdmitCorr
        implicit none
        integer(kind=c_int), value,intent(in) :: sgt_d0
        real(kind=c_double), dimension(sgt_d0),intent(in) :: sgt
        integer(kind=c_int), value,intent(in) :: sgg_d0
        real(kind=c_double), dimension(sgg_d0),intent(in) :: sgg
        integer(kind=c_int), value,intent(in) :: stt_d0
        real(kind=c_double), dimension(stt_d0),intent(in) :: stt
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), value,intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), optional,intent(in) :: mtdef
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,mtdef=mtdef&
                                ,taper_wt=taper_wt,exitstatus=exitstatus)
    end subroutine cSHBiasAdmitCorr

    subroutine cSHMTDebias(mtdebias,mtdebias_d0,mtdebias_d1,mtspectra,mtspectra_d0&
                                   ,mtspectra_d1,lmax,tapers,tapers_d0,tapers_d1,lwin&
                                   ,k,nl,lmid,lmid_d0,n,taper_wt,taper_wt_d0,exitstatus)  bind(c&
                                   , name="SHMTDebias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTDebias
        implicit none
        integer(kind=c_int), value,intent(in) :: mtdebias_d0
        integer(kind=c_int), value,intent(in) :: mtdebias_d1
        real(kind=c_double), dimension(mtdebias_d0,mtdebias_d1),intent(out) :: mtdebias
        integer(kind=c_int), value,intent(in) :: lmid_d0
        real(kind=c_double), dimension(lmid_d0),intent(out) :: lmid
        integer(kind=c_int), value,intent(in) :: mtspectra_d0
        integer(kind=c_int), value,intent(in) :: mtspectra_d1
        real(kind=c_double), dimension(mtspectra_d0,mtspectra_d1),intent(in) :: mtspectra
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: nl
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,taper_wt=taper_wt&
                                ,exitstatus=exitstatus)
    end subroutine cSHMTDebias

    subroutine cSHMTVarOpt(l,tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                            ,lwin,kmax,Sff,Sff_d0,var_opt,var_opt_d0,var_unit,var_unit_d0&
                            ,weight_opt,weight_opt_d0,weight_opt_d1,unweighted_covar&
                            ,unweighted_covar_d0,unweighted_covar_d1,nocross,exitstatus)  bind(c&
                            , name="SHMTVarOpt")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVarOpt
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        integer(kind=c_int), value,intent(in) :: var_opt_d0
        real(kind=c_double), dimension(var_opt_d0),intent(out) :: var_opt
        integer(kind=c_int), value,intent(in) :: var_unit_d0
        real(kind=c_double), dimension(var_unit_d0),intent(out) :: var_unit
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: kmax
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: weight_opt_d0
        integer(kind=c_int), value,intent(in) :: weight_opt_d1
        real(kind=c_double), optional,dimension(weight_opt_d0,weight_opt_d1),intent(out) :: weight_opt
        integer(kind=c_int), value,intent(in) :: unweighted_covar_d0
        integer(kind=c_int), value,intent(in) :: unweighted_covar_d1
        real(kind=c_double), optional,dimension(unweighted_covar_d0,unweighted_covar_d1)&
                           ,intent(out) :: unweighted_covar
        integer(kind=c_int), optional,intent(in) :: nocross
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit,weight_opt=weight_opt&
                         ,unweighted_covar=unweighted_covar,nocross=nocross,exitstatus=exitstatus)
    end subroutine cSHMTVarOpt

    subroutine cSHMTVar(l,tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0,lwin&
                         ,kmax,Sff,Sff_d0,variance,taper_wt,taper_wt_d0,unweighted_covar&
                         ,unweighted_covar_d0,unweighted_covar_d1,nocross,exitstatus)  bind(c&
                         , name="SHMTVar")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVar
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: kmax
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), value,intent(in) :: unweighted_covar_d0
        integer(kind=c_int), value,intent(in) :: unweighted_covar_d1
        real(kind=c_double), optional,dimension(unweighted_covar_d0,unweighted_covar_d1)&
                           ,intent(out) :: unweighted_covar
        integer(kind=c_int), optional,intent(in) :: nocross
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,taper_wt=taper_wt&
                      ,unweighted_covar=unweighted_covar,nocross=nocross,exitstatus=exitstatus)
    end subroutine cSHMTVar

    function cSHSjkPG(incspectra,incspectra_d0,l,m,mprime,hj_real,hj_real_d0,hk_real&
                                ,hk_real_d0,mj,mk,lwin,hkcc)  bind(c, name="SHSjkPG")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSjkPG
        implicit none
        complex(kind=c_double_complex) :: cSHSjkPG
        integer(kind=c_int), value,intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), value,intent(in) :: hj_real_d0
        real(kind=c_double), dimension(hj_real_d0),intent(in) :: hj_real
        integer(kind=c_int), value,intent(in) :: hk_real_d0
        real(kind=c_double), dimension(hk_real_d0),intent(in) :: hk_real
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: m
        integer(kind=c_int), value,intent(in) :: mprime
        integer(kind=c_int), value,intent(in) :: mj
        integer(kind=c_int), value,intent(in) :: mk
        integer(kind=c_int), value,intent(in) :: hkcc
        cSHSjkPG=SHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,hkcc)
    end function cSHSjkPG

    subroutine cSHReturnTapersMap(tapers,tapers_d0,tapers_d1,eigenvalues,eigenvalues_d0&
                                        ,dh_mask,dh_mask_d0,dh_mask_d1,n_dh,lmax,sampling&
                                        ,ntapers,degrees,degrees_d0,exitstatus)  bind(c&
                                        , name="SHReturnTapersMap")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersMap
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), value,intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), value,intent(in) :: dh_mask_d0
        integer(kind=c_int), value,intent(in) :: dh_mask_d1
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), value,intent(in) :: n_dh
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: ntapers
        integer(kind=c_int), value,intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax,sampling,ntapers=ntapers&
                                     ,degrees=degrees,exitstatus=exitstatus)
    end subroutine cSHReturnTapersMap

    subroutine cSHBiasKMask(tapers,tapers_d0,tapers_d1,lwin,k,incspectra,incspectra_d0&
                                  ,ldata,outcspectra,outcspectra_d0,taper_wt,taper_wt_d0&
                                  ,save_cg,exitstatus)  bind(c, name="SHBiasKMask")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasKMask
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), value,intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: ldata
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt=taper_wt&
                               ,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cSHBiasKMask

    subroutine cSHMultiTaperMaskSE(mtse,mtse_d0,sd,sd_d0,sh,sh_d0,sh_d1,sh_d2,lmax&
                                       ,tapers,tapers_d0,tapers_d1,lmaxt,k,taper_wt&
                                       ,taper_wt_d0,norm,csphase,exitstatus)  bind(c&
                                       , name="SHMultiTaperMaskSE")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskSE
        implicit none
        integer(kind=c_int), value,intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), value,intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), value,intent(in) :: sh_d0
        integer(kind=c_int), value,intent(in) :: sh_d1
        integer(kind=c_int), value,intent(in) :: sh_d2
        real(kind=c_double), dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: lmaxt
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,taper_wt=taper_wt,norm=norm&
                                    ,csphase=csphase,exitstatus=exitstatus)
    end subroutine cSHMultiTaperMaskSE

    subroutine cSHMultiTaperMaskCSE(mtse,mtse_d0,sd,sd_d0,sh1,sh1_d0,sh1_d1,sh1_d2&
                                        ,lmax1,sh2,sh2_d0,sh2_d1,sh2_d2,lmax2,tapers&
                                        ,tapers_d0,tapers_d1,lmaxt,k,taper_wt,taper_wt_d0&
                                        ,norm,csphase,exitstatus)  bind(c, name="SHMultiTaperMaskCSE")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskCSE
        implicit none
        integer(kind=c_int), value,intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), value,intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), value,intent(in) :: sh1_d0
        integer(kind=c_int), value,intent(in) :: sh1_d1
        integer(kind=c_int), value,intent(in) :: sh1_d2
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(kind=c_int), value,intent(in) :: sh2_d0
        integer(kind=c_int), value,intent(in) :: sh2_d1
        integer(kind=c_int), value,intent(in) :: sh2_d2
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: lmax1
        integer(kind=c_int), value,intent(in) :: lmax2
        integer(kind=c_int), value,intent(in) :: lmaxt
        integer(kind=c_int), value,intent(in) :: k
        integer(kind=c_int), value,intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,lmaxt,k,taper_wt=taper_wt&
                                     ,norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine cSHMultiTaperMaskCSE

    subroutine cComputeDMap(Dij,Dij_d0,Dij_d1,dh_mask,dh_mask_d0,dh_mask_d1,n_dh,lmax&
                               ,sampling,degrees,degrees_d0,exitstatus)  bind(c, name="ComputeDMap")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDMap
        implicit none
        integer(kind=c_int), value,intent(in) :: Dij_d0
        integer(kind=c_int), value,intent(in) :: Dij_d1
        real(kind=c_double), dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer(kind=c_int), value,intent(in) :: dh_mask_d0
        integer(kind=c_int), value,intent(in) :: dh_mask_d1
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), value,intent(in) :: n_dh
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), value,intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling=sampling,degrees=degrees,exitstatus=exitstatus)
    end subroutine cComputeDMap

    subroutine cCurve2Mask(dhgrid,dhgrid_d0,dhgrid_d1,n,sampling,profile,profile_d0&
                                 ,profile_d1,nprofile,NP,extend,exitstatus)  bind(c&
                                 , name="Curve2Mask")
        use, intrinsic :: iso_c_binding
        use shtools, only: Curve2Mask
        implicit none
        integer(kind=c_int), value,intent(in) :: dhgrid_d0
        integer(kind=c_int), value,intent(in) :: dhgrid_d1
        integer(kind=c_int), dimension(dhgrid_d0,dhgrid_d1),intent(out) :: dhgrid
        integer(kind=c_int), value,intent(in) :: profile_d0
        integer(kind=c_int), value,intent(in) :: profile_d1
        real(kind=c_double), dimension(profile_d0,profile_d1),intent(in) :: profile
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), value,intent(in) :: sampling
        integer(kind=c_int), value,intent(in) :: nprofile
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        integer(kind=c_int) :: NP
        call Curve2Mask(dhgrid,n,sampling,profile,nprofile,NP,extend=extend,exitstatus=exitstatus)
    end subroutine cCurve2Mask

    subroutine cSHBias(Shh,Shh_d0,lwin,incspectra,incspectra_d0,ldata,outcspectra&
                          ,outcspectra_d0,save_cg,exitstatus)  bind(c, name="SHBias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBias
        implicit none
        integer(kind=c_int), value,intent(in) :: Shh_d0
        real(kind=c_double), dimension(Shh_d0),intent(in) :: Shh
        integer(kind=c_int), value,intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), value,intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), value,intent(in) :: lwin
        integer(kind=c_int), value,intent(in) :: ldata
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cSHBias

    subroutine cSphericalCapCoef(coef,coef_d0,theta,lmax,exitstatus)  bind(c, name="SphericalCapCoef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SphericalCapCoef
        implicit none
        integer(kind=c_int), value,intent(in) :: coef_d0
        real(kind=c_double), dimension(coef_d0),intent(out) :: coef
        real(kind=c_double), value,intent(in) :: theta
        integer(kind=c_int), optional,intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SphericalCapCoef(coef,theta,lmax=lmax,exitstatus=exitstatus)
    end subroutine cSphericalCapCoef

    subroutine cMakeGravGridDH(cilm,cilm_dim,lmax,gm,r0,a,f,rad,rad_d0,rad_d1,theta&
                                   ,theta_d0,theta_d1,phi,phi_d0,phi_d1,total,total_d0&
                                   ,total_d1,n,sampling,lmax_calc,omega,normal_gravity&
                                   ,pot,pot_d0,pot_d1,extend,exitstatus)  bind(c, name="MakeGravGridDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGridDH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: gm
        real(kind=c_double), value,intent(in) :: r0
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: f
        integer(kind=c_int), value,intent(in) :: rad_d0
        integer(kind=c_int), value,intent(in) :: rad_d1
        real(kind=c_double), dimension(rad_d0,rad_d1),intent(out) :: rad
        integer(kind=c_int), value,intent(in) :: theta_d0
        integer(kind=c_int), value,intent(in) :: theta_d1
        real(kind=c_double), dimension(theta_d0,theta_d1),intent(out) :: theta
        integer(kind=c_int), value,intent(in) :: phi_d0
        integer(kind=c_int), value,intent(in) :: phi_d1
        real(kind=c_double), dimension(phi_d0,phi_d1),intent(out) :: phi
        integer(kind=c_int), value,intent(in) :: total_d0
        integer(kind=c_int), value,intent(in) :: total_d1
        real(kind=c_double), dimension(total_d0,total_d1),intent(out) :: total
        real(kind=c_double), optional,intent(in) :: omega
        integer(kind=c_int), value,intent(in) :: pot_d0
        integer(kind=c_int), value,intent(in) :: pot_d1
        real(kind=c_double), optional,dimension(pot_d0,pot_d1),intent(out) :: pot
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: normal_gravity
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n,sampling=sampling&
                                ,lmax_calc=lmax_calc,omega=omega,normal_gravity=normal_gravity&
                                ,pot=pot,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGravGridDH

    subroutine cMakeGravGradGridDH(cilm,cilm_dim,lmax,gm,r0,a,f,vxx,vxx_d0,vxx_d1&
                                       ,vyy,vyy_d0,vyy_d1,vzz,vzz_d0,vzz_d1,vxy,vxy_d0&
                                       ,vxy_d1,vxz,vxz_d0,vxz_d1,vyz,vyz_d0,vyz_d1&
                                       ,n,sampling,lmax_calc,extend,exitstatus)  bind(c&
                                       , name="MakeGravGradGridDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGradGridDH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: gm
        real(kind=c_double), value,intent(in) :: r0
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: f
        integer(kind=c_int), value,intent(in) :: vxx_d0
        integer(kind=c_int), value,intent(in) :: vxx_d1
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        integer(kind=c_int), value,intent(in) :: vyy_d0
        integer(kind=c_int), value,intent(in) :: vyy_d1
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        integer(kind=c_int), value,intent(in) :: vzz_d0
        integer(kind=c_int), value,intent(in) :: vzz_d1
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        integer(kind=c_int), value,intent(in) :: vxy_d0
        integer(kind=c_int), value,intent(in) :: vxy_d1
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        integer(kind=c_int), value,intent(in) :: vxz_d0
        integer(kind=c_int), value,intent(in) :: vxz_d1
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        integer(kind=c_int), value,intent(in) :: vyz_d0
        integer(kind=c_int), value,intent(in) :: vyz_d1
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling=sampling&
                                    ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGravGradGridDH

    subroutine cMakeMagGradGridDH(cilm,cilm_dim,lmax,r0,a,f,vxx,vxx_d0,vxx_d1,vyy&
                                      ,vyy_d0,vyy_d1,vzz,vzz_d0,vzz_d1,vxy,vxy_d0&
                                      ,vxy_d1,vxz,vxz_d0,vxz_d1,vyz,vyz_d0,vyz_d1&
                                      ,n,sampling,lmax_calc,extend,exitstatus)  bind(c&
                                      , name="MakeMagGradGridDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGradGridDH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: r0
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: f
        integer(kind=c_int), value,intent(in) :: vxx_d0
        integer(kind=c_int), value,intent(in) :: vxx_d1
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        integer(kind=c_int), value,intent(in) :: vyy_d0
        integer(kind=c_int), value,intent(in) :: vyy_d1
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        integer(kind=c_int), value,intent(in) :: vzz_d0
        integer(kind=c_int), value,intent(in) :: vzz_d1
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        integer(kind=c_int), value,intent(in) :: vxy_d0
        integer(kind=c_int), value,intent(in) :: vxy_d1
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        integer(kind=c_int), value,intent(in) :: vxz_d0
        integer(kind=c_int), value,intent(in) :: vxz_d1
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        integer(kind=c_int), value,intent(in) :: vyz_d0
        integer(kind=c_int), value,intent(in) :: vyz_d1
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeMagGradGridDH(cilm,lmax,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling=sampling&
                                   ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeMagGradGridDH

    subroutine cMakeGeoidGrid(geoid,geoid_d0,geoid_d1,cilm,cilm_dim,lmax,r0pot,GM&
                                   ,PotRef,omega,r,gridtype,order,nlat,nlong,interval&
                                   ,lmax_calc,a,f,extend,exitstatus)  bind(c, name="MakeGeoidGrid")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGeoidGrid
        implicit none
        integer(kind=c_int), value,intent(in) :: geoid_d0
        integer(kind=c_int), value,intent(in) :: geoid_d1
        real(kind=c_double), dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: r0pot
        real(kind=c_double), value,intent(in) :: GM
        real(kind=c_double), value,intent(in) :: r
        real(kind=c_double), value,intent(in) :: PotRef
        real(kind=c_double), value,intent(in) :: omega
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: order
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        real(kind=c_double), optional,intent(in) :: interval
        real(kind=c_double), optional,intent(in) :: a
        real(kind=c_double), optional,intent(in) :: f
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,gridtype,order&
                                ,nlat,nlong,interval=interval,lmax_calc=lmax_calc&
                                ,a=a,f=f,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeGeoidGrid

    subroutine cCilmPlus(cilm,cilm_dim,gridin,gridin_d0,gridin_d1,lmax,nmax,mass,d&
                             ,rho,gridtype,w,w_d0,zero,zero_d0,plx,plx_d0,plx_d1,n&
                             ,dref,exitstatus)  bind(c, name="CilmPlus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlus
        implicit none
        integer(kind=c_int), value,intent(in) :: gridin_d0
        integer(kind=c_int), value,intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), value,intent(in) :: mass
        real(kind=c_double), value,intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero,plx=plx&
                          ,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cCilmPlus

    subroutine cCilmMinus(cilm,cilm_dim,gridin,gridin_d0,gridin_d1,lmax,nmax,mass&
                              ,d,rho,gridtype,w,w_d0,zero,zero_d0,plx,plx_d0,plx_d1&
                              ,n,dref,exitstatus)  bind(c, name="CilmMinus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinus
        implicit none
        integer(kind=c_int), value,intent(in) :: gridin_d0
        integer(kind=c_int), value,intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), value,intent(in) :: mass
        real(kind=c_double), value,intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero,plx=plx&
                           ,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cCilmMinus

    subroutine cCilmPlusRhoH(cilm,cilm_dim,gridin,gridin_d0,gridin_d1,lmax,nmax,mass&
                                 ,d,rho,rho_d0,rho_d1,gridtype,w,w_d0,zero,zero_d0&
                                 ,plx,plx_d0,plx_d1,n,dref,exitstatus)  bind(c, name="CilmPlusRhoH")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlusRhoH
        implicit none
        integer(kind=c_int), value,intent(in) :: gridin_d0
        integer(kind=c_int), value,intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), value,intent(in) :: mass
        integer(kind=c_int), value,intent(in) :: rho_d0
        integer(kind=c_int), value,intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero&
                              ,plx=plx,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cCilmPlusRhoH

    subroutine cCilmMinusRhoH(cilm,cilm_dim,gridin,gridin_d0,gridin_d1,lmax,nmax,mass&
                                  ,d,rho,rho_d0,rho_d1,gridtype,w,w_d0,zero,zero_d0&
                                  ,plx,plx_d0,plx_d1,n,dref,exitstatus)  bind(c, name="CilmMinusRhoH")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinusRhoH
        implicit none
        integer(kind=c_int), value,intent(in) :: gridin_d0
        integer(kind=c_int), value,intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), value,intent(in) :: mass
        integer(kind=c_int), value,intent(in) :: rho_d0
        integer(kind=c_int), value,intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero&
                               ,plx=plx,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cCilmMinusRhoH

    subroutine cBAtoHilm(cilm,cilm_dim,ba,ba_d0,ba_d1,ba_d2,gridglq,gridglq_d0,gridglq_d1&
                             ,lmax,nmax,mass,r0,rho,gridtype,w,w_d0,plx,plx_d0,plx_d1&
                             ,zero,zero_d0,filter_type,filter_deg,lmax_calc,exitstatus)  bind(c&
                             , name="BAtoHilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilm
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: ba_d0
        integer(kind=c_int), value,intent(in) :: ba_d1
        integer(kind=c_int), value,intent(in) :: ba_d2
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        integer(kind=c_int), value,intent(in) :: gridglq_d0
        integer(kind=c_int), value,intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), value,intent(in) :: mass
        real(kind=c_double), value,intent(in) :: r0
        real(kind=c_double), value,intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: filter_type
        integer(kind=c_int), optional,intent(in) :: filter_deg
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call BAtoHilm(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w=w,plx=plx,zero=zero&
                          ,filter_type=filter_type,filter_deg=filter_deg,lmax_calc=lmax_calc&
                          ,exitstatus=exitstatus)
    end subroutine cBAtoHilm

    subroutine cBAtoHilmRhoH(cilm,cilm_dim,ba,ba_d0,ba_d1,ba_d2,gridglq,gridglq_d0&
                                 ,gridglq_d1,lmax,nmax,mass,r0,rho,rho_d0,rho_d1,gridtype&
                                 ,w,w_d0,plx,plx_d0,plx_d1,zero,zero_d0,filter_type&
                                 ,filter_deg,lmax_calc,exitstatus)  bind(c, name="BAtoHilmRhoH")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilmRhoH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(out) :: cilm
        integer(kind=c_int), value,intent(in) :: ba_d0
        integer(kind=c_int), value,intent(in) :: ba_d1
        integer(kind=c_int), value,intent(in) :: ba_d2
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        integer(kind=c_int), value,intent(in) :: gridglq_d0
        integer(kind=c_int), value,intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), value,intent(in) :: mass
        real(kind=c_double), value,intent(in) :: r0
        integer(kind=c_int), value,intent(in) :: rho_d0
        integer(kind=c_int), value,intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), value,intent(in) :: plx_d0
        integer(kind=c_int), value,intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), value,intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: filter_type
        integer(kind=c_int), optional,intent(in) :: filter_deg
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call BAtoHilmRhoH(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w=w,plx=plx&
                              ,zero=zero,filter_type=filter_type,filter_deg=filter_deg&
                              ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cBAtoHilmRhoH

    function cDownContFilterMA(l,half,r,d)  bind(c, name="DownContFilterMA")
        use, intrinsic :: iso_c_binding
        use shtools, only: DownContFilterMA
        implicit none
        real(kind=c_double) :: cDownContFilterMA
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: half
        real(kind=c_double), value,intent(in) :: r
        real(kind=c_double), value,intent(in) :: d
        cDownContFilterMA=DownContFilterMA(l,half,r,d)
    end function cDownContFilterMA

    function cDownContFilterMC(l,half,r,d)  bind(c, name="DownContFilterMC")
        use, intrinsic :: iso_c_binding
        use shtools, only: DownContFilterMC
        implicit none
        real(kind=c_double) :: cDownContFilterMC
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: half
        real(kind=c_double), value,intent(in) :: r
        real(kind=c_double), value,intent(in) :: d
        cDownContFilterMC=DownContFilterMC(l,half,r,d)
    end function cDownContFilterMC

    function cNormalGravity(geocentric_lat,gm,omega,a,b)  bind(c, name="NormalGravity")
        use, intrinsic :: iso_c_binding
        use shtools, only: NormalGravity
        implicit none
        real(kind=c_double) :: cNormalGravity
        real(kind=c_double), value,intent(in) :: geocentric_lat
        real(kind=c_double), value,intent(in) :: gm
        real(kind=c_double), value,intent(in) :: omega
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: b
        cNormalGravity=NormalGravity(geocentric_lat,gm,omega,a,b)
    end function cNormalGravity

    subroutine cMakeMagGridDH(cilm,cilm_dim,lmax,r0,a,f,rad_grid,rad_grid_d0,rad_grid_d1&
                                  ,theta_grid,theta_grid_d0,theta_grid_d1,phi_grid&
                                  ,phi_grid_d0,phi_grid_d1,total_grid,total_grid_d0&
                                  ,total_grid_d1,n,sampling,lmax_calc,pot_grid,pot_grid_d0&
                                  ,pot_grid_d1,extend,exitstatus)  bind(c, name="MakeMagGridDH")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGridDH
        implicit none
        integer(kind=c_int), value,intent(in) :: cilm_dim
        real(kind=c_double), dimension(2,cilm_dim,cilm_dim),intent(in) :: cilm
        real(kind=c_double), value,intent(in) :: r0
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: f
        integer(kind=c_int), value,intent(in) :: rad_grid_d0
        integer(kind=c_int), value,intent(in) :: rad_grid_d1
        real(kind=c_double), dimension(rad_grid_d0,rad_grid_d1),intent(out) :: rad_grid
        integer(kind=c_int), value,intent(in) :: theta_grid_d0
        integer(kind=c_int), value,intent(in) :: theta_grid_d1
        real(kind=c_double), dimension(theta_grid_d0,theta_grid_d1),intent(out) :: theta_grid
        integer(kind=c_int), value,intent(in) :: phi_grid_d0
        integer(kind=c_int), value,intent(in) :: phi_grid_d1
        real(kind=c_double), dimension(phi_grid_d0,phi_grid_d1),intent(out) :: phi_grid
        integer(kind=c_int), value,intent(in) :: total_grid_d0
        integer(kind=c_int), value,intent(in) :: total_grid_d1
        real(kind=c_double), dimension(total_grid_d0,total_grid_d1),intent(out) :: total_grid
        integer(kind=c_int), value,intent(in) :: pot_grid_d0
        integer(kind=c_int), value,intent(in) :: pot_grid_d1
        real(kind=c_double), optional,dimension(pot_grid_d0,pot_grid_d1),intent(out) :: pot_grid
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid,total_grid&
                               ,n,sampling=sampling,lmax_calc=lmax_calc,pot_grid=pot_grid&
                               ,extend=extend,exitstatus=exitstatus)
    end subroutine cMakeMagGridDH

    subroutine cSHMagPowerSpectrum(c,c_d0,c_d1,c_d2,a,r,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                    , name="SHMagPowerSpectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerSpectrum
        implicit none
        integer(kind=c_int), value,intent(in) :: c_d0
        integer(kind=c_int), value,intent(in) :: c_d1
        integer(kind=c_int), value,intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: r
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMagPowerSpectrum(c,a,r,lmax,spectra,exitstatus=exitstatus)
    end subroutine cSHMagPowerSpectrum

    function cSHMagPowerL(c,c_d0,c_d1,c_d2,a,r,l)  bind(c, name="SHMagPowerL")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerL
        implicit none
        real(kind=c_double) :: cSHMagPowerL
        integer(kind=c_int), value,intent(in) :: c_d0
        integer(kind=c_int), value,intent(in) :: c_d1
        integer(kind=c_int), value,intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), value,intent(in) :: a
        real(kind=c_double), value,intent(in) :: r
        integer(kind=c_int), value,intent(in) :: l
        cSHMagPowerL=SHMagPowerL(c,a,r,l)
    end function cSHMagPowerL

    subroutine cMakeCircleCoord(coord,coord_d0,coord_d1,lat,lon,theta0,cinterval,cnum&
                                     ,exitstatus)  bind(c, name="MakeCircleCoord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeCircleCoord
        implicit none
        real(kind=c_double), value,intent(in) :: lat
        real(kind=c_double), value,intent(in) :: lon
        real(kind=c_double), value,intent(in) :: theta0
        integer(kind=c_int), value,intent(in) :: coord_d0
        integer(kind=c_int), value,intent(in) :: coord_d1
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), optional,intent(in) :: cinterval
        integer(kind=c_int), optional,intent(out) :: cnum
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval=cinterval,cnum=cnum,exitstatus=exitstatus)
    end subroutine cMakeCircleCoord

    subroutine cMakeEllipseCoord(coord,coord_d0,coord_d1,lat,lon,dec,A_theta,B_theta&
                                      ,cinterval,cnum,exitstatus)  bind(c, name="MakeEllipseCoord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeEllipseCoord
        implicit none
        real(kind=c_double), value,intent(in) :: lat
        real(kind=c_double), value,intent(in) :: lon
        real(kind=c_double), value,intent(in) :: A_theta
        real(kind=c_double), value,intent(in) :: B_theta
        real(kind=c_double), value,intent(in) :: dec
        integer(kind=c_int), value,intent(in) :: coord_d0
        integer(kind=c_int), value,intent(in) :: coord_d1
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), optional,intent(in) :: cinterval
        integer(kind=c_int), optional,intent(out) :: cnum
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,cinterval=cinterval&
                                   ,cnum=cnum,exitstatus=exitstatus)
    end subroutine cMakeEllipseCoord

    subroutine cWigner3j(w3j,w3j_d0,jmin,jmax,j2,j3,m1,m2,m3,exitstatus)  bind(c, name="Wigner3j")
        use, intrinsic :: iso_c_binding
        use shtools, only: Wigner3j
        implicit none
        integer(kind=c_int), value,intent(in) :: j2
        integer(kind=c_int), value,intent(in) :: j3
        integer(kind=c_int), value,intent(in) :: m1
        integer(kind=c_int), value,intent(in) :: m2
        integer(kind=c_int), value,intent(in) :: m3
        integer(kind=c_int), intent(out) :: jmin
        integer(kind=c_int), intent(out) :: jmax
        integer(kind=c_int), value,intent(in) :: w3j_d0
        real(kind=c_double), dimension(w3j_d0),intent(out) :: w3j
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,exitstatus=exitstatus)
    end subroutine cWigner3j

    function cRandomN(idum)  bind(c, name="RandomN")
        use, intrinsic :: iso_c_binding
        use shtools, only: RandomN
        implicit none
        real(kind=c_double) :: cRandomN
        integer(kind=c_int), intent(inout) :: idum
        cRandomN=RandomN(idum)
    end function cRandomN

    function cRandomGaussian(idum)  bind(c, name="RandomGaussian")
        use, intrinsic :: iso_c_binding
        use shtools, only: RandomGaussian
        implicit none
        real(kind=c_double) :: cRandomGaussian
        integer(kind=c_int), intent(inout) :: idum
        cRandomGaussian=RandomGaussian(idum)
    end function cRandomGaussian

    subroutine cPreGLQ(x1,x2,n,zero,zero_d0,w,w_d0,exitstatus)  bind(c, name="PreGLQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: PreGLQ
        implicit none
        real(kind=c_double), value,intent(in) :: x1
        real(kind=c_double), value,intent(in) :: x2
        integer(kind=c_int), value,intent(in) :: zero_d0
        real(kind=c_double), dimension(zero_d0),intent(out) :: zero
        integer(kind=c_int), value,intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(out) :: w
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PreGLQ(x1,x2,n,zero,w,exitstatus=exitstatus)
    end subroutine cPreGLQ

    function cNGLQ(degree)  bind(c, name="NGLQ")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQ
        implicit none
        integer(kind=c_int) :: cNGLQ
        integer(kind=c_int), value,intent(in) :: degree
        cNGLQ=NGLQ(degree)
    end function cNGLQ

    function cNGLQSH(degree)  bind(c, name="NGLQSH")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQSH
        implicit none
        integer(kind=c_int) :: cNGLQSH
        integer(kind=c_int), value,intent(in) :: degree
        cNGLQSH=NGLQSH(degree)
    end function cNGLQSH

    function cNGLQSHN(degree,n)  bind(c, name="NGLQSHN")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQSHN
        implicit none
        integer(kind=c_int) :: cNGLQSHN
        integer(kind=c_int), value,intent(in) :: degree
        integer(kind=c_int), value,intent(in) :: n
        cNGLQSHN=NGLQSHN(degree,n)
    end function cNGLQSHN

    subroutine cDHaj(n,aj,aj_d0,extend,exitstatus)  bind(c, name="DHaj")
        use, intrinsic :: iso_c_binding
        use shtools, only: DHaj
        implicit none
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), value,intent(in) :: aj_d0
        real(kind=c_double), dimension(aj_d0),intent(out) :: aj
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call DHaj(n,aj,extend=extend,exitstatus=exitstatus)
    end subroutine cDHaj

    function cYilmIndexVector(i,l,m)  bind(c, name="YilmIndexVector")
        use, intrinsic :: iso_c_binding
        use shtools, only: YilmIndexVector
        implicit none
        integer(kind=c_int) :: cYilmIndexVector
        integer(kind=c_int), value,intent(in) :: i
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: m
        cYilmIndexVector=YilmIndexVector(i,l,m)
    end function cYilmIndexVector

    subroutine cEigValVecSym(ain,ain_d0,ain_d1,n,eig,eig_d0,evec,evec_d0,evec_d1,ul&
                                ,K,exitstatus)  bind(c, name="EigValVecSym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSym
        implicit none
        integer(kind=c_int), value,intent(in) :: ain_d0
        integer(kind=c_int), value,intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), value,intent(in) :: eig_d0
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        integer(kind=c_int), value,intent(in) :: evec_d0
        integer(kind=c_int), value,intent(in) :: evec_d1
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), optional,intent(in) :: ul
        integer(kind=c_int), optional,intent(in) :: K
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call EigValVecSym(ain,n,eig,evec,ul=ul,K=K,exitstatus=exitstatus)
    end subroutine cEigValVecSym

    subroutine cEigValVecSymTri(ain,ain_d0,ain_d1,n,eig,eig_d0,evec,evec_d0,evec_d1&
                                   ,ul,exitstatus)  bind(c, name="EigValVecSymTri")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSymTri
        implicit none
        integer(kind=c_int), value,intent(in) :: ain_d0
        integer(kind=c_int), value,intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), value,intent(in) :: eig_d0
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        integer(kind=c_int), value,intent(in) :: evec_d0
        integer(kind=c_int), value,intent(in) :: evec_d1
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), optional,intent(in) :: ul
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call EigValVecSymTri(ain,n,eig,evec,ul=ul,exitstatus=exitstatus)
    end subroutine cEigValVecSymTri

    subroutine cEigValSym(ain,ain_d0,ain_d1,n,eval,eval_d0,ul)  bind(c, name="EigValSym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValSym
        implicit none
        integer(kind=c_int), value,intent(in) :: ain_d0
        integer(kind=c_int), value,intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), value,intent(in) :: n
        integer(kind=c_int), value,intent(in) :: eval_d0
        real(kind=c_double), dimension(eval_d0),intent(out) :: eval
        character(kind=c_char), optional,intent(in) :: ul
        call EigValSym(ain,n,eval,ul=ul)
    end subroutine cEigValSym

    subroutine cSHRotateTapers(tapersrot,tapersrot_d0,tapersrot_d1,tapers,tapers_d0&
                                        ,tapers_d1,taper_order,taper_order_d0,lmax&
                                        ,nrot,x,x_d0,dj,dj_d0,dj_d1,dj_d2,exitstatus)  bind(c&
                                        , name="SHRotateTapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateTapers
        implicit none
        integer(kind=c_int), value,intent(in) :: tapers_d0
        integer(kind=c_int), value,intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), value,intent(in) :: x_d0
        real(kind=c_double), dimension(x_d0),intent(in) :: x
        integer(kind=c_int), value,intent(in) :: dj_d0
        integer(kind=c_int), value,intent(in) :: dj_d1
        integer(kind=c_int), value,intent(in) :: dj_d2
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer(kind=c_int), value,intent(in) :: tapersrot_d0
        integer(kind=c_int), value,intent(in) :: tapersrot_d1
        real(kind=c_double), dimension(tapersrot_d0,tapersrot_d1),intent(out) :: tapersrot
        integer(kind=c_int), value,intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nrot
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateTapers(tapersrot,tapers,taper_order,lmax,nrot,x,dj,exitstatus=exitstatus)
    end subroutine cSHRotateTapers

    subroutine cSlepianCoeffs(falpha,falpha_d0,galpha,galpha_d0,galpha_d1,flm,flm_d0&
                                    ,flm_d1,flm_d2,lmax,nmax,exitstatus)  bind(c, name="SlepianCoeffs")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffs
        implicit none
        integer(kind=c_int), value,intent(in) :: falpha_d0
        real(kind=c_double), dimension(falpha_d0),intent(out) :: falpha
        integer(kind=c_int), value,intent(in) :: galpha_d0
        integer(kind=c_int), value,intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), value,intent(in) :: flm_d0
        integer(kind=c_int), value,intent(in) :: flm_d1
        integer(kind=c_int), value,intent(in) :: flm_d2
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(in) :: flm
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SlepianCoeffs(falpha,galpha,flm,lmax,nmax,exitstatus=exitstatus)
    end subroutine cSlepianCoeffs

    subroutine cSlepianCoeffsToSH(flm,flm_d0,flm_d1,flm_d2,falpha,falpha_d0,galpha&
                                     ,galpha_d0,galpha_d1,lmax,nmax,exitstatus)  bind(c&
                                     , name="SlepianCoeffsToSH")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffsToSH
        implicit none
        integer(kind=c_int), value,intent(in) :: flm_d0
        integer(kind=c_int), value,intent(in) :: flm_d1
        integer(kind=c_int), value,intent(in) :: flm_d2
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(out) :: flm
        integer(kind=c_int), value,intent(in) :: falpha_d0
        real(kind=c_double), dimension(falpha_d0),intent(in) :: falpha
        integer(kind=c_int), value,intent(in) :: galpha_d0
        integer(kind=c_int), value,intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SlepianCoeffsToSH(flm,falpha,galpha,lmax,nmax,exitstatus=exitstatus)
    end subroutine cSlepianCoeffsToSH

    subroutine cSHSCouplingMatrix(kij,kij_d0,kij_d1,galpha,galpha_d0,galpha_d1,lmax&
                                     ,nmax,exitstatus)  bind(c, name="SHSCouplingMatrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrix
        implicit none
        integer(kind=c_int), value,intent(in) :: kij_d0
        integer(kind=c_int), value,intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        integer(kind=c_int), value,intent(in) :: galpha_d0
        integer(kind=c_int), value,intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSCouplingMatrix(kij,galpha,lmax,nmax,exitstatus=exitstatus)
    end subroutine cSHSCouplingMatrix

    subroutine cSHSlepianVar(l,galpha,galpha_d0,galpha_d1,galpha_order,galpha_order_d0&
                              ,lmax,kmax,Sff,Sff_d0,variance,exitstatus)  bind(c, name="SHSlepianVar")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSlepianVar
        implicit none
        integer(kind=c_int), value,intent(in) :: galpha_d0
        integer(kind=c_int), value,intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), value,intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), value,intent(in) :: l
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: kmax
        integer(kind=c_int), value,intent(in) :: galpha_order_d0
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSlepianVar(l,galpha,galpha_order,lmax,kmax,Sff,variance,exitstatus=exitstatus)
    end subroutine cSHSlepianVar

    subroutine cSHSCouplingMatrixCap(kij,kij_d0,kij_d1,galpha,galpha_d0,galpha_d1&
                                        ,galpha_order,galpha_order_d0,lmax,nmax,exitstatus)  bind(c&
                                        , name="SHSCouplingMatrixCap")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrixCap
        implicit none
        integer(kind=c_int), value,intent(in) :: kij_d0
        integer(kind=c_int), value,intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        integer(kind=c_int), value,intent(in) :: galpha_d0
        integer(kind=c_int), value,intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), value,intent(in) :: galpha_order_d0
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), value,intent(in) :: lmax
        integer(kind=c_int), value,intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSCouplingMatrixCap(kij,galpha,galpha_order,lmax,nmax,exitstatus=exitstatus)
    end subroutine cSHSCouplingMatrixCap

