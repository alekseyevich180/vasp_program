# 1 "wannier_interpol.F"
# 1 "./symbol.inc" 1 
!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define 1             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 59


!
!   charge density: half grid mode Z direction
!










# 101


!
!   charge density real
!




# 118


# 130

!
!   wavefunctions: full grid mode
!

# 167

!
!   wavefunctions complex
!




























!
!   common definitions
!







# 211


!
!   mpi parallel macros
!














# 248

# 254

# 259

# 266

# 273



# 282







# 297







# 306


# 319










# 336

# 2 "wannier_interpol.F" 2 

! Auxiliary functions for debugging
! Subroutines type overloaded
! No special facilities for parallel runs!
module util
  use prec
  
  interface mat_dump
     module procedure mat_dump_cq, mat_dump_rq,mat_dump_cqs, mat_dump_rqs
  end interface mat_dump

  interface mat_diag
     module procedure mat_diag_cq
  end interface mat_diag

  interface is_hermitian
     module procedure is_hermitian_cq, is_hermitian_rq
  end interface is_hermitian

  interface is_symmetric
     module procedure is_hermitian_rq
  end interface is_symmetric

contains

  subroutine mat_dump_cq(m,txt,nc,nr)
    implicit none
    complex(q) :: m(:,:)
    character(len=*), optional :: txt
    integer, optional :: nc,nr
    integer :: i1,n1,n2

    n1=size(m,1)
    n2=size(m,2)
    if(present(nc))then
       n2=nc
       n1=nc
    end if
    if(present(nr))then
       n1=nr
    end if
    
    write(*,*) "------------------------------------------------------------------"
    if(present(txt))then
       write(*,'(a)')  txt
    end if
    do i1=1,n1
       write(*,'(100e12.3)') real(m(i1,1:n2),q)
    end do
    write(*,*)
    do i1=1,n1
       write(*,'(100e12.3)') aimag(m(i1,1:n2))
    end do

  end subroutine mat_dump_cq

  subroutine mat_dump_cqs(m,txt,nc,nr)
    implicit none
    complex(qs) :: m(:,:)
    character(len=*), optional :: txt
    integer, optional :: nc,nr
    integer :: i1,n1,n2

    n1=size(m,1)
    n2=size(m,2)
    if(present(nc))then
       n2=nc
       n1=nc
    end if
    if(present(nr))then
       n1=nr
    end if
    
    write(*,*) "------------------------------------------------------------------"
    if(present(txt))then
       write(*,'(a)')  txt
    end if
    do i1=1,n1
       write(*,'(100e12.3)') real(m(i1,1:n2),qs)
    end do
    write(*,*)
    do i1=1,n1
       write(*,'(100e12.3)') aimag(m(i1,1:n2))
    end do

  end subroutine mat_dump_cqs


  subroutine mat_dump_rq(m,txt,nc,nr)
    implicit none
    real(q) :: m(:,:)
    character(len=*), optional :: txt
    integer, optional :: nc,nr
    integer :: i1,n1,n2

    n1=size(m,1)
    n2=size(m,2)
    if(present(nc))then
       n2=nc
       n1=nc
    end if
    if(present(nr))then
       n1=nr
    end if
    
    write(*,*) "------------------------------------------------------------------"
    if(present(txt))then
       write(*,'(a)')  txt
    end if
    do i1=1,n1
       write(*,'(100e12.3)') m(i1,1:n2)
    end do

  end subroutine mat_dump_rq

  subroutine mat_dump_rqs(m,txt,nc,nr)
    implicit none
    real(qs) :: m(:,:)
    character(len=*), optional :: txt
    integer, optional :: nc,nr
    integer :: i1,n1,n2

    n1=size(m,1)
    n2=size(m,2)
    if(present(nc))then
       n2=nc
       n1=nc
    end if
    if(present(nr))then
       n1=nr
    end if
    
    write(*,*) "------------------------------------------------------------------"
    if(present(txt))then
       write(*,'(a)')  txt
    end if
    do i1=1,n1
       write(*,'(100e12.3)') m(i1,1:n2)
    end do

  end subroutine mat_dump_rqs

  function is_hermitian_cq(m,tol)
    implicit none
    logical :: is_hermitian_cq
    complex(q) :: m(:,:)
    integer :: n1,n2,i,j
    real(q), optional :: tol
    real(q) :: tol0=1e-10_q
    n1=size(m,1)
    n2=size(m,2)
    if(n1/=n2)then
       write(*,*) 'IS_HERMITIAN: matrix not quadratic'
       stop
    end if
    if(present(tol)) tol0=tol
    is_hermitian_cq=.true.
    do i=1,n1
       do j=1,n2
          if( abs(m(i,j)-conjg(m(j,i)))>tol0 )then
             is_hermitian_cq=.false.
             return
          end if
       end do
    end do
  end function is_hermitian_cq


  function is_hermitian_rq(m,tol)
    implicit none
    logical :: is_hermitian_rq
    real(q) :: m(:,:)
    integer :: n1,n2,i,j
    real(q), optional :: tol
    real(q) :: tol0=1e-8_q
    n1=size(m,1)
    n2=size(m,2)
    if(n1/=n2)then
       write(*,*) 'IS_HERMITIAN: matrix not quadratic'
       stop
    end if
    if(present(tol)) tol0=tol
    is_hermitian_rq=.true.
    do i=1,n1
       do j=1,n2
          if( abs(m(i,j)-m(j,i))>tol0 )then
             is_hermitian_rq=.false.
             return
          end if
       end do
    end do
  end function is_hermitian_rq


  subroutine mat_diag_cq(m,typ,ev,dump)
    implicit none
    complex(q), intent(in) :: m(:,:)
    character(len=1) :: typ
    real(q), intent(inout), optional :: ev(:)
    logical, optional :: dump
!
    real(q),allocatable :: w(:),rwork(:)
    complex(q),allocatable :: mcopy(:,:),work(:)
    integer :: n,lwork,info,j
    logical :: ldump
!
    n=size(m,2)
    lwork=max(1,2*n)
    allocate(w(n),rwork(3*n),work(lwork))
    allocate(mcopy(size(m,1),size(m,2)))
    mcopy=m
    ldump=.false.
    if(present(dump))ldump=dump

    select case (typ)
    case ('h','H')
! hermitian
       call zheev('N','U',n,mcopy,size(mcopy,1),w,work,lwork,rwork,info)
       if(info/=0)then
          write(*,*) 'MAT_DIAG: '
          stop
       endif
       if(ldump)then
          write(*,*) 'MAT_DIAG: Eigenvalues'
          do j=1,n
             write(*,*) w(j)
          end do
       end if
    case default
! general complex
! TODO
       write(*,*) 'MAT_DIAG: general complex matrix not implemented!'
       stop
    end select
    if(present(ev)) ev=w
    deallocate(w,rwork,work,mcopy)
  end subroutine mat_diag_cq
  
end module util


module interpolate
  use prec
  implicit none

  interface inter_1d_linear
     module procedure inter_1d_linear_rq, inter_1d_linear_cq, inter_1d_linear_rq_v, inter_1d_linear_cq_v
  end interface inter_1d_linear

  interface inter_1d_cubic
     module procedure inter_1d_cubic_rq, inter_1d_cubic_cq, inter_1d_cubic_rq_v, inter_1d_cubic_cq_v
  end interface inter_1d_cubic

  interface inter_2d_cubic
     module procedure inter_2d_cubic_rq
  end interface inter_2d_cubic

  interface inter_3d_cubic
     module procedure inter_3d_cubic_rq, inter_3d_cubic_cq
  end interface inter_3d_cubic


contains

!=================================================
! LINEAR INTERPOLATION of a function at point r in [0,1]
! given function values
! f = a*r + b
!
! f0=f(0)
! f1=f(1)
!-------------------------------------------------
  subroutine inter_1d_linear_rq(r,f0,f1,f)
    implicit none
    real(q), intent(in) :: r,f0,f1
    real(q), intent(out) :: f
!
    f=f0+r*(f1-f0)
  end subroutine inter_1d_linear_rq

  subroutine inter_1d_linear_cq(r,f0,f1,f)
    implicit none
    real(q), intent(in) :: r
    complex(q), intent(in) :: f0,f1
    complex(q), intent(out) :: f
!
    f=f0+r*(f1-f0)
  end subroutine inter_1d_linear_cq
!-------------------------------------------------
! vector versions
! Usage:
! 1. calculate weights w(2)
! 2. interpolate
! (Advantage: weights can be used multiple times)
!-------------------------------------------------
  subroutine inter_1d_linear_weights(r,w)
    implicit none
    real(q), intent(in) :: r
    real(q), intent(out) :: w(2)
!
! weights for function interpolation
    w(1)=1-r     ! f(0)
    w(2)=r       ! f(1)
  end subroutine inter_1d_linear_weights

  subroutine inter_1d_linear_rq_v(w,f0,f1,f)
    implicit none
    real(q), intent(in) :: w(2),f0(:),f1(:)
    real(q), intent(out) :: f(:)
!
    f=w(1)*f0+w(2)*f1
  end subroutine inter_1d_linear_rq_v

  subroutine inter_1d_linear_cq_v(w,f0,f1,f)
    implicit none
    real(q), intent(in) :: w(2)
    complex(q), intent(in) :: f0(:),f1(:)
    complex(q), intent(out) :: f(:)
!
    f=w(1)*f0+w(2)*f1
  end subroutine inter_1d_linear_cq_v
!=================================================
  

!=================================================
! QUADRATIC INTERPOLATION
! NOT IMPLEMENTED
!=================================================


!=================================================
! CUBIC INTERPOLATION of a function at point r in [0,1]
! given function and derivative values
! f = a*r^3 + b*r^2 + c*r + d
! g = 3*a*r^2 + 2*b*r + c
!
! f0=f(0)
! f1=f(1)
! g0=df/dr(0)=a*df/dx(x0)
! g1=df/dr(1)=a*df/dx(x0+a)
! if x=x0+a*r (i.e. scaled coordinate with scaling factor a)
! then df/dr = df/dx*dx/dr = a*df/dx
!-------------------------------------------------
  subroutine inter_1d_cubic_rq(r,f0,f1,g0,g1,f,g)
    implicit none
    real(q), intent(in) :: r,f0,f1,g0,g1
    real(q), intent(out) :: f
    real(q), intent(out), optional :: g
    real(q) a,b,c,d
!
    a=2*f0-2*f1+g0+g1
    b=-2*g0-g1+3*f1-3*f0
    c=g0
    d=f0
!
    f=d + c*r + b*r**2 + a*r**3
    if(present(g)) g=c + 2*b*r + 3*a*r**2
  end subroutine inter_1d_cubic_rq

  subroutine inter_1d_cubic_cq(r,f0,f1,g0,g1,f,g)
    implicit none
    real(q), intent(in) :: r
    complex(q), intent(in) :: f0,f1,g0,g1
    complex(q), intent(out) :: f
    complex(q), intent(out), optional :: g
    complex(q) a,b,c,d
!
    a=2*f0-2*f1+g0+g1
    b=-2*g0-g1+3*f1-3*f0
    c=g0
    d=f0
!
    f=d + c*r + b*r**2 + a*r**3
    if(present(g)) g=c + 2*b*r + 3*a*r**2
  end subroutine inter_1d_cubic_cq
!-------------------------------------------------
! vector versions
!-------------------------------------------------
  subroutine inter_1d_cubic_weights(r,w,v)
    implicit none
    real(q), intent(in) :: r
    real(q), intent(out) :: w(4)
    real(q), intent(out), optional :: v(4)
!
! weights for function interpolation
    w(1)=1-3*r**2+2*r**3      ! f(0)
    w(2)=3*r**2-2*r**3        ! f(1)
    w(3)=r-2*r**2+r**3        ! f'(0)
    w(4)=-r**2+r**3           ! f'(1)
! weights for derivative interpolation
    if(present(v))then
       v(1)=-6*r+6*r**2     ! f(0)
       v(2)=6*r-6*r**2      ! f(1)
       v(3)=1-4*r+3*r**2    ! f'(0)
       v(4)=-2*r+3*r**2     ! f'(1)
    end if
  end subroutine inter_1d_cubic_weights

  subroutine inter_1d_cubic_rq_v(w,f0,f1,g0,g1,f,v,g)
    implicit none
    real(q), intent(in) :: w(4),f0(:),f1(:),g0(:),g1(:)
    real(q), intent(out) :: f(:)
    real(q), intent(in), optional :: v(4)
    real(q), intent(out), optional :: g(:)
    real(q) a,b,c,d
!
    f=w(1)*f0+w(2)*f1+w(3)*g0+w(4)*g1
    if(present(g)) g=v(1)*f0+v(2)*f1+v(3)*g0+v(4)*g1
  end subroutine inter_1d_cubic_rq_v

  subroutine inter_1d_cubic_cq_v(w,f0,f1,g0,g1,f,v,g)
    implicit none
    real(q), intent(in) :: w(4)
    complex(q), intent(in) :: f0(:),f1(:),g0(:),g1(:)
    complex(q), intent(out) :: f(:)
    real(q), intent(in), optional :: v(4)
    complex(q), intent(out), optional :: g(:)
!
    f=w(1)*f0+w(2)*f1+w(3)*g0+w(4)*g1
    if(present(g)) g=v(1)*f0+v(2)*f1+v(3)*g0+v(4)*g1
  end subroutine inter_1d_cubic_cq_v
!============================================


!============================================
! 2D INTERPOLATION
!-------------------------------------------------
  subroutine inter_2d_cubic_rq(r,f,gx,gy,v)
    implicit none
    real(q), intent(in) :: r(2),f(2,2),gx(2,2),gy(2,2)
    real(q), intent(out) :: v
    real(q) :: f1(2),gy1(2)
!
! interpolate in x direction
    call inter_1d_cubic(r(1),f(1,1),f(2,1),gx(1,1),gx(2,1),f1(1))
    call inter_1d_linear(r(1),gy(1,1),gy(2,1),gy1(1))
    call inter_1d_cubic(r(1),f(1,2),f(2,2),gx(1,2),gx(2,2),f1(2))
    call inter_1d_linear(r(1),gy(1,2),gy(2,2),gy1(2))
!
! interpolate in y direction
    call inter_1d_cubic(r(2),f1(1),f1(2),gy1(1),gy1(2),v)
  end subroutine inter_2d_cubic_rq
!============================================


!============================================
! 3D INTERPOLATION
!-------------------------------------------------
  subroutine inter_3d_cubic_rq(r,f,gx,gy,gz,v)
    implicit none
    real(q), intent(in) :: r(3),f(2,2,2),gx(2,2,2),gy(2,2,2),gz(2,2,2)
    real(q), intent(out) :: v
    real(q) :: f1(2,2),gy1(2,2),gz1(2,2),f2(2),gz2(2)
!
! interpolate in x direction on 4 edges
! edge 1,1
    call inter_1d_cubic(r(1),f(1,1,1),f(2,1,1),gx(1,1,1),gx(2,1,1),f1(1,1))
    call inter_1d_linear(r(1),gy(1,1,1),gy(2,1,1),gy1(1,1))
    call inter_1d_linear(r(1),gz(1,1,1),gz(2,1,1),gz1(1,1))
! edge 1,2
    call inter_1d_cubic(r(1),f(1,1,2),f(2,1,2),gx(1,1,2),gx(2,1,2),f1(1,2))
    call inter_1d_linear(r(1),gy(1,1,2),gy(2,1,2),gy1(1,2))
    call inter_1d_linear(r(1),gz(1,1,2),gz(2,1,2),gz1(1,2))
! edge 2,1
    call inter_1d_cubic(r(1),f(1,2,1),f(2,2,1),gx(1,2,1),gx(2,2,1),f1(2,1))
    call inter_1d_linear(r(1),gy(1,2,1),gy(2,2,1),gy1(2,1))
    call inter_1d_linear(r(1),gz(1,2,1),gz(2,2,1),gz1(2,1))
! edge 2,2
    call inter_1d_cubic(r(1),f(1,2,2),f(2,2,2),gx(1,2,2),gx(2,2,2),f1(2,2))
    call inter_1d_linear(r(1),gy(1,2,2),gy(2,2,2),gy1(2,2))
    call inter_1d_linear(r(1),gz(1,2,2),gz(2,2,2),gz1(2,2))
    
! interpolate in y direction on 2 edges
! edge x,1
    call inter_1d_cubic(r(2),f1(1,1),f1(2,1),gy1(1,1),gy1(2,1),f2(1))
    call inter_1d_linear(r(2),gz1(1,1),gz1(2,1),gz2(1))    
! edge x,2
    call inter_1d_cubic(r(2),f1(1,2),f1(2,2),gy1(1,2),gy1(2,2),f2(2))
    call inter_1d_linear(r(2),gz1(1,2),gz1(2,2),gz2(2))    

! interpolate in z direction on 1 line
    call inter_1d_cubic(r(3),f2(1),f2(2),gz2(1),gz2(2),v)

  end subroutine inter_3d_cubic_rq

  subroutine inter_3d_cubic_cq(r,f,gx,gy,gz,v)
    implicit none
    real(q), intent(in) :: r(3)
    complex(q), intent(in) :: f(2,2,2),gx(2,2,2),gy(2,2,2),gz(2,2,2)
    complex(q), intent(out) :: v
    complex(q) :: f1(2,2),gy1(2,2),gz1(2,2),f2(2),gz2(2)
!
! interpolate in x direction on 4 edges
! edge 1,1
    call inter_1d_cubic(r(1),f(1,1,1),f(2,1,1),gx(1,1,1),gx(2,1,1),f1(1,1))
    call inter_1d_linear(r(1),gy(1,1,1),gy(2,1,1),gy1(1,1))
    call inter_1d_linear(r(1),gz(1,1,1),gz(2,1,1),gz1(1,1))
! edge 1,2
    call inter_1d_cubic(r(1),f(1,1,2),f(2,1,2),gx(1,1,2),gx(2,1,2),f1(1,2))
    call inter_1d_linear(r(1),gy(1,1,2),gy(2,1,2),gy1(1,2))
    call inter_1d_linear(r(1),gz(1,1,2),gz(2,1,2),gz1(1,2))
! edge 2,1
    call inter_1d_cubic(r(1),f(1,2,1),f(2,2,1),gx(1,2,1),gx(2,2,1),f1(2,1))
    call inter_1d_linear(r(1),gy(1,2,1),gy(2,2,1),gy1(2,1))
    call inter_1d_linear(r(1),gz(1,2,1),gz(2,2,1),gz1(2,1))
! edge 2,2
    call inter_1d_cubic(r(1),f(1,2,2),f(2,2,2),gx(1,2,2),gx(2,2,2),f1(2,2))
    call inter_1d_linear(r(1),gy(1,2,2),gy(2,2,2),gy1(2,2))
    call inter_1d_linear(r(1),gz(1,2,2),gz(2,2,2),gz1(2,2))
    
! interpolate in y direction on 2 edges
! edge x,1
    call inter_1d_cubic(r(2),f1(1,1),f1(2,1),gy1(1,1),gy1(2,1),f2(1))
    call inter_1d_linear(r(2),gz1(1,1),gz1(2,1),gz2(1))    
! edge x,2
    call inter_1d_cubic(r(2),f1(1,2),f1(2,2),gy1(1,2),gy1(2,2),f2(2))
    call inter_1d_linear(r(2),gz1(1,2),gz1(2,2),gz2(2))    

! interpolate in z direction on 1 line
    call inter_1d_cubic(r(3),f2(1),f2(2),gz2(1),gz2(2),v)

  end subroutine inter_3d_cubic_cq
!============================================

end module interpolate



MODULE wannier_interpolation
CONTAINS
  SUBROUTINE wannier_interpolation_dummy
    WRITE(*,*)'Im a DEC compiler so I need this line'
  END SUBROUTINE wannier_interpolation_dummy
END MODULE wannier_interpolation
# 3671

