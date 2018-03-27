! test_mod_ls65.f90
!
! gfortran -Wall -O -o test_mod_ls65.ex mod_least65.f90 test_mod_ls65.f90 -llapack -lblas
! time ./test_mod_ls65.ex < tsls65.dat > tsls65.out
! ifort -mkl -O -o test_mod_ls65.ex mod_least65.f90 test_mod_ls65.f90
!
!-----
MODULE least_var
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4), parameter :: set_NSMAX=2048
  real(DP) :: Q(set_NSMAX)
!
END MODULE least_var
!
!-----
program ts_modls65
!
  use mod_least, only : NPR,NS,NP,NPFIT,IPFIT,PA,IPRSW,IDRSW,OBS,WSQRT,LEAST,least_alloc,least_dealloc
  use least_var, only : Q,set_NSMAX
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  real(DP) :: EPSR,SUMS,SIGST
  integer(4) :: I,MAXX,ICON
!
  NPR=6
  IDRSW=1  ! NUMERICAL DERIVATIVE
!
  READ(*,*) NS,NP,NPFIT
  WRITE(*,'(a,3I5)')'  NS,NP,NPFIT=',NS,NP,NPFIT
  if(NS > set_NSMAX) then
    write(*,*)' invalid input'  ;  stop
  endif
  call least_alloc
  READ(*,*) (IPFIT(I),I=1,NPFIT)
  WRITE(*,'(a)')'  IPFIT'
  WRITE(*,'(5I3,2X,5I3)')(IPFIT(I),I=1,NPFIT)
  READ(*,*) (PA(I),I=1,NP)
  WRITE(*,'(a)')'  PA'
  WRITE(*,'(2G15.7)')(PA(I),I=1,2)
  WRITE(*,'(3G15.7)')(PA(I),I=3,NP)
  READ(*,*) EPSR,MAXX,IPRSW
  WRITE(*,'(a,g15.7,I10)')'  EPSR,MAXX=',EPSR,MAXX
  WRITE(*,'(a,I10)')'  IPRSW=',IPRSW
  WRITE(*,'(a)')'  IDRSW=1   NUMERICAL DERIVATIVE'
  do i=1,NS
    READ(*,*)Q(I),OBS(I),WSQRT(I)
  end do
!
  CALL LEAST(EPSR,MAXX,SUMS,SIGST,ICON)
  WRITE(*,'(''    LEAST ICON='',I7)')ICON
  call least_dealloc
!
end program ts_modls65
!
!-----
subroutine FMODL
!
  use mod_least, only : PA,CALC,NS
  use least_var, only : Q
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
  real(DP) :: xx
  integer(4) :: I,j
!
      DO I=1,NS
        xx=PA(1)+PA(2)*Q(I)
        do j=0,18
         xx=xx+PA(3+j*3)/( ( ( Q(I)-PA(4+j*3) )/PA(5+j*3) )**2 + 1.0d0 )
        end do
        CALC(I)=xx
      end do
!
end subroutine FMODL
!
!-----
subroutine DFMODL
  IMPLICIT NONE
!
  continue
  return
!
end subroutine DFMODL
!

