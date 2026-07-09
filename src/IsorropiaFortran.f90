module IsorropiaFortran
  use, intrinsic :: iso_c_binding
  implicit none

  ! Modern Fortran types mapping exactly to standard Isorropia C layout structures
  type, bind(C) :: IsorropiaInput
    real(c_double) :: w(8)
    real(c_double) :: org(3)
    real(c_double) :: waer(8)
    real(c_double) :: temp
    real(c_double) :: rh
    integer(c_int) :: iprob
    integer(c_int) :: nadj
  end type IsorropiaInput

  type, bind(C) :: IsorropiaState
    real(c_double) :: temp
    real(c_double) :: rh
    real(c_double) :: w(8)
    real(c_double) :: waer(8)
    real(c_double) :: org(3)

    real(c_double) :: molal(10)
    real(c_double) :: molalr(23)
    real(c_double) :: gama(23)
    real(c_double) :: zz(23)
    real(c_double) :: z(10)
    real(c_double) :: gamou(23)
    real(c_double) :: gamin(23)
    real(c_double) :: m0(23)
    real(c_double) :: gasaq(3)
    integer(c_int) :: actmod
    real(c_double) :: epsact
    real(c_double) :: coh
    real(c_double) :: chno3
    real(c_double) :: chcl
    real(c_double) :: water
    real(c_double) :: ionic
    real(c_double) :: watcmp(24)
    integer(c_int) :: frst
    integer(c_int) :: calain
    integer(c_int) :: calaou
    integer(c_int) :: dryf

    real(c_double) :: ch2so4, cnh42s4, cnh4hs4, cnacl, cna2so4, cnano3, cnh4no3, cnh4cl, cnahso4, clc
    real(c_double) :: ccaso4, ccano32, ccacl2, ck2so4, ckhso4, ckno3, ckcl, cmgso4, cmgno32, cmgcl2
    real(c_double) :: gnh3, ghno3, ghcl

    integer(c_int) :: num_errors
  end type IsorropiaState

  ! Define standard interlanguage binding interfaces to link directly to flat C symbols
  interface
    subroutine isorropia_solve_c(input, state) bind(C, name="isorropia_solve_c")
      import :: IsorropiaInput, IsorropiaState
      type(IsorropiaInput), intent(in)  :: input
      type(IsorropiaState), intent(out) :: state
    end subroutine isorropia_solve_c
  end interface

contains

  ! Native high-level Fortran subroutine wrapping C solves for host integration (like CATChem)
  subroutine isorropia_solve_f(input, state)
    type(IsorropiaInput), intent(in)   :: input
    type(IsorropiaState), intent(out)  :: state
    call isorropia_solve_c(input, state)
  end subroutine isorropia_solve_f

end module IsorropiaFortran
