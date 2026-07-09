program test_linkage
  use IsorropiaFortran
  implicit none

  type(IsorropiaInput) :: input
  type(IsorropiaState) :: state

  print *, "=== Isorropia Fortran-to-C++ Linkage Test ==="

  ! Initialize inputs
  input%w = 0.0d0
  input%w(2) = 1.0d0 ! H2SO4 component (1-indexed inside Fortran array)
  input%w(3) = 2.0d0 ! NH3 component
  input%w(4) = 1.0d0 ! HNO3 component

  input%org = 0.0d0
  input%org(1) = 10.0d0 ! org (concentration)
  input%org(2) = 0.15d0 ! korg (hygroscopicity)
  input%org(3) = 1000.0d0 ! density (kg/m3)

  input%rh = 0.80d0
  input%temp = 298.15d0
  input%iprob = 0
  input%nadj = 1

  ! Call native standard Fortran solver wrapper
  call isorropia_solve_f(input, state)

  print *, "Simulation completed successfully."
  print *, "Calculated Dynamic Aerosol Water:", state%water
  print *, "Calculated Gaseous Ammonia (NH3) :", state%gnh3
  print *, "Calculated Gaseous Nitric Acid  :", state%ghno3
  print *, "Diagnostic Error Counts         :", state%num_errors

  ! Mathematically verify the inter-language link and calculations
  if (state%water > 8.0d0 .and. state%gnh3 > 1.5d0 .and. state%ghno3 > 0.7d0 .and. state%num_errors == 0) then
     print *, "✅ FORTRAN-TO-C++ BINDING LINKAGE TEST PASSED SUCCESSFULLY!"
  else
     print *, "❌ FORTRAN-TO-C++ BINDING LINKAGE TEST FAILED."
     call exit(1)
  end if

end program test_linkage
