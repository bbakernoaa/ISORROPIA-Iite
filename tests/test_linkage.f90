program test_linkage
  use IsorropiaFortran
  implicit none

  type(IsorropiaInput) :: input
  type(IsorropiaState) :: state

  print *, "=== Isorropia Fortran-to-C++ Linkage Test ==="

  ! Initialize inputs (scaled to standard internal physical units)
  input%w = 0.0d0
  input%w(2) = (1.0d0 / 98.0d0) * 1.0d-6 ! H2SO4 component (1-indexed inside Fortran array) -> mol/m3
  input%w(3) = (2.0d0 / 17.0d0) * 1.0d-6 ! NH3 component -> mol/m3
  input%w(4) = (1.0d0 / 63.0d0) * 1.0d-6 ! HNO3 component -> mol/m3

  input%org = 0.0d0
  input%org(1) = 10.0d0 * 1.0d-9 ! org (concentration ug/m3 to kg/m3)
  input%org(2) = 0.15d0          ! korg (hygroscopicity)
  input%org(3) = 1000.0d0        ! density (kg/m3)

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

  ! Mathematically verify the inter-language link and calculations using standard physical scales
  if (state%water > 7.0d-9 .and. state%gnh3 > 9.0d-8 .and. state%ghno3 > 1.0d-8 .and. state%num_errors == 0) then
     print *, "✅ FORTRAN-TO-C++ BINDING LINKAGE TEST PASSED SUCCESSFULLY!"
  else
     print *, "❌ FORTRAN-TO-C++ BINDING LINKAGE TEST FAILED."
     call exit(1)
  end if

end program test_linkage
