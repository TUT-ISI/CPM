FUNCTION DIFFC(Temp,Press,I)
  use precisi
  implicit none
!!$C
!!$C This function is used to calculate the diffusion coefficient from the
!!$C vapor in the gas as a function of temperature and pressure after expans.
!!$C I is species number.
!!$C
  REAL(doubp), intent(in) :: Temp, Press
  REAL(doubp) :: DIFFC
  integer :: I
!!$  REAL(dp) ::   MYY(6), DN(6)
  real(doubp), dimension(1:6), parameter :: &
       & MYY = (/1.6658_doubp, 1.75_doubp , 1.75_doubp , 1.75_doubp , 1.75_doubp , 0.0_doubp  /)
  real(doubp), dimension(1:6), parameter :: &
       & DN = (/22.0748e-6_doubp , 12.98e-6_doubp , 9.380307e-6_doubp , 20.0e-6_doubp, &
       & 14.76518e-6_doubp , 0.0_doubp  /)

  DIFFC = DN(I) * ( (TEMP/273.15_doubp )**MYY(I) ) * 101325.0_doubp  / Press
END FUNCTION DIFFC




