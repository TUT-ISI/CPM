FUNCTION SURFACE_TENSION(T,wtH2SO4,wtHNO3)

!!$C
  use precisi
  implicit none
  real(doubp) :: SURFACE_TENSION
  REAL(doubp), intent(in) :: T, wtH2SO4, wtHNO3
  real(doubp) :: g253, SUM293, g293, SUM253, s
  REAL(doubp), dimension(0:5,0:6,2) :: Y
  integer I, J
!!$C
!!$C     T=253K
!!$C
  Y(0,0:5,1) = (/ 3.601e7_doubp,  5.894e5_doubp,  5.397e4_doubp, &
       &       -2.994e3_doubp,6.919e1_doubp, -5.648e-1_doubp /)
  Y(1,0:5,1) = (/ 8.127e5_doubp,  0.0_doubp,     0.0_doubp, &
       &       0.0_doubp, 1.001e-2_doubp, 0.0_doubp /)
  Y(2,0:5,1) = (/ 2.194e4_doubp,  0.0_doubp,     0.0_doubp, &
       &       -9.681e-2_doubp, 0.0_doubp,     0.0_doubp /)
  Y(3,0:5,1) = (/ -4.554e2_doubp,  0.0_doubp,     8.724e-2_doubp, &
       &       0.0_doubp, 0.0_doubp,     0.0_doubp /)
  Y(4,0:5,1) = (/ 7.115_doubp, -3.648e-2_doubp, 0.0_doubp, &
       &       0.0_doubp, 0.0_doubp,     0.0_doubp /)
  Y(5,0:5,1) = (/ -4.483e-2_doubp, 0.0_doubp,     0.0_doubp, &
       &       0.0_doubp,  0.0_doubp,     0.0_doubp /)
!!$C
  Y(0,6,1) = 6.726e2_doubp
  Y(1,6,1) = 9.692_doubp
  Y(2,6,1) = 8.276_doubp
!!$C
!!$C     T=293 K
!!$C
  Y(0,0:5,2) = (/ 7.201e7_doubp,  3.893e6_doubp,  9.736e4_doubp, &
       &       -1.832e3_doubp, 1.282e1_doubp,  1.076e-1_doubp /)
  Y(1,0:5,2) = (/ 1.073e6_doubp,  0.0_doubp,     0.0_doubp, &
       &       0.0_doubp, -2.811e-1_doubp, 0.0_doubp /)
  Y(2,0:5,2) = (/ 8.580e3_doubp,  0.0_doubp,     0.0_doubp, &
       &       3.358e-1_doubp, 0.0_doubp,     0.0_doubp /)
  Y(3,0:5,2) = (/ -1.036e2_doubp,  0.0_doubp,    -1.866e-1_doubp, &
       &       0.0_doubp, 0.0_doubp,     0.0_doubp /)
  Y(4,0:5,2) = (/ 2.270_doubp,  4.895e-3_doubp, 0.0_doubp, &
       &       0.0_doubp, 0.0_doubp,     0.0_doubp /)
  Y(5,0:5,2) = (/ -2.333e-2_doubp, 0.0_doubp,     0.0_doubp, &
       &       0.0_doubp, 0.0_doubp,     0.0_doubp /)
!!$C
  Y(0,6,2) = 9.949e2_doubp
  Y(1,6,2) = 7.321_doubp
  Y(2,6,2) = 2.917e1_doubp
!!$C
  SUM253=0.0_doubp
!!$C
  DO I=0,5
     DO J=0,5-I
        SUM253=SUM253+Y(I,J,1)*(wtH2SO4**i)*(wtHNO3**j)
     END DO
  END DO
  g253=SUM253/(Y(0,6,1)+Y(1,6,1)*wtH2SO4+Y(2,6,1)*wtHNO3)**2
!!$C
  SUM293=0.0_doubp
!!$C
  DO I=0,5
     DO J=0,5-I
        SUM293=SUM293+Y(I,J,2)*(wtH2SO4**i)*(wtHNO3**j)
     END DO
  END DO
!!$C
  g293=SUM293/(Y(0,6,2)+Y(1,6,2)*wtH2SO4+Y(2,6,2)*wtHNO3)**2
!!$C
  s=(g293-g253)/40.0_doubp
  SURFACE_TENSION = s*(T-253.0_doubp)+g253
  SURFACE_TENSION = 1.0e-3_doubp*SURFACE_TENSION
END FUNCTION SURFACE_TENSION
