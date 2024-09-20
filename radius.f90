SUBROUTINE RADIUS(IBin)
  !C
  USE headfile
  use precisi
  implicit none
  real(doubp) :: Volume
  integer, intent(in) :: IBin
  integer :: ISPec
  !C
  Volume = 0.0_doubp
  DO ISpec = 1, NALiquids + NASolids
     Volume = Volume + C(InC(ISpec,IBin)) / C(Ibin)*WtMol(ISpec) & 
          &  / Dens(ISpec) * 1.e-06_doubp
  END DO
  !C
  Rp(IBin) = (Volume / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!  IF(IBIN.EQ.20) write(*,*)  Rp(IBin),C(InC(1,IBin)),CBin(Ibin)
END SUBROUTINE RADIUS





SUBROUTINE RADIUS_Y2(IBin)
  !C
  USE headfile
  use precisi
  implicit none
  real(doubp) :: Volume
  integer, intent(in) :: IBin
  integer :: ISPec
  !C

  Volume = 0.0_doubp

  IF (C_Y2(Ibin) > 1.0e-40) THEN
     DO ISpec = 1, NALiquids + NASolids
        Volume = Volume + C_Y2(InCL(ISpec,IBin)) / C_Y2(Ibin)*WtMol(ISpec) & 
             &  / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_CL(IBin) = (Volume / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
  ELSE
     R_CL(IBin) = 0.5_doubp*r_Y2(IBin) + 0.5_doubp*r_Y2(IBin+1)
  ENDIF
!  IF(IBIN.EQ.20) write(*,*)  Rp(IBin),C(InC(1,IBin)),CBin(Ibin)
END SUBROUTINE RADIUS_Y2


SUBROUTINE RADIUS_DRY(IBin)
  !C
  USE headfile
  use precisi
  implicit none
  real(doubp) :: Volume
  integer, intent(in) :: IBin
  integer :: ISPec
  !C
  Volume = 0.0_doubp
  DO ISpec = 2, NALiquids + NASolids
     Volume = Volume + C(InC(ISpec,IBin)) / C(Ibin)*WtMol(ISpec) & 
          &  / Dens(ISpec) * 1.e-06_doubp
  END DO
  !C
  R_dry(IBin) = (Volume / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!  IF(IBIN.EQ.20) write(*,*)  Rp(IBin),C(InC(1,IBin)),CBin(Ibin)
END SUBROUTINE RADIUS_DRY



SUBROUTINE RADIUS_Y2dry(IBin)
  !C
  USE headfile
  use precisi
  implicit none
  real(doubp) :: Volume
  integer, intent(in) :: IBin
  integer :: ISPec
  !C

  Volume = 0.0_doubp

  IF (C_Y2(Ibin) > 1.0e-8) THEN
     DO ISpec = 2, NALiquids + NASolids
        Volume = Volume + C_Y2(InCL(ISpec,IBin)) / C_Y2(Ibin)*WtMol(ISpec) & 
             &  / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_CLdry(IBin) = (Volume / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
  ELSE
     R_CLdry(IBin) = 0.5_doubp*r_Y2(IBin) + 0.5_doubp*r_Y2(IBin+1)
  ENDIF
!  IF(IBIN.EQ.20) write(*,*)  Rp(IBin),C(InC(1,IBin)),CBin(Ibin)
END SUBROUTINE RADIUS_Y2dry
