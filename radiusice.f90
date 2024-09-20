SUBROUTINE RADIUSICE(IBin)
!C
  use headfile
  use precisi
  implicit none
  real(doubp) :: Volumeice, rho_ice
  integer, intent(in) :: IBin
  integer :: ISpec
!C
  Volumeice = 0.0_doubp
  !C     First, the volume is calculated for ice
  Volumeice=C_i(InCi(IH2Ol,IBin))/C_i(IBin)*WtMol(IH2Ol)/(rho_ice(T)*1.0e3_doubp)
!C     and here for other species (how reliable are the densities used here
  !C     with these cold temperatures? need to be improved)!
  DO ISpec = 2, NALiquids + NASolids
     Volumeice = Volumeice + C_i(InCi(ISpec,IBin)) / C_i(IBin)* &
          &      WtMol(ISpec) / Dens(ISpec) * 1.e-06_doubp
  END DO
!C
  R_ice(IBin) = (Volumeice / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!  R_ice(IBin)=0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!  write(*,*) R_ice(IBin), 0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!C
END SUBROUTINE RADIUSICE

SUBROUTINE RADIUSICE2(IBin)
!C
  use headfile
  use precisi
  implicit none
  real(doubp) :: Volumeice, rho_ice
  integer, intent(in) :: IBin
  integer :: ISpec
!C
  IF(C_Y1(IBin)>1.0e-8) THEN
     Volumeice = 0.0_doubp
     !C     First, the volume is calculated for ice
     Volumeice=C_Y1(InCi(IH2Ol,IBin))/C_Y1(IBin)*WtMol(IH2Ol)/(rho_ice(T)*1.0e3_doubp)
     !C     and here for other species (how reliable are the densities used here
     !C     with these cold temperatures? need to be improved)!
     DO ISpec = 2, NALiquids + NASolids
        Volumeice = Volumeice + C_Y1(InCi(ISpec,IBin)) / C_Y1(IBin)* &
             &      WtMol(ISpec) / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_ice(IBin) = (Volumeice / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!     IF(IBIN.EQ.20) write(*,*)  R_ice(IBin),C_Y1(InCi(IH2Ol,IBin)),C_Y1(20)
  ELSE
     R_ice(IBin)=R_Y1(IBIN)*.5+R_Y1(IBin+1)*.5
  ENDIF
!  R_ice(IBin)=0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!  write(*,*) R_ice(IBin), 0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!C
END SUBROUTINE RADIUSICE2

SUBROUTINE RADIUS_ICE_DRY(IBin)
!C
  use headfile
  use precisi
  implicit none
  real(doubp) :: Volumeice, rho_ice
  integer, intent(in) :: IBin
  integer :: ISpec
!C
  IF(C_Y1(IBin)>1.0e-8) THEN
     Volumeice = 0.0_doubp

     DO ISpec = 2, NALiquids + NASolids
        Volumeice = Volumeice + C_Y1(InCi(ISpec,IBin)) / C_Y1(IBin)* &
             &      WtMol(ISpec) / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_ice_dry(IBin) = (Volumeice / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!     IF(IBIN.EQ.20) write(*,*)  R_ice(IBin),C_Y1(InCi(IH2Ol,IBin)),C_Y1(20)
  ELSE
     R_ice_dry(IBin)=R_Y1(IBIN)*.5+R_Y1(IBin+1)*.5
  ENDIF
!  R_ice(IBin)=0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!  write(*,*) R_ice(IBin), 0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!C
END SUBROUTINE RADIUS_ICE_DRY

SUBROUTINE RADIUSICE3(IBin)
!C
  use headfile
  use precisi
  implicit none
  real(doubp) :: Volumeice, rho_ice
  integer, intent(in) :: IBin
  integer :: ISpec
!C
  IF(C_Y3(IBin)>1.0e-8) THEN
     Volumeice = 0.0_doubp
     !C     First, the volume is calculated for ice
     Volumeice=C_Y3(InC3(IH2Ol,IBin))/C_Y3(IBin)*WtMol(IH2Ol)/(rho_ice(T)*1.0e3_doubp)
     !C     and here for other species (how reliable are the densities used here
     !C     with these cold temperatures? need to be improved)!
     DO ISpec = 2, NALiquids + NASolids
        Volumeice = Volumeice + C_Y3(InC3(ISpec,IBin)) / C_Y3(IBin)* &
             &      WtMol(ISpec) / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_C3(IBin) = (Volumeice / (4.0_doubp/3.0_doubp*Pi))**(1./3.)
!     IF(IBIN.EQ.20) write(*,*)  R_ice(IBin),C_Y1(InCi(IH2Ol,IBin)),C_Y1(20)
  ELSE
     R_C3(IBin)=R_Y3(IBIN)*.5+R_Y3(IBin+1)*.5
  ENDIF
!  R_ice(IBin)=0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!  write(*,*) R_ice(IBin), 0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!C
END SUBROUTINE RADIUSICE3

SUBROUTINE RADIUSICE4(IBin)
!C
  use headfile
  use precisi
  implicit none
  real(doubp) :: Volumeice, rho_ice
  integer, intent(in) :: IBin
  integer :: ISpec
!C
  IF(C_Y4(IBin)>1.0e-8) THEN
     Volumeice = 0.0_doubp
     !C     First, the volume is calculated for ice
     Volumeice=C_Y4(InC3(IH2Ol,IBin))/C_Y4(IBin)*WtMol(IH2Ol)/(rho_ice(T)*1.0e3_doubp)
     !C     and here for other species (how reliable are the densities used here
     !C     with these cold temperatures? need to be improved)!
     DO ISpec = 2, NALiquids + NASolids
        Volumeice = Volumeice + C_Y4(InC4(ISpec,IBin)) / C_Y4(IBin)* &
             &      WtMol(ISpec) / Dens(ISpec) * 1.e-06_doubp
     END DO
     !C
     R_C4(IBin) = (Volumeice / (4.0_doubp/3.0_doubp*Pi))**(1.0_doubp/3.0_doubp)
!     IF(IBIN.EQ.20) write(*,*)  R_ice(IBin),C_Y1(InCi(IH2Ol,IBin)),C_Y1(20)
  ELSE
     R_C4(IBin)=R_Y4(IBIN)*.5_doubp+R_Y3(IBin+1)*.5_doubp
  ENDIF
!  R_ice(IBin)=0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!  write(*,*) R_ice(IBin), 0.5_doubp*1.0e-2_doubp*(so_mass(IBin)/0.0219_doubp)**(1./2.6)
!C
END SUBROUTINE RADIUSICE4
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$C      Subroutine radiusLiq2Ice(IBin)
!!$CC
!!$C      INCLUDE 'headfile.f'
!!$CC
!!$C      double precision Volumeice
!!$C      Volumeice = 0.d0
!!$CC     First, the volume is calculated for ice
!!$C      Volumeice=C(InC(1,IBin))/CBin(IBin)*WtMol(1)/(rho_ice(T)*1.d3)
!!$CC     and here for other species (how reliable are the densities used here
!!$CC     with these cold temperatures? need to be improved)
!!$C      DO 100 ISpec = 2, NALiquids + NASolids
!!$C         Volumeice = Volumeice + C(InC(ISpec,IBin)) / CBin(IBin)*
!!$C     >        WtMol(ISpec) / Dens(ISpec) * 1.d-06
!!$C 100  CONTINUE
!!$CC
!!$C      R_ice_temp = (Volumeice/ (4./3.*Pi))**(1./3.)
!!$CC     
!!$C      end subroutine radiusLiq2Ice
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$C      Subroutine radiuscombice(IBin,JBin)
!!$CC
!!$C      INCLUDE 'headfile.f'
!!$CC
!!$C      double precision Volumeice
!!$C      Volumeice = 0.d0
!!$CC     First, the volume is calculated for ice
!!$C      Volumeice=(C_i(InCi(1,JBin))/C_i(JBin)+C(InC(1,IBin))/CBin(IBin))
!!$C     >     *WtMol(1)/(rho_ice(T)*1.d3)
!!$CC     and here for other species (how reliable are the densities used here
!!$CC     with these cold temperatures? need to be improved)
!!$C      DO 100 ISpec = 2, NALiquids + NASolids
!!$C         Volumeice = Volumeice + (C_i(InCi(ISpec,JBin)) / C_i(JBin) +
!!$C     >        C(InC(ISpec,IBin)) / CBin(IBin))*
!!$C     >        WtMol(ISpec) / Dens(ISpec) * 1.d-06
!!$C 100  CONTINUE
!!$CC
!!$C      R_ice(JBin) = (Volumeice/ (4./3.*Pi))**(1./3.)
!!$CC     
!!$C      end subroutine radiuscombice
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$C      Subroutine radiuscombexice(IBin,JBin,ration)
!!$CC
!!$C      INCLUDE 'headfile.f'
!!$CC     
!!$C      double precision Volumeice, ration
!!$C      Volumeice = 0.d0
!!$CC     First, the volume is calculated for ice
!!$C      Volumeice=(C_i(InCi(1,JBin))/C_i(JBin)+C(InC(1,IBin))/CBin(IBin)*
!!$C     >     ration)*WtMol(1)/(rho_ice(T)*1.d3)
!!$CC     and here for other species (how reliable are the densities used here
!!$CC     with these cold temperatures? need to be improved)
!!$C      DO 100 ISpec = 2, NALiquids + NASolids
!!$C         Volumeice = Volumeice + (C_i(InCi(ISpec,JBin)) / C_i(JBin) +
!!$C     >        C(InC(ISpec,IBin)) / CBin(IBin)*ration)*
!!$C     >        WtMol(ISpec) / Dens(ISpec) * 1.d-06
!!$C 100  CONTINUE
!!$CC
!!$C      R_ice(JBin) = (Volumeice/ (4./3.*Pi))**(1./3.)
!!$CC     
!!$C      end subroutine radiuscombexice
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$C      Subroutine radiuscomb2ice(IBin,JBin)
!!$CC
!!$C      INCLUDE 'headfile.f'
!!$CC
!!$C      double precision Volumeice
!!$C      Volumeice = 0.d0
!!$CC     First, the volume is calculated for ice
!!$C      Volumeice=(C_i(InCi(1,JBin))/C_i(JBin)+C_i(InCi(1,IBin))/C_i(IBin)
!!$C     >     )*WtMol(1)/(rho_ice(T)*1.d3)
!!$CC     and here for other species (how reliable are the densities used here
!!$CC     with these cold temperatures? need to be improved)
!!$C      DO 100 ISpec = 2, NALiquids + NASolids
!!$C         Volumeice = Volumeice + (C_i(InCi(ISpec,JBin)) / C_i(JBin) +
!!$C     >        C_i(InCi(ISpec,IBin)) / C_i(IBin))*
!!$C     >        WtMol(ISpec) / Dens(ISpec) * 1.d-06
!!$C 100  CONTINUE
!!$CC
!!$C      R_ice(IBin) = (Volumeice/ (4./3.*Pi))**(1./3.)
!!$CC     
!!$C      end subroutine radiuscomb2ice
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
