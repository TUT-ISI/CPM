SUBROUTINE SIZEDIST
!!$C
!!$C     Calculate the number size distribution and
!!$C     
!!$C
  USE headfile
  use precisi
  implicit none
  real(doubp), dimension(NBins) :: Vlo, Vhi, dD, Vtemp
  real(doubp) :: Dlo, Dhi, Vrat, VFac,Vrat2,vrat3,c_kpl, sum_d 
  integer :: IBin, ISpec, IBini, I, NBinLo, NBinHi,II
  integer :: LRAD

!!$C     Initialization for ice (InCi)
  DO IBini     = 1,NY1
     DO ISpec  = 1, NALiquids + NASolids
        InCi(ISpec,IBini) = (NFBins*(ISpec - 1)) + IBini + NFBins
        C_Y1(Inci(ISpec,IBini)) = 0.0_doubp
     END DO
  END DO
  !!     Initialization for CL (InCl)
  DO IBini     = 1,NY2
     DO ISpec  = 1, NALiquids + NASolids
        InCl(ISpec,IBini) = (NY2*(ISpec - 1)) + IBini + NY2
        C_Y2(Incl(ISpec,IBini)) = 0.0_doubp
     END DO
  END DO
  
  DO IBini     = 1,NY3
     DO ISpec  = 1, NALiquids + NASolids
        InC3(ISpec,IBini) = (NY3*(ISpec - 1)) + IBini + NY3
     END DO
  END DO
  
  
  DO IBini     = 1,NY4
     DO ISpec  = 1, NALiquids + NASolids
        InC4(ISpec,IBini) = (NY4*(ISpec - 1)) + IBini + NY4
     END DO
  END DO


!!$C
!!$C     If modes are inseparable:
!!$C
  
  IF(IMModes == 1 .OR. NAModes == 1) THEN
!!$C
!!$C     Initialize vectors and matrices and
!!$C     Set index for where a given species 
!!$C     ISPEC in size bin IBIN can be found 
!!$C     in vector C
!!$C
     DO IBin     = 1,NABins
        CBin(IBin)  = 0.0_doubp
        Vtemp(IBin) = 0.0_doubp
        DO ISpec  = 1, NALiquids + NASolids
           InC(ISpec,IBin) = (NABins*(ISpec - 1)) + IBin + NABins
           C(Inc(ISpec,IBin)) = 0.0_doubp
        END DO
     END DO
!!$C
     ! Here also for gases (gases has non-bin-structure so the value is in last
     ! "slot" (NABins))
      DO ISpec = NALiquids + NASolids + 1, NAGases + NALiquids + NASolids
        InC(ISpec,NABins) = (NABins + NABins*NALiquids + NABins*NASolids + &
             & ISpec - NALiquids - NASolids)
     end DO
             

!!$C
     IF(GMD(NAModes) <= GMD(1) .AND. NAModes /= 1) THEN
        WRITE(*,*) 'Geometric mean diameters must be in'
        WRITE(*,*) 'ascending order in "aerosol.dat"'
        STOP
     END IF
!!$C
!!$C     Create sections for volume ratio size distribution:
!!$C
!!$C     Dlo  = diameter of the smallest size bin (m)
!!$C     Dhi  = diameter of the largest size bin (m)
!!$C     Vrat = volume ratio
!!$C     Vlo  = min volume in the size bin (m^3)
!!$C     Vhi  = max volume in the size bin (m^3)
!!$C     dD   = width of size bin (m)
!!$C
     Dlo    = GMD(1) / (1.5*STD(1))**2
     Dhi    = GMD(NAModes) * (1.5 * STD(NAModes))**2
     
     
     Vrat   = (Dhi/Dlo)**(3./REAL(NABins-1))
     Vlo(1) = (1./3.*PI*Dlo**3) / (1+Vrat)
!!$C     
     DO IMode     = 1, NAModes
!!$C     
        DO IBin    = 1, NABins
!!$C     
           IModeBin(Ibin) = Imode
           Vhi(IBin)  = Vrat*Vlo(IBin)
           dD(IBin)   = (6./PI)**(1./3.) * (Vhi(IBin)**(1./3.) - Vlo(IBin)**(1./3.))
           Dp(IBin)   = (6./2.*(Vhi(IBin) + Vlo(IBin))/pi)**(1./3.)
           Rp(IBin)   = Dp(IBin) / 2.
           CBin(IBin) = CBin(IBin) + CMode(IMode)*dD(IBin) / & 
                &              (Dp(IBin)*SQRT(2.*PI)*LOG(STD(IMode))) * & 
                &              EXP(-(log(Dp(IBin)/GMD(IMode)))**2 / &
                &              (2.*(log(STD(IMode)))**2))
           C(IBin)    = CBin(IBin)
           Ctot       = Ctot + C(IBin)
           IF(IBin < NABins) Vlo(IBin + 1) = Vhi(IBin)
!!$C
!!$C     Convert weight fraction to moles. The value is to be corrected later.
!!$C     Also, the volume of particles in size bin is calculated.
!!$C
           DO ISpec = 1, NALiquids + NASolids
              C(InC(ISpec,IBin)) = C(InC(ISpec,IBin)) + WtFrac(ISpec,IMode) / &
                   & WtMol(ISpec)
              Vtemp(IBin) = Vtemp(IBin) + C(InC(ISpec,IBin)) * WtMol(ISpec) / &
                   & Dens(ISpec) * 1.0e-06_doubp
           END DO
!!$C
        END DO
!!$C     
     END DO
!!$C
!!$C     Convert moles to molar concentrations and normalize the 
!!$C     concentrations so that particles size equals Rp
!!$C     
     DO IBin  = 1, NABins
        VFac      = (4./3.*PI*Rp(IBin)**3) / Vtemp(IBin)
        DO ISpec = 1, NALiquids + NASolids
           C(InC(ISpec,IBin)) = C(InC(ISpec,IBin))*VFac*CBin(IBin)
!!$C               WRITE(*,*) FormSpec(ISpec)
!!$C               WRITE(*,*) C(InC(ISpec,IBin))
!!$C               PAUSE
        END DO
     END DO
!!$C
!!$C     ELSE If modes are separable:
!!$C
  ELSE
!!$C
     NBinLo           = 1
     NBinHi           = 0
     NABins           = 0
!!$C         
!!$C     Calculate the number of size bins used in a simulation
!!$C
     DO IMode     = 1, NAModes
        FirstBinMode(IMode) = NABins+1
        NABins        = NABins + NBinsPerMode(IMode) 
     END DO
!!$C
     DO IMode     = 1, NAModes
!!$C     
!!$C     NBinHi = top limit of bin number
!!$C
        NBinHi        = NBinHi + NBinsPerMode(IMode)
!!$C
        Dlo    = GMD(IMode) / (1.5_doubp * STD(IMode))**2
        Dlo=max(Dlo,1e-8)
        Dhi    = GMD(IMode) * (1.5_doubp * STD(IMode))**2
        Vrat   = (Dhi/Dlo)**(3./REAL(NBinsPerMode(IMode)-1))
        Vlo(NBinLo) = (1./3.*PI*Dlo**3) / (1+Vrat)
!!$C
        DO IBin    = NBinLo, NBinHi
!!$C     
           IModeBin(Ibin) = Imode
           Vhi(IBin)   = Vrat*Vlo(IBin)
           dD(IBin)    = (6./PI)**(1./3.) * (Vhi(IBin)**(1./3.) - Vlo(IBin)**(1./3.))
           Dp(IBin)    = (6./2.*(Vhi(IBin) + Vlo(IBin))/pi)**(1./3.)
!           Dp(IBin)    = (6.*Vhi(IBin)/pi)**(1./3.)*.5 + (6.*Vlo(IBin)/pi)**(1./3.)*.5
           Rp(IBin)    = Dp(IBin) / 2.
           c_kpl = 0.0
           sum_d = 0.0
           DO  II = 1,50
              
              c_kpl = c_kpl+ CMode(IMode)*dD(IBin)*.02/  &
                   (((6.*(Vlo(IBin))/pi)**(1./3.) +(II-.5)*dD(IBin)*.02)*SQRT(2.*PI)*LOG(STD(IMode)))* &
                   EXP(-(log(((6.*(Vlo(IBin))/pi)**(1./3.) +(II-.5)*dD(IBin)*.02)/GMD(IMode))**2 / &
                &        (2.*(log(STD(IMode)))**2)))
              sum_d = sum_d + CMode(IMode)*dD(IBin)*.02/  &
                   (((6.*(Vlo(IBin))/pi)**(1./3.) +(II-.5)*dD(IBin)*.02)*SQRT(2.*PI)*LOG(STD(IMode)))* &
                   EXP(-(log(((6.*(Vlo(IBin))/pi)**(1./3.) +(II-.5)*dD(IBin)*.02)/GMD(IMode))**2 / &
                &        (2.*(log(STD(IMode)))**2)))*((6.*(Vlo(IBin))/pi)**(1./3.) +(II-.5)*dD(IBin)*.02)**3.

           ENDDO


           

!           CBin(IBin)  = CMode(IMode)*dD(IBin) / &
!                &        (Dp(IBin)*SQRT(2.*PI)*LOG(STD(IMode))) * &
!                &        EXP(-(log(Dp(IBin)/GMD(IMode)))**2 / &
!                &        (2.*(log(STD(IMode)))**2))

           CBin(IBin) = c_kpl
!           write(*,*)  CBin(IBin)
           C(IBin)     = c_kpl
           Dp(IBin)    = (sum_d /c_kpl)**(1./3.)
           Rp(IBin)    = Dp(IBin) / 2.

           Ctot        = Ctot + C(IBin)
           IF(IBin < NBinHi) Vlo(IBin + 1) = Vhi(IBin)
!!$C
!!$C     Convert weight fraction to moles. The value is to be corrected later.
!!$C     Also, the volume of particles in size bin is calculated.
!!$C
           Vtemp(IBin) = 0.0_doubp
           DO ISpec = 1, NALiquids + NASolids
              InC(ISpec,IBin) = (NABins*(ISpec - 1)) + IBin + NABins
              C(InC(ISpec,IBin)) = WtFrac(ISpec,IMode) / WtMol(ISpec) 
              Vtemp(IBin) = Vtemp(IBin) + C(InC(ISpec,IBin)) * WtMol(ISpec) / & 
                   &        Dens(ISpec) * 1.0e-06_doubp
           END DO
!!$C
! Here also for gases (gases has non-bin-structure so the value is in last
! "slot" (NABins))
           DO ISpec = NALiquids + NASolids + 1, NAGases + NALiquids + NASolids
              InC(ISpec,NABins) = (NABins + NABins*NALiquids + NABins*NASolids + &
                   & ISpec - NALiquids - NASolids)
           end DO
             
         
!!$C
!!$C     Convert moles to molar concentrations and normalize the 
!!$C     concentrations so that particles size equals Rp
!!$C     
           VFac         = (4./3.*Pi*Rp(IBin)**3) / Vtemp(IBin)
!!$C               WRITE(*,*) Rp(IBin)
           DO ISpec = 1, NALiquids + NASolids
              C(InC(ISpec,IBin))=C(InC(ISpec,IBin))*VFac*CBin(IBin)
           END DO
        END DO
        NBinLo        = NBinHi + 1
!        DO IBini     = 1,NFBins
!           DO ISpec  = 1, NALiquids + NASolids
!              InCi(ISpec,IBini) = (NFBins*(ISpec - 1)) + IBini + NFBins
!           END DO
!        END DO
!
 !        DO IBini     = 1,NY2
 !          DO ISpec  = 1, NALiquids + NASolids
 !             InCl(ISpec,IBini) = (NY2*(ISpec - 1)) + IBini + NY2
 !          END DO
 !       END DO
        

     END DO
  END IF
!!$C     
  DO I = 1, NABins
     CFrac(I) = C(I) / CTot
  END DO
!!$C
  CSumOld = CTot


!!$C ICE BINS: Distribution for ice bins is made here

  R_Y1(1)=1.0d-8
!  Vrat=(1.e-7/1.e-8)**(1.0/(4.001))!1.4_doubp**(1.0_doubp/3.0_doubp)
!  Vrat2=(3.e-5/1.e-7)**(1.0/(20.001))
!  Vrat3=(2.e-3/3.e-5)**(1.0/(26.001))
  Vrat3=(2.e-3/1.e-8)**(1.0/(NY1+1.e-5))
!  DO IBin = 2, 4+1
!     R_Y1(IBin)= Vrat**(1.-IBIN*.0)* R_Y1(IBin-1)
!  END DO
!  DO IBin = 6, 24+1
!     R_Y1(IBin)= Vrat2**(1.-IBIN*.0)* R_Y1(IBin-1)
!  END DO
!  DO IBin = 26, NY1+1
!     R_Y1(IBin)= Vrat3**(1.-IBIN*.0)* R_Y1(IBin-1)
!  END DO
  DO IBin = 2, NY1+1
     R_Y1(IBin)= Vrat3**(1.-IBIN*.0)* R_Y1(IBin-1)
  END DO

!!$C CLOUD DROPLETS: Distribution for cloud drops bins is made here
  r_Y2(1)=1.e-6
  Vrat=(2.e-3/( r_Y2(1)))**(1.0/(NY2+1e-5))
  DO IBin = 2, NY2+1
     r_Y2(IBin)= Vrat*r_Y2(IBin-1)
  END DO

!!$C GRAUBEL: Distribution for graubel bins is made here
  r_Y3(1)=1.e-5
  Vrat=(5.e-3/( r_Y3(1)))**(1.0/(NY3+1e-5))
  DO IBin = 2, NY3+1
     r_Y3(IBin)= Vrat*r_Y3(IBin-1)
  END DO
  
  !!$C SNOW: Distribution for snow bins is made here
  r_Y4(1)=1.e-5
  Vrat=(5.e-3/( r_Y4(1)))**(1.0/(NY4+1e-5))
  DO IBin = 2, NY3+1
     r_Y4(IBin)= Vrat*r_Y4(IBin-1)
  END DO

  LRAD=NY1 + NY1 * NAerT
  C_Y1(1:LRAD)=0.0
  LRAD=NY2 + NY2 * NAerT
  C_Y2(1:LRAD)=0.0
  LRAD=NY3 + NY3 * NAerT
  C_Y3(1:LRAD)=0.0
  LRAD=NY4 + NY4 * NAerT
  C_Y4(1:LRAD)=0.0
END SUBROUTINE SIZEDIST
