SUBROUTINE EQLB_cle(x)

!
! This subroutine is used to call liquid phase thermodynamics
! when Clegg is used for thermodynamics
!

  use headfile
  use precisi
  implicit none
  integer,PARAMETER :: NCmax=3, NAmax=5, NNmax=1
  real(doubp), dimension(NBins) ::  CSVI,CHNO3,CNH3
  real(doubp), dimension(-NAmax:NCmax) ::  MOLAL(-NAmax:NCmax)
  integer :: L3, ISPEC
  real(doubp) :: CSVITEMP, CHNO3EMP, CNH3TEMP, TEMPERATURE,PH2O
  real(doubp) :: PHNO3, PNH3,SURFACE_TENSION,aw,pwe,pH2Oeq,sum_ions
  real(doubp) :: x

!  write(*,*) IF_CHANGE(1:NABins)
  DO L3  = Lfirst, LAst
!     IF(IF_CHANGE(L3) == 1) THEN
     CSVITEMP = 0.0_doubp
     ISpec = ISO42l
     IF(ISO42l /= 0) CSVITEMP = CSVITEMP + C(InC(ISpec,L3))
     CSVI(L3) = CSVITEMP/(C(InC(IH2Ol,L3))*WtGas(IH2Og)*1.e-3_doubp)
     CHNO3EMP = 0.0_doubp
     ISpec    = INO3l
     IF(INO3l /= 0) CHNO3EMP = CHNO3EMP + C(InC(ISpec,L3))
     CHNO3(L3) = CHNO3EMP/(C(InC(IH2Ol,L3))*WtGas(IH2Og)*1.e-3_doubp)

     CNH3TEMP = 0.0_doubp
     ISpec    = INH3l
     IF(INH3l /= 0) CNH3TEMP = CNH3TEMP + C(InC(ISpec,L3))
     CNH3(L3) = CNH3TEMP/(C(InC(IH2Ol,L3))*WtGas(IH2Og)*1.e-3_doubp)
!     END IF
  end DO

!  write(*,*) IF_CHANGE(1:NABins)
!  write(*,*) Lfirst, last
!  pause
  DO L3 = Lfirst, last
    
!     IF(IF_CHANGE(L3) == 1) THEN
     CALL RADIUS(L3)
     MOLAL(-NAmax:NCmax)=0.0_doubp
     MOLAL(-2) = CSVI(L3)
     MOLAL(-3) = CHNO3(L3)
     MOLAL(2)  = CNH3(L3)
     MOLAL(1)  = 2.0_doubp*CSVI(L3) + CHNO3(L3) - CNH3(L3)
     
     temperature = T
!     if(x > 800)  temperature = T-min(0.1,Rp(L3)**2*1e9)
!     write(*,*) Rp(L3)
     sum_ions=0.0
     DO ISpec=2,NALiquids
        sum_ions= sum_ions+C(InC(Ispec,L3))
 !       write(*,*) namespec(Ispec),C(InC(Ispec,L3))
     END DO
     pwe         = pH2Oeq(temperature) * 101325.0_doubp
!
! If contact freezing is on, it is supposed that particles in mode I_MODE_CONTACT
! does not contain water, so exact calculation of liquid phase thermodynamics not 
! needed.
!
! Also the calculation of thermodynamics for very dilute droplets is not needed.
!
     IF(IF_FREEZ_CONTACT == 1) THEN
        IF (IModeBin(L3) == I_MODE_CONTACT) THEN
           aw =C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
           pH2O = pwe*C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
        ELSEIF(IF_Thermo(L3) == 1) THEN
 
!           CALL THERMO_cle(temperature,MOLAL,pH2O,pHNO3,pNH3,aw)
        ELSE
           aw =C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
           pH2O = pwe*C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
           
        ENDIF
     ELSEIF(IF_Thermo(L3) == 1) THEN
        
!        CALL THERMO_cle(temperature,MOLAL,pH2O,pHNO3,pNH3,aw)
     ELSE
        aw =C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
        pH2O = pwe*C(InC(1,L3))/(C(InC(1,L3))+sum_ions)
        
     ENDIF

     Water_activity(L3)=aw

!
! Add Kelvin effect to pressures     
!     
     call mass_frac(L3)
     PSurf(IH2Og,L3) = pH2O &
          & * EXP(2.0_doubp*surface_tension(T,wtSO4(L3)*100.0_doubp, &
          & wtHNO3(L3)*100.0_doubp)*WtGas(IH2Og)/(Rp(L3)*R*T*1.e6_doubp))
!     C         PSurf(IH2Og,L3) = pH2O
!C     >        * EXP(2.*0.072*WtGas(IH2Og)/(Rp(L3)*R*T*1.d6))
!         WRITE(*,*) PSurf(IH2Og,1),'-'
     IF(IHNO3g /= 0) PSurf(IHNO3g,L3) = pHNO3 &
          & * EXP(2.0_doubp*surface_tension(T,wtSO4(L3)*100.0_doubp, &
          & wtHNO3(L3)*100.0_doubp)*WtGas(IHNO3g)/(Rp(L3)*R*T*1.e6_doubp))
!     C         IF(IHNO3g.NE.0) PSurf(IHNO3g,L3) = pHNO3
!C     >        * EXP(2.*0.072*WtGas(IHNO3g)/(Rp(L3)*R*T*1.d6))
!C         WRITE(*,*) PSurf(INH3g,L3)
     IF(INH3g /= 0) PSurf(INH3g,L3) = pNH3 &
          & * EXP(2.0_doubp*0.072_doubp*WtGas(INH3g)/(Rp(L3)*R*T*1.e6_doubp))
!     END IF
  end DO
END SUBROUTINE EQLB_cle
