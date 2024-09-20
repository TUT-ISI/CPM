FUNCTION pH2Oeq(T)
  use precisi
  implicit none
!!$C     
!!$C     Below 2oC an expression for pH2O0 based on Goff-Gratch is used
!!$C     but modified so that the heat capacities of (supercooled) pure
!!$C     water are those implied by an extrapolation of Hill's equation
!!$C     of state. Above 2oC, an empirical fit to pH2O0 given in the CRC 
!!$C     Handbook is used, valid to 100oC. At the 'join', the fit gives 
!!$C     6.96289E-3 atm, and the lower temperature expression yields 
!!$C     6.96157E-3 atm, a difference of 0.019%.
!!$C     
!!$C     
  real(doubp), intent(in) :: T
  real(doubp) :: LNK1, LNK2, R, Agas, Bgas, Cgas, Dgas, Egas
  real(doubp) :: Aliq, Bliq, Cliq, Dliq, Eliq, pH2Oeq, delHr
  real(doubp) :: a, DUM, T1, T2, delA, delB, delC, delD, delE
  integer :: METHOD
  LNK1 = 1.809429005_doubp
  delHr = 45.077441e3_doubp
  R = 8.3144_doubp 
  Agas = 33.269811_doubp
  Bgas = 0.00113261_doubp
  Cgas = -1.09982e-5_doubp
  Dgas = 3.573575e-8_doubp
  Egas = 0.0_doubp
!!$  C     Hill-based liquid phase heat capacities:
  Aliq = 295.1612_doubp
  Bliq = -1.540498_doubp
  Cliq = 2.7023e-3_doubp
  Dliq = 0.0_doubp
  Eliq = 0.0_doubp
!!$  C     
  METHOD = 1
!!$  C     
  IF(METHOD == 1) THEN
!!$     C     
     IF(T <= 275.1) THEN
!!$        C     ..below 2 degrees C, used the Goff-Gratch based expression:
!!$        C     
        delA=Agas - Aliq
        delB=Bgas - Bliq
        delC=Cgas - Cliq
        delD=Dgas - Dliq
        delE=Egas - Eliq
        
        T2=T
        T1=273.15 
        DUM=delA/R*LOG(T2/T1) + delB/(2.*R)*(T2-T1) &
             &          + delC/(6.*R)*(T2**2-T1**2) &
             &          + delD/(12.*R)*(T2**3-T1**3) &
             &          + delE/(2.*R)*(1.D0/T2**2-1.D0/T1**2) &
             &          + (1./R)*(-delHr+delA*T1+delB/2.*T1**2+delC/3.*T1**3 &
             &          + delD/4.*T1**4 - delE/T1)*(1./T2 - 1./T1)
!!$C     
        LNK2 = DUM + LNK1
        pH2Oeq = EXP(LNK2) * 0.000986923_doubp
!!$C     ^- converts to atm.
!!$C     
     ELSEIF(T > 275.1 .AND. T <= 275.2) THEN
         delA=Agas - Aliq
        delB=Bgas - Bliq
        delC=Cgas - Cliq
        delD=Dgas - Dliq
        delE=Egas - Eliq
        
        T2=T
        T1=273.15 
        DUM=delA/R*LOG(T2/T1) + delB/(2.*R)*(T2-T1) &
             &          + delC/(6.*R)*(T2**2-T1**2) &
             &          + delD/(12.*R)*(T2**3-T1**3) &
             &          + delE/(2.*R)*(1.D0/T2**2-1.D0/T1**2) &
             &          + (1./R)*(-delHr+delA*T1+delB/2.*T1**2+delC/3.*T1**3 &
             &          + delD/4.*T1**4 - delE/T1)*(1./T2 - 1./T1)
!!$C     
        LNK2 = DUM + LNK1
        pH2Oeq = (EXP(LNK2) * 0.000986923_doubp)*(275.2_doubp-T)/0.1_doubp + &
             &  (EXP( 23.54872_doubp  - 6459.987931_doubp /T - 0.022791752_doubp *T & 
             & + 1.6290826e-5_doubp*T**2 ))*(T-275.1_doubp)/0.1_doubp

     ELSE
!!$C     ..else use empirical fit to CRC tabulated vapour pressures:
!!$C     
        DUM=23.54872_doubp  - 6459.987931_doubp /T - 0.022791752_doubp *T & 
             & + 1.6290826e-5_doubp*T**2
        pH2Oeq = EXP(DUM)
     END IF
  END IF
!!$C
!!$C
!!$C
  IF(METHOD == 2) THEN
     a           = 1.0-(373.15_doubp/T)
     pH2Oeq      = exp(13.3185_doubp*a-1.976_doubp*a**2- &
          &        0.6445_doubp*a**3-0.1299_doubp*a**4)
  END IF
END FUNCTION pH2Oeq

