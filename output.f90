SUBROUTINE OUTPUT(time,Y)
!!$C
  USE headfile
  use precisi
  implicit none
!!$C      
  REAL(doubp), intent(in) ::  time
  REAL(doubp), DIMENSION(NODEs), intent(in) ::  y
  REAL(doubp) :: RHI, temperature, pw, pwe, phno3, phcl
  REAL(doubp) :: pH2Oeq, so_mass_sum,LWC,r_eff,r2,r3,N100
  integer :: J, IFam, L3, JJ, I, ISpec
!!$C
!!
!! Firs put values from Y back to , C_Y1 and C_Y2
!! 
  IF (time /= 0) THEN
     DO J=1,NYLiquids
        C(IY2C(J))=y(J)
     END DO
!!$C     
     IF (ice_mod == 1) THEN
        DO J=1,NY1
           C_Y1(InCi(IH2Ol,J)) = Y(NYConcs + J)
        END DO
     END IF

     IF (cloud_mod == 1) THEN
        DO J=1,NY2
           C_Y2(InCl(IH2Ol,J)) = Y(NYConcs +NY1+ J)
        END DO
     END IF
!!$C
!! 
!!    
     J = NYLiquids
     DO IFam = 1,NAFl
        IF (NAFMbr(IFam) >= 2) THEN
           DO L3 = 1,NABins
             
              J = J+1
              TotF(IFam) = 0.0e0_doubp
              DO JJ = 1, NFMbr(IFam)
                 IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
                    TotF(IFam) = TotF(IFam) + C(InC(IFMbr(IFam,JJ),L3))
                 END IF
              END DO
              DO JJ = 1, NFMbr(IFam)
                 IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
                    C(InC(IFMbr(IFam,JJ),L3)) = C(InC(IFMbr(IFam,JJ),L3)) / & 
                         & TotF(IFam)*y(J) 
                 END IF
              END DO
           END DO
        END IF
     END DO
!!$C     
     
     DO J = 1, NAGases
        C(NAero+J)=y(J + NYAero)
     END DO
!!$C     
     IF (cloud_mod == 1) c_cl_tot = y(NODEs-6)

     IF (ice_mod == 1) THEN
        C_ice_tot                = y(NODEs-5)
     END IF
     Adia_fex                    = y(NODEs-4)
     z                           = y(NODEs-3)
     P                           = y(NODEs-2)
     Ctot                        = y(NODEs-1)
     T                           = y(NODEs)
!!$C     
!!$C     pw  = WATER VAPOR PRESSURE
!!$C     pwe = WATER SATURATION VAPOR PRESSURE, Seinfeld & Pandis 1998 (page 49)
!!$C     
     pw          = C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp
 
     !  pwe         = pH2Oeq(T) * 101325.0
     temperature=T
     select case(Ithermo_mod)
     case(1) !Luo's thermo
        temperature=T
        CALL s_luo(temperature,0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
        pwe=pwe*100.0_doubp
     case(2,3,4) !jacobson
        pwe         = pH2Oeq(T) * 101325.0_doubp
     end select
     RH          = pw / pwe
!     write(*,*) RH
     CALL sat_vap_press_ice
     RHI=pw/e_sat_i*100.0_doubp
!!$C
     N_drops=0.0
     N100=0.0
     LWC=0.0
     r2=0.0
     r3=0.0
     DO I = 1, NABins
        CALL RADIUS(I)
        CALL RADIUS_DRY(I)
        IF(Rp(I)>2.0e-6)  N_drops=N_drops+C(I)
        IF(Rp(I)>1.0e-6) LWC=LWC+4./3.*Pi*Rp(I)**3*C(I)*1e9
        IF(Rp(I)>1.0e-6) r2=r2+Pi*Rp(I)**2*C(I)
        IF(Rp(I)>1.0e-6) r3=r3+Pi*Rp(I)**3*C(I)
        IF(R_dry(I)>3.5e-8) N100=N100+C(I)
     END DO
     IF(CLOUD_MOD==1) THEN  
        N_drops=N_drops+c_cl_tot
        DO I = 1, NY2
           CALL RADIUS_Y2(I)
           LWC=LWC+4./3.*Pi*R_CL(I)**3*C_Y2(I)*1e9
        END DO
    END IF
!!$C
     IF (ice_mod == 1) THEN
        DO  I = 1,NY1
           so_mass_sum = 0.0_doubp
           IF(C_Y1(I) > 0.0) THEN
              DO ISpec = 1, NALiquids + NASolids
                 so_mass_sum = so_mass_sum + C_Y1(InCi(ISpec,I))*WtMol(ISpec)/C_Y1(I)
              END DO
           ENDIF
           so_mass(I) = so_mass_sum
           CALL RADIUSICE2(I)
        END DO
     END IF
!!$C

     

     
  ELSE 
     pw          = C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp
   
     !  pwe         = pH2Oeq(T) * 101325.0
     select case(Ithermo_mod)
     case(1) !Luo's thermo
        temperature=T
        CALL s_luo(temperature,0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
        pwe=pwe*100.0_doubp
     case(2,3,4) !jacobson
        pwe         = pH2Oeq(T) * 101325.0_doubp
     end select
     RH          = pw / pwe
     CALL sat_vap_press_ice
     RHI=pw/e_sat_i*100.0_doubp

  END IF
!!$C

  IF (ice_mod == 0) THEN
     IF (time == 0) WRITE(*,310)'time (s)','RH','T (K)','RHi'
     WRITE(*,300) time, RH, T, z
  ELSE
     IF (time == 0) WRITE(*,330)'time (s)','RH','T (K)', &
          & 'RHi','NFABins','N_ice (1/cm3)','N_CL_BINS','N_cl(1/cm3)'
     WRITE(*,320) time, RH, T, RHI,INT(NFABins2),'/',INT(NFBins),C_ice_tot,N_CL2,'/',NY2 &
          ,w/0.8 &
          ,N_drops,N100
!     write(*,*) 
  END IF
  !  END IF
!!$C     
300 FORMAT(F8.2,TR6,F8.6,TR4,F7.3,TR2,F10.2)
310 FORMAT(TR2,A8,TR10,A2,TR6,A5,TR7,A5) 
320 FORMAT(F8.2,TR4,F8.6,TR4,F7.3,TR2,F10.2,TR4,I3,A1,I3,TR5,F8.4,TR5,I3,A1,I3,TR3,F8.2,TR3,F9.2,TR3,F9.2)
330 FORMAT(TR2,A8,TR6,A2,TR6,A5,TR7,A5,TR4,A7,TR1,A15,TR1,A10,TR1,A12)
!!$C     
  IF (time == 0) CALL out_dat
!!$C     
459 FORMAT(8001(E18.10))
!  write(*,*) Rp(20)/R_dry(20)
  WRITE(KOUT,459)time,(Rp(I),I=1,NABins),RH,T, &
       &     (C(NAero+I),I=1,NAGases), C(InC(IH2Ol,1:NABins)), &
       &     C(1:NABins),(R_dry(I),I=1,NABins),LWC,N_drops/rho,r3/r2
!!$C
!FORMAT(8001(E14.6))
  IF (ice_mod == 1) THEN
     WRITE(IOUT,459)time,(R_Y1(I),I=1,NY1+1),(C_Y1(i),i=1,NY1), &
          &     (r_Y2(I),I=1,NY2+1),(C_Y2(i),i=1,NY2),RHI,C_ice_tot,sum(c(NABins+1:2*NABins)),C(InC(IH2Ogi,NABins)), &
          &  sum(c_Y1(NY1+1:int(2)*NY1)),sum(c_Y2(NY2+1:int(2)*NY2))
  END IF



END SUBROUTINE OUTPUT

!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
SUBROUTINE out_dat
!!$C
  USE headfile
!!$C
  INTEGER :: RLEN
  CHARACTER(len=7) :: hno3_fr, hno3_ad
  CHARACTER(len=10) :: nuc_method
  CHARACTER(len=8) :: thermo_model
!!$C
  select case(Ithermo_mod)
  case(1) ! Luo's thermodynamics
     thermo_model = 'Luo'
  case(2) !Jacobson
     thermo_model = 'Jacobson'
  case(3)
     thermo_model = 'Joku'
  case(4)
     thermo_model = 'Clegg'
  end select
! 
! time, RH, T, R(nbin), C_gas(ngas), C_w(nbin), C_hno3(nbin)
  RLEN=3*NABins+NAGases+3
! 
  select case (ice_mod)
  case(0) !ice model off
!!$  IF (ice_mod == 0) THEN
!!$C     output initial conditions (number of size bins, aerosol 
!!$C     number concentrations in size bins and dry radius)
     IF (IMModes == 1) THEN
!!$C     modes are inseparable
        WRITE(KOUT,FMT=390)thermo_model,1,NABins
     ELSE
        WRITE(KOUT,FMT=390)thermo_model,NAModes, NBinsPerMode(1:NAModes)
     END IF
!!$C     
     WRITE(KOUT,FMT=400)NAGases,RLEN,TotSec,w,CBin(1:NABins)
!!$C     
     WRITE(KOUT,FMT=410) Rp(1:NABins)
     WRITE(KOUT,'(A1)')'%' ! the last comment line
!!$C     
390  FORMAT('% output.dat'/,'%'/, &
          &             '% Thermodynamical model =     ', A8/, &
          &             '% Number of separable modes = ',I4/, &
          &             '% Number of bins in modes =   ',10(I4,TR1))
400  FORMAT('% Number of gases =     ',I4/, &
          &             '% Record length =       ',I4/,  &! for one time step
          &             '% Simulation time =     ',F15.2/, &
          &             '% Vertical speed =      ',F7.3/, &
          &             '% Number concentrations (#/cm-3) for size bins', &
          &             200(/,'%',5(TR2,E11.5)))
410  FORMAT('% Dry radius (m) for size bins', &
          &             200(/,'%',5(TR2,E11.5)))
!!$  ELSE
  case(1) !ice model on
!!$C     let's check some initial variables which ice model uses
     select case(rate_method)
     case(1) !Classic icenucleation
        nuc_method='Classic' 
     case(2) !Koop's parameterization
        nuc_method='Koop   '
     end select

!!$     IF (rate_method == 1) THEN
!!$        nuc_method='Classic' 
!!$     ELSE IF (rate_method == 2) THEN
!!$        nuc_method='Koop   '
!!$     END IF
!!$C     
     select case(method_HNO3)
     case(1) !hno3: freezes, adsorption : off
        hno3_fr='yes'
        hno3_ad='off'
     case(2)!hno3: gas phase, adsorption : off
        hno3_fr='off'
        hno3_ad='off'
     case(3) !hno3: freezes, adsorption : on
        hno3_fr='yes'
        hno3_ad='yes'
     case(4) !hno3: gas phase, adsorption : on
        hno3_fr='off'
        hno3_ad='yes'
     end select
!!$     IF (method_HNO3 == 1) THEN
!!$        hno3_fr='yes'
!!$        hno3_ad='off'
!!$     ELSE IF (method_HNO3 == 2) THEN
!!$        hno3_fr='off'
!!$        hno3_ad='off'
!!$     ELSE IF (method_HNO3 == 3) THEN
!!$        hno3_fr='yes'
!!$        hno3_ad='yes'
!!$     ELSE IF (method_HNO3 == 4) THEN
!!$        hno3_fr='off'
!!$        hno3_ad='yes'
!!$     END IF
     IF (IMModes == 1) THEN
!!$C     modes are inseparable
        WRITE(KOUT,FMT=391)thermo_model,1,NABins
     ELSE
        WRITE(KOUT,FMT=391)thermo_model,NAModes, NBinsPerMode(1:NAModes)
     END IF
!!$C     
     WRITE(KOUT,FMT=401)NAGases,RLEN,TotSec,w, &
          &             nuc_method,NFBins,hno3_fr,hno3_ad,CBin(1:NABins)
!!$C     
     WRITE(KOUT,FMT=411) Rp(1:NABins)
     WRITE(KOUT,'(A1)')'%' ! the last comment line
!!$C     
391  FORMAT('% output.dat'/,'%'/, &
          &             '% Thermodynamical model =     ', A8/, &
          &             '% Number of separable modes = ',I4/, &
          &             '% Number of bins in modes =   ',10(I4,TR1))
401  FORMAT('% Number of gases =     ',I4/, &
          &             '% Record length =       ',I4/,  &! for one time step
          &             '% Simulation time =     ',F15.2/, &
          &             '% Vertical speed =      ',F7.3/, &
          &             '% Ice nucleation =      ',A7/, &
          &             '% Number of ice bins =  ',I4/, &
          &             '% HNO3 freez =          ',A3/, &
          &             '% HNO3 adsorption =     ',A3/, &
          &             '% Number concentrations (#/cm-3) for size bins', &
          &             200(/,'%',5(TR2,E11.5)))
411  FORMAT('% Dry radius (m) for size bins', &
          &             200(/,'%',5(TR2,E11.5)))
  end select
!!$C  END IF
!!$C
END SUBROUTINE out_dat
