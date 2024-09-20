
MODULE funs

contains

  SUBROUTINE fex(n,x,y,dydx)
!!$     
    USE headfile
    use precisi
    use aikaa_mittaavat_parametrit
    use jacobi
    implicit none
!!$     
    INTEGER, intent(in) :: n
    REAL(doubp), intent(in), DIMENSION(n) :: y
    REAL(doubp), intent(in) :: x
    REAL(doubp), intent(out), DIMENSION(n) :: dydx
    REAL(doubp), DIMENSION(NBins + NBins * NAerT + NGases) :: dCdx
    REAL(doubp) :: PSurfice, CSurfice, Sum_ice, Tot_Wat_ice, helpa
    REAL(doubp) :: sum_ads, K_eq, CpTot, Condice
    REAL(doubp) :: Cond, Sum1, Tot_Wat, Heat, CSurf, PMol, rho_ice
    REAL(doubp) :: RMol, CAQ, TotC, Rate, pw, pwe!, Tnow, Told_fex
    REAL(doubp) :: phno3,phcl, temperature, pH2Oeq, surfE_ice,T_tmp,sum_ions
    REAL(doubp) :: R2T3CK
    REAL(doubp), DIMENSION(NABins,NAGases) :: cond_tmp
    REAL(doubp), DIMENSION(n) :: Y2          ! All positive Y values
    INTEGER :: II, IIOld, Ji, I_idx, J_J, Idx, Jgas_i, I_I, JGas, JGas1
    INTEGER :: IC, NR, Iice, ISpec, I, LWat, JJ, L3, IFam, J, Ibin
    REAL(doubp) ::SURFACE_TENSION,TR,help_term
!    INTEGER :: LFirst, Last

!!$ *****************************************************************
!!$ *   SUBROUTINE FEX DEFINES THE DIFFERENTIAL EQUATIONS:          *
!!$ *   Dy/DX = DISSOLUTION + !!$HEMICAL PRODUCTION - CHEMICAL LOSS *
!!$ *   Dy/DX = TEMPERATURE CHANGE                                  *
!!$ *****************************************************************
!!$

!!$
!!    write(*,*) x,n
    dydx(1:n)=0.0_doubp
!    IF_Change(1:NaBins)=0
!!$
!!$     
    LFirst = 0
    Last = 0
    IIOld = 1000000
!!$
!!$CCC  HERE SOME CHECKING IS MADE FOR NEGATIVE VALUES IN CONCENTRATIONS
!!
!!  
! This is here as 0 values are causing problems. Thing need to be solved....
!
!    IF(N_C4==0) Y(N-8)=max(1.0d-20,Y(N-8))
!    IF(N_C3==0) Y(N-7)=max(1.0d-20,Y(N-7))
!    IF(cloud_mod==0) Y(N-6)=max(1.0d-20,Y(N-6))
    Y2=Y
    do  II=1,N-NDyna
       if (Y(II) < 0.0_doubp) then
          
          IF(II.NE.NODES) THEN 
             write(*,*) II, 'SOMETHING NEGATIVE IN Y, subroutine fex',Yold(II),NODES,x
!             write(*,*) Y(II-1),Y(II),Y(II+1)
             !             pause
             !             Y(II)=Yold(II)*0.01
             Y2(II)=Yold(II)*0.01
            pause
          END IF
          
          
       else
          Y2(II)=Y(II)
       endif
    end do
!
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Here it is decided which concentrations changes.
! Call to thermodynamic is needed only if temperature or some liquid
! phase concentration changes. At the moment this is quite unclear 
! and needs to be rewritten and optimized. Currrently might be wrong in same cases.
! 
! 
!
    DO II = 1,NYAero + NAGases
       IF(yOld(II) /= y(II)) THEN
          IF (II < IIOld) THEN
             LFirst = II-NABins*((II-1)/NABins)
             IIOld = II
            
          END IF
          IF (II > IIOld)  Last = II - NABins*((II-1)/NABins)
!          IF_CHANGE(II) = 1 
       END IF
    END DO
!!$
    IF (ice_mod == 1 .AND. cloud_mod ==1) THEN
       DO II = NYAero+NAGases+NFABins2+N_CL2+1, NYAero + NAGases + NFABins2 +N_CL2+ NDyna
          IF(yOld(II) /= y(II)) THEN
             LFirst = 1
             Last = NABins
!             IF_CHANGE(1:NABins) = 1 
          END IF
       END DO
    ELSEIF (ice_mod == 1) THEN
       DO II = (NYAero + NAGases + NFABins2 + 1),NYAero + NAGases + NFABins2 + NDyna
          IF(yOld(II) /= y(II)) THEN
             IF (II < IIOld) THEN
                LFirst = II-NABins*((II-1)/NABins)
                IIOld = II
             END IF
             IF (II > IIOld)  Last = II - NABins*((II-1)/NABins)
          END IF
       END DO
    ELSE
       DO II = (NYAero + NAGases + 1),NYAero + NAGases + NDyna
          IF(yOld(II) /= y(II)) THEN
             IF (II < IIOld) THEN
                LFirst = II-NABins*((II-1)/NABins)
                IIOld = II
             END IF
             IF (II > IIOld)  Last = II - NABins*((II-1)/NABins)
          END IF
       END DO
    END IF
!!$
!!$     
    DO II = NYAero+1,NYAero + NAGases
       IF(yOld(II) /= y(II)) THEN
          LFirst = 1
          Last = NABins
!          IF_CHANGE(1:NABins) = 1 
       END IF
    END DO

    
     DO II = NODEs-Ndyna+1,NODEs
          IF(yOld(II) /= y(II)) THEN
             LFirst = 1
             Last = NABins
          END IF
       END DO


!!$      
    IF (Last == 0) Last = LFirst
    IF (LFirst > Last) THEN
       LFirst = 1
       Last = NABins
!       IF_CHANGE(1:NABins) = 1 
    END IF
!!$  
    IF (NODEs /= NODEs_old) THEN
       LFirst = 1
       Last = NABins
!       IF_CHANGE(1:NABins) = 1 
    END IF
!

!    write(*,*)  LFirst,Last

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  PLACE CONCENTRATIONS FROM Y TO C AND C_Y1 AND ALSO C_Y2

    NODES_old=NODEs
!!$ 
    DO J=1,NYLiquids
        C(IY2C(J))=y2(J)
    END DO
!!$     
    DO J = 1, NAGases
       C(NAero+J)=y2(J + NYAero)
    END DO
!!$
    IF (ice_mod == 1) THEN
       DO J=1,NFABins2
          C_Y1(InCi(IH2Ol,Y1C(J)-NYConcs)) = Y2(NYConcs + J)
       END DO
    END IF
!!$    
    IF (cloud_mod == 1) THEN
       DO J=1,N_CL2
          C_Y2(InCl(IH2Ol,Y2C(J)-NYConcs-NY1)) = Y2(NYConcs+NFABins2+J)
       END DO
    END IF

!    IF (N_C3 > 1) THEN
!       DO J=1,NY3
!          C_Y3(InC3(IH2Ol,Y3C(J)-NYConcs-NY1-NY2)) = Y(NYConcs+NFABins2+NY2+J)
!       END DO
!    END IF
!
!    IF (N_C4 > 1) THEN
!       DO J=1,NY4
!          C_Y4(InC4(IH2Ol,Y4C(J)-NYConcs-NY1-NY2-NY3)) = Y(NYConcs+NFABins2+NY2+NY3+J)
!       END DO
!    END IF

!!$     
!!$  CONCENTRATIONS FROM FAMILIES (EG SULFATE HAS THREE FORMS: H2SO4(L), HSO4-, and SO4-- WITH ONLY
!!$  TOTAL SUM SOLVED WITH DIFFERENTIAL EQUATIONS AND LIQUID PHASE EQUILIBRIUM IN THERMODYNAMICS).
!!   
!!   I'M NOT SURE THIS SYSTEM IS THOROUGHLY TESTED, AND WITH ONLY WATER NOT NEEDED --Sami--
!!
!    J = NYLiquids
!    DO IFam = 1,NAFl
!       IF (NAFMbr(IFam) >= 2) THEN
!          DO L3 = 1,NABins
!             J = J+1
!             TotF(IFam) = 0.0_doubp
!             DO JJ = 1, NFMbr(IFam)
!                IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
!                   TotF(IFam) = TotF(IFam) + C(InC(IFMbr(IFam,JJ),L3))
!                END IF
!             END DO
!             DO JJ = 1, NFMbr(IFam)
!                IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
!                   C(InC(IFMbr(IFam,JJ),L3)) = C(InC(IFMbr(IFam,JJ),L3)) / &
!                        & TotF(IFam)*y(J) 
!                END IF
!             END DO
!          END DO
!       END IF
!    END DO
!!$     
!!$
    IF (ice_mod == 1) THEN
       C_ice_tot   = y2(NODEs-5)
    END IF
    IF (cloud_mod == 1 ) c_cl_tot    = y2(NODEs-6)
    Adia_fex    = y2(NODEs-4)
    z           = y2(NODEs-3)
    P           = y2(NODEs-2)
    Ctot        = y2(NODEs-1)
    T           = y2(NODEs)
    TR=T*R
!!$     
!!$     pw  = WATER VAPOR PRESSURE
!!$     pwe = WATER SATURATION VAPOR PRESSURE, Seinfeld & Pandis 1998 (page 49)
!!$     
    pw          = C(InC(IH2Ogi,NABins)) * TR * 1.0e6_doubp
    select case(Ithermo_mod)
    case(1) !Luo's thermo
       temperature=T
       CALL s_luo(temperature,0.0_doubp,0.0_doubp,0.0_doubp, &
            & 0.0_doubp,pwe,phno3,phcl)
       pwe=pwe*100.0_doubp
    case(2,3,4) !jacobson
       pwe         = pH2Oeq(T) * 101325.0_doubp
    end select
    !     write(*,*) pw
    RH          = pw / pwe
!!$     
!! NUMBER CONCENTRATIONS IN LIQUID BINS
!!$
!    DO L3   = 1, NABins
!       C(L3)    = CTot * CFrac(L3)
    C(1:NaBins)    = CTot * CFrac(1:NaBins)
!    END DO
!!$
!  Number concentration in ice bins
!!$
    IF (ice_mod == 1) THEN
       IF (NFABins2 >= 1) THEN
!          DO I=1,NY1
!             C_i(I)=C_ice_tot*C_iceFrac(I)
!             C_Y1(I)=C_ice_tot*C_iceFrac(I)             
          C_Y1(1:NY1)=C_ice_tot*C_iceFrac(1:NY1)
!          END DO
       END IF
    END IF
!!$
!  Number concentration in cl bins
!!$
    IF (cloud_mod == 1) THEN
       C_Y2(1:NY2)=C_cl_tot*C_clFrac(1:NY2)  
    END IF
!!$
!  ION CONCENTRATIONS FOR CONCENTRATIONS NOT SOLVED
!!$     
    IF (Adia_fex - Adia_fex_old /= 0.0) THEN
       DO ISpec = 1, NALiquids + NASolids
          IF (IfInY(ISpec) == 0 .AND. InFamily(ISpec) == 0) THEN
             DO I = 1,NABins
                C(InC(ISpec,I)) = C(InC(ISpec,I))* &
                     & Adia_fex/Adia_fex_old
             END DO
          END IF
       END DO
       IF (ice_mod == 1) THEN
!!$     For ice
          DO Iice=1,NY1
             C_Y1(Iice) = C_Y1(Iice)*Adia_fex/Adia_fex_old
             DO ISpec=2,NALiquids + NASolids !muutos
                C_Y1(InCi(ISPec,Iice))=C_Y1(InCi(ISPec,Iice))* &
                     & Adia_fex/Adia_fex_old
             END DO
          END DO
       END IF
       
       IF (cloud_mod == 1) THEN
!!$     For cl
          DO ISpec = 1, NALiquids + NASolids
             IF (IfInY(ISpec) == 0 .AND. InFamily(ISpec) == 0) THEN
                DO I = 1,NY2
                   C_Y2(InCl(ISpec,I)) = C_Y2(InCl(ISpec,I))* &
                        & Adia_fex/Adia_fex_old
                END DO
             END IF
          END DO
       END IF

!!!!! NOTE!!!!!
! here should be similar system added for C_Y3 and C_Y4.   
!         -ADDED, but not tested-
!!!!!!!!!!!!!!!

!       if (N_C3>0) then
!          DO ISpec = 1, NALiquids + NASolids
!             IF (IfInY(ISpec) == 0 .AND. InFamily(ISpec) == 0) THEN
!                DO I = 1,NY3
!                   C_Y3(InC3(ISpec,I)) = C_Y3(InC3(ISpec,I))* &
!                        & Adia_fex/Adia_fex_old
!                END DO
!             END IF
!          END DO
!       endif
       
!       if (N_C4>0) then
!          DO ISpec = 1, NALiquids + NASolids
!             IF (IfInY(ISpec) == 0 .AND. InFamily(ISpec) == 0) THEN
!                DO I = 1,NY4
!                   C_Y4(InC4(ISpec,I)) = C_Y4(InC4(ISpec,I))* &
!                        & Adia_fex/Adia_fex_old
!                END DO
!             END IF
!          END DO
!       endif


    END IF




!!$     
!!$
    CSumOld = CTot
    Adia_fex_old=Adia_fex
!!$     
!  IF(Last /= 0)  CALL CMODEL_EQUILIBRIUM ! CALC_EQUILIBRIUM(.FALSE.)
!
!  CALL THERMODYNAMICS
!

!    CALL CPU_TIME(TIME_IN)


!    select case(Ithermo_mod)
!    case(1) !Luo thermo
!       IF(Last /= 0) CALL EQLB !(LFirst,Last)
!    case(2) !Jacobson thermo
!       IF(Last /= 0) CALL AERPROC2 !(LFirst,Last)
!    case(3) !
!       IF(Last /= 0) CALL CMODEL_EQUILIBRIUM !(LFirst,Last)
!    case(4) !
    IF(Last /= 0) CALL EQLB_cle(x) !(LFirst,Last)
!    end select
!    CALL CPU_TIME(TIME_OUT)
!    TOTAL_EQMODEL_TIME=TOTAL_EQMODEL_TIME+(TIME_OUT-TIME_IN)
!    EQ_CALCULATIONS=EQ_CALCULATIONS+1

    CALL DIFFCOEF
!!$     
!!$     Calculation of altitude and pressure added here because chance 
!!$     of total pressure needed later... 
!!$     
!!$     dz/dt...
!!$     
!    IF (x.GT.800.0) w = 0.0
    dydx(NODEs-3) = w    
   
!!$     
!!$     dP/dt...
!!$      
    dydx(NODEs-2)  = - P * WtAir * 1.0e-3_doubp * G /(TR)*dydx(NODEs - 3)
!!$     
!!$
!!$
!!$                  START CONDENSATION          
!!$
!!$   Calculate differential equations with 
!!$   condensation
!!$   
!!$   H_Term takes into account temperature differecne between gas and droplet.
!!$   Calculation is only approximative...
!!$   
!!$   
    !    CSurf=0.0_doubp ! just for initializing

    R2T3CK=R**2 * T**3 * CKappa * 1.0e-2_doubp
    DO I = 1,NABins 
       HTerm(I) = 0.e0_doubp
!       write(*,*) I,IF_Cond(I)
!       pause
       IF(IF_Cond(I) ==1) THEN
          DO JGas  = 1,NAGases
             JGas1     = JGas + NAero
             CSurf     = PSurf(JGas,I) / (TR) * 1.0e-6_doubp
             HTerm(I)  = HTerm(I)+HeatL(JGas)*WtGas(JGas)*D(JGas) &
                  &    * TR * (C(JGas1) - CSurf)/ &
                  &    (R2T3CK * BT(JGas,I) &
                  &    + WtGas(JGas) * WtGas(JGas) * HeatL(JGas) * &
                  &    HeatL(JGas) * D(JGas) * CSurf * TR)
!!$   
          END DO
       ENDIF
    END DO
!!$C
    IF (iffreez == 1) THEN
       help_term=HeatLice*WtGas(1)* D(1)* TR
!       CALL sat_vap_press_ice
       PSurfice = e_sat_i !!*1.d2
       CSurfice  = PSurfice / (TR) * 1.0e-6_doubp
       DO J_J = 1,NY1
          IF (Y1C(J_J) > 0) THEN
             I_I=Y1C(J_J)-NYConcs
             HTerm_ice(I_I)=0.e0_doubp
             DO JGas = 1,1
                Jgas_i=JGas + NAero 
                HTerm_ice(I_I) = HTerm_ice(I_I)+help_term* &
                     &         (C(JGas_i) - CSurfice)/ &
                     &         (R2T3CK * BTice(JGas,I_I) &
                     &         + WtGas(JGas) * WtGas(JGas) * HeatLice * &
                     &         HeatLice * D(JGas) * CSurf * TR)
             END DO
          END IF
       END DO
    END IF
!!$C
    Heat             = 0.0_doubp
    kaasut:  DO JGas  = 1, NAGases
!!$   
       IF(Diss(JGas) == 'none') CYCLE kaasut
       JGas1         = JGas + NAero
!!$   
       Tot_Wat       = 0.0_doubp
       Sum1          = 0.0_doubp
       Sum_ice       = 0.0_doubp
!!$   
!!$   Check if the gas belongs to a family
!!$   
       IF(InFamily(NALiquids + NASolids + JGas) == 1) THEN
          DO J   = 1,NAFl
             IF(IFMbr(J,1) == (NALiquids + NASolids + JGas)) THEN
                IFam=J
             END IF
          END DO
       END IF
!!$   
       DO I      = 1,NABins 
          Cond=0.0
          IF(IF_Cond(I) == 1) THEN
!!$   
!!$   CSurf   = GAS PHASE CONCENTRATION AT THE DROPLET SURFACE #/CM**3
!!$   Tot_Wat = TOTAL AMOUNT OF LIQUID WATER MOL/CM**3
!!$   
             CSurf      = PSurf(JGas,I) / (TR) * 1.e-6_doubp

             Idx        = IC2Y(IdxBin(JGas,I))
!!$   
             IF(InFamily(NALiquids + NASolids + JGas) == 1) THEN
                Idx     = NYLiquids + NABins * (IFam-1) + I
             END IF

!!$   CONDENSATION RATE
!!$   
             Cond = FourPi*C(I)*Rp(I)*1.e+02_doubp*D(JGas)*Bm(JGas,I)* &
                  & (C(JGas1) - CSurf*(1.0_doubp + WtGas(JGas)*HeatL(JGas)*HTerm(I)))

            
!!$   SUM1 IS EQUIVALENT TO SUM IN EQ (19.48), Jacobson 2000
!!$   
             Sum1 = Sum1 + Cond
!!$   
!!$   dydx is equivalent to eq (19.40), Jacobson 2000
!!$   dydx(Idx) already includes chemical loss and production
!!$   
!!!!$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$$
             dydx(Idx)  = Cond + dydx(Idx)
!!!!$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!$$
          END IF
          cond_tmp(I,JGas)=Cond
       END DO

!  
!  CONDENSATION INTO CLOUD DROPLETS IN C_Y2
!
       IF(CLOUD_MOD == 1) THEN
          
           DO J_J = 1,NY2
              sum_ions=0.0
              IF (Y2C(J_J) > 0) THEN
                 I_I=Y2C(J_J)-NYConcs-NY1
!
! THIS IS SIMPLIFIED TOO MUCH, CORRECT IF NEEDED
!               
                 DO ISpec=2,NALiquids
                    sum_ions= sum_ions+C_Y2(InCl(Ispec,I_I))
                 END DO



!                 CALL RADIUS_Y2(I_I)
                 Water_activity(NaBins+I_I)=C_Y2(InCl(1,I_I))/(C_Y2(InCl(1,I_I))+sum_ions)
!                      &   /EXP(2.0_doubp*surface_tension(T,0.0, &
!                      & 0.0)*WtGas(IH2Og)/(R_CL(I_I)*R*T*1.e6_doubp))

                 CSurf = pwe*C_Y2(InCl(1,I_I))/(C_Y2(InCl(1,I_I))+sum_ions)/ (TR) * 1.e-6_doubp* &
                      &  EXP(2.0_doubp*surface_tension(T,0.0, &
                      & 0.0)*WtGas(IH2Og)/(R_CL(I_I)*R*T*1.e6_doubp))


                 HTermCL(I_I)  = HeatL(JGas)*WtGas(JGas)*D(JGas) &
                      &    * TR * (C(JGas1) - CSurf)/ &
                      &    (R2T3CK * BT(JGas,NaBins+I_I) &
                      &    + WtGas(JGas) * WtGas(JGas) * HeatL(JGas) * &
                      &    HeatL(JGas) * D(JGas) * CSurf * TR) 
                 

                 Cond = FourPi*C_Y2(I_I)*R_CL(I_I)*1.e+02_doubp*D(JGas)*Bm(JGas,NaBins+I_I)* &
                      & (C(JGas1) - CSurf*(1.0_doubp + WtGas(JGas)*HeatL(JGas)*HTermCL(I_I)))


                 I_idx = NYConcs+NFABINS2+J_J

                 dydx(I_idx)  = Cond + dydx(I_Idx)

                 Sum1 = Sum1 + Cond
              ENDIF
           ENDDO
       ENDIF

       
       if (nfabins2>0) then
          surfE_ice=-0.25_doubp*T+168.2875_doubp !ergs/cm2
          DO J_J = 1,NY1
             IF (Y1C(J_J) > 0) THEN
                I_I=Y1C(J_J)-NYConcs
              
!                CALL radiusice2(I_I)

!!!!!  NOTE HERE KELVIN EFFECT IS COMMENTED 
!!!!!  IT SHOULD NOT BE              

                PSurfice = e_sat_i!*exp(4.0_doubp*WtMol(IH2Ol)*surfE_ice*1.0e-7_doubp/&
!                     & (R*T*rho_ice(T)*1.0e-3_doubp*R_ice(I_I)*100.0_doubp))
                CSurfice  = PSurfice / (TR) * 1.e-6_doubp
                Condice = FourPi*C_Y1(I_I)*R_ice(I_I)*1.e+02_doubp*D(JGas)* &
                     & BmIce(JGas,I_I)*(C(JGas1) - CSurfice &
                     & *(1.0_doubp + WtGas(JGas)*HeatLice*HTerm_ice(I_I)))
                Sum_ice = Sum_ice + Condice
                I_idx = NYConcs+J_J

                dydx(I_idx) = Condice + dydx(I_idx)
             ENDIF
          END DO
       endif

!        if (N_C3>0) then
!           surfE_ice=-0.25_doubp*T+168.2875_doubp !ergs/cm2
!
! NOTE: CONDENSATION EQUATIONS & OTHER NEED TO CORRECTED
!
!          DO J_J = 1,NY3
!             IF (Y3C(J_J) > 0) THEN
!                I_I=Y3C(J_J)-NYConcs-NY1-NY2
!                CALL sat_vap_press_ice            
!                CALL radiusice3(I_I)

!!!!!  NOTE HERE KELVIN EFFECT IS COMMENTED 
!!!!!  IT SHOULD NOT BE              

!                PSurfice = e_sat_i!*exp(4.0_doubp*WtMol(IH2Ol)*surfE_ice*1.0e-7_doubp/&
!                     & (R*T*rho_ice(T)*1.0e-3_doubp*R_ice(I_I)*100.0_doubp))
!                CSurfice  = PSurfice / (TR) * 1.e-6_doubp
!                Condice = FourPi*C_Y3(I_I)*R_C3(I_I)*1.e+02_doubp*D(JGas)* &
!                     & BmIce(JGas,I_I)*(C(JGas1) - CSurfice &
!                     & *(1.0_doubp + WtGas(JGas)*HeatLice*HTerm_ice(I_I)))
!                Sum_ice = Sum_ice + Condice
!                I_idx = NYConcs+NFABINS2+N_CL2+J_J
!
!                dydx(I_idx) = Condice + dydx(I_idx)
!             ENDIF
!          END DO
!       endif
       

!       if (N_C4>0) then
!          surfE_ice=-0.25_doubp*T+168.2875_doubp !ergs/cm2
          
!
! NOTE: CONDENSATION EQUATIONS & OTHER NEED TO CORRECTED
!
          
!          DO J_J = 1,NY4
!             IF (Y4C(J_J) > 0) THEN
!                I_I=Y4C(J_J)-NYConcs-NY1-NY2-NY3
!                CALL sat_vap_press_ice
!                CALL radiusice4(I_I)

!!!!!  NOTE HERE KELVIN EFFECT IS COMMENTED 
!!!!!  IT SHOULD NOT BE              
!
!                PSurfice = e_sat_i!*exp(4.0_doubp*WtMol(IH2Ol)*surfE_ice*1.0e-7_doubp/&
!                     & (R*T*rho_ice(T)*1.0e-3_doubp*R_ice(I_I)*100.0_doubp))
!                CSurfice  = PSurfice / (TR) * 1.e-6_doubp
!                Condice = FourPi*C_Y4(I_I)*R_C4(I_I)*1.e+02_doubp*D(JGas)* &
!                     & BmIce(JGas,I_I)*(C(JGas1) - CSurfice &
!                     & *(1.0_doubp + WtGas(JGas)*HeatLice*HTerm_ice(I_I)))
!                Sum_ice = Sum_ice + Condice
!                I_idx = NYConcs+NFABINS2+N_CL2+N_C3+J_J
!
!                dydx(I_idx) = Condice + dydx(I_idx)
!             ENDIF
!          END DO
!       endif


!!!!!!!!!!! 881     continue
!!$   
!!$   
!!$   DIFFERENTIAL EQUATION FOR GAS PHASE CONCENTRATION 
!!$   EQ (17.55), Jacobson 2000
!!$   

       
       IF (ice_mod == 1) THEN
          IF (JGas == 1) THEN
             dydx(NYAero+JGas)=dydx(NYAero+JGas) - Sum1 - Sum_ice
             Heat = Heat + HeatL(JGas) * Sum1 * WtGas(JGas) + &
                  & HeatLice * Sum_ice * WtGas(JGas)
          ELSE
             dydx(NYAero+JGas)=dydx(NYAero+JGas) - Sum1 - sum_ads
             Heat = Heat + HeatL(JGas) * Sum1 * WtGas(JGas)
          END IF
       ELSE
          dydx(NYAero+JGas)=dydx(NYAero+JGas) - Sum1
          Heat = Heat + HeatL(JGas) * Sum1 * WtGas(JGas)
       END IF
!!$
!!$   CONSTANT PRODUCTION OF GAS PHASE SPECIES (PPB/DAY)
!!$   
!!$   VMIXR = GASPROD(JGas)
!!$     c         IF (X.GT.((finmin-10.)*60.)) VMIXR = 0.0
!!$   dydx(NYAero+JGas) = dydx(NYAero+JGas)+
!!$   1        VMIXR/(24.*60.*60.)*101325./(R*T)/1.d6
!!$   
!!$   Calculation for heat produced by condensation
!!$   
!!$       if (Jgas .eq. 1 .and. ice_mod .eq. 1) then
!!$          Heat = Heat + HeatL(JGas) * Sum1 * WtGas(JGas) + 
!!$   >           HeatLice * Sum_ice * WtGas(JGas)
!!$       else
!!$          Heat = Heat + HeatL(JGas) * Sum1 * WtGas(JGas)
!!$       end if
!!$   
    END DO kaasut
!!$   
!!$   DIFFERENTIAL EQUATION FOR TEMPERATURE
!!$   NOTE: Specific heats are C_p values...
!!$   
!!$   Cptot = (m(dry air)*Cp(dry air) + sum(m(other gases)*Cp(other)) +
!!$   sum(m(liquids)*Cp(liquids))) / mtot
!!$   
!!$   Here Cptot is only approximated
!!$    GOTO 596
!!$   


! 
! Calculate the dropelt temperature for freezing module. This is 
! causing promlems so it is not included now.
!

    IF (IfFreez == 1) THEN
       bins: do IBin=1,NABins
          T_tmp = 0.0_doubp
          gases: do JGas=1,NAGases
             T_tmp = T_tmp + (HeatL(JGas)*cond_tmp(IBin,JGas)/C(IBin) * &
                  & WtMol(JGas))/(FourPi*BT(JGas,IBin)*Rp(IBin) * &
                  & CKappa)
          end do gases

          T_drop(IBin) = T_tmp+T
       end do bins
 
! Add temperature calculation like above    
       IF(cloud_mod ==1) THEN
          T_drop(NaBins+1:NABins+NY2) = T
       ENDIF
    ENDIF

!    IF(IF_NUC_FEX==1) THEN
!       CALL freez_prob(1)
!       do IBin=1,NABins
!          IF(prob(IBin)>0) dydx(NOdes-Ndyna)=dydx(NOdes-Ndyna)+4./3.*pi*((Rp(IBin)**3))*z_rate(IBin)*C(Ibin)
!       end do
!       IF (cloud_mod == 1) THEN
!          drops: do IBin = 1,NY2
!             IF(C_Y2(Ibin)<1e-10) cycle drops
!             dydx(NOdes-Ndyna)=dydx(NOdes-Ndyna)+4./3.*pi*((R_Y2(IBin)**3))*z_rate(Ibin+NaBins)*C_Y2(Ibin)
!          end do drops
!       ENDIF
!    ENDIF

    
 

    !    if (w /= 0.0_doubp) then
    Tot_wat=SUM(C(Inc(IH2Ol,1:NABins)))
    IF(NFBins>0) Tot_wat_ice=SUM(C_Y1(InCi(IH2Ol,1:NFBins)))

!    cp_ice= -2.0572_doubp+0.14644_doubp*T+0.06163_doubp*T*exp(-(T/125.1_doubp)**2)

! TOTAL HEAT CAPACITY

    Cptot =  (C(InC(IH2Ogi,NABins)) * WtGas(IH2Og) * 1.95_doubp + &! water gas
         &   Rho * CpAir * 1.0e-3_doubp + &! air
         &   Tot_Wat*WtMol(IH2Ol)*CpWat*1.0e-3_doubp + &! liquid water
         &   Tot_Wat_ice*WtMol(IH2Ol)*CpIce*1.0e-3_doubp)/ &! ice
         &   (Rho +C(InC(IH2Ogi,NABins))*WtGas(IH2Og) + Tot_Wat*WtMol(IH2Ol) + & 
         &   Tot_Wat_ice*WtMol(IH2Ol))
!!$   
!! CHANGE IN TEMPERATURE 
!! CAN BE SOME PREDESCRIPED OR ADIABATIC
!!
    IF (IFADIA == 0) THEN
       dydx(NODEs) = -540.0_doubp/(4000.0_doubp)
    ELSE
       dydx(NODEs) = -g / Cptot*dydx(NODEs-3)*1.e-3_doubp + &
            &        1.0_doubp/( &
            &        C(InC(IH2Ogi,NABins)) * WtGas(IH2Og) * 1.95_doubp + &! water gas
            &        Rho * CpAir*1.0e-3_doubp + &! air
            &        Tot_Wat * WtMol(IH2Ol)*CpWat*1.0e-3_doubp + &! liquid water
            &        Tot_Wat_ice * WtMol(IH2Ol)*CpIce*1.0e-3_doubp &! ice
            &        )*Heat 
!       IF (x.gt.4500)        dydx(NODEs)=0.1
    END IF

!!$
!!$
    IF (IDRIVER_EQ == 1) dydx(NODEs)=0.0_doubp
!!$
!!$   dC/dt  = CHANGE IN NUMBER CONCENTRATION DUE TO THE CHANGE
!!$   IN PRESSURE AND TEMPERATURE
!!$   
    help_term=(1.0_doubp/Y2(NODEs-2)*dydx(NODEs-2) &
         &  -1.0_doubp/Y2(NODEs)*dydx(NODEs))
    dydx(NODEs-1) = Ctot * help_term
!!$   Same for ice
    IF (ice_mod == 1) THEN
       dydx(NODEs-5) = C_ice_tot *  help_term
       if (NFABins2 == 0) dydx(NODEs-5) = 0.0_doubp
    END IF
!!$   Adiabatic variable
    dydx(NODEs-4) = Adia_fex *  help_term
    if (w == 0) dydx(NODEs-4) = 0.0
!!$   Adiabatic change in cloud droplets in C_Y2
    IF (cloud_mod ==1) THEN
       dydx(NODEs-6) = dydx(NODEs-6)+y2(NODEs-6) * help_term    
       if(c_cl_tot == 0)  dydx(NODEs-6) = 0.0_doubp
    ENDIF
!!$   
!!$   CALCULATE THE EFFECT OF PRESSURE AND TEMPERATURE CHANGE
!!$   ON LIQUID PHASE CONCENTRATIONS !!!!!!!!!!!!!
!!$   

    DO J=1,NYAero
       dydx(J)=  dydx(J)+ Y2(J)*help_term
    END DO

    IF (ice_mod == 1) THEN
       DO Ji=1,NFABins2
          dydx(NYAero+NAGases+Ji) = dydx(NYAero+NAGases+Ji)+ &
               & Y2(NYAero+NAGases+Ji)*help_term
       END DO
    END IF
    
    IF (CLOUD_MOD == 1) THEN
       DO Ji=1,N_CL2
          dydx(NYAero+NAGases+NFABins2+Ji) = dydx(NYAero+NAGases+NFABins2+Ji)+ &
               & Y2(NYAero+NAGases+NFABins2+Ji)*help_term
          
       END DO
    END IF

!!$   CALCULATE THE EFFECT OF PRESSURE AND TEMPERATURE CHANGE
!!$   ON GAS PHASE CONCENTRATIONS
    DO JGas= 1,NAGases
       J       = NYAero+JGas
       dydx(J) = dydx(J) &
            &    + Y2(J) *help_term
    END DO

   
    IF(IF_NUC_FEX==1) THEN
          dydx(NOdes-Ndyna)=dydx(NOdes-Ndyna)+ Y2(NOdes-Ndyna) &
               &   * help_term
    ENDIF 

!!$596 aaak = 1.
!!$   
!    DO I = 1,NODEs
       yOld(1:NODEs) = y2(1:NODEs)
!    END DO
!      write(*,*)  dydx(Nodes),Y(Nodes)

  END SUBROUTINE fex

!!$
!!$  SUBROUTINE jex(n,x,y,ml,mu,dfdy,nrowpd)
!!!$!!$C
!!$  END SUBROUTINE jex


!!$  SUBROUTINE JEX(NEQ,T,Y,IA,JA,NZ,P)
!!$    !       Load the Jacobian as a sparse matrix.
!!$    use precisi
!!$    IMPLICIT NONE
!!$    INTEGER NEQ, IA, JA, NZ, ML, MU, COL, ROW, I, NROWPD
!!$    !       INTEGER NZSAVE
!!$    REAL(doubp) :: T, Y, P
!!$    DIMENSION Y(*), IA(*), JA(*), P(*)
!!$    IF (NZ<=0) THEN
!!$       NZ = NEQ*NEQ
!!$       RETURN
!!$    END IF
!!$    ML = NEQ - 1
!!$    MU = NEQ - 1
!!$    NROWPD = NEQ
!!$    CALL JACD(NEQ,T,Y,ML,MU,P,NROWPD)
!!$    IA(1) = 1
!!$    DO I = 1, NEQ
!!$       IA(I+1) = IA(I) + NEQ
!!$    END DO
!!$    I = 0
!!$    DO COL = 1, NEQ
!!$       I = NEQ*(COL-1)
!!$       DO ROW = 1, NEQ
!!$          I = I + 1
!!$          JA(I) = ROW
!!$       END DO
!!$    END DO
!!$    RETURN
!!$  END SUBROUTINE JEX

!!$  SUBROUTINE JACD(NEQ,T,Y,ML,MU,PD,NROWPD)
!!$    use precisi
!!$    !       Load the Jacobian as a dense matrix.
!!$    IMPLICIT NONE
!!$    INTEGER NEQ, ML, MU, NROWPD, I, J
!!$    REAL(doubp) ::  T, Y, PD,DY
!!$    DIMENSION Y(NEQ), PD(NROWPD,NEQ)
!!$    CALL PDERV(T,Y)
!!$    DO J = 1, NEQ
!!$       DO I = 1, NEQ
!!$          PD(I,J) = DY(I+(J-1)*NROWPD)
!!$       END DO
!!$    END DO
!!$    RETURN
!!$  END SUBROUTINE JACD










END MODULE funs
