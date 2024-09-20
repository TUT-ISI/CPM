SUBROUTINE DRIVER(NEq)
  
  USE HEADFILE
  use precisi
  use dvode_f90_m
  use funs, only : fex!, jex
  use aikaa_mittaavat_parametrit
  implicit none
!!$
!!$ *****************************************************************
!!$ * SUBROUTINE DRIVER IS A DRIVER ROUTINE FOR ODE-SOLVER LSODE    *
!!$ *****************************************************************
!!$j
!!$  EXTERNAL Fex, Jex
  integer, intent(inout) :: NEq
  real(doubp), dimension(22 +  9*NEq + NEq**2) :: RWork
  real(doubp), dimension(NEq) :: ATol
  real(doubp), dimension(NEq + NY1 + NY2) :: y,Y2
  real(doubp) :: Time, Tout, pw, pwe, pH2Oeq, phcl,Adia_fex_old2
  real(doubp) :: phno3,c_tot_old,c_cl_old,SURFACE_TENSION
  real(doubp) :: temperature,c_ice_tot_old,t_old,set_vel,c_cl_tot_old
  real(doubp) :: c_Y3_tot_old,c_Y4_tot_old
  real(doubp),dimension(NABins+NY1 + NY2) ::set_vel2
  INTEGER, dimension(20 + NEq) :: IWork
  integer :: IBin, J, IDriv, NEQ2!, JOUT_init
  integer :: ITask, IState2, I, LIW, LRW ,ISpec,II,II2,i_out2,L3,JJ,IFAM,I_VEL
  integer :: J_BOX,NFABINS_OLD,N_CL2_OLD,II3,II4
  TYPE (VODE_OPTS) :: OPTIONS
!!$
!!$     neq    = number of first order ode-s.
!!$     lrw    = declared length of rwork (in user-s dimension).
!!$     liw    = declared length of iwork (in user-s dimension).
!!$
!!  write(*,*) 'in driver'
  IDriv=0
  LRW    = 22 + 9 * NEq + NEq**2
  LIW    = 20 + NEq
!!$
  P = P - C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$     can be found in liq_ice.f
  call initialize
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  Y(1:NEq + NY1 + NY2)=0.0
  Y2(1:NEq + NY1 + NY2)=0.0
  N_C3_OLD = 0
  N_C4_OLD = 0
     
  CALL C2Y(y,NEQ)
!!$
  i_out2=0
  DO  IBin = 1,NABins
     Rp_old(IBin) = Rp(IBin)
  END DO
!!$
!     IFADIA = 0
  IFADIA = 1                ! OBS!
!!$
  IF(IFADIA == 0) THEN
     w = 0.e0_doubp
  ENDIF
!!$
!  
! Put some values to Y
!
  if (ice_mod == 1) then
     y(NEq - 5)    = C_ice_tot
  end if
  if (cloud_mod ==1) y(NEq - 6)    = c_cl_tot
  y(NEq - 4)       = Adia_fex    
  y(NEq - 3)       = z
  y(NEq - 2)       = P
  y(NEq - 1)       = Ctot
  y(NEq)           = T

  Adia_fex_old2=Adia_fex

  DO  IBin = 1,NABins
     IF_Thermo(IBin)=0
  END DO

!!$
!!$     time   = the initial value of time
!!$     tOut   = first point where output is desired
!!$ 

  time   = 0.0_doubp
  tOut   = time+DTOUT    !1.d-5
!!$     
!!$     
  DO I=1,LRW
     RWork(I) = 0.e0_doubp
  END DO
!!$     
  DO I=1,LIW
     IWork(I)=0
  END DO
!!$
!!$     RWork(5) = the step size to be attempted on the first step.
!!$     RWork(6) = the maximum absolute step size allowed.
!!$     RWork(7) = the minimum absolute step size allowed
!!$
  RWork(5) = 1.0e-8
  RWork(6) = StepMax!1.e-0
  RWork(7) = 1.0e-20
  IWork(6) = 1000000
!!$
!!$     ITol    = 1 or 2 according as atol (below) is a scalar or array.
!!$     RTol    = relative tolerance parameter (scalar).
!!$     ATol    = absolute tolerance parameter (scalar or array).
!!$               the estimated local error in y(i) will be controlled so as
!!$               to be roughly less (in magnitude) than
!!$                  ewt(i) = RTol*abs(y(i)) + ATol     if ITol = 1, or
!!$                  ewt(i) = RTol*abs(y(i)) + ATol(i)  if ITol = 2.
!!$               thus the local error test passes if, in each component,
!!$               either the absolute error is less than atol (or ATol(i)),
!!$               or the relative error is less than rtol.
!!$               use RTol = 0.0 for pure absolute error control, and
!!$               use ATol = 0.0 (or ATol(i) = 0.0) for pure relative error
!!$               control.  caution.. actual (global) errors may exceed these
!!$               local tolerances, so choose them conservatively.
!!$     ITask   = 1 for normal computation of output values of y at t = tOut.
!!$     IOpt    = 0 to indicate no optional inputs used.
!!$     MF      = method flag, 22 for stiff method, internally generated 
!!$               full jacobian.
!!$
  IState2 = 1
!  ITol    = 2
!  RTol    = RelTol
  ATol(1:NEq) = AbsTol 
  ITask   = 1
!  IOpt    = 1
!  MF      = 22
! Option Types

! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! RELERR                 - real(wp) scalar
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! NEVENTS                - integer
! Options:
! ABSERR                 = Absolute error tolerance
! ABSERR_VECTOR          = Vector of absolute error tolerances
! RELERR                 = Scalar relative error tolerance
! NEVENTS                = Number of event functions (requires
!                          user-supplied GFUN)
! DENSE_J                = Use dense linear algebra if .TRUE.
! BANDED_J               = Use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH      = Lower bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH      = Upper bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
! SPARSE_J               = Use sparse linear algebra if .TRUE.
! USER_SUPPLIED_JACOBIAN = Exact Jacobian option
!                          (requires user-supplied JAC;
!                          ignored for SPARSE_J=.TRUE.)


  OPTIONS = SET_OPTS(DENSE_J=.TRUE. ,&
       USER_SUPPLIED_JACOBIAN=.FALSE.,   &
       RELERR=RelTol,ABSERR_VECTOR=ATOL(1:NEq),  &
       MXSTEP=1000000,H0=RWork(5),HMAX=StepMax)!,        &
!       CONSTANT_JACOBIAN=.FALSE.,            &
!       USER_SUPPLIED_SPARSITY=.TRUE.,    &
!       JACOBIAN_BY_JACSP=.TRUE.,                &
!       MA28_RPS=.TRUE.)
!  OPTIONS = SET_OPTS(SPARSE_J=.TRUE., ABSERR=AbsTol, RELERR=RTOL)
!!$     
!
  if (indx_ret == 0) then
     CALL OUTPUT(time,Y)
  else 
     indx_ret=0
     y(1:NEq)=Y_init(1:NEq)
  end if

  IF_Cond(1:NaBins)=1

!
!  No condensation on particles is I_MODE_CONTACT 
!
  IF(IF_FREEZ_CONTACT==1) IF_Cond(FirstBinMode(I_MODE_CONTACT) : FirstBinMode(I_MODE_CONTACT)+NBinsPerMode(I_MODE_CONTACT)-1) = 0
  IF(IF_FREEZ_DEPOS==1) IF_Cond(FirstBinMode(I_MODE_DEPOS) : FirstBinMode(I_MODE_DEPOS)+NBinsPerMode(I_MODE_DEPOS)-1) = 0
!
!
!!$  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!   ACTUAL DIFFERENTIAL EQUATION SOLVIN LOOP STARTS HERE
!
!!$  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!
  DO !JOUT = JOUT_init, (TotSec/DTOUT+1)*10000
!
!     
!! 
!! Here new system is added to move Y2 to solver, which contains only those
!! bins with something to condense on.
!!

     c_tot_old = y(NODEs-1)
     c_cl_tot_old=y(NODEs-6)
! NYCONS=water concentrations from 1..NABins and gas phase water     
     Y2(1:NYCONCS)= Y(1:NYCONCS)


!
! Ice particle concentrations with high enough water content and number concentration
!     
     II=1
     Y1C(1:NY1)=0
     NYCONCS2=NYCONCS
     DO I=NYCONCS+1,NYCONCS+NY1
        IF((Y(I)>1.0e-25.AND.C_Y1(I-NYCONCS).GT.1e-6)) THEN
           Y2(NYCONCS+II)=Y(I)
           Y1C(II)=I
           II=II+1
           NYCONCS2=NYCONCS2+1
        ENDIF
        IF(C_Y1(I-NYCONCS).LT.1e-6) THEN
           Y(I)=C_Y1(NY1+I-NYCONCS)
        ENDIF
     ENDDO

     II=II-1
 
!
! Cloud droplets from CY2 with high enough concentration
!
     II2=1
     N_CL2=0
     Y2C(1:NY2)=0
     IF(CLOUD_MOD == 1) THEN
        DO I=NYCONCS+NY1+1,NEQ-NDYNA
           IF(Y(I)>1.0e-25.AND.C_Y2(I-NYCONCS-NY1).GT.1e-10) THEN
              Y2(NYCONCS2+1)=Y(I)
              Y2C(II2)=I
              II2=II2+1
              NYCONCS2=NYCONCS2+1
              N_CL2=N_CL2+1

           ELSE
              Y(I)=C_Y2(NY2+I-NYCONCS-NY1)
           ENDIF
        ENDDO
     ENDIF
     II2=II2-1
!
! graupel particles from CY3 with high enough concentration
!
     II3=1
     N_C3=0
!     Y3C(1:NY3)=0
!     DO I=NYCONCS+NY1+NY2+1,NEQ-NDYNA-NY4
!        IF(Y(I)/=0.0.AND.C_Y3(I-NYCONCS-NY1-NY2).GT.1e-8) THEN
!           Y2(NYCONCS2+1)=Y(I)
!           Y3C(II3)=I
!           II3=II3+1
!           NYCONCS2=NYCONCS2+1
!           N_C3=N_C3+1
!           
!        ELSE
!           Y(I)=C_Y3(NY3+I-NYCONCS-NY1-NY2)
!        ENDIF
!     ENDDO
     II3=II3-1
!
! snow particles from CY3 with high enough concentration
!
     II4=1
     N_C4=0
!     Y4C(1:NY4)=0
!     DO I=NYCONCS+NY1+NY2+NY3+1,NEQ-NDYNA
!        IF(Y(I)/=0.0.AND.C_Y4(I-NYCONCS-NY1-NY2-NY3).GT.1e-8) THEN
!           Y2(NYCONCS2+1)=Y(I)
!           Y4C(II4)=I
!           II4=II4+1
!           NYCONCS2=NYCONCS2+1
!           N_C4=N_C4+1
!           
!        ELSE
!           Y(I)=C_Y4(NY4+I-NYCONCS-NY1-NY2-NY3)
!        ENDIF
!     ENDDO    
     II4=II4-1

     IInuc=0
     IF(IF_NUC_FEX ==1 ) THEN
        
        Y2(NYCONCS+II+II2+II3+II4+1) = c_ice2
        IInuc=IInuc+1
     ENDIF
     
!
! Values like total numbr concentrations and temperature/pressure/altitude from Y to Y2
!     

     Y2(NYCONCS+II+II2+1 + IInuc : II+II2+NDYNA+NYCONCS + IInuc)=Y(NEQ-NDYNA+1:NEQ)

     NEQ2=NYCONCS+II+II2+II3+II4+NDYNA + IInuc
     NFABINS2=II
     N_CL2=II2
 

     IF((NFABINS2.NE.NFABINS_OLD.OR.(C_ice_tot_old/c_ice_tot-1)**2>1.0e-6)) IState2 = 1
     IF((N_CL2.NE.N_CL2_OLD.OR.(C_cl_tot_old/c_cl_tot-1)**2>1.0e-8)) IState2 = 1

! Uncomment these as soon as C_Y3 and C_Y4 are working

!     IF((N_C3.NE.N_C3_OLD.OR.(C_Y3_tot_old/c_Y3_tot-1)**2>1.0e-8)) IState2 = 1
!     IF((N_C4.NE.N_C4_OLD.OR.(C_Y3_tot_old/c_Y4_tot-1)**2>1.0e-8)) IState2 = 1

!     write(*,*) ISTATE2,(C_ice_tot_old/c_ice_tot-1)**2
!     write(*,*) NFABINS2,NFABINS_OLD,(C_ice_tot_old/c_ice_tot-1)**2
     C_ice_tot_old=c_ice_tot
     C_cl_tot_old=c_cl_tot
     C_Y3_tot_old=c_Y3_tot
     C_Y4_tot_old=c_Y4_tot
     NODEs=neq2
     t_old=T
!!
!! Call stiff differential equation solver
!! Notice it is Y2 that goes in
!!


!     write(*,*) tout
     CALL DVODE_F90(fex,NEq2,Y2(1:NEq2),time,tout,ITASK,ISTATE2,OPTIONS)

!     CALL CPU_TIME(TIME_OUT)
!!
!!  Move values back from Y2 to Y
!!  
     NODEs=neq
     Y(1:NYCONCS)= Y2(1:NYCONCS)
     II=NYCONCS+1
     DO I=1,NY1
        IF(Y1C(I)/=0) THEN
          Y(Y1C(I))=Y2(II)
          II=II+1
        ENDIF
     ENDDO

     IF(CLOUD_MOD==1) THEN
        DO I=1,NY2
           IF(Y2C(I)/=0) THEN
              Y(Y2C(I))=Y2(II)
              II=II+1
           ENDIF
        ENDDO
     ENDIF


     Y(NEQ-NDYNA+1:NEQ)=Y2(NEQ2-NDYNA+1:NEQ2)

     IF(IF_NUC_FEX ==1 ) THEN
        
        c_ice2=Y2(NEQ2-NDYNA)
        
     ENDIF
     
     DO J=1,NYLiquids
        C(IY2C(J))=y(J)
     END DO
     DO J = 1, NAGases
        C(NAero+J)=y(J + NYAero)
     END DO
     IF (istate2 < 0) THEN
        write(*,90) istate2
        stop
     end IF

!!$     
!!$     CALL VISIB
!!$     

!
!  Calculate relative humidity
!
     pw          = C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp

     select case(Ithermo_mod)
     case(1) !Luo's thermo
        temperature=T
        CALL s_luo(temperature,0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
        pwe=pwe/100.0_doubp
     case(2,3,4) !jacobson
        pwe         = pH2Oeq(T) * 101325.0_doubp
     end select

     RH          = pw / pwe
    
 !    RH_OLD = RH


!
! Correct the number concentrations due to adiabatic expansion.
! 
     Ctot = y(NODEs-1)
     if (ice_mod == 1) then
        C_ice_tot   = y(NODEs-5)
        if (C_ice_tot>1.0e-20) then
           DO I=1,NFBins
              C_Y1(I)=C_ice_tot*C_iceFrac(I)
           END DO
        end if
     end if

     do I=1,NABins
        C(I)=Ctot*CFrac(I)
     end do
     C_cl_tot = y(NODEs-6)
     IF(CLOUD_MOD == 1.AND.c_cl_old>0.0) C_Y2(1:NY2) =  C_clfrac(1:NY2)*c_cl_tot
!     IF(c_y3_tot >0) C_Y3(1:NY3) =  C_Y3frac(1:NY3)*c_y3_tot
!     IF(c_y4_tot >0) C_Y4(1:NY4) =  C_Y4frac(1:NY4)*c_y4_tot
!     
! But the values from Y back to C, C_Y! and C_Y2     
!
     DO J=1,NYLiquids
        C(IY2C(J))=y(J)
     END DO
 

!!$C     
     IF (ice_mod == 1) THEN
        DO J=1,NY1
           C_Y1(InCi(IH2Ol,J)) = Y(NYConcs + J)
        END DO
     END IF
!!C
     IF (cloud_mod == 1) THEN
        DO J=1,NY2
           C_Y2(InCl(IH2Ol,J)) = Y(NYConcs +NY1+ J)
        END DO
     END IF
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
!!     
     DO J = 1, NAGases
        C(NAero+J)=y(J + NYAero)
     END DO
!!$C     
!     IF (cloud_mod == 1) c_cl_tot = y(NODEs-6)
!
!     IF (ice_mod == 1) THEN
!        C_ice_tot                = y(NODEs-5)
!     END IF
!     Adia_fex                    = y(NODEs-4)
!     z                           = y(NODEs-3)
!     P                           = y(NODEs-2)
!     Ctot                        = y(NODEs-1)
!     T                           = y(NODEs)
!!$C     
!!$C     pw  = WATER VAPOR PRESSURE
!!$C     pwe = WATER SATURATION VAPOR PRESSURE, Seinfeld & Pandis 1998 (page 49)
!!$C     
     pw          = C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp
     !  pwe         = pH2Oeq(T) * 101325.0
     temperature=T
    


!     IF(I_BOX>10) write(*,*) R_cl(27:35),'0'
 
!!
!! Move values from Y to C, CY1 and CY2
!!
   
     CALL init_Y2C(Y,NEQ)
!!

     dtemp=(T-T_old)/max(1.0e-3,DtOut_used)   ! Used in liq_ice to calculate deposition nucleation

!!
!!  The call for coagulation could be placed to somewhere else
!!

!     ISTATE2=1
!     CALL CPU_TIME(TIME_OUT)
     TOTAL_EQMODEL_TIME=TOTAL_EQMODEL_TIME+(TIME_OUT-TIME_IN)
     EQ_CALCULATIONS=EQ_CALCULATIONS+1

!
! This is here to decrease calculation time. Thermodynamic is not needed for those bins
! condensation not in. Also thermodynamics can be turned of for dilute cloud droplets if 
! there is no chemistry
!
! Uncomment line "IF_Thermo(IBin)=0" and model thermodynamics is not called to any bins
!
     DO  IBin = 1,NABins
        Rp_old(IBin) = Rp(IBin)
!        IF_Thermo(IBin)=0
        IF( water_activity(Ibin)> 1.0-1.0e-4.AND. RH>0.995)  IF_Thermo(IBin)=0
        IF( IF_Cond(IBin)==0)  IF_Thermo(IBin)=0
!        IF_Thermo(IBin)=0
     END DO

!
! Values back from C to Y. I'm not sure if this is needed!!!!
!
!     CALL C2Y(y,NEQ)
   
!
! Fix bins in dists C_Y1, C_Y2 (and also C_Y3, C_Y4)
!

     CALL FIX_BINS(Y,NEQ)
  
     if (ice_mod == 1) then
        IDriv=NFABins

     
        IF(IfFreez==1) CALL liq_ice(Y,NEq)

        DO II = 1,NY2
           CALL RADIUS_Y2(II)
        END DO
        !     CALL CPU_TIME(TIME_IN)

        DO II = 1,NaBins
           CALL RADIUS_dry(II)
        END DO
!
! Because of coagulation outside the solver, the Istate is 1 after every call
! Not smart....
!

        IF(IfCoag == 1) THEN
           CALL COAGULATION(Y,NEQ,tout-DtOut_used,tout)   
           istate2 = 1
        ENDIF

!
! Bins are also fixxed after coagulation/coalescence
!
        CALL FIX_BINS(Y,NEQ)


        DO II = 1,NaBins
           CALL RADIUS_dry(II)
        END DO

!
!   Call for sedimentation. It is only ON, if IF_BOXES=1 
!
!
!        write(*,*) time,I_out2,int(t_file/DtOut_used)
        IF(IF_BOXES==1.AND.I_BOX>0) THEN 
           CALL SEDIMENTATION2(NEQ,Y)
         END IF

        IF( IfCloud == 1 ) THEN
           DO II = 1,NY2
              CALL RADIUS_Y2(II)
           END DO
        ENDIF
        
        IF( IfCoag == 1 ) THEN
           DO II = 1,NY1
              CALL RADIUSICE2(II)
           END DO
        ENDIF
        
        NFABINS_OLD = NFABINS2
        N_CL2_OLD=N_CL2
        N_C3_OLD = N_C3
        N_C4_OLD = N_C4


        DO II=1,NY1
           IF(tout>180000.0.OR.C_Y1(II)<1.0e-40) THEN
              C_Y1(II) = 0.0
              do  ISpec=1,NALiquids + NASolids
                 C_Y1(InCi(ISpec,II)) = 0.0
              end do
              Y(NYCONCS+II)= C_Y1(NY1+II)
           END IF
        END DO
           
        
        DO II=1,NY2
           IF(tout>180000.0.OR.C_Y2(II)<1.0e-40) THEN
              C_Y2(II) = 0.0
              do  ISpec=1,NALiquids + NASolids
                 C_Y2(InCl(ISpec,II)) = 0.0
              end do
              Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
           END IF
        END DO
           
        


        IF(I_out2==min(2,int(t_file/DtOut_used))) CALL OUTPUT(time,Y(1:NEq))

        adia_fex_old2=adia_fex_old

        IF(I_out2==int(t_file/DtOut_used)) THEN
           I_step=I_step+1
!           write(*,*) I_BOX,I_step
        END IF

        IF(I_out2==int(t_file/DtOut_used)) I_out2=0
        I_OUT2=I_OUT2+1


!!
        if (NFABins == 0.AND.IFfreez == 1) call freez_prob(2)
        wat_act_old(1:NABins) = Water_activity(1:NABins)

        if (ice_step_start == 1) sum_prob=sum(prob(1:NABins))
        IF (sum_prob > 1.0e-7_doubp) THEN
!           IF(NFABINS2>0) IState2 = 1
           tOut = tOut + DtOutIce
           DtOut_used = DtOutIce
        else
!           IF(NFABINS2>0) IState2 = 1
           tOut = tOut + DtOut
           DtOut_used = DtOut
        END IF
!        if( tOut.GT.1382) pause
     else
        tOut = tOut + DtOut
     end if

!!$
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$         if (NFABins .gt. 0) return
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     if (ice_mod == 1) then
        do I=1,NY1
           C_iceFrac(I)=C_Y1(I)/C_ice_tot
        end do
     end if
     if (cloud_mod == 1 .AND. C_cl_tot>0.0) then
        do I=1,NY2
           C_clFrac(I)=C_Y2(I)/C_cl_tot
        end do
     endif
     DO I = 1, NABins
        CFrac(I) = C(I) / CTot
     END DO


!     if (IDriv /= NFABins .and. ice_mod == 1) then 
!        I=NBins + NBins * NAerT + NGases+NFBins
!        YOld(1:I) = 0.0_doubp
!        NODEs=NEq
!        indx_ret=1
!        time_help=time
!        tOut_help=tOut
 !       JOUT_help=JOUT+1
!        NEq_help=NEq
!        write(*,*) 'here4',NEQ,NFABins,IDriv
        
!        Y_init(1:(NEq))=Y(1:(NEq))
!             write(*,*) 'here5'
!        return
!     end if
     if (time >= TotSec) return
  END DO
!!$
!  StepSize=RWork(11)
!!$     
!!$      write(*,*) time, stepsize,iwork(16),aerintv
!!$     
  DO I=1, NEq
     IF(y(I) < 0.0_doubp) THEN
        WRITE(*,*) I,y(I)
        write(*,*) 'Y neg in driver'
        stop
        !        pause 'hah'
!        C(I)=0.e0
     ENDIF
  END DO
!!$     
 
  
  if (IState2 >= 0) THEN !go to 80
!!$     
     write(*,60) iwork(11),iwork(12),iwork(13)
!!$
     return
!!$20 format(7h at t =,e12.4,6h   y =,3e14.6)
60   format(/12h no. steps =,i4,11h  no. f-s =,i4,11h  no. j-s =,i4)
!!$
!!$     ODE-SOLVER FAILED, STOP PROGRAM
!!$
  ELSE
!!$80   write(*,90) istate2
90   format(///22h error halt.. istate =,i3)
!!$     
END if
!!$666 FORMAT(99(E11.4,1X))
!!$     
  stop
end SUBROUTINE DRIVER



SUBROUTINE SEDIMENTATION(NEQ,Y)
  
  !
  ! This subroutine can be used to calculate the sedimentation of large particles
  ! if model is run in columm mode. i.e. box after box.
  ! At the moment this is quite simplified. 
  ! 
  !
  ! CY1MULTI contains concentrations in ice particles for all previous time steps 
  ! CY2multi contains concentrations is cloud droplets
  ! set_vel is the function to calculate gravitational settling velocity of particle
  ! 

  USE HEADFILE
  use precisi
  integer, intent(in) :: NEq
  real(doubp), intent(inout), dimension(NEQ) :: y
  integer :: ISpec,II,L3,I_VEL
  integer :: J_BOX
  real(doubp) ::set_vel
  real(doubp),dimension(NABins+NY1 + NY2) ::set_vel2


  
  CY1multi(I_BOX,I_step,1)=time
  CY1multi(I_BOX,I_step,2)=z 
  CY1multi(I_BOX,I_step,3:2+NY1 + NY1 * (NALiquids + NASolids))=  &
       &   C_Y1(1:NY1 + NY1 * (NALiquids + NASolids))
  CY2multi(I_BOX,I_step,3:2+NY2 + NY2 * (NALiquids + NASolids))=  &
       &   C_Y2(1:NY2 + NY2 * (NALiquids + NASolids))
     



  
  if(I_step>I_step_max) I_step_max= I_step
  
! Calculate settling velocity
  DO II=1,NY1-1
     set_vel2(NaBins+II)=set_vel(R_ice(II))
  END DO
  DO II=1,NY2-1
     set_vel2(NaBins+NY1+II)=set_vel(R_cl(II))
  END DO
    
!    DO II=1,NY3-1
!     set_vel2(NaBins+NY1+NY2+II)=1.0e-20
!  END DO
!  DO II=1,NY3-1
!     set_vel2(NaBins+NY1+NY2+NY3+II)=1.0e-20
!  END DO
  
!
! First ice particles
!
  DO II=1,NY1-1
     I_VEL=NaBins+II
     
     IF(I_BOX==1) THEN   ! If thi is first box, then only sedimentation downwards
        IF(set_vel2(I_VEL)<w) THEN
           C_Y1(II) = C_Y1(II) -C_Y1(II)*(set_vel2(I_VEL)/w)
           do  ISpec=1,NALiquids + NASolids
              C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                   &  -C_Y1(InCi(ISpec,II))*(set_vel2(I_VEL)/w)
           end do
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ELSE
           C_Y1(II) =0.0
           do  ISpec=1,NALiquids + NASolids
              C_Y1(InCi(ISpec,II)) = 0.0
           end do
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ENDIF
     ENDIF
        
        
        
     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/w)) THEN
            
! If set_vel less than w, then gradient between this and upper box. Else boxes for gradiens are 
! calculated based on set_vel/w
! 
        IF(set_vel2(I_VEL)<w) THEN
           
           C_Y1(II) = C_Y1(II) + max(-C_Y1(II)*.9999, &
                &  (CY1multi(I_BOX-1,I_step+1,II+2)-C_Y1(II)) * (set_vel2(I_VEL)/w))
           
           do  ISpec=1,NALiquids + NASolids
              
              C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                   &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
                   &  (CY1multi(I_BOX-1,I_step+1,InCi(ISpec,II)+2) &
                   &  -C_Y1(InCi(ISpec,II))) * (set_vel2(I_VEL)/w))
              
           end do
              
           Y(NYCONCS+II)= C_Y1(NY1+II)
           
        ELSE
           J_BOX=max(1,I_BOX-1-int(set_vel2(I_VEL)/w))
           

           C_Y1(II) =CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),II+2) &
                &  + max(-CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),II+2)*.9999,  & 
                &  (CY1multi(J_BOX,I_step+1+int(set_vel2(I_VEL)/w),II+2)  &
                &  -CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w) &
                &  ,II+2)) * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))

           do  ISpec=1,NALiquids + NASolids
                       
!                       C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
!                            &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
              C_Y1(InCi(ISpec,II)) = CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)  & 
                   &  + max(-CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)*.9999,     &
                   &  (CY1multi(J_BOX,I_step+1+int(set_vel2(I_VEL)/w),  &
                   &  InCi(ISpec,II)+2)-CY1multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)) &
                   &  * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
              
           end do
           
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ENDIF
     ENDIF
  END DO
  
              
  C_ice_tot=sum(C_Y1(1:NY1))
  y(NODEs-5)=C_ice_tot
  
  
  
     
  DO II=5,NY2-1
     I_VEL=NaBins+NY1+II
     
     IF(I_BOX==1) THEN
        IF(set_vel2(I_VEL)<w) THEN
           C_Y2(II) = C_Y2(II) -C_Y2(II)*(set_vel2(I_VEL)/w)
           do  ISpec=1,NALiquids + NASolids
              C_Y2(InCl(ISpec,II)) = C_Y2(InCl(ISpec,II))  & 
                   &  -C_Y2(InCl(ISpec,II))*(set_vel2(I_VEL)/w)
           end do
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
        ELSE
           C_Y2(II) =0.0
           do  ISpec=1,NALiquids + NASolids
              C_Y2(InCl(ISpec,II)) = 0.0
           end do
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
        ENDIF
     ENDIF
                 
                 
        
     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/w)) THEN
        
        IF(set_vel2(I_VEL)<w) THEN
           if(I_BOX>3.and.CY2multi(I_BOX-1,I_step+1,II+2)>C_Y2(II)) then
              write(*,*)  C_Y2(II) ,(CY2multi(I_BOX-1,I_step+1,II+2)-C_Y2(II)) * (set_vel2(I_VEL)/w),II
           endif
           C_Y2(II) = C_Y2(II) + max(-C_Y2(II)*.9999, &
                &  (CY2multi(I_BOX-1,I_step+1,II+2)-C_Y2(II)) * (set_vel2(I_VEL)/w))
           
           do  ISpec=1,NALiquids + NASolids
              
              C_Y2(InCl(ISpec,II)) = C_Y2(InCl(ISpec,II))  & 
                   &  + max(-C_Y2(InCl(ISpec,II))*.9999,     &
                   &  (CY2multi(I_BOX-1,I_step+1,InCl(ISpec,II)+2) &
                   &  -C_Y2(InCl(ISpec,II))) * (set_vel2(I_VEL)/w))
              
           end do
           
             Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
           
        ELSE
           J_BOX=max(1,I_BOX-1-int(set_vel2(I_VEL)/w))
              
              !                    C_Y1(II) = C_Y1(II) + max(-C_Y1(II)*.9999,  & 
           C_Y2(II) =CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),II+2) &
                &  + max(-CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),II+2)*.9999,  & 
                &  (CY2multi(J_BOX,I_step+1+int(set_vel2(I_VEL)/w),II+2)  &
                &  -CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w) &
                &  ,II+2)) * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
           do  ISpec=1,NALiquids + NASolids
              
!                       C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                 !                            &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
              C_Y2(InCl(ISpec,II)) = CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCl(ISpec,II)+2)  & 
                   &  + max(-CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCl(ISpec,II)+2)*.9999,     &
                   &  (CY2multi(J_BOX,I_step+1+int(set_vel2(I_VEL)/w),  &
                   &  InCl(ISpec,II)+2)-CY2multi(J_BOX+1,I_step+int(set_vel2(I_VEL)/w),InCl(ISpec,II)+2)) &
                   &  * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
           end do
              
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
        ENDIF
     ENDIF
  END DO
  
  C_cl_tot=sum(C_Y2(1:NY2))
  y(NODEs-6)=c_cl_tot
  
  
END SUBROUTINE SEDIMENTATION





SUBROUTINE SEDIMENTATION2(NEQ,Y)
  
  !
  ! This subroutine can be used to calculate the sedimentation of large particles
  ! if model is run in columm mode. i.e. box after box.
  ! At the moment this is quite simplified. 
  ! 
  !
  ! CY1MULTI contains concentrations in ice particles for all previous time steps 
  ! CY2multi contains concentrations is cloud droplets
  ! set_vel is the function to calculate gravitational settling velocity of particle
  ! 

  USE HEADFILE
  use precisi
  integer, intent(in) :: NEq
  real(doubp), intent(inout), dimension(NEQ) :: y
  integer :: ISpec,II,L3,I_VEL
  integer :: J_BOX, J_STEP
  real(doubp) ::set_vel
  real(doubp),dimension(NABins+NY1 + NY2+NY3+NY4) ::set_vel2


  eetta = 1.8325d-5*(416.16/(T+120.0))*(T/296.16)**1.5 ![kg m-1 s-1]
  GASPEED2 = SQRT(8. * R * T/(Pi *0.029))      
  CY1multi(1,I_step,1)=time
  CY1multi(1,I_step,2)=z 
  CY1multi(1,I_step,3:2+NY1 + NY1 * (NALiquids + NASolids))=  &
       &   C_Y1(1:NY1 + NY1 * (NALiquids + NASolids))
  CY2multi(1,I_step,3:2+NY2 + NY2 * (NALiquids + NASolids))=  &
       &   C_Y2(1:NY2 + NY2 * (NALiquids + NASolids))
     
  

  
  if(I_step>I_step_max) I_step_max= I_step
  
! Calculate settling velocity
  DO II=1,NY1-1
     
     set_vel2(NaBins+II)=set_vel(R_ice(II))

  END DO
  DO II=1,NY2-1
     set_vel2(NaBins+NY1+II)=set_vel(R_cl(II))
  END DO
          

!
! First ice particles
!
  DO II=1,NY1-1
     I_VEL=NaBins+II
     
     IF(I_BOX==1) THEN   ! If thi is first box, then only sedimentation downwards
        IF(set_vel2(I_VEL)<w) THEN
           C_Y1(II) = C_Y1(II) -C_Y1(II)*(set_vel2(I_VEL)/w)
           do  ISpec=1,NALiquids + NASolids
              C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                   &  -C_Y1(InCi(ISpec,II))*(set_vel2(I_VEL)/w)
           end do
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ELSE
           C_Y1(II) =0.0
           do  ISpec=1,NALiquids + NASolids
              C_Y1(InCi(ISpec,II)) = 0.0
           end do
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ENDIF
     ENDIF
        
        

     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/w)) THEN
            
! If set_vel less than w, then gradient between this and upper box. Else boxes for gradiens are 
! calculated based on set_vel/w
! 
        IF(set_vel2(I_VEL)<w*10) THEN
!        IF(set_vel2(I_VEL)<.5*w) THEN
           C_Y1(II) = C_Y1(II) + max(-C_Y1(II)*.9999, &
                &  (CY1multi(2,I_step+1,II+2)-C_Y1(II)) * (set_vel2(I_VEL)/w))
!
!          
!
!
           do  ISpec=1,NALiquids + NASolids
              
              C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                   &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
                   &  (CY1multi(2,I_step+1,InCi(ISpec,II)+2) &
                   &  -C_Y1(InCi(ISpec,II))) * (set_vel2(I_VEL)/w))
              
           end do

              
           Y(NYCONCS+II)= C_Y1(NY1+II)

           
        ELSE
           J_BOX=max(1,1+int(set_vel2(I_VEL)/(w*10.0)))

!           J_BOX=max(2,2+int(set_vel2(I_VEL)/w-.5))
           

           C_Y1(II) =CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2) &
                &  + max(-CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2)*.9999,  & 
                &  (CY1multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),II+2)  &
                &  -CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w) &
                &  ,II+2)) * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!            C_Y1(II) =CY1multi(J_BOX,I_step+J_BOX-1,II+2)

           do  ISpec=1,NALiquids + NASolids
                       
!                       C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
!                            &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
              C_Y1(InCi(ISpec,II)) = CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)  & 
                   &  + max(-CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)*.9999,     &
                   &  (CY1multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),  &
                   &  InCi(ISpec,II)+2)-CY1multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InCi(ISpec,II)+2)) &
                   &  * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!              C_Y1(InCi(ISpec,II)) = CY1multi(J_BOX,I_step+J_BOX-1,InCi(ISpec,II)+2)
           end do
           
           Y(NYCONCS+II)= C_Y1(NY1+II)
        ENDIF
     ENDIF
  END DO
  
              
  C_ice_tot=sum(C_Y1(1:NY1))
  y(NODEs-5)=C_ice_tot
  


     
  DO II=1,NY2-1
     I_VEL=NaBins+NY1+II
     
     J_STEP=I_step+floor(set_vel2(I_VEL)/(w*DtOut))
     J_BOX=max(1,1+floor(set_vel2(I_VEL)/(w*DtOut)))

     IF(I_BOX==1) THEN
        IF(set_vel2(I_VEL)<w*DtOut) THEN
           C_Y2(II) = C_Y2(II) -C_Y2(II)*(set_vel2(I_VEL)/(w*DtOut))
           do  ISpec=1,NALiquids + NASolids
              C_Y2(InCl(ISpec,II)) = C_Y2(InCl(ISpec,II))  & 
                   &  -C_Y2(InCl(ISpec,II))*(set_vel2(I_VEL)/(w*DtOut))
           end do
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
        ELSE
           C_Y2(II) =0.0
           do  ISpec=1,NALiquids + NASolids
              C_Y2(InCl(ISpec,II)) = 0.0
           end do
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
        ENDIF
     ENDIF
                 
     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/(w*DtOut))) THEN
        
!        IF(set_vel2(I_VEL)<w) THEN
!
!           C_Y2(II) = C_Y2(II) + max(-C_Y2(II)*.9999, &
!                &  (CY2multi(2,I_step+1,II+2)-C_Y2(II)) * (set_vel2(I_VEL)/w))
!
!           do  ISpec=1,NALiquids + NASolids
!              
!              C_Y2(InCl(ISpec,II)) = C_Y2(InCl(ISpec,II))  & 
!                   &  + max(-C_Y2(InCl(ISpec,II))*.9999,     &
!                   &  (CY2multi(2,I_step+1,InCl(ISpec,II)+2) &
!                   &  -C_Y2(InCl(ISpec,II))) * (set_vel2(I_VEL)/w))
!              
!           end do
!           
!           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
           
!        ELSE 
          
           C_Y2(II) =CY2multi(J_BOX,J_step,II+2) &
                &  + max(-CY2multi(J_BOX,J_step,II+2)*.9999,  & 
                &  (CY2multi(J_BOX+1,J_step+1,II+2)  &
                &  -CY2multi(J_BOX,J_step &
                &  ,II+2)) * (set_vel2(I_VEL)/(w*DtOut)-floor(set_vel2(I_VEL)/(w*DtOut))))
          
           do  ISpec=1,NALiquids + NASolids
              
              C_Y2(InCl(ISpec,II)) = CY2multi(J_BOX,J_step,InCl(ISpec,II)+2)  & 
                   &  + max(-CY2multi(J_BOX,J_step,InCl(ISpec,II)+2)*.9999,     &
                   &  (CY2multi(J_BOX+1,J_step+1,  &
                   &  InCl(ISpec,II)+2)-CY2multi(J_BOX,J_step,InCl(ISpec,II)+2)) &
                   &  * (set_vel2(I_VEL)/(w*DtOut)-floor(set_vel2(I_VEL)/(w*DtOut))))
           end do
              
           Y(NYCONCS+NY1+II)= C_Y2(NY2+II)
!        ENDIF
     ENDIF
  END DO
  
  C_cl_tot=sum(C_Y2(1:NY2))
  y(NODEs-6)=c_cl_tot



! Sedimentation for particles in C_Y3

!   DO II=1,NY3-1
!     I_VEL=NaBins+NY1+NY2+II
!     
!     IF(I_BOX==1) THEN
!        IF(set_vel2(I_VEL)<w) THEN
!           C_Y3(II) = C_Y3(II) -C_Y3(II)*(set_vel2(I_VEL)/w)
!           do  ISpec=1,NALiquids + NASolids
!              C_Y3(InC3(ISpec,II)) = C_Y3(InC3(ISpec,II))  & 
!                   &  -C_Y3(InC3(ISpec,II))*(set_vel2(I_VEL)/w)
!           end do
!           Y(NYCONCS+NY1+NY2+II)= C_Y3(NY3+II)
!        ELSE
!           C_Y3(II) =0.0
!           do  ISpec=1,NALiquids + NASolids
!              C_Y3(InC3(ISpec,II)) = 0.0
!           end do
!           Y(NYCONCS+NY1+NY2+II)= C_Y3(NY3+II)
!        ENDIF
!     ENDIF
                 
        
!     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/w)) THEN
!        
!        IF(set_vel2(I_VEL)<w) THEN
!
!           C_Y3(II) = C_Y3(II) + max(-C_Y3(II)*.9999, &
!                &  (CY3multi(2,I_step+1,II+2)-C_Y3(II)) * (set_vel2(I_VEL)/w))
!
!           do  ISpec=1,NALiquids + NASolids
!              
!              C_Y3(InC3(ISpec,II)) = C_Y3(InC3(ISpec,II))  & 
!                   &  + max(-C_Y3(InC3(ISpec,II))*.9999,     &
!                   &  (CY3multi(2,I_step+1,InC3(ISpec,II)+2) &
!                   &  -C_Y3(InC3(ISpec,II))) * (set_vel2(I_VEL)/w))
!              
!           end do
!           
!           Y(NYCONCS+NY1+NY2+II)= C_Y3(NY3+II)
!           
!        ELSE
!           J_BOX=max(1,1+int(set_vel2(I_VEL)/w))
!             
!          
              !                    C_Y1(II) = C_Y1(II) + max(-C_Y1(II)*.9999,  & 
!           C_Y3(II) =CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2) &
!                &  + max(-CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2)*.9999,  & 
!                &  (CY3multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),II+2)  &
!                &  -CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w) &
!                &  ,II+2)) * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!          
!           do  ISpec=1,NALiquids + NASolids
              
!                       C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
                 !                            &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
!              C_Y3(InC3(ISpec,II)) = CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC3(ISpec,II)+2)  & 
!                   &  + max(-CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC3(ISpec,II)+2)*.9999,     &
!                   &  (CY3multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),  &
!                   &  InC3(ISpec,II)+2)-CY3multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC3(ISpec,II)+2)) &
!                   &  * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!           end do
!              
!           Y(NYCONCS+NY1+NY2+II)= C_Y3(NY3+II)
!        ENDIF
!     ENDIF
!  END DO
  
!  C_cl_tot=sum(C_Y2(1:NY2))
!  y(NODEs-6)=c_cl_tot


! Sedimentation for particles in C_Y4
  
!   DO II=1,NY4-1
!     I_VEL=NaBins+NY1+NY2+NY3+II
!     
!     IF(I_BOX==1) THEN
!        IF(set_vel2(I_VEL)<w) THEN
!           C_Y4(II) = C_Y4(II) -C_Y4(II)*(set_vel2(I_VEL)/w)
!           do  ISpec=1,NALiquids + NASolids
!              C_Y4(InC4(ISpec,II)) = C_Y4(InC4(ISpec,II))  & 
!                   &  -C_Y4(InC4(ISpec,II))*(set_vel2(I_VEL)/w)
!           end do
!           Y(NYCONCS+NY1+NY2+NY3+II)= C_Y4(NY4+II)
!        ELSE
!           C_Y4(II) =0.0
!           do  ISpec=1,NALiquids + NASolids
!              C_Y4(InC3(ISpec,II)) = 0.0
!           end do
!           Y(NYCONCS+NY1+NY2+NY3+II)= C_Y4(NY4+II)
!        ENDIF
!     ENDIF
                 
                 
        
!     IF(I_BOX>1.AND.I_step<I_STEP_MAX-1-int(set_vel2(I_VEL)/w)) THEN
!        
!        IF(set_vel2(I_VEL)<w) THEN
!
!           C_Y4(II) = C_Y4(II) + max(-C_Y4(II)*.9999, &
!                &  (CY4multi(2,I_step+1,II+2)-C_Y4(II)) * (set_vel2(I_VEL)/w))
!
!           do  ISpec=1,NALiquids + NASolids
!              
!              C_Y4(InC4(ISpec,II)) = C_Y4(InC4(ISpec,II))  & 
!                   &  + max(-C_Y3(InC4(ISpec,II))*.9999,     &
!                   &  (CY4multi(2,I_step+1,InC4(ISpec,II)+2) &
!                   &  -C_Y4(InC4(ISpec,II))) * (set_vel2(I_VEL)/w))
!              
!           end do
!           
!           Y(NYCONCS+NY1+NY2+NY3+II)= C_Y4(NY4+II)
!           
!        ELSE
!           J_BOX=max(1,1+int(set_vel2(I_VEL)/w))
             
          
!              !                    C_Y1(II) = C_Y1(II) + max(-C_Y1(II)*.9999,  & 
!           C_Y4(II) =CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2) &
!                &  + max(-CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),II+2)*.9999,  & 
!                &  (CY4multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),II+2)  &
!                &  -CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w) &
!                &  ,II+2)) * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!          
!           do  ISpec=1,NALiquids + NASolids
              
!!                       C_Y1(InCi(ISpec,II)) = C_Y1(InCi(ISpec,II))  & 
!                 !                            &  + max(-C_Y1(InCi(ISpec,II))*.9999,     &
!              C_Y4(InC4(ISpec,II)) = CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC4(ISpec,II)+2)  & 
!                   &  + max(-CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC4(ISpec,II)+2)*.9999,     &
!                   &  (CY4multi(J_BOX+1,I_step+1+int(set_vel2(I_VEL)/w),  &
!                   &  InC4(ISpec,II)+2)-CY4multi(J_BOX,I_step+int(set_vel2(I_VEL)/w),InC4(ISpec,II)+2)) &
!                   &  * (set_vel2(I_VEL)/w-int(set_vel2(I_VEL)/w)))
!           end do
!              
!           Y(NYCONCS+NY1+NY2+NY3+II)= C_Y4(NY4+II)
!        ENDIF
!     ENDIF
!  END DO
  
!  C_cl_tot=sum(C_Y2(1:NY2))
!  y(NODEs-6)=c_cl_tot
  


END SUBROUTINE SEDIMENTATION2











SUBROUTINE FIX_BINS(Y,n)
  USE headfile
  use precisi
  implicit none
  real(doubp) :: mass_tot, mass_tot_old,so_mass_sum
  integer, intent(inout) :: n
  real(doubp), dimension(n), intent(inout) :: Y
  integer :: I,ISpec,J,j_in_y2,j_step,JJ
  real(doubp) :: T_C, Tot_Wat, Cptot,Tot_wat_ice,frac,frac2
  real(doubp) :: rho_ice
!!$    
!!$  ********************************************************************
!!$  *   SUBROUTINE FIX_BINS IS TO MOVE PARTICLES FROM BINS TO ANOTHER  *
!!$  *   WHEN MOVING CENTER SYSTEM IS USED                              *
!!$  ********************************************************************
!!$
!!$


  mass_tot_old=0.0_doubp
  
  do I=1,NY1
     mass_tot_old=mass_tot_old+so_mass(I)*C_Y1(I)
     CALL RADIUSICE2(I)
  end do
   do I=1,NY2

     CALL RADIUS_Y2(I)
  end do
  !!
  !! First lets move the ice from one bin to the next one if it has
  !! grown due to condensation over the maximum size.
  !! Now only possible to move max 5 bins upwards. 
  !! 
  !! 
  !! Also added possibility for melting, i.e. to move ice particles from Y1 to Y2
  !! 
  !!
  



  do JJ=1,10
     do I=NY1-1,1,-1
        IF (C_Y1(I) > 0.0) THEN
           !        write(*,*) C_Y1(NY1+I)/C_Y1(I)
           if(R_ice(I)>R_Y1(I+1)) then
              C_Y1(I+1)= C_Y1(I+1)+C_Y1(I)
              C_Y1(I) = 0.0
              Y(NYConcs + I+1)= Y(NYConcs + I+1)+Y(NYConcs + I)
              Y(NYConcs + I)=0.0
              do  ISpec=1,NALiquids + NASolids
                 C_Y1(InCi(ISpec,I+1)) = C_Y1(InCi(ISpec,I+1)) + C_Y1(InCi(ISpec,I))
                 C_Y1(InCi(ISpec,I))=0.0
              end do
              F_num(I+1)=1
              F_num(I)=0
              CALL RADIUSICE2(I+1)
           elseif(R_ice(I)<R_Y1(I).AND.I>1) then
              IF(C_Y1(NY1+I)/C_Y1(I) < 1.0e-20) then       
                 
                 call RADIUS_ice_dry(I)
                 DO J=1,NBinsPerMode(1)
                    CALL RADIUS_DRY(J)
                    IF(R_dry(J).GT.R_ice_dry(I).OR.J.EQ.NBinsPerMode(1)) THEN
                       frac=max(0.0,(R_dry(J)**3-R_ice_dry(I)**3)/(R_dry(J)**3-R_dry(J-1)**3))
                       C(J-1)= C(J-1)+C_Y1(I)*frac
                       C(J)= C(J)+C_Y1(I)*(1.0-frac)
                       C_Y1(I) = 0.0
                       do  ISpec=1,NALiquids + NASolids
                          C(InC(ISpec,J-1)) = C(InC(ISpec,J-1)) + C_Y1(InCi(ISpec,I))*frac
                          C(InC(ISpec,J)) = C(InC(ISpec,J)) + C_Y1(InCi(ISpec,I))*(1.0-frac)
                          C_Y1(InCi(ISpec,I))=0.0
                       end do
                       Ctot=sum(C(1:NaBins))
                    
                       GOTO 10
                    ENDIF
                    
                 END DO
!        C(NaBins)= C(NaBins)+C_Y2(1)
!        C_Y2(1) = 0.0
!        do  ISpec=1,NALiquids + NASolids
!           C(InC(ISpec,NaBins)) = C(InC(ISpec,NaBins)) + C_Y2(InCl(ISpec,1))
!           C_Y2(InCl(ISpec,1))=0.0
!        end do
10               Ctot=sum(C(1:NaBins))
                 c_ice_tot=sum(C_Y1(1:NY1))
                 Y(NODES-1)= Ctot
                 Y(NODES-5) = C_ice_tot
                 

              else
                 C_Y1(I-1)= C_Y1(I-1)+C_Y1(I)
                 C_Y1(I) = 0.0
                 Y(NYConcs + I-1)= Y(NYConcs + I-1)+Y(NYConcs + I)
                 Y(NYConcs + I)=0.0
                 do  ISpec=1,NALiquids + NASolids
                    C_Y1(InCi(ISpec,I-1)) = C_Y1(InCi(ISpec,I-1)) + C_Y1(InCi(ISpec,I))
                    C_Y1(InCi(ISpec,I))=0.0
                 end do
                 F_num(I-1)=1
                 F_num(I)=0
                 CALL RADIUSICE2(I-1)
                 CALL RADIUSICE2(I)
              endif
           elseif(T>274.15.AND.cloud_mod == 1 ) then           !Move particles from Y1 to Y2        
              j_in_y2=NY2/2
              finder: do j_step=1,NY2/2+1 
                 if (R_Y2(j_in_y2) < R_ice(I)*(rho_ice(T)*.001)**(1./3.) .and. & 
                      R_Y2(j_in_y2-1) > R_ice(I)*(rho_ice(T)*.001)**(1./3.) ) then
                    j_in_y2 = j_in_y2-1
                    exit finder
                 endif
                 
                 if (R_Y2(j_in_y2) > R_ice(I)*(rho_ice(T)*.001)**(1./3.)) then
                    j_in_y2=max(2,j_in_y2-1)
                 else
                    j_in_y2=min(j_in_y2+1,NY2-1)
                 endif
              end do finder
              
              j_in_y2=max(1,j_in_y2)
              C_Y2(j_in_y2)= C_Y2(j_in_y2)+C_Y1(I)*min(1.0,DtOut/500.)
              C_Y1(I) = C_Y1(I)-C_Y1(I)*min(1.0,DtOut/500.)
              Y(NYConcs + NY1+j_in_y2)= Y(NYConcs+ NY1 + j_in_y2)+Y(NYConcs+ I)*min(1.0,DtOut/500.)
              Y(NYConcs + I)=   Y(NYConcs + I)-Y(NYConcs + I)*min(1.0,DtOut/500.)
              do  ISpec=1,NALiquids + NASolids
                 C_Y2(InCl(ISpec,j_in_y2)) = C_Y2(InCl(ISpec,j_in_y2)) + C_Y1(InCi(ISpec,I))*min(1.0,DtOut/500.)
                 C_Y1(InCi(ISpec,I))=C_Y1(InCi(ISpec,I))-C_Y1(InCi(ISpec,I))*min(1.0,DtOut/500.)
              end do
!              write(*,*)  C_Y2(j_in_y2),C_Y2(j_in_y2+NY2),j_in_y2
!              write(*,*)   C_Y1(I), C_Y1(I+NY1),I
!              pause
!              F_num(I)=0
              CALL RADIUS_Y2(j_in_y2)
           endif
           
        ENDIF
     enddo
     
  enddo
!
! IF ICE PARTICLE HAS GROWN OUT OF GRID IT IS REMOVED
!
  if(R_ice(NY1)>R_Y1(NY1+1)) then
     C_Y1(NY1) = 0.0
     Y(NYConcs + NY1)=0.0
     do  ISpec=1,NALiquids + NASolids
        C_Y1(InCi(ISpec,NY1))=0.0
     end do
     F_num(NY1)=0
     CALL RADIUSICE2(NY1)
  ENDIF

  
!
! New total ice mass
!
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
  do I=1,NY1
     mass_tot=mass_tot+so_mass(I)*C_Y1(I)
  end do
  
  IF(mass_tot.NE.mass_tot_old) THEN
     T_C=T-T_ntp
     HeatL_liq2ice=(4.1752_doubp*T_C+0.5_doubp*(-11.319e-3_doubp)*T_C**2+(1./3.)* &
          &     (-9.7215e-5_doubp)*T_C**3+(1./4.)*(1.8315e-5_doubp)*T_C**4+(1./5.)* &
          &     (1.1354e-6_doubp)*T_C**5-2.1046_doubp*T_C-0.0037_doubp*T_C**2+&
          &     334.0_doubp)*(mass_tot-mass_tot_old)
!!$
     Tot_wat=sum(C(Inc(IH2Ol,1:NABins)))
     IF(cloud_mod ==1)  Tot_wat=sum(C(Inc(IH2Ol,1:NABins)))+sum(C_Y2(Incl(IH2Ol,1:NY2)))
     Tot_wat_ice=sum(C_Y1(Inci(IH2Ol,1:NFBins)))
!!$C
     Cptot = (C(InC(IH2Ogi,NABins)) * WtGas(IH2Og) * 1.95_doubp + &! water gas
          &   Rho * CpAir * 1.0e-3_doubp + &! air
          &   Tot_Wat*WtMol(IH2Ol)*CpWat*1.0e-3_doubp + &! liquid water
          &   Tot_Wat_ice*WtMol(IH2Ol)*CpIce*1.0e-3_doubp)!/ ! ice
!!$
     T=T+HeatL_liq2ice/Cptot
!     write(*,*) HeatL_liq2ice/Cptot,mass_tot-mass_tot_old
!     pause
     Y(NODEs)= Y(NODEs)+HeatL_liq2ice/Cptot

  ENDIF
 

 !! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC !!
!      THIS PART BETTER TO MOVE INTO OWN SUBROUTINE 
!
! Move cloud droplets in C_Y2 from one bin to an other if needed.
! This is done based on wet radius
! At the moment possible to move to next or to bin after that. 
! Time step too long if more needed
!
  if(cloud_mod == 1) then
     do J=1,3
        do I=NY2-1,2,-1
           IF (C_Y2(I) > 0.0) THEN
           
!           IF(R_CL(I)>R_Y2(I+1).aND.I_BOX>10) THEN
!              write(*,*) I,R_CL(I),R_Y2(I+1),R_Y2(I+2)
!              pause
!           ENDIF
!              write(*,*) R_CL(I),R_Y2(I),R_Y2(I+1)
              if(R_CL(I)>R_Y2(I+1)) then
                 C_Y2(I+1)= C_Y2(I+1)+C_Y2(I)
                 C_Y2(I) = 0.0
                 Y(NYConcs + NY1+I+1)= Y(NYConcs+ NY1 + I+1)+Y(NYConcs + NY1+ I)
                 Y(NYConcs + I+ NY1)=0.0
                 do  ISpec=1,NALiquids + NASolids
                    C_Y2(InCl(ISpec,I+1)) = C_Y2(InCl(ISpec,I+1)) + C_Y2(InCl(ISpec,I))
                    C_Y2(InCl(ISpec,I))=0.0
                 end do
                 CALL RADIUS_Y2(I+1)
                 CALL RADIUS_Y2(I)
!                 write(*,*) 'here'
              endif
              
!              if(R_CL(I)>R_Y2(I)) then
                
!                 frac=min(2.0,(0.5*R_Y2(I+1)**3+0.5*R_Y2(I)**3)/R_CL(I)**3)
!                 frac2=min(0.1,(R_CL(I)-R_Y2(I))/(R_Y2(I+1)-R_Y2(I)))
!                 frac=R_Y2(I+1)**2/R_Y2(I)**2*frac2
!                 write(*,*)  I,frac2, frac,C_Y2(I),(R_CL(I)-R_Y2(I))/(R_Y2(I+1)-R_Y2(I))
!                 
!                 write(*,*) I,frac2,frac
!                 C_Y2(I+1)= C_Y2(I+1)+frac2*C_Y2(I)
!                 C_Y2(I) = (1.-frac2)*C_Y2(I)
!                 Y(NYConcs + NY1+I+1)= Y(NYConcs+ NY1 + I+1)+frac*Y(NYConcs + NY1+ I)
!                 Y(NYConcs + I+ NY1)=(1.0-frac)*Y(NYConcs + NY1+ I)
!                 do  ISpec=1,NALiquids + NASolids
!                    C_Y2(InCl(ISpec,I+1)) = C_Y2(InCl(ISpec,I+1)) + frac*C_Y2(InCl(ISpec,I))
!                    C_Y2(InCl(ISpec,I))=(1-frac)*C_Y2(InCl(ISpec,I))
!                 end do
!                 CALL RADIUS_Y2(I+1)
!                 CALL RADIUS_Y2(I)
!              
!              endif
              

           ENDIF
        end do
     end do


     IF (C_Y2(NY2) > 0.0) THEN
        if(R_CL(NY2)>R_Y2(NY2+1)) then
           C_Y2(NY2) = 0.0
           Y(NYConcs + NY2+ NY1)=0.0
           do  ISpec=1,NALiquids + NASolids
              C_Y2(InCl(ISpec,NY2))=0.0
           end do
        endif
        
     ENDIF


!
! Move evaporating cloud droplets from C_Y2 to C. 
!
!!!!!!!!!!!!!! REDO THIS, PROBABLY BETTER TO DISTRIBUTE BETWEEN TO TWO BINS !!!!!!!!!!!!!!
!

!
! Move evaporating cloud droplets in C_Y2 to smaller bin. 
! Cloud droplet evaporation can be fast, so at the moment moving 3 bins 
! downwards possible (outermost loop).

     do j=1,7
        do I=2,NY2
           IF (C_Y2(I) > 0.0) THEN
              if(R_CL(I)<0.5*R_Y2(I-1)+0.5*R_Y2(I)) then
                 C_Y2(I-1)= C_Y2(I-1)+C_Y2(I)
                 C_Y2(I) = 0.0
                 Y(NYConcs + NY1 +I-1)= Y(NYConcs + NY1 + I-1)+Y(NYConcs +NY1+ I)
                 Y(NYConcs + NY1+I)=0.0
                 do  ISpec=1,NALiquids + NASolids
                    C_Y2(InCl(ISpec,I-1)) = C_Y2(InCl(ISpec,I-1)) + C_Y2(InCl(ISpec,I))
                    C_Y2(InCl(ISpec,I))=0.0
                 end do
                 CALL RADIUS_Y2(I-1)
                 CALL RADIUS_Y2(I)
              endif
           ENDIF
        enddo
     enddo
     c_cl_tot=sum(C_Y2(1:NY2))
     Y(NODES-6) = C_cl_tot
  endif
  

  
     if(R_CL(1)<R_Y2(1)) then
        call RADIUS_Y2dry(1)
        
        DO I=1,NBinsPerMode(1)
           CALL RADIUS_DRY(I)
           IF(R_dry(I).GT.R_CLdry(1).OR.I.EQ.NBinsPerMode(1)) THEN
              frac=max(0.0,(R_dry(I)**3-R_CLdry(1)**3)/(R_dry(I)**3-R_dry(I-1)**3))
              C(I-1)= C(I-1)+C_Y2(1)*frac
              C(I)= C(I)+C_Y2(1)*(1.0-frac)
              C_Y2(1) = 0.0
              do  ISpec=1,NALiquids + NASolids
                 C(InC(ISpec,I-1)) = C(InC(ISpec,I-1)) + C_Y2(InCl(ISpec,1))*frac
                 C(InC(ISpec,I)) = C(InC(ISpec,I)) + C_Y2(InCl(ISpec,1))*(1.0-frac)
                 C_Y2(InCl(ISpec,1))=0.0
              end do
              Ctot=sum(C(1:NaBins))
           
              GOTO 100
           ENDIF
           
        END DO

!        C(NaBins)= C(NaBins)+C_Y2(1)
!        C_Y2(1) = 0.0
!        do  ISpec=1,NALiquids + NASolids
!           C(InC(ISpec,NaBins)) = C(InC(ISpec,NaBins)) + C_Y2(InCl(ISpec,1))
!           C_Y2(InCl(ISpec,1))=0.0
!        end do

        
100  endif

     DO J = 1, NY2
        CALL RADIUS_Y2(J)
     END DO

     Ctot=sum(C(1:NaBins))
     c_cl_tot=sum(C_Y2(1:NY2))
     c_ice_tot=sum(C_Y1(1:NY1))
     Y(NODES-5) = C_ice_tot
     Y(NODES-1)= Ctot
     Y(NODES-6) = C_cl_tot
  


ENDSUBROUTINE FIX_BINS
