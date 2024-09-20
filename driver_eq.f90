SUBROUTINE DRIVER_EQ(NEq)

  USE headfile
  use precisi
  use dvode_f90_m
  use funs, only : fex!, jex
  implicit none
!!$C
!!$C ******************************************************************************************
!!$C * SUBROUTINE DRIVER_EQ IS A (KINETIC) EQUILIBRIUM DRIVER ROUTINE FOR ODE-SOLVER LSODE    *
!!$C ******************************************************************************************
!!$C
!  EXTERNAL Fex, Jex
!!$  integer, intent(inout) :: NEq
!!$  real(doubp), dimension(22 +  9*NEq + NEq**2):: RWork
!!$  real(doubp), dimension(NEq) :: y, ATol
!!$  real(doubp) :: Time, Tout, StepSize, pw, pwe, pH2Oeq, phcl
!!$  real(doubp) :: RH_OLD, dt, RTol, w_init, temperature, phno3
!!$  INTEGER, dimension(20 + NEq) :: IWork
!!$  integer :: i_1, i_2, i_3, JJ, L3, IFam, IBin, J, JOUT
!!$  integer :: MF, IOpt, ITask, ITol, IState2, I, LIW, LRW 
  integer, intent(inout) :: NEq
  real(doubp), dimension(22 +  9*NEq + NEq**2) :: RWork
  real(doubp), dimension(NEq) :: ATol
  real(doubp), dimension(NEq + NFBins) :: Y
  real(doubp) :: Time, Tout, pw, pwe, pH2Oeq, phcl
  real(doubp) :: phno3
  real(doubp) :: temperature
  INTEGER, dimension(20 + NEq) :: IWork
  integer :: IBin, J, IDriv!, JOUT_init
  integer :: ITask, IState2, I, LIW, LRW ,i_1, i_2, i_3
  integer :: IFam, JJ, L3
  TYPE (VODE_OPTS) :: OPTIONS
  i_1=0
  i_2=0
  i_3=0
!!$C     neq    = number of first order ode-s.
!!$C     lrw    = declared length of rwork (in user-s dimension).
!!$C     liw    = declared length of iwork (in user-s dimension).
!!$C
  LRW    = 22 + 9 * NEq + NEq**2
  LIW    = 20 + NEq
  P = P - C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  CALL C2Y(y,NEQ)
!!$
  IFADIA = 0
!  w_init=w
!!$
  IF(IFADIA.EQ.0) THEN
     w = 0.e0_doubp
  ENDIF
!!$
  if (ice_mod == 1) then
     y(NEq - 6)    = ads_ice
     y(NEq - 5)    = C_ice_tot
  end if
  y(NEq - 4)       = Adia_fex    
  y(NEq - 3)       = z
  y(NEq - 2)       = P
  y(NEq - 1)       = Ctot
  y(NEq)           = T
!!$
!!$     time   = the initial value of time
!!$     tOut   = first point where output is desired
!!$ 
  time   = 0.0_doubp
  tOut   = time+DTOUT_eq    !1.d-5
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
  RWork(5) = 0.1e-15_doubp
  RWork(6) = StepMax_eq!1.e-0
  RWork(7) = 1.0e-20_doubp
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
!!$      ATol(1:NEq) = 0
  ATol(1:NEq) = AbsTol
!!$      ATol    = AbsTol
  ITask   = 1
!  IOpt    = 1
!  MF      = 22
!!$     
!!  dt = 0.0001_doubp
  OPTIONS = SET_OPTS(DENSE_J=.TRUE. ,&
       USER_SUPPLIED_JACOBIAN=.FALSE.,   &
       RELERR=RelTol,ABSERR_VECTOR=ATOL(1:NEq),  &
       MXSTEP=1000000,H0=RWork(5),HMAX=StepMax_eq)!,   
!!$
!!$      WRITE(*,*) C(1:NAero)
!!$      PAUSE 'driver'
!!$      if (indx_ret .eq. 0) then
!!$         JOUT_init=1
!!$         CALL OUTPUT(time,Y)
!!$      else 
!!$         JOUT_init=JOUT_help
!!$         indx_ret=0
!!$         y(1:NEq)=Y_init(1:NEq)
!!$c      end if
!!$
!!$      DO 40 JOUT = JOUT_init, TotSec/DTOUT+1
  DO !JOUT = 1, DRIV_EQ_TIME/DTOUT+1
      if (i_1 == 0) then
        if (time > (DRIV_EQ_TIME/4.0_doubp)) then 
           write(*,*) '25% done...'
           i_1=1
        end if
     end if
     if (i_2 == 0 .and. i_1 == 1 ) then
        if (time > (DRIV_EQ_TIME/2.0_doubp)) then 
           write(*,*) '50% done...'
           i_2=1
        end if
     end if
     if (i_3 == 0 .and. i_2 == 1) then
        if (time > (DRIV_EQ_TIME*(3./4.))) then 
           write(*,*) '75% done...'
           i_3=1
        end if
     end if
     CALL DVODE_F90(fex,NEq,Y(1:NEq),time,tout,ITASK,ISTATE2,OPTIONS)
!     CALL DLSODE(fex,NEq,y,time,tOut,ITol,RTol,ATol,ITask, &
!          & IState2,IOpt,RWork,LRW,IWork,LIW,jex,MF)
!!$     
     DO J=1,NYLiquids
        C(IY2C(J))=y(J)
     END DO
     IF (istate2 < 0) THEN
        write(*,90) istate2
        stop
     end IF

!!$         tOut = tOut + DtOut
!!$         WRITE(*,*) time, RH, z
!!$                     
!!$     CALL VISIB
!!$     
!!$     
     pw          = C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp
     select case(Ithermo_mod)
     case(1) !Luo's thermo
        temperature=T
        CALL s_luo(temperature,0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
        pwe=pwe/100.0_doubp
     case(2,3) !jacobson
        pwe         = pH2Oeq(T) * 101325.0_doubp
     end select
     
     RH          = pw / pwe
!!!     RH_OLD = RH
!!$
     DO IBin = 1,NABins
        Rp_old(IBin) = Rp(IBin)
     END DO
     Ctot = y(NODEs-1)
!!$
     do I=1,NABins
        C(I)=Ctot*CFrac(I)
     END do

     tOut = tOut + DTOUT_eq

!!$
!!$         CALL OUTPUT(time,Y)
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC from output.f CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     DO J=1,NYLiquids
        C(IY2C(J))=y(J)
     END DO
!!$     
     J = NYLiquids
     DO IFam = 1,NAFl
        IF (NAFMbr(IFam) >= 2) THEN
           DO L3 = 1,NABins
              J = J+1
              TotF(IFam) = 0.0e0
              DO JJ = 1, NFMbr(IFam)
                 IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
                    TotF(IFam) = TotF(IFam) + C(InC(IFMbr(IFam,JJ),L3))
                 END IF
              END DO
              DO JJ = 1, NFMbr(IFam)
                 IF(State(IFMbr(IFam,JJ)).NE.'G') THEN
                    C(InC(IFMbr(IFam,JJ),L3)) = C(InC(IFMbr(IFam,JJ),L3)) / & 
                         & TotF(IFam)*y(J) 
                 END IF
              END DO
           END DO
        END IF
     END DO
!!$     
     DO J = 1, NAGases
        C(NAero+J)=y(J + NYAero)
     END DO
!!$     
     z           = y(NODEs-3)
     P           = y(NODEs-2)
     Ctot        = y(NODEs-1)
     T           = y(NODEs)
!!$      write(*,*) y
!!$     
!!$     pw  = WATER VAPOR PRESSURE
!!$     pwe = WATER SATURATION VAPOR PRESSURE, Seinfeld & Pandis 1998 (page 49)
!!$     
!     pw          = C(NAero + IH2Og) * R * T * 1.0e6
!     temperature = T
!     CALL s_luo(temperature,0.,0.,0.,0.,pwe,phno3,phcl)
!!!     pwe         = pH2Oeq(T) * 101325.0
!!$     WRITE(*,*) pwe
!!$     PAUSE
!     RH          = pw / pwe
!!$     
     DO I = 1, NABins
        CALL RADIUS(I)
     END DO
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCend output.f CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     DO I = 1, NABins
        CFrac(I) = C(I) / CTot
     END DO
     if (time > DRIV_EQ_TIME) return
  END DO

!  StepSize=RWork(11)
!!$     
!!$      write(*,*) time, stepsize,iwork(16),aerintv
!!$     
  DO I=1, NEq
     IF(y(I) < 0.e0) THEN
        WRITE(*,*) I,y(I)
        write(*,*) 'y neg in driver_eq'
        stop
!        pause 'hah'
!        C(I)=0.0 
     ENDIF
  END DO

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
80   write(*,90) istate2
90   format(///22h error halt.. istate =,i3)
!!$     
  END if
!!$666 FORMAT(99(E11.4,1X))
!!$     
  stop
end SUBROUTINE DRIVER_EQ
