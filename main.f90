PROGRAM CMODEL
!!$C     
!
!  SAMI ROMAKKANIEMI 25.9.2007
!
!  THIS IS AN ADIABATIC AIR PARCEL MODEL TO SIMULATE CLOUD DROPLET FORMATION 
!  AND ICE PARTICLE FORMATION.
!
!
!  INCLUDED IN THE MODEL: CONDENSATION OF WATER, OTHER GASES NEED SOME WORK
!                       : COAGULATION/COALESCENCE BETWEEN DIFFERENT PARTICLES AND DROPS
!                       : ICE PARTICLE FORMATION, HETEROGENEOUS AND HOMOGENEOUS
!
!  INPUT FILES: aerosol.dat     AEROSOL SIZE DISTRIBUTION
!             : environ.dat     LENGTH AND INITIAL VALUES OF SIMULATION
!             : freez.dat       FREEZING OPTIONS
!             : liquids.dat     LIQUID PROPERTIES
!             : solids.dat      SOLID PROPERTIES
!             : gases.dat       GAS PROPERTIES
!             : chemset.dat     POSSIBLE CHEMICAL REACTIONS, NOT WORKING AT THE MOMENT
!             : thermoset.dat   IS THIS STILL NEEDED????
!
!
!  CODE USES DLSODE TO SOLVE CONDENSATION OF WATER. COAGULATION/COALESCENCE AND FREEZING ARE SOLVED
!  OUTSIDE DLSODE FOR EVERY TIMESTEP (DTOUT). THIS CAUSES UNACCURACY IF TOO LONG TIME STEP IS USED.
!
!  COMMON VARIABLES IN headfile.f90
!
!  FIVE DIFFERENT AEROSOL PARTICLE SIZE GRIDS (concentration vectors): 
!          C FOR AEROSOL PARTICLES, FULLY MOVING
!          C_Y1 FOR ICE PARTICLES, MOVING CENTER
!          C_Y2 FOR CLOUD DROPLETS (FORMED THROUGH COALESCENCE), MOVING CENTER
!          C_Y3 FOR SOMETHING
!          C_Y4 FOR SOMETHING
!
!  Initially all particles are in C, they are moved to C_Y1 due to freezing or to C_Y2
!  due to coalescence between cloud droplets. Cloud droplets are moved back   
!  from C_Y2 to C due to evaporation. At the moment ice particle melting
!  happens instantly at some temeperature
!
!  Sedimentation of particles in C_Y1 and C_Y2 can be estimated, and it is possible to run 
!  code in columm mode. However, this part is far from ready. Some issues with mass conservation...  
!
!  At the moment only liquid phase thermodynamic working is AIM2, but model is very slow with that.
!  For testing purpose it is better to use only mole fractions.  This is easiest to do by setting
!  variable IF_Thermo(IBin) to 0 for all size bins. That is now done in driver.f90.
!
!
!  Work is in progress to describe all most important variables in headfile. Size of the 
!  variables CY1multi and CY2multi are quite big, and might cause problems.
!
  USE headfile
  use precisi
  use aikaa_mittaavat_parametrit
  implicit none
  real(doubp) :: alku_e, loppu_e, w_init, v_apu, pwe
  real(doubp) :: phno3, phcl, pH2Oeq,CNH42SO4,sum_ions,C_temp
  integer :: I, J, NEq,index_ice, L3,LRAD
  integer :: IFEQUISOLV,Ispec

!  write(*,*) '! COAGULATION KERNELS NEED TO BE TESTED !'
!  write(*,*)
!  write(*,*) '!!  COAGULATION SHOULD BE FIXXED SO THAT AEROSOL SIZE IS NOT THE ONE AT THE END OT TIME STEP BUT ONE IN THE MIDDLE !!'
!  write(*,*)
!  write(*,*) '!!  COAGULATION KERNELS AND SEDIMENTATION VELOCITY FOR ICE PARTICLES NEED TO BE WRITTEN  !!'
! write(*,*)
 write(*,*) '!! KELVIN EFFECT FOR ICE PARTICLES IS TURNED OFF. FOR SOME REASON THERE WAS EVAPORATION OF ICE PARTICLES!!'
 write(*,*)
 write(*,*) '!! NOTE, AT THE MOMENT CALL TO THERMODYNAMICS IS COMMENTED AND MOLE FRACTIONS ONLY USED!!'
 write(*,*)
 write(*,*) '!! SOME PROBLEMS WITH SIZE DISTRIBUTIONS!! FOR MULTIPLE MODES THE NUMBER AND MEAN SIZE IN A BIN NOW ITERATED'
 write(*,*)
 write(*,*) '!! COAGULATION KERNEL LIMITED'
 write(*,*)
 write(*,*) '!! IN THE CASE OF DRIZZLE FORMATION CARE SHOULD BE TAKEN TO KEEP TIME STEP SHORT ENOUGH.'  
 write(*,*) '!!INSIDE COAGULATION.F90 POSSIBLE SHORTEN TIME STEP USED TO ESTIMATE COAGULATION. '
 write(*,*)
! write(*,*) '!! FORMAT OF OUTPUT FIELS COULD BE BETTER. THE SHORTER TIME STEP FOR ICE NOT WORKING PROPERLY'
! write(*,*)
 write(*,*) '!! ISTATE2=1 NOT WORKING WHEN ICE FORMED IN CASE OF NO COAGULATION'
 write(*,*)
 

 CALL FileSet(0)     !Open the input files
 CALL ReadFile(0)       !Read in those files
 WRITE(*,*) n_box 
 CLOSE(KENV)

 I_step_max=0

 DO I_BOX=1,N_box**IF_BOXES
    DtOut_used = DtOut
  CALL CPU_TIME(BEGIN_TIME)
  EQ_CALCULATIONS=1 ! first equilibrium
!!$C
  ALN10     = LOG(10.0_doubp)
  ALN101    = 3.0_doubp / ALN10
  ALOGEXP     = LOG10(2.718281828_doubp)
!!$C
!!$C
  I_step=1
  CALL FileSet(1)     !Open the input files
  CALL ReadFile(1)    !Read in those files
  call out_dat_print  !Print some values on screen
  CALL Liquids        !Read data from liquids.dat
  CALL Solids         !Read data from solids.dat
  CALL Gases          !Read data from gases.dat
  CALL Families       !Set up families if needed
  CALL ChemSet        !Read data for chemistry 
  CALL SizeDist       !Make aerosol size distribution
  CALL SetIndex       !Make index matrices
  CALL FileSet(2)     !Open the output files


  f_act_max=0.0
  c_cl_tot=0.0
  C_ice_tot=0.0

  IF_NUC_FEX = 0
  c_ice2=1.0e-6
  IF_Change(1:NaBins)=1
!  RH=min(RH,85.0+I_Box/10.0)
!  write(*,*) RH
  !  totsec=min(20*I_Box+2000.1,FINSEC)
  totsec=150/w
  DO  LRAD = 1,NABins
     CALL RADIUS(LRAD)
  END DO

  
  CALL CPU_TIME(TIME_IN)


!
! Here thermodynamics is called once just to check everything is fine
!
! Only Clegg works currently
!

  LFirst=1
  Last=NABins

  write(*,*) Ithermo_mod
  select case(Ithermo_mod)
  case(1) ! Luo's thermodynamics

     CALL s_luo(T,0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
     C(InC(IH2Ogi,NABins)) = RH * pwe  * 100.0_doubp / (R * T) * 1.e-8_doubp
     DO  I=1,NABins
        v_apu=0.0_doubp
        do J=1,NALiquids-1
           v_apu=v_apu+C(InC(IH2Ol,I)+NABIns*(J))
        end do
        C(InC(IH2Ol,I))=v_apu/(1.0_doubp-RH/100.0_doubp)
     end DO    
!     CALL EQLB !(LFirst,Last)

!  case(2) !Jacobson thermo 
!     C(InC(IH2Ogi,NABins)) = RH * pH2Oeq(T) * 101325.0_doubp &
!          & / (R * T) * 1.e-8_doubp
!     IFEQUISOLV=1
!     IF(IFEQUISOLV == 1) CALL Equilset !!!!!!!!
!     DO J = 1,NEEQUAT
!        MSTT  = STATE(IEQSPEC(J,1))
!        IF(IEQSPEC(J,2) == 0 .and. MSTT == 'G') THEN !!!!!! KORJAA!! BUG?: f -> f90???
!!$C     
!!$C     NEQGAS = INDEX OF THE EQUILIBRIUM REACTION USED TO CALCULATE
!!$C              THE CONCENTRATION AT THE DROPLET SURFACE
!!$C
           NEQGAS(IEQSPEC(J,1)-NAERTY) = J
!!$C
!        END IF
!        IF(FORMSPEC(IEQSPEC(J,1)) == 'NH3' .AND. FORMSPEC(IEQSPEC(J,2)) &
!             & == 'HNO3') THEN
!           NEQGAS(IEQSPEC(J,1)-NAERTY) = J
!        END IF
!        IF(FORMSPEC(IEQSPEC(J,1)) == 'NH3' .AND. FORMSPEC(IEQSPEC(J,2)) &
!             & == 'H+') THEN
!           NEQGAS(IEQSPEC(J,1)-NAERTY) = J
!        END IF
!     END DO
!     CALL AERPROC
!     DO LRAD = 1,NABins
!        CALL RADIUS(LRAD)
!     END DO
!     CALL AERPROC2 !(LFirst,Last)

  case(3) !Raatikainen

     C(InC(IH2Ogi,NABins)) = RH * pH2Oeq(T) * 101325.0_doubp &
     & / (R * T) * 1.e-8_doubp
!     CALL CMODEL_EQUILIBRIUM

  case(4) !Clegg

     ! For thermodynamics, different solids must be divided into ions.

     pwe         = pH2Oeq(T) * 101325.0_doubp
     C(NAERO + IH2Og) = RH * pwe / (R * T) * 1.e-8_doubp
     DO L3 = 1,NABins
        IF(INH42SO4s /= 0) THEN
           IF(C(InC(INH42SO4s,L3)) /= 0.0_doubp) THEN
              CNH42SO4 = C(InC(INH42SO4s,L3))
              C(InC(INH3l,L3)) = C(InC(INH3l,L3)) + CNH42SO4 * 2.0_doubp
              C(InC(ISO42l,L3)) = C(InC(ISO42l,L3)) + CNH42SO4
              C(InC(INH42SO4s,L3)) = 0.0_doubp
           END IF
        END IF
        IF(IH2SO4l /= 0) THEN
           IF(C(InC(IH2SO4l,L3)) /= 0.0_doubp) THEN
              C_temp = C(InC(IH2SO4l,L3))
              C(InC(IHPlus,L3)) = C(InC(IHPlus,L3)) +  C_temp * 2.0_doubp
              C(InC(ISO42l,L3)) = C(InC(ISO42l,L3)) + C_temp 
              C(InC(IH2SO4l,L3)) = 0.0_doubp
           END IF
        END IF
        IF(IHNO3l /= 0) THEN
           IF(C(InC(IHNO3l,L3)) /= 0.0_doubp) THEN
              C_temp = C(InC(IHNO3l,L3))
              C(InC(IHPlus,L3)) = C(InC(IHPlus,L3)) +  C_temp
              C(InC(INO3l,L3)) = C(InC(INO3l,L3)) + C_temp 
              C(InC(IHNO3l,L3)) = 0.0_doubp
           END IF
        END IF
         IF(INH4NO3s /= 0) THEN
           IF(C(InC(INH4NO3s,L3)) /= 0.0_doubp) THEN
              C_temp = C(InC(INH4NO3s,L3))
              C(InC(INO3l,L3)) = C(InC(INO3l,L3)) +  C_temp
              C(InC(INH3l,L3)) = C(InC(INH3l,L3)) + C_temp 
              C(InC(INH4NO3s,L3)) = 0.0_doubp
           END IF
        END IF
!        
!        
     end DO
!
!     CALL EQLB_cle(0.0)
  case default ! strange input
     write(*,*) 'check your Ithermo_mod input (environ.dat)'
     stop
  end select
  
!
! Only quess for the initial amount of water
!
  IF(Ithermo_mod == 4 .AND. IF_FREEZ_CONTACT ==1) THEN
     DO L3 = 1,NaBins
        sum_ions=0.0
        DO ISpec=2,NALiquids
           sum_ions= sum_ions+C(InC(Ispec,L3))
        END DO
        IF(IModeBin(L3) /= I_MODE_CONTACT) C(InC(1,L3))=0.01*RH*sum_ions/(1.0-RH*0.01)
     END DO
  ENDIF


  CALL CPU_TIME(TIME_OUT)

  TOTAL_INITIALIZE_TIME=TIME_OUT-BEGIN_TIME
  !
  FIRST_EQ_TIME=TIME_OUT-TIME_IN

  Adia_fex=100.0_doubp
  Adia_fex_old=100.0_doubp
!!$C
  NEq = NODEs
  NODEs_old=0
!!$C
!
! Initial size distribution can be equilibrated through condensation.
! Not needed if only dondensable gas is water.
!
  ind_driver=1
  if (IDRIVER_EQ .eq. 1) then
     write(*,*) 'Starting driver_equilibrium...'
     index_ice=ice_mod
     w_init=w
     ice_mod=0
     call CPU_TIME(alku_e)
     call driver_eq(NEq)
     call CPU_TIME(loppu_e)
     w=w_init
     IDRIVER_EQ=0
     ice_mod=index_ice
     print 133, loppu_e - alku_e
133  format(' driver_equilibrium done, time used = ', F8.3, ' s')
  end if
!!$C

  NEQ=NYConcs + NDyna+NY1+NY2*CLOUD_MOD
  NOdes=NEQ
!  write(*,*) 'here1'
!
!  CALL DRIVER, WHICH INCLUDES CALL TO ODE-SOLVER
!
  CALL CPU_TIME(TIME_IN)
  CALL DRIVER(NEq)
  CALL CPU_TIME(TIME_OUT)
  NEQ=NYConcs + NDyna+NY1 +NY2
  NOdes=NEQ
  write(*,*) 'Back from driver'

  if (ice_mod /= 0 .and. NFABins > 0) then 
     ! ice formation is on and there has been some freezing
     ! it's possible that there hasn't been any freezing and 
     ! we end the simulation (or ice form. is off)
     if (indx_ret == 1) then 
        ! we got the "signal" (indx_ret = 1) that some freezing
        ! has occured and now the cycling of driver starts...
        driverloop: do 
           WRITE(*,*) 'FREEZING'
           indx_ret=2
           ind_driver=2
           CALL Driver(NEq)
           if (indx_ret == 0) exit driverloop
        end do driverloop
     end if
  end if
!!$
  if (ice_mod == 1) CALL CHECK_ICE_BINS
  CALL INFO_PRINT
  
  

  DO I=min(I_BOX,25),1,-1
     CY1multi(I+1,1:N_step,1:2+NY1 + NY1 * (NALiquids + NASolids))=CY1multi(I,1:N_step,1:2+NY1 + NY1 * (NALiquids + NASolids))
     CY2multi(I+1,1:N_step,1:2+NY2 + NY2 * (NALiquids + NASolids))=CY2multi(I,1:N_step,1:2+NY2 + NY2 * (NALiquids + NASolids))
     CY3multi(I+1,1:N_step,1:2+NY3 + NY3 * (NALiquids + NASolids))=CY3multi(I,1:N_step,1:2+NY3 + NY3 * (NALiquids + NASolids))
     CY4multi(I+1,1:N_step,1:2+NY4 + NY4 * (NALiquids + NASolids))=CY4multi(I,1:N_step,1:2+NY4 + NY4 * (NALiquids + NASolids))
  END DO

  CLOSE(KGAS)
  CLOSE(KLIQ)
  CLOSE(KSLD)
  CLOSE(KSPC)
  CLOSE(KCHM)
  CLOSE(KENV)
  CLOSE(KFRE)
  CLOSE(KTRM)
  CLOSE(KEQQ)
  CLOSE(KACT)


  CLOSE(KOUT)
  CLOSE(IOUT)
  CLOSE(JAOU) 
END DO
!!$C
END PROGRAM CMODEL
!!$C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INFO_PRINT
  use aikaa_mittaavat_parametrit
  use headfile, ONLY : Ithermo_mod
  implicit none
  REAL :: END_TIME, TOTAL_TIME, TOT_CALC
  !
  ! Time in the end
  CALL CPU_TIME(END_TIME) 
  ! Total calculation time
  TOTAL_TIME=END_TIME-BEGIN_TIME 
  !
  ! Total time for calculations (no first equilibrium)
  TOT_CALC=TOTAL_TIME-TOTAL_INITIALIZE_TIME
  ! This includes time for the equilibrium model (TOTAL_EQMODEL_TIME) and
  ! for non-equilibrium calculations (TOT_CALC-TOTAL_EQMODEL_TIME)
  ! 
  ! Data for the first equilbrium:
  !	TOTAL_INITIALIZE_TIME	Initializing time for eq. and non-eq. models
  !	FIRST_EQ_TIME			Time (initialzing+first eq.) for eq. model only
  !
  WRITE(*,FMT='(/,A,/A)')	'CPU time (s) used during calculations:', &
       '--------------------------------------' 
  ! Model initializing time 
  WRITE(*,'(A,F10.4)')'Total initializing time: ',TOTAL_INITIALIZE_TIME 
  WRITE(*,'(A,TR6,F10.4)')'  initializing cloud model:',  &
       TOTAL_INITIALIZE_TIME-FIRST_EQ_TIME 
  WRITE(*,'(A,TR6,F10.4)')'  first equilibrium:       ',FIRST_EQ_TIME 
  ! Other calculations
  WRITE(*,'(A,F12.4)')'Total calculation time:',TOT_CALC
  WRITE(*,'(A,TR6,F12.4)')'  cloud model:           ',  &
       TOT_CALC-TOTAL_EQMODEL_TIME 
  WRITE(*,'(A,TR6,F12.4)')'  equilibrium model:     ',TOTAL_EQMODEL_TIME
  ! If data is available, print number of function calls and iterations
  WRITE(*,'(TR4,A,I8)') 'function calls:   ',EQ_CALCULATIONS
  ! Total time
  WRITE(*,'(A,/,A5,TR16,F20.2,/)') &
  '----------------------------------------------', &
  'Total',TOTAL_TIME
  !
END SUBROUTINE INFO_PRINT


SUBROUTINE CHECK_ICE_BINS
  use headfile
  implicit none

  if (screen(1) == 1 .and. screen(2) == 1) then
     write(*,*) 'subroutines new_bin_right and combine_closest_bins ',&
          & 'has been used -> indexes of ice size classes are "messed up"!!!'
  else if (screen(1) ==1 .and. screen(2) == 0) then
     write(*,*) 'subroutine new_bin_right ',&
          & 'has been used -> indexes of ice size classes are "messed up"!!!'
  else if (screen(1) == 0 .and. screen(2) == 1) then
     write(*,*) 'subroutine combine_closest_bins ',&
          & 'has been used -> indexes of ice size classes are "messed up"!!!'
  end if

END SUBROUTINE CHECK_ICE_BINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE out_dat_print
!!$C
  USE headfile
!!$C
  CHARACTER(len=7) :: hno3_fr, hno3_ad
  CHARACTER(len=10) :: nuc_method
  CHARACTER(len=8) :: thermo_model
!!$C
  select case(Ithermo_mod)
  case(1) ! Luo's thermodynamics
     thermo_model = 'Luo'
  case(2) !Jacobson
     thermo_model = 'Jacobson'
  case(3) !
     thermo_model = 'Joku'
  case(4) !
     thermo_model = 'Clegg'
  end select

  select case (ice_mod)
  case(0) !ice model off
!!$  IF (ice_mod == 0) THEN
     write(*,*) ' ' 
     write(*,*) 'Cloud model starting with the initial values:'
!!$C     output initial conditions (number of size bins, aerosol 
!!$C     number concentrations in size bins and dry radius)
     IF (IMModes == 1) THEN
!!$C     modes are inseparable
        WRITE(*,FMT=390)thermo_model,1,NABins
     ELSE
        WRITE(*,FMT=390)thermo_model,NAModes, NBinsPerMode(1:NAModes)
     END IF
!!$C     
     WRITE(*,FMT=400)NAGases,NABins+NAGases+3,TotSec,w!,CBin(1:NABins)
     WRITE(*,FMT=395)
!!$C     
     !     WRITE(*,FMT=410) Rp(1:NABins)
     !     WRITE(*,'(A1)')'%' ! the last comment line
!!$C     
395  FORMAT('--------------------------------------------------')
390  FORMAT('--------------------------------------------------'/, &
          &             '  Thermodynamical model =      ', A8/, &
          &             '  Number of separable modes = ',I4/, &
          &             '  Number of bins in modes =   ',10(I4,TR1))
400  FORMAT('  Number of gases =     ',I4/, &
          &             '  Record length =       ',I4/,  &! for one time step
          &             '  Simulation time =     ',F8.2/, &
          &             '  Vertical speed =      ',F7.3)!,5(TR2,E11.5)))!, &
     !         &             '  Number concentrations (#/cm-3) for size bins', &
     !          &             200(/,' ',5(TR2,E11.5)))
     !410  FORMAT('  Dry radius (m) for size bins', &
     !          &             200(/,' ',5(TR2,E11.5)))
!!$  ELSE
  case(1) !ice model on
!!$C     let's check some initial variables which ice model uses
     select case(rate_method)
     case(1) !Classic icenucleation
        nuc_method='Classic' 
     case(2) !Koop's parameterization
        nuc_method='Koop   '
     end select

!!$C 
     write(*,*) ' ' 
     write(*,*) 'Cloud model starting with the initial values:'
     IF (IMModes == 1) THEN
!!$C     modes are inseparable
        WRITE(*,FMT=391)thermo_model,1,NABins
     ELSE
        WRITE(*,FMT=391)thermo_model,NAModes, NBinsPerMode(1:NAModes)
     END IF
!!$C     
     WRITE(*,FMT=401)NAGases,NABins+NAGases+3,TotSec,w, &
          &             nuc_method,NFBins,hno3_fr,hno3_ad!,CBin(1:NABins)
     write(*,FMT=396)
!!$C     
     !     WRITE(*,FMT=411) Rp(1:NABins)
     !     WRITE(*,'(A1)')' ' ! the last comment line
!!$C   
396  FORMAT('--------------------------------------------------')
391  FORMAT('--------------------------------------------------'/, &
          &             '  Thermodynamical model =     ', A8/, &
          &             '  Number of separable modes = ',I4/, &
          &             '  Number of bins in modes =   ',10(I4,TR1))
401  FORMAT('  Number of gases =     ',I4/, &
          &             '  Record length =       ',I4/,  &! for one time step
          &             '  Simulation time =     ',F8.2/, &
          &             '  Vertical speed =      ',F7.3/, &
          &             '  Ice nucleation =      ',A7/, &
          &             '  Number of ice bins =  ',I4/, &
          &             '  HNO3 freez =          ',A3/, &
          &             '  HNO3 adsorption =     ',A3)!, 5(TR2,E11.5)))!, &
     !         &             '  Number concentrations (#/cm-3) for size bins', &
     !         &             200(/,' ',5(TR2,E11.5)))
     !411  FORMAT('  Dry radius (m) for size bins', &
     !          &             200(/,' ',5(TR2,E11.5)))
  end select
!!$C  END IF
!!$C
END SUBROUTINE out_dat_print


