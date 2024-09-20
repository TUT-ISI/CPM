SUBROUTINE ReadFile(I)
  use headfile
  use precisi
  implicit none
  integer, intent(in) :: I
  real(doubp) :: Temperature, Pressure, RelativeHumidity, HeightOfCloudBase
  real(doubp) :: VerticalSpeed, RelativeTolerance, AbsoluteTolerance
  real(doubp) :: Factor, Wt1, Wt2, Wt3, Wt4, Wt5, Wt6, VolMixR
  real(doubp) :: WtTot1, WtTot2, WtTot3, WtTot4, WtTot5, WtTot6
  integer :: IGas, ISPec, NeverMindThis
!  integer :: 
!  integer :: 
!!$
!!$     Read initial conditions for environmental
!!$     variables from 'environ.dat'
!!$
  NAMELIST /Conditions/ &
       &     Temperature, &
       &     Pressure, &
       &     RelativeHumidity, &
       &     HeightOfCloudBase, &
       &     VerticalSpeed, &
       &     i_out, &
       &     IfCloud, &
       &     IfCoag, &
       &     IfFreez
!!$
  NAMELIST /Timing/ &
       &     FinHour, &
       &     FinMin,    & 
       &     FinSec,  &
       &     DtOut, &
       &     DtOutIce, &
       &     t_file
!!$
  NAMELIST /ODESolver/ &
       &     RelativeTolerance, &
       &     AbsoluteTolerance, &
       &     StepMax
!!$
  NAMELIST /THERMO_MOD/ &
       &     Ithermo_mod
!!$
  NAMELIST /DRIVER_EQUILIBRIUM/ &
       &     IDRIVER_EQ, &
       &     DRIV_EQ_TIME, &
       &     DTOUT_eq, &
       &     StepMax_eq
  NAMELIST /MULTIBOX/ &
       &     CLOUD_MOD,&
       &     rmin_cloud, &
       &     IF_BOXES, &
       &     N_BOX
!!$ These are for ice model
  NAMELIST /Icemodo/ &
       &     ice_mod
  NAMELIST /Nucleation/ &
       &     rate_method
  NAMELIST /Probability/ &
       &     upper_limit, &
       &     ower_limit
  NAMELIST /IceStepStart/ &
       &     ice_step_start
  NAMELIST /Divider/ &
       &     dr_c, &
       &     freez_cons_lo, &
       &     C_min_limit, &
       &     y_min_limit
  NAMELIST /Satvappress/ &
       &     e_sat_method
  NAMELIST /Heterogeneous/ &
       &     IF_FREEZ_CONTACT, &
       &     I_MODE_CONTACT, &
       &     IF_FREEZ_DEPOS, &
       &     I_MODE_DEPOS, &
       &     IF_FREEZ_IMMERSION, &
       &     I_MODE_IMMERSION, &
       &     r_min_hm  , &
       &     I_HM_METHOD, &
       &     I_HM_NEW
!!$
  DO
     READ(KENV,*) Heading 
     IF(Heading == 'BEGIN') EXIT
  END DO
!!$ 
  READ(KENV,Conditions) 
!!$ 
  T      = Temperature 
  P      = Pressure * 1.e-2_doubp ! convert from Pa to mbar 
  RH     = RelativeHumidity 
  z      = HeightOfCloudBase 
  w      = VerticalSpeed 
!!$ 
  READ(KENV,Timing) 
  TotSec = FinHour * 3600.0_doubp + FinMin * 60.0_doubp + FinSec 
  READ(KENV,ODESolver)       
  READ(KENV,DRIVER_EQUILIBRIUM) 
  READ(KENV,MULTIBOX) 
  IF(I==0) return
  READ(KTRM,THERMO_MOD) 
!!$ 
  RelTol = RelativeTolerance  
  AbsTol = AbsoluteTolerance 
!!$
!!$     Read size distribution data from 'aerosol.dat'
!!$
  DO
     READ(KSPC,*) Heading
     IF(Heading == 'Size') EXIT
  END DO
!!$
  NAModes  = 0
  READ(KSPC,*) Heading     

  RLOOP: DO
     NAModes  = NAModes + 1
     READ(KSPC,1000) Heading, GMD(NAModes)
     GMD(NAModes) = GMD(NAModes) * 1.e-9_doubp
     READ(KSPC,1000) Heading, STD(NAModes)
     READ(KSPC,1000) Heading, CMode(NAModes)
     READ(KSPC,1500) Heading, NBinsPerMode(NAModes)
     READ(KSPC,1000) Heading
     READ(KSPC,1000) Heading
!!$
     IF(Heading /= ' '.AND. NAModes < NModes) CYCLE RLOOP
     EXIT RLOOP
  END DO RLOOP
  DO
     READ(KSPC,1500) Heading, NABins
     IF(NABins /= 0) EXIT
  END DO
  READ(KSPC,1750) Heading, JST
  IF(JST == 'manual') THEN
     READ(KSPC,1500) Heading, NAModes
  ELSE
     READ(KSPC,1500) Heading, NeverMindThis
  END IF
  READ(KSPC,1750) Heading, JST
  IF(JST == 'yes') THEN
     IMMODES = 0
  ELSE
     IMMODES = 1
  END IF
  IF(IMMODES == 0 .AND. NAModes > 1) THEN
     NABins = 0
     DO IMode = 1, NAModes
        NABins = NABins + NBinsPerMode(IMode)
     END DO
  END IF
!!$
  IF(NAModes == 1) NABins = NBinsPerMode(1)
!!$ 
!!$     Read data for liquids
!!$
  ISpec  = 0
  WtTot1 = 0.0_doubp
  WtTot2 = 0.0_doubp
  WtTot3 = 0.0_doubp
  WtTot4 = 0.0_doubp
  WtTot5 = 0.0_doubp
  WtTot6 = 0.0_doubp
  NALiquids = 0
  DO
     READ(KSPC,*) Heading
     IF(Heading == 'Liquid') EXIT
  END DO
  READ(KSPC,*)
  LOOP400: DO
     READ(KSPC,3000) JST, Formula, Wt1, Wt2, Wt3, Wt4, Wt5, Wt6
!!$
     IF(JST == 'A') THEN
        NALiquids       = NALiquids + 1
        ISpec           = ISpec + 1
        FormSpec(ISpec) = Formula
        NameSpec(ISpec) = 'none'
!!$
        WtFrac(ISpec,1) = Wt1
        WtFrac(ISpec,2) = Wt2
        WtFrac(ISpec,3) = Wt3
        WtFrac(ISpec,4) = Wt4
        WtFrac(ISpec,5) = Wt5 
        WtFrac(ISpec,6) = Wt6 
        WtTot1          = WtTot1 + Wt1
        WtTot2          = WtTot2 + Wt2
        WtTot3          = WtTot3 + Wt3
        WtTot4          = WtTot4 + Wt4
        WtTot5          = WtTot5 + Wt5 
        WtTot6          = WtTot6 + Wt6 
!!$
        NRLo(ISpec)     = ISpec*NABins+1
     END IF
!!$
     IF(JST == 'A' .OR. JST == 'D') CYCLE LOOP400
     EXIT LOOP400
  END DO LOOP400
!!$
!!$     Read data for solids
!!$
  NASolids = 0
  DO
     READ(KSPC,*) Heading
     IF(Heading == 'Solid') EXIT
  END DO
!!$
  LOOP600: DO
     READ(KSPC,3000) JST, Formula, Wt1, Wt2, Wt3, Wt4, Wt5, Wt6
     IF(JST == 'A') THEN
        NASolids        = NASolids + 1
        ISpec           = ISpec + 1
        FormSpec(ISpec) = Formula
        NameSpec(ISpec) = 'none'
        WtFrac(ISpec,1) = Wt1
        WtFrac(ISpec,2) = Wt2
        WtFrac(ISpec,3) = Wt3
        WtFrac(ISpec,4) = Wt4 
        WtFrac(ISpec,5) = Wt5 
        WtFrac(ISpec,6) = Wt6 
        WtTot1          = WtTot1 + Wt1
        WtTot2          = WtTot2 + Wt2
        WtTot3          = WtTot3 + Wt3
        WtTot4          = WtTot4 + Wt4 
        WtTot5          = WtTot5 + Wt5 
        WtTot6          = WtTot6 + Wt6
        NRLo(ISpec)     = ISpec*NABins+1
     END IF
!!$
     IF(JST == 'A' .OR. JST == 'D') CYCLE LOOP600
     EXIT LOOP600
  END DO LOOP600
!!$
  DO
     READ(KSPC,*) Heading
     IF(Heading == 'Gases') EXIT
  END DO
!!$
!!$     Water
!!$
  READ(KSPC,*) JST, Formula
  ISpec              = ISpec + 1
  NAGases            = 1
  FormSpec(ISpec)    = Formula
!!$
!!$     Read state, chemical formula and volume mixing ratio
!!$     of gas phase compounds
!!$
  LOOP200: DO
     READ(KSPC,2000) JST,Formula,VolMixR,unit
     IF(JST == 'A') THEN
        NAGases         = NAGases + 1
        ISpec           = ISpec + 1
        FormSpec(ISpec) = Formula
        NameSpec(ISpec) = 'none'
        IF(unit == 'ppm') THEN
           Factor       = 1.e-6_doubp
        ELSE IF(unit == 'ppb') THEN
           Factor       = 1.e-9_doubp
        ELSE IF(unit =='ppt') THEN
           Factor       = 1.e-12_doubp
        ELSE
           WRITE(*,*) 'ERROR: Unknown unit:',unit
           WRITE(*,*) '       for volume mixing ratio of ',Formula
           STOP
        END IF
!!$
!!$     Gas phase concentration (mol/cm3)
!!$
        IGas = NABins * (1 + NALiquids + NASolids) + NAGases
        C(IGas) = VolMixR * Factor * P * 100.0_doubp/ (R * T) &
             &    * 1.e-6_doubp
!!$
     END IF
!!$         
     IF(JST == 'A' .OR. JST == 'D') CYCLE LOOP200
     EXIT LOOP200
  END DO LOOP200

!!$     Normalize weight fractions
!!$
  DO ISpec = 1, NALiquids + NASolids
     WtFrac(ISpec,1) = WtFrac(ISpec,1) / WtTot1
     WtFrac(ISpec,2) = WtFrac(ISpec,2) / WtTot2
     WtFrac(ISpec,3) = WtFrac(ISpec,3) / WtTot3
     WtFrac(ISpec,4) = WtFrac(ISpec,4) / WtTot4
     WtFrac(ISpec,5) = WtFrac(ISpec,5) / WtTot5 
     WtFrac(ISpec,6) = WtFrac(ISpec,6) / WtTot6
  END DO
!!$     
1000 FORMAT(A36,1PE10.3)
1500 FORMAT(A36,I6)
1750 FORMAT(A36,A6)
2000 FORMAT(A3,A10,1PE9.3,A3)
3000 FORMAT(A3,A17,6(1PE9.3))
!!$     Read freez.dat for freezing parameters and methods

  DO
     READ(KFRE,*) Heading 
     IF(Heading == 'BEGIN') EXIT
  END DO

  READ(KFRE,Icemodo)
  READ(KFRE,Nucleation)
  READ(KFRE,Probability)
  READ(KFRE,IceStepStart)
  READ(KFRE,Divider)
  READ(KFRE,Satvappress)
  READ(KFRE,Heterogeneous) 
END SUBROUTINE ReadFile
