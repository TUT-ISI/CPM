MODULE headfile
  USE precisi
  IMPLICIT NONE
  !  C-----------------------------------------------------------------------
  !  C     VARIABLE MODULE FOR THE CLOUD MODEL
  !  C-----------------------------------------------------------------------
!!$  C
!!$  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!!$  C
!!$  C     Constants:
!!$  C     R      = Gas constant (J K-1 mol-1)
!!$  C     Pi     = PiS
!!$  C     WtAir  = Molecular weight of air (g / mol)
!!$  C     g      = Gravitational coefficient (m/s2)
!!$  C     CWat   = Molality of water (mol / kg)
!!$  C     CpWat  = Heat capacity of water at constant pressure at 15
!!$  C     CpAir  = Heat capacity of dry air at constant pressure
!!$  C     h_planck = Planck's constant (m2 kg / s)
!!$  C
  REAL(doubp), PARAMETER :: R        = 8.314_doubp
  REAL(doubp), PARAMETER :: Pi       = 3.141592653589793_doubp
  REAL(doubp), PARAMETER :: FourPi   = 4.0_doubp*3.141592653589793_doubp
  REAL(doubp), PARAMETER :: WtAir    = 28.966_doubp
  REAL(doubp), PARAMETER :: g        = 9.806_doubp
  REAL(doubp), PARAMETER :: Avog     = 6.02252E+23_doubp
  REAL(doubp), PARAMETER :: Boltz    = R / Avog
  REAL(doubp), PARAMETER :: CWat     = 55.508681E+00_doubp
  REAL(doubp), PARAMETER :: CpWat    = 4.1858E+03_doubp
  REAL(doubp), PARAMETER :: CpIce    = 2.060E+03_doubp
  REAL(doubp), PARAMETER :: CpAir    = 1.00467E+03_doubp
  REAL(doubp), PARAMETER :: EAir     = 97.0_doubp
  REAL(doubp), PARAMETER :: Sigmair  = 3.617_doubp
  REAL(doubp), PARAMETER :: T_ntp    = 273.15_doubp
  REAL(doubp), PARAMETER :: h_planck = 6.626068E-34_doubp
!!$  C
!!$  C     NReact = Number of possible chemical reactions
!!$  C
  INTEGER, PARAMETER :: NReact   = 20
!!$  C
!!$  C     NGases    = Max number of gases
!!$  C     NLiquids  = Max number of liquid phase species
!!$  C     NSolids   = Max number of solids
!!$C     NAerT     = Max number of particle phase species
!!$C     NSpecT    = Max numser of species
!!$C     NModes    = Max number of modes
!!$C     NBins     = Max number of size bins
!!$C     NFamilies = Max number of families
!!$C     NDyna     = Number of variables describing dynamics of air parcel, and also total concentrations
!!$C
  INTEGER, PARAMETER :: NGases         = 7
  INTEGER, PARAMETER :: NLiquids       = 15
  INTEGER, PARAMETER :: NSolids        = 16
  INTEGER, PARAMETER :: NAerT          = NLiquids + NSolids
  INTEGER, PARAMETER :: NSpecT         = NGases + NLiquids + NSolids
  INTEGER, PARAMETER :: NModes         = 6
  INTEGER, PARAMETER :: NBins          = 400
  INTEGER, PARAMETER :: NFamilies      = 11
  INTEGER, PARAMETER :: NFamilyMembers = 10
  INTEGER, PARAMETER :: NDyna          = 7
  INTEGER, PARAMETER :: NFBins         = 50
  ! for Jacobson --->
  INTEGER, PARAMETER :: Mxaqeeq        = 50
  INTEGER, PARAMETER :: NMPartic       = 5
  INTEGER, PARAMETER :: IAERTY         = 50
  INTEGER, PARAMETER :: MXEQUIL        = 113
  INTEGER, PARAMETER :: MAXSALT        = 50
  INTEGER, PARAMETER :: MXPOLY         = 14
  INTEGER, PARAMETER :: NMDEAD         = 150
  INTEGER, PARAMETER :: IMISC          = 100
  INTEGER, PARAMETER :: IAERG          = IAERTY + NGases
  INTEGER, PARAMETER :: MXAQCONS       =  100
  INTEGER, PARAMETER :: MXMSULF        =   79
  INTEGER, PARAMETER :: MXTSULF        =  144
  INTEGER, PARAMETER :: KCOUNT         = 1431
  INTEGER, PARAMETER :: MXITOT         = 8*MAXSALT
  INTEGER, PARAMETER :: MXINDX         = MAXSALT
  REAL(doubp), PARAMETER  :: HINUMB    = 1.0E+38_doubp
  REAL(doubp), PARAMETER  :: SMAL3     = 1.0E-37_doubp
  REAL(doubp), PARAMETER  :: SMAL2     = 1.0E-37_doubp
  REAL(doubp), PARAMETER  :: EPS1      = 1.0_doubp - 1.2e-16_doubp
  ! --->
  !C    
  CHARACTER(len=80),SAVE :: Heading, Comment, Name
  !C
  CHARACTER(len=18),SAVE, DIMENSION(3,NReact) :: ReacSpec, ProdSpec
  CHARACTER(len=18),SAVE, DIMENSION(NSpecT) :: FormSpec, NameSpec, State
  CHARACTER(len=18),SAVE, DIMENSION(NGases) :: Diss
  CHARACTER(len=18),SAVE :: Formula, DissolvesTo, JST, Unit
  !C
  CHARACTER(len=18),SAVE, DIMENSION(NFamilies,NFamilyMembers) ::FMBR
  CHARACTER(len=18),SAVE, DIMENSION(NFamilyMembers) ::FInp

  !for Jacobson --->
  CHARACTER(len=14),save :: ISTT, MSTT, JST2
  CHARACTER(len=14),save, dimension(IMISC) :: SINP
  CHARACTER(len=14), save, dimension(0:MXEQUIL) :: NAMEQUIL
  CHARACTER(len=14),save,dimension(IAERG) :: MSTATE          
  CHARACTER(len=14),save,dimension(IMISC) :: XINP  
  CHARACTER(len=14), save, dimension(IMISC) :: RINP
  CHARACTER(len=14), save,dimension(0:NGASES) :: NAMEGAS
  CHARACTER(len=14),save,dimension(MXEQUIL) :: ISTATE
  CHARACTER(len=14),save,dimension(NMDEAD) :: NAMD

  CHARACTER(len=4), save, dimension(MXAQCONS) :: EQID
  CHARACTER(len=4), save :: KOLD, DINP, DINPLAST
  CHARACTER(len=4), save, dimension(8,MXAQCONS) :: EQIDR
  ! --->

!!$  CHARACTER(len=18), save :: &
!!$       &  State(NSpecT),            DissolvesTo, &
!!$       &  Heading,                  Comment, &
!!$       &  Name,                     Formula, &
!!$       &  JST,                      unit
!!$C
!!$C     File units or input/output files
!!$C     KCHM ~ 'chemset.dat'
!!$C     KGAS ~ 'gases.dat'
!!$C     KLIQ ~ 'liquids.dat'
!!$C     KSLD ~ 'solids.dat'
!!$C     KSPC ~ 'aerosol.dat'
!!$C     KENV ~ 'environ.dat'
!!$C     KTRM ~ 'thermoset.dat'
!!$C     KOUT ~ 'output.dat'
!!$C
!  COMMON /FILES/
  INTEGER, SAVE :: KCHM, KGAS, KLIQ, KSLD, KSPC, KENV, KTRM
  INTEGER, SAVE :: KFRE, KOUT, IOUT, JAOU, KACT, KEQQ
!!$C
!!$C     C   = Concentrations of gases, liquids and solids (mol/cm3)
!!$C     InC(I,J) = Index corresponding to where species I in size bin J
!!$C                can be found in vector C
!!$C     
!  COMMON /CONCENTRATIONS/  
  REAL(doubp),SAVE,DIMENSION(NBins + NBins * NAerT + NGases) :: C
  INTEGER,SAVE,DIMENSION(NAerT + NGases, NBins) :: InC
  INTEGER,SAVE,DIMENSION(NAerT + NGases, NFBins) :: InCi
  INTEGER,SAVE,DIMENSION(NGases) :: IGas2Liq
  REAL(doubp),SAVE,DIMENSION(NBins + NBins * NAerT + NGases+NFBins) :: YOld
  REAL(doubp),SAVE,DIMENSION(NBins) :: CFrac
  REAL(doubp),SAVE :: CSum, CSumOld
!!$C
!!$C     Number of active species
!!$C
!  COMMON /ACTIVE/
  INTEGER, SAVE :: NAGases, NALiquids, NASolids, NAerty
!!$C
!!$C     WtMol = Molar weight of species
!!$C     WtGas = Molar weight of gases
!!$C     Dens  = Density
!!$C
!  COMMON /PROPERTIES/
  REAL(doubp),SAVE,DIMENSION(NSpecT) :: WtMol, Dens
  REAL(doubp),SAVE,DIMENSION(NGases) :: WtGas
  REAL(doubp),SAVE,DIMENSION(NBins) :: wtH2O, wtHNO3, wtSO4
  INTEGER, SAVE, DIMENSION(NSpecT) :: IZLiq
!!$C
!!$C     ALPHA = sticking coefficient
!!$C     CpV   = heat capacity, vapor
!!$C     HpV   = Helmholz free energy, vapor
!!$C     HpL   = Helmholz free energy, liquid
!!$C     Sigma = Characteristic length
!!$C     E     = Characteristic energy
!!$C     Vb    = Volume at boiling point
!!$C     Tb    = Boiling temperature
!!$C     DipM  = Dipole moment
!!$C     Diss  = Formula of compound formed when gas dissolves in the liquid phase
!!$C
!  COMMON /GASPROPERTIES/
  REAL(doubp),SAVE,DIMENSION(NGases) :: Alpha, CpV, HpV, HpL, Sigma
  REAL(doubp),SAVE,DIMENSION(NGases) :: E, Vb, Tb, DipM, H0	! H0 is not needed
!!$C
!!$C     Environmental variables
!!$C     T      = Temperature (K)
!!$C     P      = Pressure (mbar) (during the time dependent calculations,
!!$C                             partial pressure of water is taken out
!!$C                             from P)
!!$C                                 
!!$C     TotP   = Pressure including gaseous water
!!$C     RH     = Relative humidity (%)
!!$C     z      = Vertical position of air parcel (m)
!!$C     w      = Vertical speed of air parcel (m/s)
!!$C     Rho    = Air density 
!!$C     CKappa = Thermal conductivity of air (J m-1 s-1 K-1)
!!$C
!  COMMON /ENVIRON/
  REAL(doubp),SAVE :: T, P, TotP, RH, z, w, Rho, CKappa
  REAL(doubp),SAVE :: DRIV_EQ_TIME,StepMax_eq, DTOUT_eq
  INTEGER, SAVE :: IDRIVER_EQ, Ithermo_mod
!!$C
!!$C     Size distribution data
!!$C     CBin         = Number concentration of one size bin (cm-3)
!!$C     Dp           = Particle diameter (m)
!!$C     Rp           = Particle radius (m)
!!$C     R_rdy        = Particle dry radius (m)
!!$C     GMD          = Geometric mean diameter (m)
!!$C     STD          = Geometric standard deviation 
!!$C     Ctot         = Total number concentration in all size bins (cm-3)
!!$C     CMode        = Number concentration of a mode (cm-3)
!!$C     WtFrac       = Weight fraction of a compound in a mode
!!$C     NBinsPerMode = Number of size bins per on mode,
!!$C                    when separable modes are used
!!$C     FirstBinMode = Index of the first bin in mode
!!$C     IModeBin     = Mode of bin
!!$C     NABins       = Number of size bins used in a simulation
!!$C     IMode        = Index of a mode
!!$C     NAModes      = Number of modes used in a simulation
!!$C 
!  COMMON /DISTR/     
  REAL(doubp),SAVE,DIMENSION(NBins * NModes) :: CBin, Dp, Rp, Rp_old,R_dry
  REAL(doubp),SAVE,DIMENSION(NModes) :: STD, GMD,CMode
  INTEGER,SAVE,DIMENSION(NModes) ::NBinsPerMode,FirstBinMode
  INTEGER,SAVE,DIMENSION(NBins*NModes) ::  IModeBin
  REAL(doubp),SAVE,DIMENSION(NSpecT,NModes) ::  WtFrac
  REAL(doubp),SAVE :: Ctot
  INTEGER,SAVE :: NABins, IMode, NAModes
  REAL(doubp),SAVE, dimension(NBins * NModes) :: C_OLD
!!$C
!!$C     Variables for irreversible chemical reactions:
!!$C
!!$C     ReacMol      = Stoichiometric numbers of the reacting species
!!$C     ProdMol      = Stoichiometric numbers of the reaction products
!!$C     A            = Reaction coefficient at reference temperature 
!!$C                    (see 'chemset.dat')
!!$C     B            = Parameter for temperature dependence of reaction
!!$C                    coefficient (see 'chemset.dat')
!!$C     BroadF       = Parameter for temperature dependence of reaction
!!$C                    coefficient (see 'chemset.dat')
!!$C     NReacType    = type of reaction (1 for liquid, 2 for gas)
!!$C     NReacEq      = Number of chemical reactions where species is a 
!!$C                    reactive species
!!$C     NProdEq      = Number of chemical reactions where species is 
!!$C                    produced
!!$C     IReacSpec    = Index of reacting species
!!$C     IProdSpec    = Index of reaction products
!!$C     NReactions   = Number of active chemical reactions
!!$C
!  COMMON /CHEM/
  REAL(doubp),SAVE,DIMENSION(3,NReact) :: ReacMol, ProdMol
  REAL(doubp),SAVE,DIMENSION(NReact) :: A, B, BroadF
  INTEGER,SAVE,DIMENSION(NSpecT) :: NReacEq, NProdEq, NReacType
  INTEGER,SAVE,DIMENSION(NSpecT,NReact) :: IReacSpec, IProdSpec
  INTEGER,SAVE :: NReactions
!!$C
!!$C     NAFL   = number of active families
!!$C     NFMBR  = number of family members
!!$C     NAFMBR = number of active family members
!!$C     IFMBR  = index of family member
!!$C
!  COMMON /NFML/
  REAL(doubp),SAVE,DIMENSION(NFamilies) :: TotF
  INTEGER,SAVE,DIMENSION(NFamilies,NFamilyMembers) :: IFMBR
  INTEGER,SAVE,DIMENSION(NSpecT) :: InFamily
  INTEGER,SAVE,DIMENSION(NFamilyMembers) :: NFMBR, NAFMBR
  INTEGER,SAVE :: NaFl
!!$C
!!$C     FMbr     = name of family members
!!$C     ReacSpec = Name of reacting species
!!$C     ProdSpec = Name of reaction products
!!$C     FormSpec = Formula of species
!!$C     NameSpec = Name of species
!!$C
!!$  COMMON /CHMCHAR/
!!$  &  FMbr(NFamilies,NFamilyMembers), 
!!$  &  FInp(NFamilyMembers),
!!$  &  ReacSpec(3,NReact),       ProdSpec(3,NReact),
!!$  &  FormSpec(NSpecT),         NameSpec(NSpecT)
!!$C
!!$C     FinHour  = Hours of simulation time
!!$C     FinMin   = Minutes of simulation time
!!$C     FinSec   = Seconds of simulation time
!!$C
!  COMMON /TIMESET/
  REAL(doubp),SAVE :: FinHour, FinMin, FinSec, TotSec, DtOut, DtOutIce,t_file
!!$C
!!$C     RelTol   = Relative tolerance for the ODE solver
!!$C     AbsTol   = Absolute tolerance for the ODE solver
!!$C     StepMax  = Maximum absolute step size allowed.
!!$C
!  COMMON /ODE/
  REAL(doubp),SAVE :: RelTol, AbsTol, StepMax
!!$C     
!!$C     IY2C     = Index in Y corresponding to index in C
!!$C     IC2Y     = Index in C corresponding to index in Y
!!$C     IMModes  = 1, if modes are separable
!!$C                0, if modes are inseparable
!!$C     If_Cond  = If condensation on for bin
!!$C
!  COMMON /CONTROL/
  INTEGER,SAVE,DIMENSION(NGases, NBins) :: IdxBin
  INTEGER,SAVE,DIMENSION(NBins + NBins * NAerT + NGases) :: IY2C, IC2Y
  INTEGER,SAVE,DIMENSION(NLiquids + NSolids) :: IfInY
  INTEGER,SAVE,DIMENSION(NBins) :: If_Cond,If_Thermo,IF_Change
  INTEGER,SAVE :: ICheck, IMModes, IfAdia
  INTEGER,SAVE,DIMENSION(NAert) :: NRLo
!!$C
!!$C     IH2Ol    = Index for liquid water
!!$C     IH2Og    = Index for gaseous water
!!$C     IHPlus   = Index for H+ ion
!!$C
!  COMMON /COMPOUNDS/
  INTEGER, SAVE :: IH2Ol, IH2Og, IHPlus, IOHMinus
  INTEGER, SAVE :: IHNO3g, INH3g, IHClg, INO3l, ISO2g
  INTEGER, SAVE :: ICO2g, IH2O2g, IO3g, IHNO3l
  INTEGER, SAVE :: INH3l, IHCll, ISO2l, ICO2l
  INTEGER, SAVE :: IH2O2l, IO3l, ISO42l, IORGL
  INTEGER, SAVE :: IORGS, IHNO3gi, INH3gi, IHClgi
  INTEGER, SAVE :: ISO2gi, ICO2gi, IH2O2gi, IO3gi
  INTEGER, SAVE :: IH2Ogi, IHSO4l, INH42SO4s,IH2SO4l,INH4NO3s
!!$C
!!$C       |-NYNumConc-|                             number concentrations in each bin
!!$C       |--------NYLiquids--------|              + concentrations of liquids in bins
!!$C       |-----------NYAero---------------|          + concentrations of grouped species
!!$C       |--------------NYConcs-------------|    + concentrations of gas phase compounds 
!!$C       |----------------NYIce-----------------|      + IWC (if ice mode on)
!!$C     Y(                                          )
!!$C
!!$C     NODEs = Number of ordinary differential equations to be solved
!!$C     NAero = Number of bins reserved for particle phase species
!!$C
!  COMMON /COUNTERS/
  INTEGER, SAVE :: NYNumConcs, NYLiquids, NYAero
  INTEGER, SAVE :: NYConcs, NODEs,NAero, NODEs_old, NYIce, Last, LFirst

!!$ C
  ! for Jacobson --->
  ! Following common blocks are needed if EQUISOLV II is used as thermodynamics
  ! But don't even try, it probably doesn't work.....
  REAL(doubp), SAVE :: RHLASTK, RELHUM3K
  
  !      COMMON /EQUISOLV/
  REAL(doubp), save, dimension(NBins) :: VOLUMET
  integer, save, dimension(3) :: NKNONSOL
  integer, save,dimension(2) :: NSOLEQUAT
  integer, save :: NONSOLID, NTOTION, NTOTI2, NTOTBINY, NEEQUAT, LCO2, LCO32

  !      COMMON /IMAQEEQ/
  integer, save, dimension(MXAQEEQ) :: NUMSPEQ, NACTIV, NUMWATS, NKTYPE
  integer, save, dimension(MXAQEEQ) :: NEQUSE, KCVAL,NUMTOPS
  integer, save, dimension(MXAQEEQ) :: NUMBOTS, LIQONLY, NKNEWS, KSALTEQ
  integer, save, dimension(MXAQEEQ) :: IEQTYP

  !      COMMON /DMAQEEQ/
  REAL(doubp), save, dimension(MXAQEEQ) :: AEQK, BEQK, CEQK, AEQT, BEQT, CEQT
  REAL(doubp), save, dimension(MXAQEEQ) :: DEQK, DEQT, STARTHI, STARTLO
  REAL(doubp), save, dimension(MXAQEEQ) :: GEXP, HEXP   

  !      COMMON /DMXAQEEQ2/
  REAL(doubp), save, dimension(MXAQEEQ,NMPARTIC) :: AKOEF3, SUMOL10
  
  ! COMMON /IMXAQEEQ2/
  integer, save, dimension(MXAQEEQ,NMPARTIC) :: KOEF,KOEF3
  integer, save, dimension(MXAQEEQ,NMPARTIC) :: IEQSPEC, IEQSPEC3,NACTSALT
  integer, save, dimension(MXAQEEQ,4,2) :: NACTSPEC2
  integer, save, dimension(MXAQEEQ,IAERTY) :: IFUSENK
  integer, save, dimension(MXAQEEQ,3) :: NKNEWR !!! alunperin NKEWR((MXAQEEQ,2)
  integer, save, dimension(MXAQEEQ,2) :: ISOLEQAT

!      COMMON /IMXEQUIL/
  integer, save, dimension(MXEQUIL) :: NEQUIMAP, JNION, IEQSOLID, NGASEQN
  integer, save, dimension(MXEQUIL) :: JUNIQCAT, JUNIQAN, NUMCATS, NUMANION
  integer, save, dimension(MXEQUIL) :: ISBACT, ISWACT, JUNIQAN2, NUMANI2
  integer, save, dimension(MXEQUIL) :: JNION2

  !      COMMON /IMAXSALT/
  integer, save, dimension(MAXSALT) :: IWATSALT, JBINSALT, ITPOLY, INDXBIN
  integer, save, dimension(MAXSALT) :: IPOLYTYP, IWATLIQ, NWWATLIQ, NUMCOEFS

  !      COMMON /DMAXSALT/
  REAL(doubp), save, dimension(MAXSALT) :: RELHLO, RTMXMOL, RHWAT1, RHWAT2

  !      COMMON /DIAERTY/
  REAL(doubp), save, dimension(IAERTY) :: CIONST, BMASS, RHCRYST, RHDTEMP
  REAL(doubp), save, dimension(IAERTY) :: RHD298, FRACMASS, KSALTDIS

  !      COMMON /IIAERG/
  integer, save, dimension(IAERG) :: IFASOLID, IFEQSOLV

  !      COMMON /DMSALT2/
  REAL(doubp), save, dimension(MAXSALT,MXPOLY,MXPOLY) :: WATT, BINT


  !      COMMON /DMXEQUIL/
  REAL(doubp), save, dimension(MXEQUIL) :: WCAT, WANI, G12CAT, G12ANI, Z12
  REAL(doubp), save, dimension(MXEQUIL) :: Z12CON, ZSUM22, BMOLARY

  !      COMMON /DMXPOLY/
  REAL(doubp), save, dimension(MXPOLY) :: BWT

  !      COMMON /DMISC/
  REAL(doubp), save, dimension(20) :: PINP
  REAL(doubp), save, dimension(10) :: EINP
  REAL(doubp), save, dimension(11) :: APHFAC
  REAL(doubp), save, dimension(4) :: APHXC
  REAL(doubp), save, dimension(11) :: TFACA 
  REAL(doubp), save :: APHI 

!      COMMON /KMISC/
  integer, save, dimension(400) :: IUVAL1
  integer, save, dimension(12) :: LDMONTH
  integer, save, dimension(20) :: KINP

  !      COMMON /DMXAQC/
  REAL(doubp), save, dimension(MXAQCONS,4) :: EQKCT
  
  !      COMMON /IMXAQC/
  integer, save, dimension(8,MXAQCONS) :: ICOEF
  

  !  COMMON /DMXMSULF/
  REAL(doubp),save,dimension(MXMSULF) :: ELECV  
  REAL(doubp),save,dimension(MXTSULF,MXMSULF*2) :: VALTABL

  !      COMMON /SOLVEQUI/
  REAL(doubp),save,dimension(1,NBins + NBins * NAerT + NGases) :: CBLK
  REAL(doubp),save,dimension(NBins + NBins * NAerT + NGases) :: CINIT, CORIG
  REAL(doubp),save,dimension(NBins+1,NBins + NBins * NAerT + NGases) :: CEST
  REAL(doubp),save :: TNOTT, TNOTT1, TNOTT2, TNOTT3, T3, T3K, T3VANHA, TDIFRHD
  REAL(doubp),save :: aphil
  REAL(doubp),save,dimension(Mxaqeeq) :: RRATE
  REAL(doubp),save,dimension(MXAQEEQ,NBins) :: SRATE
  REAL(doubp),save,dimension(NGases) :: BKELV
  REAL(doubp),save,dimension(MXITOT) :: g0tot,taulk 
  REAL(doubp),save,dimension(NBins,0:Ngases) :: akelfac
  REAL(doubp),save,dimension(NBins+1) :: GLOSS
  REAL(doubp),save,dimension(IAERTY) :: PPD
  REAL(doubp),save,dimension(iaerg,NBins) :: Bgammat
  integer, save :: Nsapoly, ITLOSULF
  integer, save, dimension(MXINDX, NBins) :: INDX
    
  !      COMMON /bromley/
  REAL(doubp),save,dimension(NBins) :: summl, summl2, smlzz, hipsu_B, osmcoe, watakti
  REAL(doubp),save,dimension(iaerg) :: Big_B
  REAL(doubp),save :: surftden, aln10, aln101, conswat, alogexp, told
  REAL(doubp),save,dimension(NBins) :: warning
  integer,save,dimension(NGases) :: neqgas
  REAL(doubp),save,dimension(Ngases,NBins) :: hdash
  integer, save :: NEG1
  ! ---> END OF EQUISOLV II


  !      COMMON /Icebin
! FOLLOWING VARIABLES ADDED BY SAMI
! NY1 = Number of ice bins in moving center distribution
! DY_1 = Bin limits in moving center bin (actually radius, fix the name) 
! C_Y1 = Concentrations in ice bins, similar to C
! Y1C = Index matrix for moving concentrations from Y1 to C
! B_coagul = Matrix including coagulation coefficients
! G_ind = target bins for colliding particles
!
  INTEGER, PARAMETER :: NY1         = 50
  INTEGER :: NFABINS2
  REAL(doubp),save,dimension(NY1) :: Y1
  REAL(doubp),save,dimension(NY1) :: R_ice_dry
  REAL(doubp),save,dimension(NY1+1) :: R_Y1
  REAL(doubp),SAVE,DIMENSION(NY1 + NY1 * NAerT) :: C_Y1
  INTEGER,SAVE,DIMENSION(NY1) :: Y1C
 


  !      COMMON /clouddrop
! NY2 = Number of cloud droplet bins in moving center distribution
! rY_2 = Bin limits in moving center grid for cloud droplets
! C_Y2 = Concentrations in cloud droplet bins, similar to C
! Y2C = Index matrix for moving concentrations from Y2 to C
! InCl(I,J) = Index corresponding to where species I in size bin J
!                can be found in vector C_Y2
! R_CL = Radius of cloud droplets
! R_CLdry = Dry radius of cloud droplets
  INTEGER, PARAMETER :: NY2         = 100
  REAL(doubp),save,dimension(NY2+1) :: r_Y2
  REAL(doubp),SAVE,DIMENSION(NY2 + NY2 * NAerT) :: C_Y2 
  INTEGER, save :: MIN_CL,IFCLOUD,NYCONCS2,NYCONCS_FEX,CLOUD_MOD
  INTEGER,SAVE,DIMENSION(NY2) :: Y2C
  INTEGER,SAVE,DIMENSION(NAerT + NGases, NY2) :: InCl
  REAL(doubp), SAVE, DIMENSION(NY2) :: R_CL,R_CLdry,HtermCL, C_clFrac
  INTEGER,SAVE :: N_CL2,MIN_COAG,MAX_COAG,IFCOAG
  REAL(doubp),SAVE :: rmin_cloud, c_cl_tot,N_drops,c_ice2
  INTEGER,SAVE :: IInuc,IF_NUC_FEX

  !      COMMON /Iceparticles
! NY3 = Number of graubel bins in moving center distribution
! IFCOAGL3 = matrix for coagulation including information if bin coagulates
! G_ind = target bins for colliding particles
  INTEGER, PARAMETER :: NY3         = 2
  INTEGER, PARAMETER :: NY4         = 2
  INTEGER,SAVE,DIMENSION(NBINS+NY1+NY2+NY3+NY4) :: IFCOAGL3
  REAL(doubp), SAVE,DIMENSION(NBins+NY1+NY2+NY3+NY4,NBins+NY1+NY2+NY3+NY4) :: G_ind
  REAL(doubp),save,dimension(NY3+1) :: r_Y3
  REAL(doubp),save,dimension(NY4+1) :: r_Y4
  REAL(doubp),SAVE,DIMENSION(NY3 + NY3 * NAerT) :: C_Y3
  REAL(doubp),SAVE,DIMENSION(NY4 + NY4 * NAerT) :: C_Y4
  INTEGER,SAVE,DIMENSION(NY3) :: Y3C
  INTEGER,SAVE,DIMENSION(NY4) :: Y4C
  INTEGER,SAVE,DIMENSION(NAerT + NGases, NY3) :: InC3
  INTEGER,SAVE,DIMENSION(NAerT + NGases, NY4) :: InC4
  REAL(doubp),SAVE :: c_Y3_tot,C_Y4_tot
  REAL(doubp), SAVE, DIMENSION(NY3) :: R_C3,R_C3dry,HtermC3, C_Y3Frac
  REAL(doubp), SAVE, DIMENSION(NY4) :: R_C4,R_C4dry,HtermC4, C_Y4Frac
  INTEGER,SAVE :: N_C3, N_C3_OLD
  INTEGER,SAVE :: N_C4, N_C4_OLD
  REAL(doubp), SAVE,DIMENSION(NBins+NY1+NY2+NY3+NY4,NBins+NY1+NY2+NY3+NY4) :: B_coagul
  !      COMMON /heterogeneous
! I_HM = 1, If Hallet Mossop for used, =0, if not
! I_HM_METHOD = Method to calculate Hallet Mossop
! I_HM_NEW = Ice bin to place HM formed particles 
! I_MODE_CONTACT = Number of lognormal mode where particles for contact freezing are
! I_MODE_DEPOS =  Number of lognormal mode where particles for deposition freezing are
! I_MODE_IMMERSION = Number of lognormal mode where particles for immersion freezing are
! IF_FREEZ_CONTACT = Switch to turn on contact freezing (0 or 1)
! IF_FREEZ_DEPOS = Switch to turn on deposition freezing (0 or 1)
! IF_FREEZ_IMMERSION =  Switch to turn on immersion freezing (0 or 1)
  REAL(doubp),save,dimension(NY1) :: dmdt    
  REAL(doubp),SAVE :: r_min_hm
  INTEGER,SAVE :: I_HM,I_HM_METHOD,I_HM_NEW
  INTEGER,SAVE :: I_MODE_CONTACT
  INTEGER,SAVE :: I_MODE_IMMERSION
  INTEGER,SAVE :: IF_FREEZ_CONTACT
  INTEGER,SAVE :: IF_FREEZ_DEPOS,IF_FREEZ_IMMERSION
  INTEGER,SAVE :: I_MODE_DEPOS
  INTEGER,SAVE :: IfFreez
  REAL(doubp),SAVE :: a_kiehl,B_kiehl,dtemp
  REAL(doubp),SAVE,dimension(NY2+NBins) ::  J_immersion
  REAL(doubp),SAVE ::f_act_dep,f_act_max

  !      COMMON /multibox
! IF_BOXES = Switch for column mode (0 or 1)
! N_BOX  = Number of consecutive runs in columm mode
! N_boxes = Maximum number of updrafts kept in memory
! N_step = Maximum number of vertical levels
  INTEGER, PARAMETER :: N_boxes         = 40
  INTEGER, PARAMETER :: N_step         = 500
  REAL(doubp),save,dimension(N_boxes,N_step,NY1 + NY1 * (NAerT-5) + NGases+10) :: CY1multi
  REAL(doubp),save,dimension(N_boxes,N_step,NY2 + NY2 * (NAerT-5) + NGases+10) :: CY2multi
  REAL(doubp),save,dimension(N_boxes,N_step,NY3 + NY3 * (NAerT-5) + NGases+10) :: CY3multi
  REAL(doubp),save,dimension(N_boxes,N_step,NY4 + NY4 * (NAerT-5) + NGases+10) :: CY4multi
  INTEGER,SAVE :: I_step,I_box,I_step_max
  INTEGER,SAVE :: IF_BOXES
  INTEGER,SAVE :: N_BOX
 !      COMMON /coalescence
  REAL(doubp),save :: eetta,GASPEED2
!!$C
!!$C     BT = Thermal correction factor
!!$C     Bm = Mass transfer correction factor
!!$C     D  = Diffusion coefficient
!!$C
!  COMMON /DIFFUSION/
  REAL(doubp), SAVE, DIMENSION(NGases,NBins+NY1+NY2) :: BT, Bm
  REAL(doubp), SAVE, DIMENSION(NBins) :: Hterm
  REAL(doubp), SAVE, DIMENSION(NGases) :: HeatL, D
!  COMMON /FREEZ/
  REAL(doubp), SAVE, DIMENSION(NBins+NY1+NY2) :: z_rate, T_drop, prob, F_HNO3 
  REAL(doubp), SAVE, DIMENSION(NBins) :: F_H2SO4, liq_mass, a_w, wat_act_old
  REAL(doubp), SAVE, DIMENSION(NFBins) :: T_ice, R_ice, so_mass, F_num, C_iceFrac, HTerm_ice
  REAL(doubp), SAVE, DIMENSION(NFBins) :: thetaHNO3, C_ice_HNO3
  REAL(doubp), SAVE, DIMENSION(NGases,NBins) :: BmIce, BTIce
  REAL(doubp), SAVE :: dr_c, HeatLice, upper_limit, ower_limit, e_sat_i
  REAL(doubp), SAVE :: freez_cons_lo, div_limit, C_ice_tot, so_mass_temp, time_help, tOut_help 
  REAL(doubp), SAVE :: NEq_help, y_min_limit, sum_prob, HeatL_liq2ice, Adia_fex
  REAL(doubp), SAVE :: theta_old, Adia_fex_old, ads_ice, C_min_limit
  REAL(doubp), SAVE :: DtOut_used
  INTEGER, SAVE :: rate_method, ice_mod, NFABins, ind_driver 
  INTEGER, SAVE :: e_sat_method, ice2y, istart_y, ind_sma, q_help, indx_ret 
  INTEGER, SAVE :: i_screen, ice_step_start, i_out, idiss_ads
  INTEGER, SAVE, DIMENSION(2) :: screen !CBini jäi pois
  REAL(doubp), SAVE, DIMENSION(NSpecT+NFBins+NBins) :: Y_init
  REAL(doubp), SAVE, DIMENSION(NFBins+NFBins*NLiquids+NFBins*NSolids) :: C_i
!  COMMON /KELVIN/				! CKelv is not needed 
  REAL(doubp), SAVE, DIMENSION(NGases, NBins+NY2) :: PSurf, CKelv
  REAL(doubp), SAVE, DIMENSION(NBins+NY2) :: Water_activity
  !C

END MODULE headfile


MODULE aikaa_mittaavat_parametrit
  IMPLICIT NONE
  ! Optional module for timing calculations. This module is used in files
  ! main.f90 and funs.f90.
  !
  ! BEGIN_TIME
  !		CPU time when calculations started
  ! TOTAL_INITIALIZE_TIME
  !		Time from the beginning to the end of first equilibrium calculation
  ! FIRST_EQ_TIME
  !		Time used in calculating the first equilibrium. This includes 
  !		intialization of the equilibrium model
  ! TOTAL_EQMODEL_TIME
  !		Time used to calculate equilibriums during time dependent calculations.
  !		The first equilibrium is not included
  ! EQ_CALCULATIONS
  !		Number of equilibrium model calls 
  REAL, SAVE :: BEGIN_TIME, TOTAL_INITIALIZE_TIME, FIRST_EQ_TIME, &
       TOTAL_EQMODEL_TIME
  INTEGER, SAVE :: EQ_CALCULATIONS
  ! Additional variables, needed for calculating time difference
  REAL :: TIME_IN, TIME_OUT
  !
END MODULE aikaa_mittaavat_parametrit


MODULE jacobi
  IMPLICIT NONE

  INTEGER :: IFJACOP,NJACOP

END MODULE jacobi
