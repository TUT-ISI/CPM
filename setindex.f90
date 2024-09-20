SUBROUTINE SetIndex
!!$
  USE headfile
  use precisi
  implicit none
  
  integer :: ISpec, IndexY, IndexC, IBin, IGas!, NYNumConc 
  integer :: IFam, I, L3, IfCond, JGas, JLiq
!!$
  NAero = NABins + NALiquids * NABins + NASolids  * NABins
  IGas2Liq(1:NGases) = 0
!!$
!!$     Make vectors for moving values between y and C
!!$
  IBin = 1
!!$
  IndexY         = 0
  IndexC         = 0
!!$
!!$     Number concentrations
!!$
!!$      DO 50 I        = 1, NABins
!!$         IY2C(I)     = I
!!$         I!!$2Y(I)     = I
!!$         IndexY      = IndexY + 1
!!$         IndexC      = IndexC + 1
!!$ 50   CONTINUE
!!$
!  NYNumConc = IndexY
!!$     
!!$     Liquid phase species
!!$
  DO ISpec      = 1, NALiquids + NASolids
!!$     
     IfInY(ISpec)   = 0
     IfCond         = 0
!!$
     DO  IGas = 1, NAGases
        IF(Diss(IGas) == FormSpec(ISpec)) THEN
           IfCond   = 1
        END IF
     END DO
     DO IBin    = 1, NABins   
        IndexC      = InC(ISpec,IBin)
!!$
!!$     Exclude species grouped in families for now
!!$     
 !       write(*,*) InFamily(ISpec), FormSpec(ISpec), IfCond
        IF(InFamily(ISpec) /= 1 &
!!$     >           .AND.FormSpec(ISpec).EQ.'H2O(l)'
             &           .AND.FormSpec(ISpec) /= 'H+' &
             &           .AND.FormSpec(ISpec) /= 'OH-' &
             &           .AND.State(ISpec) /= 'S' &
             &           .AND.IfCond == 1) THEN
           IndexY       = IndexY + 1
           IY2C(IndexY) = IndexC
           IC2Y(IndexC) = IndexY
           IfInY(ISpec) = 1
        END IF
     END DO
  END DO
!!$
  NYLiquids = IndexY
!!$
!!$     Species grouped in families
!!$
  DO IFam              = 1, NAFl
     IF(NAFMbr(IFam) > 0) THEN
        DO IBin     = 1, NABins
           IndexY       = IndexY + 1
!!$               IndexC       = IndexC + 1
!!$               IY2C(IndexY) = IndexC
!!$               IC2Y(IndexC) = IndexY
        END DO
     END IF
  END DO
!!$
  NYAero = IndexY
!!$     
!!$     Gas phase species
!!$
  DO IGas              = 1, NAGases
     IndexY                = IndexY + 1
     IndexC                = IndexC + 1
     IY2C(IndexY)          = IndexC
     IC2Y(IndexC)          = IndexY
     DO I              = 1, NALiquids
        IF(FormSpec(I) == Diss(IGas)) THEN
           DO L3          = 1, NABins
              IdxBin(IGas,L3) = InC(I,L3)
           END DO
        END IF
     END DO
  END DO
!!$
  NYConcs = IndexY
  NYIce   = NYConcs
  NODEs   = NYConcs + NDyna
!!$
!!$     Identify some key species
!!$
  DO I                 = 1, NALiquids + NASolids + NAGases
     IF(FormSpec(I) == 'H2O(l)')  IH2Ol  = I
     IF(FormSpec(I) == 'H2O')     IH2Og  = I -NALiquids-NASolids
     IF(FormSpec(I) == 'H2O')     IH2Ogi = I !-NALiquids-NASolids these are for InC indexes
     IF(FormSpec(I) == 'H+')      IHPlus = I
     IF(FormSpec(I) == 'OH-')     IOHMinus = I
     IF(FormSpec(I) == 'HNO3')    IHNO3g = I -NALiquids-NASolids
     IF(FormSpec(I) == 'HNO3')    IHNO3gi= I !-NALiquids-NASolids
     IF(FormSpec(I) == 'NH3')     INH3g  = I -NALiquids-NASolids
     IF(FormSpec(I) == 'NH3')     INH3gi = I !-NALiquids-NASolids
     IF(FormSpec(I) == 'HCl')     IHClg  = I -NALiquids-NASolids
     IF(FormSpec(I) == 'HCl')     IHClgi = I !-NALiquids-NASolids
     IF(FormSpec(I) == 'SO2')     ISO2g  = I -NALiquids-NASolids
     IF(FormSpec(I) == 'SO2')     ISO2gi = I !-NALiquids-NASolids
     IF(FormSpec(I) == 'CO2')     ICO2g  = I -NALiquids-NASolids
     IF(FormSpec(I) == 'CO2')     ICO2gi = I !-NALiquids-NASolids
     IF(FormSpec(I) == 'H2O2')    IH2O2g = I -NALiquids-NASolids
     IF(FormSpec(I) == 'H2O2')    IH2O2gi= I !-NALiquids-NASolids
     IF(FormSpec(I) == 'O3')      IO3g   = I -NALiquids-NASolids
     IF(FormSpec(I) == 'O3')      IO3gi  = I !-NALiquids-NASolids
     IF(FormSpec(I) == 'NO3-')    INO3l = I
     IF(FormSpec(I) == 'HNO3(l)') IHNO3l = I
     IF(FormSpec(I) == 'NH4+')    INH3l  = I
     IF(FormSpec(I) == 'Cl-')     IHCll  = I
     IF(FormSpec(I) == 'SO2(l)')  ISO2l  = I
     IF(FormSpec(I) == 'CO2(l)')  ICO2l  = I
     IF(FormSpec(I) == 'H2O2(l)') IH2O2l = I
     IF(FormSpec(I) == 'O3(l)')   IO3l   = I
     IF(FormSpec(I) == 'SO42-')   ISO42l = I
     IF(FormSpec(I) == 'HSO4-')   IHSO4l = I
     IF(FormSpec(I) == 'C10H16O3(l)') IORGL = I
     IF(FormSpec(I) == 'C10H16O3(s)') IORGS = I
     IF(FormSpec(I) == '(NH4)2SO4(s)') INH42SO4s = I
     IF(FormSpec(I) == 'NH4NO3(s)') INH4NO3s = I
     IF(FormSpec(I) == 'H2SO4(l)') IH2SO4l = I
  END DO
!!$         
  DO JGas = 1, NAGases
     IF(JGas == IH2Og)  ISpec = IH2Ol
     IF(JGas == IHNO3g) ISpec = IHNO3l
     IF(JGas == INH3g)  ISpec = INH3l
     IF(JGas == IHClg)  ISpec = IHCll
     IF(JGas == ISO2g)  ISpec = ISO2l
     IF(JGas == ICO2g)  ISpec = ICO2l
     IF(JGas == IH2O2g) ISpec = IH2O2l
     IF(JGas == IO3g)   ISpec = IO3l
!!$
!     write(*,*) I, ISpec
     !! there are no bin-structure in gas phase
!     IdxBin(I,1:NABins-1) = 0 !InC(ISpec,NABins)
     IdxBin(JGas,1:NABins) = InC(ISpec,1:NABins)
!     DO IBin = 1, NABins
!        write(*,*) InC(ISpec,IBin)
!        IdxBin(I,IBin) = InC(ISpec,IBin)
!     END DO
     ! This is for finding the index of dissolving gas in the liquid phase
     DO JLiq = 1, NALiquids
        IF(DISS(Jgas) == FormSpec(JLiq)) IGas2Liq(JGas) = JLiq
     END DO
     
     !write(*,*) formspec(1:NALIquids), DISS(2)
!!$
  END DO
!!$
END SUBROUTINE SetIndex
