SUBROUTINE ReadGas
!!$
!!$     SUBROUTINE FOR READING INPUT VALUES FOR GASES
!!$     
  USE headfile
  use precisi
  implicit none
  
  real(doubp) :: WeightMole, StickingCoefficient, HeatCapacityVapor
  real(doubp) :: StandardEntalphyV, StandardEntalphyL, CharacteristicLength
  real(doubp) :: CharacteristicEnergy, VolumeAtBoilingPt, BoilingTemperature
  real(doubp) :: DipoleMoment
  integer :: IGas, LGas
!!$
  NAMELIST /Compound/ &
  &     Name, &
  &     Formula, &
  &     WeightMole, &
  &     StickingCoefficient, &
  &     HeatCapacityVapor,  &
  &     StandardEntalphyV, &
  &     StandardEntalphyL, &
!!$     &     HenrysLawConstant,
  &     CharacteristicLength, &
  &     CharacteristicEnergy, &
  &     VolumeAtBoilingPt, &
  &     BoilingTemperature, &
  &     DipoleMoment, &
  &     DissolvesTo
!!$
  READ(KGAS,Compound)
!!$
  DO IGas = 1,NAGases
     LGas     = NALiquids + NASolids + IGas
!!$
!!$     Convert read values to a more compact form
!!$
     IF(Formula == FormSpec(LGas)) THEN
        ICheck         = ICheck + 1
        NameSpec(LGas) = Name
        WtMol(LGas)    = WeightMole
        WtGas(IGas)    = WeightMole
        Alpha(IGas)    = StickingCoefficient
        CpV(IGas)      = HeatCapacityVapor / WtGas(IGas)
        HpV(IGas)      = StandardEntalphyV / WtGas(IGas)
        HpL(IGas)      = StandardEntalphyL / WtGas(IGas)
!!$            H0(IGas)       = HenrysLawConstant
        Sigma(IGas)    = CharacteristicLength
        E(IGas)        = CharacteristicEnergy
        Vb(IGas)       = VolumeAtBoilingPt
        Tb(IGas)       = BoilingTemperature
        DipM(IGas)     = DipoleMoment
        State(LGas)    = 'G'
        Diss(IGas)     = DissolvesTo
     END IF
  END DO
!!$
END SUBROUTINE ReadGas
