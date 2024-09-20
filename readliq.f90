SUBROUTINE ReadLiq
!!$C
!!$C     SUBROUTINE FOR READING INPUT VALUES FOR GASES
!!$C     
  USE headfile
  use precisi
  implicit none

  real(doubp) :: WeightMole, Density
  integer :: IZLiquid, ILiq
!!$C
  NAMELIST /Compound/ &
  &     Name, &
  &     Formula, &
  &     WeightMole, &
  &     Density, &
  &     IZLiquid
!!$C
  READ(KLIQ,Compound)
!!$C
  DO ILiq    = 1,NALiquids
     IF(Formula == FormSpec(ILiq)) THEN
        ICheck         = ICheck + 1
        NameSpec(ILiq) = Name
        WtMol(ILiq)    = WeightMole
        Dens(ILiq)     = Density
        IZLiq(ILiq)    = IZLiquid
        State(ILiq)    = 'L'
     END IF
  END DO
!!$C
END SUBROUTINE ReadLiq
