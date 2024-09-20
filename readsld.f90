SUBROUTINE ReadSld
!!$C
!!$C     Subroutine for reading input values for solids
!!$C     
  USE headfile
  use precisi
  implicit none
  
  real(doubp) :: WeightMole, Density
  integer :: ISolid, ISld
!!$C
  NAMELIST /Compound/ &
  &     Name, &
  &     Formula, &
  &     WeightMole, &
  &     Density
!!$C
  READ(KSLD,Compound)
!!$C
  DO  ISolid = 1,NASolids
     ISld       = NALiquids + ISolid
     IF(Formula == FormSpec(ISld)) THEN
        ICheck         = ICheck + 1
        NameSpec(ISld) = Name
        WtMol(ISld)    = WeightMole
        Dens(ISld)     = Density
        State(ISld)    = 'S'
     END IF
  END DO
!!$C
END SUBROUTINE ReadSld
