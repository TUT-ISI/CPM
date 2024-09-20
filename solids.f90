SUBROUTINE Solids
!!$C
!!$C     Read data for volatile compounds from 'solids.dat'
!!$C
  use headfile
  implicit none
  
  integer ISolid, ISld
!!$C
  ICheck = 0
  DO
     READ(KSLD,*) Heading
     IF(Heading == 'BEGIN') EXIT
  end DO
!!$C
  readloop: DO
     READ(KSLD,*) Heading
     IF(Heading == 'END' .AND. ICheck /= NASolids) THEN
        DO ISolid = 1,NASolids
           ISld  = NALiquids + ISolid
           IF(NameSpec(ISld) == 'none') THEN
              WRITE(*,*) 'Unknown solid compound in aerosol.dat'
              WRITE(*,*) 'No data found for ',FormSpec(ISld)
              WRITE(*,*) 'in solids.dat'
              STOP
           END IF
        END DO
     END IF
     IF(Heading == 'END') RETURN
     CALL ReadSld
!!$C
     IF(Heading == 'END') EXIT readloop
  end DO readloop
  
!!$C
END SUBROUTINE Solids


SUBROUTINE Liquids
!!$C
!!$C     Read data for volatile compounds from 'liquids.dat'
!!$C
  USE headfile
  implicit none
  integer ILiq
!!$C
  ICheck = 0
  DO
     READ(KLIQ,*) Heading
     IF(Heading == 'BEGIN') EXIT
  end DO
  
!!$C
  readloop: DO
     READ(KLIQ,*) Heading
     IF(Heading == 'END' .AND. ICheck /= NALiquids) THEN
        DO ILiq = 1,NALiquids
           IF(NameSpec(ILiq) == 'none') THEN
              WRITE(*,*) 'Unknown liquid compound in aerosol.dat'
              WRITE(*,*) 'No data found for ',FormSpec(ILiq)
              WRITE(*,*) 'in liquids.dat'
              STOP
           END IF
        END DO
     END IF
     IF(Heading == 'END') RETURN
     CALL ReadLiq
!!$C
     IF(Heading == 'END') EXIT readloop
  end DO readloop
!!$C
END SUBROUTINE Liquids

SUBROUTINE Gases
!!$C
!!$C     Read data for volatile compounds from 'gases.dat'
!!$C
  USE headfile
  implicit none
  integer LGases, LGas
!!$C
  ICheck = 0
  DO
     READ(KGAS,*) Heading
     IF(Heading == 'BEGIN') EXIT
  end DO
  
!!$C
  readloop: DO
     READ(KGAS,*) Heading
     IF(Heading == 'END'.AND. ICheck /= NAGases) THEN
        DO LGases = 1,NAGases
           LGas  = NALiquids + NASolids + LGases
           IF(NameSpec(LGas) == 'none') THEN
              WRITE(*,*) 'Unknown gas in aerosol.dat'
              WRITE(*,*) 'No data found for ',FormSpec(LGas)
              WRITE(*,*) 'in gases.dat'
              STOP
           END IF
        END DO
     END IF
     IF(Heading == 'END') RETURN
     CALL ReadGas
     
     IF(Heading == 'END') EXIT readloop
  end DO readloop
  
END SUBROUTINE Gases
