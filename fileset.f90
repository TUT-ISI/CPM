SUBROUTINE FileSet(i)

  USE headfile
  implicit none
  
  CHARACTER(len=3) :: NUMB
  INTEGER :: i
!!$C     INPUTS

  KGAS=5
  KLIQ=7
  KSLD=8
  KSPC=9
  KCHM=10
  KFRE=11
  KACT = 12
  KENV=13
  KEQQ =17
  KTRM=18
!!$C     OUTPUTS

  KOUT = 14  
  IOUT = 15
  JAOU = 16
!!$  KGAS=5
!!$  KLIQ=7
!!$  KSLD=8
!!$  KSPC=9
!!$  KCHM=10
!!$  KFRE=11
!!$  KENV=13
!!$aa
!!!C     OUTPUTS
!!$
!!$  KOUT =14  
!!$  IOUT =15
!!$  C      DATA KACT / 15 /  
!!$C
  IF (I == 0) THEN
     OPEN(KENV,FILE='input/environ.dat')
  ELSEIF (I == 1) THEN
     OPEN(KGAS,FILE='input/gases.dat')
     OPEN(KLIQ,FILE='input/liquids.dat')
     OPEN(KSLD,FILE='input/solids.dat')
     OPEN(KSPC,FILE='input/aerosol.dat')
     OPEN(KCHM,FILE='input/chemset.dat')
     OPEN(KENV,FILE='input/environ.dat')
     OPEN(KFRE,FILE='input/freez.dat')
     OPEN(KTRM,FILE='input/thermoset.dat')
  ELSE IF (I == 2) THEN
     IF(IF_BOXES==1) i_out=I_BOX

     NUMB = char(48+(I_out/100))//char(48+((I_out-100*(i_out/100))/10))  &
          &          //char(48+(I_out-10*(I_out/10)))
!     OPEN(KOUT,FILE='output/outputi'//NUMB//'.out')
!     OPEN(IOUT,FILE='output/iceouti'//NUMB//'.out')
     OPEN(KOUT,FILE='output/output.out')
     OPEN(IOUT,FILE='output/output_ice.out')
     OPEN(JAOU,FILE = 'xx1', FORM = 'FORMATTED', STATUS = 'UNKNOWN') 
  END IF

END SUBROUTINE FileSet
