SUBROUTINE Families
  USE headfile
  use precisi
  implicit none

  INTEGER :: K, NHL, I, J
  character(len=3), dimension(2) :: check_phase
!!$C
!!$C     Species involved in dissociation and hydration reactions are 
!!$C     grouped in "families". The families are read from 'chemset.dat'
!!$C


  NAFL = 0

  DO 
     READ(KCHM,*) Heading
     IF(Heading /= 'BEGIN') EXIT
  end DO
!!$     
!!$     NAFL   = number of active families
!!$     NFMBR  = number of family members
!!$     NAFMBR = number of active family members
!!$     FMBR   = name of the family member
!!$     IFMBR  = index of family member
!!$
!!$     Read families from chemset.dat
!!$     
200 READ(KCHM,1200) Heading, (FInp(J),J=1,4)
!!$
  IF(Heading /= 'END') THEN
!!$
     NAFL = NAFL + 1
     NHL  = 0
     I    = 0
!!$
     DO K=1,4
        DO J = 1,NAGases + NALiquids + NASolids
           IF(FInp(K) == FormSpec(J)) THEN
              I = I + 1
              FMBR(NAFL,I)=FInp(K)
              IFMBR(NAFL,I)=0
              IF(FInp(I) /= ' ') then 
                 NHL=NHL+1
                 IF (NHL <= 2) check_phase(NHL) = State(J)
              end IF
           END IF
        END DO
     END DO
!!$
     !  IF (NHL == 2 .and. 
     !  write(*,*) formspec(1), namespec(1)
     !  IF(NHL <= 2) THEN !!!!!!!!! LOL
     !IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
     IF(NHL <= 1) THEN 
        NAFl = NAFl - 1 
     ELSE IF(NHL == 2 .and. check_phase(1) /= check_phase(2)) THEN
        ! here the length of fam memb is 2 but they are in different phases (no family)
        NAFl = NAFl - 1 
     ELSE
        NFMBR(NAFl) = NHL
     END IF
!!$ 
     GOTO 200
  END IF
!!$
  DO I=1,NAGases + NALiquids + NASolids
     InFamily(I) = 0
  END DO
!!$
  DO I=1,NAFL
!!$
     NHL=0
!!$
!!$ What is the limit for family, 2 species here...
     DO J=1,NFMBR(I)
!!$
!!$     Identify family members
!!$
        DO K=1,NAGases + NALiquids + NASolids
           IF(FormSpec(K) == FMBR(I,J)) THEN
              IFMBR(I,J)  = K
              NHL=NHL+1
              InFamily(K) = 1
           END IF
        END DO
!!$
     END DO
!!$
     NAFMBR(I) = NHL
  END DO
!!$     
!!$
!!! IF (NAFMBR(2) < 2) NAFL = 0
!!$     
1200 FORMAT(A9,4(A13))
!!$
END SUBROUTINE Families
