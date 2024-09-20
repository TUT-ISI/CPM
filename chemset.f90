SUBROUTINE CHEMSET
!!$C
!!$C     Read reaction equations and reaction coefficients used in
!!$C     calculations from 'chemset.dat'
!!$C
  USE headfile
  use precisi
  implicit none
  real(doubp) :: PMol, RMol
  integer :: I, J, NR, NDEAD

  DO I=1,NReact
     NProdEq(I)=0
     NReacEq(I)=0
  END DO

  NR = 0

  readloop: DO
     READ(KCHM,*) Heading
     IF(Heading == 'BEGIN') EXIT readloop
  end DO readloop
  
  READ(KCHM,*)
  
200 NR    = NR + 1
  NDead = 0
  
  READ(KCHM,2000) JST, NReacType(NR),(ReacMol(J,NR), &
       & ReacSpec(J,NR),J=1,3)
  READ(KCHM,*)
  READ(KCHM,2100) (ProdMol(J,NR), ProdSpec(J,NR),J=1,3)
  READ(KCHM,*)
  READ(KCHM,3000) A(NR),B(NR),BroadF(NR)
  READ(KCHM,*)

  IF(JST == 'D' .OR. JST == 'END') NR=NR-1
  IF(JST == 'D') GOTO 200
  IF(JST == 'END') GOTO 500

  DO I=1,3
     IReacSpec(I,NR) = 0
     RMol = ReacMol(I,NR)
     PMol = ProdMol(I,NR)
!!$C    
!!$C     ReacMol = Stoichiometric numbers of the reacting species
!!$C     ProdMol = Stoichiometric numbers of the reaction products
!!$C
     IF(RMol < 0.5 .AND. ReacSpec(I,NR) /= ' ') ReacMol(I,NR)=1
     IF(PMol < 0.5 .AND. ProdSpec(I,NR) /= ' ') ProdMol(I,NR)=1
     DO J=1,NSpecT
!!$C
!!$C     Identify reactive species and species produced 
!!$C     in the chemical reaction NR
!!$C
        IF(ReacSpec(I,NR) == FormSpec(J)) THEN
!!$C     
!!$C     NReacEq = Number of chemical reactions where J is a reactive species
!!$C     
           NReacEq(J)=NReacEq(J)+1
           IReacSpec(I,NR)=J
        ENDIF
!!$C     
!!$C     NProdEq = Number of chemical reactions where J is produced
!!$C     
        IF(ProdSpec(I,NR) == FormSpec(J)) THEN
           NProdEq(J)=NProdEq(J)+1
           IProdSpec(I,NR)=J
        END IF
     END DO
!!$C
!!$C     See if there are any dead species in the reaction
!!$C 
     IF(IReacSpec(I,NR) == 0 .AND. ReacMol(I,NR) == 1) NDEAD = 1
     IF(IProdSpec(I,NR) == 0 .AND. ProdMol(I,NR) == 1) NDEAD = 1

  END DO
!!$C     
!!$C     Skip the reactions in which there are "dead" species
!!$C
  IF(NDead /= 0) THEN 
     WRITE(*,*) 'chemset.f skipped some chemical reactions with'
     WRITE(*,*) 'inactive species in chemical reactions.'
     NR = NR - 1
  ENDIF

GOTO 200

500 CLOSE(KCHM)

NReactions = NR

2000 FORMAT(A3,I1,F3.0,A14,3(1X,F2.0,A14))
2100 FORMAT(F7.0,A14,3(1X,F2.0,A14))
3000 FORMAT(3(4X,D8.3))

END SUBROUTINE CHEMSET
