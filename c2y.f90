SUBROUTINE C2Y(Y,NEq)

  USE headfile
  USE precisi
  IMPLICIT NONE
  INTEGER, INTENT(in) :: NEq
  INTEGER I,IFam,IMbr,IBin
  REAL(doubp), INTENT(inout), DIMENSION(NEq) :: Y
  


  DO I = 1, NYLiquids
!!$     Number concentrations
     Y(I) = C(IY2C(I))
!!$
  END DO

!!$     Family species
  I = NYLiquids
  
  DO IFam = 1,NAFl
     IF(NAFMbr(IFam).GT.0) THEN
        DO IBin = 1, NABins
           I = I + 1
           Y(I) = 0.
           DO IMbr = 1, NFMbr(IFam)
              IF(IFMbr(IFam,IMbr).NE.0) THEN
                 IF(State(IFMbr(IFam,IMbr)) /= 'G') THEN
                    Y(I) = Y(I) + C(InC(IFMbr(IFam,IMbr),IBin))
                 ENDIF
              ENDIF
           END DO
        END DO
     ENDIF
  END DO

!!$     Gases

  DO I = NYAero + 1, NYConcs
     Y(I) = C(IY2C(I))
  END DO
!  write(*,*) Y
!  pause
END SUBROUTINE C2Y
