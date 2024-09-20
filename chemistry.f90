SUBROUTINE CHEMISTRY(Y,n,x_v,x)

!!C     
!!C *****************************************************************
!!C *   SUBROUTINE COAGULATION CALCULATES COAGULATION RATE:         *
!!C *                                                               *
!!C *                                                               *
!!C *****************************************************************
!!C     
  USE headfile
  use precisi
  implicit none
  INTEGER, intent(in) :: n
  REAL(doubp), intent(inout), DIMENSION(n) :: y, dydx
  REAL(doubp), intent(inout) :: x
  REAL(doubp), DIMENSION(NBins + NBins * NAerT + NGases) :: dCdx
  REAL(doubp) :: K_eq, dH_ads, dG_ads
  REAL(doubp) :: PMol
  REAL(doubp) :: RMol, CAQ, Rate, pw, pwe
  INTEGER :: II, IIOld, Ji, I_idx, J_J, Idx, Jgas_i, I_I, JGas, JGas1
  INTEGER :: IC, NR, Iice, ISpec, I, LWat, JJ, L3, IFam, J, Ibin


  IF(NReactions >= 1) THEN
       DO II = 1, NABins
!!$     
          LWat  = InC(IH2Ol,II)
!!$     
          DO NR = 1, NReactions
!!$     
!!$     Rate = RATE OF CHEMICAL REACTION
!!$     (for example d[Z]/dt = k[X][Y])
!!$     
!!$     k = A exp(B(298.15/T-1))
!!$     
             IF(NReacType(NR) == 1) THEN
                Rate=A(NR)*EXP(B(NR)*(298.15_doubp/T-1.0_doubp))
             ELSE IF(NReacType(NR) == 2) THEN
                TotC = P/1.3807e-19_doubp/T
                Rate = A(NR)*(300.0_doubp/T)**0.6_doubp*B(NR)*(300.0_doubp/T)**2.9_doubp &
                     & *TotC/(A(NR)*(300.0_doubp/T)**0.6_doubp+B(NR)*(300.0_doubp/T) &
                     & **2.9_doubp*TotC)*BROADF(NR)**(1/(1+LOG10(B(NR) &
                     & *(300.0_doubp/T)**2.9_doubp*TotC/ &
                     & (A(NR)*(300.0_doubp/T)**0.6_doubp))**2))
             ELSE
                write(*,*) 'Something is messed up in fex (funs.f90)'
                STOP
!!$ PAUSE 'JOTAIN VIKAA'
             END IF
!!$   
             DO I=1,3
                IF(IReacSpec(I,NR) /= 0) THEN
!!$   
!!$   IC  = Index in vector C
!!$   CAQ = Molality of species ISpec (mol/kg) !!!!!!!!
!!$   
                   IF(State(IReacSpec(I,NR)) == 'G' .AND. II == NABins) THEN
                      IC   = InC(IReacSpec(I,NR),II)
                      Rate = Rate*C(IC)**ReacMol(I,NR)
                   END IF
                   IF(State(IReacSpec(I,NR)) /= 'G') THEN
                      IC=InC(IReacSpec(I,NR),II)
                      CAQ=C(IC) * CWat / C(LWat)
                      Rate=Rate*CAQ**ReacMol(I,NR)
                   END IF
                END IF
             END DO
!!$   
             IF(State(IReacSpec(1,NR)) /= 'G') THEN
                Rate=Rate*C(LWat)/CWat
             END IF
!!$   
!!$   Calculate the chemical production and loss rate
!!$   for species ISpec
!!$   
             DO I=1,3
                IF(IReacSpec(I,NR) /= 0) THEN 
                   IF(State(IReacSpec(I,NR)) == 'G' .AND. II == NABins) THEN
                      ISpec=InC(IReacSpec(I,NR),II)
                      RMol=ReacMol(I,NR)
                      IF(NAMESpec(IReacSpec(I,NR)) /= 'OH') THEN
                         dCdx(ISpec)=dCdx(ISpec)-RMol*Rate
                      END IF
                   END IF
                   IF(State(IReacSpec(I,NR)) /= 'G') THEN
                      ISpec=InC(IReacSpec(I,NR),II)
                      RMol=ReacMol(I,NR)
                      dCdx(ISpec)=dCdx(ISpec)-RMol*Rate                     
                   END IF
                END IF
                IF(IPRODSpec(I,NR) /= 0) THEN
                   IF(State(IReacSpec(I,NR)) == 'G' .AND. II == NABins) THEN
                      ISpec=InC(IPRODSpec(I,NR),II)
                      PMol=ProdMol(I,NR)
                      dCdx(ISpec)=dCdx(ISpec)+PMol*Rate
                   END IF
                   IF(State(IReacSpec(I,NR)) /= 'G') THEN
                      ISpec=InC(IPRODSpec(I,NR),II)
                      PMol=ProdMol(I,NR)
                      dCdx(ISpec)=dCdx(ISpec)+PMol*Rate
                   END IF
                END IF
             END DO
!!$   
          END DO
!!$   CONTINUE DO 250 NR=1,NReactions
!!$   
          dCdx(InC(IH2Ol,II))=0.e0_doubp
!!$   
       END DO
!!$   CONTINUE DO 300 II = 1, NABins
!!$   
!!$   Move values of dCdx into dydx
!!$   
       DO I=1,NYLiquids
          dydx(I)=dCdx(IY2C(I))
       END DO
    
       DO I=1, NAGases
           dydx(NYAero+I)=dCdx(NAero+I)
       END DO
!!$   
       I = NYAero-NAFL*NABins
       DO IFam = 1,NAFL
          IF (NAFMbr(IFam) >= 2) THEN!!!!!!!!!!!ok?
             DO L3 = 1,NABins 
                I = I+1
                dydx(I) = 0.0_doubp
                DO II = 1, NFMbr(IFam)
                   IF(State(IFMbr(IFam,II)) /= 'G') THEN
                      dydx(I) = dydx(I) + dCdx(InC(IFMbr(IFam,II),L3))
                   END IF
                END DO
             END DO
          END IF
       END DO
    END IF











END SUBROUTINE CHEMISTRY
