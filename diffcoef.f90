SUBROUTINE DIFFCOEF
!!$C     
!!$C     Subroutine that calculates the diffusion coefficients for the
!!$C     gas phase species and the transition correction factor for heat 
!!$C     and mass transfer for each size bin (Seinfeld and Pandis 1998
!!$C     
  USE headfile
  use precisi
  implicit none
  REAL(doubp) :: ek_AB, sigma_AB, WT_AB, T_star, DIFFC
  REAL(doubp) :: SIGMA_D, GASPEED, ZZ, GMeanFP, CpVD, DH, VelocTH 
  REAL(doubp) :: AlphaH, TMeanFP, CoeffKn, so_mass_sum
  integer :: JGas, I, ii, jj, ISpec
!!$C     
!!$C     TotP      = Total pressure (mbar)
!!$C     CKappa    = Thermal conductivity of air (J m-1 s-1 K-1)
!!$C     Rho       = Air density (G CM-3)
!!$C     
  TotP             = P + C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  CKappa           = 0.023807_doubp + 7.1128e-5_doubp * (T - 273.16_doubp)
  Rho              = TotP / (10000.0_doubp  * R / WtAir * T)
!!$c     3        + C(JLOOP,IH2O + NAERO) / AVG * WtGas(IH2O)
!!$C     

   HeatLice=(3.14566e6_doubp-(2.36164e3_doubp)*T + 4180.0_doubp*(79.7_doubp + &
       & 0.708_doubp*(T-273.15_doubp)-(2.5e-3_doubp)*(T-273.15_doubp)**2.)) &
       & /1000.0_doubp
  
  DO JGas=1,NAGases
!     JGas1          = JGas + NAero
     HeatL(JGas)    = HPV(JGas) + (T - 298.15_doubp)*CpV(JGas) - HPL(JGas)
!!$C     
!!$C     Diffusion coefficients calculated according to Reid et. al.
!!$C     The Properties of Gases and Liquids
!!$C     
!!$C     ek_AB    (K)              (eq 11-3.12)
!!$C     sigma_AB (Å)              (eq 11-3.13)
!!$C     sigma_D  dimensionless    (eq 11-3.7)
!!$C     
!!$C     WT_AB     = 2/[(1/M_A)+(1/M_B)], where M_A, M_B = molar weights of
!!$C     A and B (G MOL-1)
!!$C     T_star    = kT/ek_AB (dimensionless)
!!$C     COEFFDF   = THE MOLECULAR DIFFUSION COEFFICIENT (CM**2 S-1)
!!$C     Reid. et. al
!!$C     
!     ek_AB         = SQRT(E(JGas)*EAIR)
!     sigma_AB      = (SIGMAIR + SIGMA(JGas))/2.0_doubp
!     WT_AB         = 2.0_doubp/(1.0_doubp/WtGas(JGas) + 1.0_doubp/WTAIR)
!     T_star        = T/ek_AB
!!$C     DELTA_A       = 1.9D3*DIPM(JGas)**2/(VB(JGas)*TB(JGas))
!!$C     DELTA_AB      = SQRT(DELTA_A*DELTAIR)
!     SIGMA_D       = 1.06036_doubp/(T_star**0.1561_doubp) & 
!          &        + 0.193_doubp/exp(T_star*0.47635_doubp) &
!          &        + 1.03587_doubp/exp(T_star*1.52996_doubp) &
!          &        + 1.76474_doubp/exp(T_star*3.89411_doubp)
!!$C     
!!$C     To use Brokaw method, uncomment next line 
!!$C     (and find a value for DELTAIR or DELTA_AB)
!!$C     4                    + 0.19*DELTA_AB**2/T_star
!!$C     
!!$C     Wilke and Lee
!!$C     
!     D(JGas) = (3.03_doubp-(0.98_doubp/SQRT(WT_AB))) * T**(3./2.) / &
!          &    (TotP * SQRT(WT_AB)*SIGMA_AB**2*SIGMA_D)
    
!!$     C     WRITE(*,*) COEFFDF(JGas)
!!$C     
     IF(JGas == IH2Og)   D(JGas) = DIFFC(T,TOTP,1) *1.e2_doubp
     IF(JGas == IHNO3g)  D(JGas) = DIFFC(T,TOTP,2) *1.e2_doubp
     IF(JGas == INH3g)   D(JGas) = DIFFC(T,TOTP,4) *1.e2_doubp
     IF(JGas == IHCLg)   D(JGas) = DIFFC(T,TOTP,5) *1.e2_doubp
!!$C     
!!$C     Chapman and Enskog
!!$C     
!!$C     COEFFDF(JGas) = 2.66*T3(JLOOP)**(3./2.) /
!!$C     1           (TOTP * SQRT(WT_AB)*SIGMA_AB**2*SIGMA_D)
!!$C     
!!$C     GASPEED   = THE MEAN SPEED OF GAS MOLECULES (CM S-2)
!!$C     GMEANFP   = THE MEAN FREE PATH OF TRACE GASES (CM)
!!$C     
     GASPEED       = SQRT(8.0_doubp * R * T/ &
          & (Pi *WtGas(JGas)*1.e-3_doubp))*1.e2_doubp
     ZZ            = WtGas(JGas) / WTAIR
     GMeanFP       = 32.0_doubp * D(JGas) / &
          &    (3.0_doubp* Pi * (1.0_doubp+ ZZ) * GASPEED)
!!$C     
!!$C     GMeanFP   = THE THERMAL MEAN FREE PATH
!!$C     THVELOC   = THERMAL VELOCITY OF AN AIR MOLECULE 
!!$C     AlphaH    = THERMAL ACCOMODATION COEFFICIENT
!!$C     Alpha     = MASS ACCOMODATION COEFFICIENT
!!$C     
!!     CpVD       = 717.63_doubp
     CpVD       = 1004.0_doubp
     DH         = CKappa / (Rho * CpVD)
     VelocTH    = SQRT(8.0_doubp * R * 1.e7_doubp * T / Pi * WtAir)
     AlphaH     = 1.0_doubp
     IF(JGas == IH2Og) AlphaH = 0.96_doubp
!     IF(JGas == IH2Og) AlphaH = 0.7_doubp
     TMeanFP    = 3.0_doubp * DH / VelocTH
     Alpha(JGas) = 1.0_doubp
     DO I=max(1,Lfirst), last
!!$C     
!!$C     CoeffKn   = THE KNUDSEN NUMBER OF THE CONDENCING GAS
!!$C     BM        = THE TRANSITION CORRECTION FACTOR FOR MASS TRANSFER
!!$C     
        CALL RADIUS(I)
        CoeffKn     = GMeanFP / (Rp(I) * 1.e2_doubp)
       
!!$C     
!!$C     Fuchs and Sutugin (1971) (taken from Seinfeld and Pandis 1998)
!!$C     
        Bm(JGas,I) = (0.75_doubp * Alpha(JGas) * (1.0_doubp + CoeffKn) / &
             &     (CoeffKn**2 + CoeffKn + &
             &     0.283_doubp * CoeffKn * Alpha(JGas) + &
             &     0.75_doubp * Alpha(JGas)))
!!$C     
!!$C     Fuchs and Sutugin (1971), Pruppacher and Klett (1997)
!!$C     (taken from Jacobson 2000)
!!$C     
!!$C     BM(JGas,I) = 1. / (1 + 
!!$C     >              ((1.33 + 0.71 / CoeffKn) / 
!!$C     >              (1 + 1. / CoeffKn) + 
!!$C     >              4. * (1 - Alpha(JGas))/(3.*Alpha(JGas)))
!!$C     >              * CoeffKn)
!!$C     
!!$C     CoeffKn   = THE KNUDSEN NUMBER FOR ENERGY
!!$C     BT        = THE TRANSITION CORRECTION FACTOR FOR HEAT TRANSFER
!!$C     
        CoeffKn    = TMeanFP / (Rp(I) * 1.e2_doubp)

        BT(JGas,I) = (0.75_doubp * AlphaH * (1.0_doubp + CoeffKn) / &
             &           (CoeffKn**2 + CoeffKn + &
             &           0.283_doubp * CoeffKn * AlphaH + &
             &           0.75_doubp * AlphaH))
!!$C           PAUSE 'diffcoef'
!!$C     
!!$C     BM(JGas,I) = 1. / (1 + 
!!$C     >              ((1.33 + 0.71 / CoeffKn) / 
!!$C     >              (1 + 1. / CoeffKn) + 
!!$C     >              4. * (1 - AlphaH)/(3.*AlphaH))
!!$C     >              * CoeffKn)
!!$C   
!!$C            WRITE(*,*) BM(JGAS,I),BT(JGAS,I),'-'
     END DO


!
! CALCULATE SAME ALSO FOR "NEW" CLOUD DROPLETS IN Y2
!     
     IF (CLOUD_MOD ==1) THEN
        Alpha(JGas) = 1.0_doubp
        DO JJ=1, NY2
           IF(C_Y2(JJ)>1.0e-10) THEN
           I=NABINS+JJ
!!$C     
!!$C     CoeffKn   = THE KNUDSEN NUMBER OF THE CONDENCING GAS
!!$C     BM        = THE TRANSITION CORRECTION FACTOR FOR MASS TRANSFER
!!$C     
           CALL RADIUS_Y2(JJ)
           CoeffKn     = GMeanFP / (R_CL(JJ) * 1.e2_doubp)
         
!!$C     
!!$C     Fuchs and Sutugin (1971) (taken from Seinfeld and Pandis 1998)
!!$C     
           Bm(JGas,I) = (0.75_doubp * Alpha(JGas) * (1.0_doubp + CoeffKn) / &
                &     (CoeffKn**2 + CoeffKn + &
                &     0.283_doubp * CoeffKn * Alpha(JGas) + &
                &     0.75_doubp * Alpha(JGas)))
!!$C     
!!$C     Fuchs and Sutugin (1971), Pruppacher and Klett (1997)
!!$C     (taken from Jacobson 2000)
!!$C     
!!$C     BM(JGas,I) = 1. / (1 + 
!!$C     >              ((1.33 + 0.71 / CoeffKn) / 
!!$C     >              (1 + 1. / CoeffKn) + 
!!$C     >              4. * (1 - Alpha(JGas))/(3.*Alpha(JGas)))
!!$C     >              * CoeffKn)
!!$C     
!!$C     CoeffKn   = THE KNUDSEN NUMBER FOR ENERGY
!!$C     BT        = THE TRANSITION CORRECTION FACTOR FOR HEAT TRANSFER
!!$C     
           CoeffKn    = TMeanFP / (R_CL(JJ) * 1.e2_doubp)

           BT(JGas,I) = (0.75_doubp * AlphaH * (1.0_doubp + CoeffKn) / &
                &           (CoeffKn**2 + CoeffKn + &
                &           0.283_doubp * CoeffKn * AlphaH + &
                &           0.75_doubp * AlphaH))
!!$C           PAUSE 'diffcoef'
!!$C     
!!$C     BM(JGas,I) = 1. / (1 + 
!!$C     >              ((1.33 + 0.71 / CoeffKn) / 
!!$C     >              (1 + 1. / CoeffKn) + 
!!$C     >              4. * (1 - AlphaH)/(3.*AlphaH))
!!$C     >              * CoeffKn)
!!$C   
!!$C            WRITE(*,*) BM(JGAS,I),BT(JGAS,I),'-'
        ENDIF
        END DO
     ENDIF



     if (JGas == 1 .AND. ice_mod == 1) then
        Alpha(JGas) = 0.1_doubp
        DO jj=1, NY1 
           if(C_Y1(jj)>0.0) then
!              ii=(NFBins-NFABins+jj)
              so_mass_sum = 0.0_doubp
              ! Let's update the solid mass in bin I
              DO ISpec = 1, NALiquids + NASolids
                 so_mass_sum = so_mass_sum+ C_Y1(InCi(ISpec,JJ))*WtMol(ISpec)/C_Y1(JJ)
              END DO
              so_mass(jj) = so_mass_sum
              call RADIUSICE2(jj)
              
!!$C     Same than above, but for ice
              CoeffKn     = GMeanFP / (R_ice(jj) * 1.e2_doubp)
              
              BmIce(JGas,jj) = (0.75_doubp * Alpha(JGas) * (1.0_doubp + CoeffKn) / &
                   &         (CoeffKn**2 + CoeffKn + &
                   &         0.283_doubp * CoeffKn * Alpha(JGas) + &
                   &         0.75_doubp * Alpha(JGas)))
              
              CoeffKn    = TMeanFP / (R_ice(jj) * 1.e2_doubp)
              BTIce(JGas,jj) = (0.75_doubp * AlphaH * (1.0_doubp + CoeffKn) / &
                   &           (CoeffKn**2 + CoeffKn + &
                   &           0.283_doubp * CoeffKn * AlphaH + &
                   &           0.75_doubp * AlphaH))
           endif
        END DO
     end if
  END DO

END SUBROUTINE DIFFCOEF
