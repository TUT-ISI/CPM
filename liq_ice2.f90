SUBROUTINE liq_ice(Y,n)
!!$    
  USE headfile
  use precisi
  implicit none
  real(doubp) :: C_ice_help, mass_tot, mass_tot_old, Tot_Wat_ice
  real(doubp) :: T_C, Tot_Wat, Cptot,so_mass_sum,frac
  integer, intent(inout) :: n
  real(doubp), dimension(n), intent(inout) :: Y
  real(doubp) :: f_act,S0,a_f
  integer I, J,ISpec,j_in_Y2,j_step,JJ
!!$    
!!$  *****************************************************************
!!$  *   SUBROUTINE LIQ_ICE SIMULATES PARTICLE FREEZING AND MAKES    *
!!$  *   THE ICE PARTICLE DISTRIBUTION FOR FURTHER ANALYSIS        *
!!$  *****************************************************************
!!$
!!$
 
  ! ICE MASS BEFORE ADDITION OF NEW PARTICLES
  mass_tot_old=0.0_doubp

  do I=1,NY1
     mass_tot_old=mass_tot_old+so_mass(I)*C_Y1(I)
  end do


  
  
  ! Amount of mass in individual ice bins


  C_ice_help = 0.0_doubp
  C_ice_help = C_ice_tot
 
  
  call ice_bin_divider(Y,n)
  


  IF(IF_FREEZ_DEPOS ==1) CALL ice_deposition(Y,n)
  IF (ice_mod == 1) THEN
     DO  I = 1,NY1
        so_mass_sum = 0.0_doubp
        IF(C_Y1(I) > 0.0) THEN
           DO ISpec = 1, NALiquids + NASolids
              so_mass_sum = so_mass_sum + C_Y1(InCi(ISpec,I))*WtMol(ISpec)/C_Y1(I)
           END DO
        ENDIF
        so_mass(I) = so_mass_sum
        CALL RADIUSICE2(I)
     END DO
  END IF
     
  
  NYIce = NYConcs + NFABins

  mass_tot=0.0_doubp
  do I=1,NY1
     mass_tot=mass_tot+so_mass(I)*C_Y1(I)
  end do
!
!!$
!
!  C_ice_tot=0.0_doubp
  C_ice_tot=sum(C_Y1(1:NY1))
  Ctot = sum(C(1:NaBins)) ! Ctot - C_ice_tot + C_ice_help
  Y(NODEs-1) =  Ctot  !Y(NODEs-1) - C_ice_tot + C_ice_help
  Y(NODEs-5) = C_ice_tot
!$
!$    Latent heat for phase transform liquid->solid (J/cm3)
!$
!$
  if (NFABins2 > 0) then
     T_C=T-T_ntp
     HeatL_liq2ice=(4.1752_doubp*T_C+0.5_doubp*(-11.319e-3_doubp)*T_C**2+(1./3.)* &
          &     (-9.7215e-5_doubp)*T_C**3+(1./4.)*(1.8315e-5_doubp)*T_C**4+(1./5.)* &
          &     (1.1354e-6_doubp)*T_C**5-2.1046_doubp*T_C-0.0037_doubp*T_C**2+&
          &     334.0_doubp)*(mass_tot-mass_tot_old)
!!$
     Tot_wat=sum(C(Inc(IH2Ol,1:NABins)))
     IF(cloud_mod ==1)  Tot_wat=sum(C(Inc(IH2Ol,1:NABins)))+sum(C_Y2(Incl(IH2Ol,1:NY2)))
     Tot_wat_ice=sum(C_Y1(Inci(IH2Ol,1:NFBins)))
!!$C
     Cptot = (C(InC(IH2Ogi,NABins)) * WtGas(IH2Og) * 1.95_doubp + &! water gas
          &   Rho * CpAir * 1.0e-3_doubp + &! air
          &   Tot_Wat*WtMol(IH2Ol)*CpWat*1.0e-3_doubp + &! liquid water
          &   Tot_Wat_ice*WtMol(IH2Ol)*CpIce*1.0e-3_doubp)!/ ! ice
!!$
     T=T+HeatL_liq2ice/Cptot
!     write(*,*) HeatL_liq2ice/Cptot,mass_tot-mass_tot_old
!     pause
     Y(NODEs)= Y(NODEs)+HeatL_liq2ice/Cptot
!!$C
  end if
!!$     
  DO J = 1, NABins
     CALL RADIUS(J)
  END DO

  DO J = 1, NY2
     CALL RADIUS_Y2(J)
  END DO

  Ctot=sum(C(1:NaBins))
  c_cl_tot=sum(C_Y2(1:NY2))
  c_ice_tot=sum(C_Y1(1:NY1))
  Y(NODES-5) = C_ice_tot
  Y(NODES-1)= Ctot
  Y(NODES-6) = C_cl_tot
  


END SUBROUTINE liq_ice

!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine ice_bin_divider(Y,n)
!!$
!!$     This program uses the freez_prob subroutine and divides the frozen liquid
!!$     phase into solid phase distribution.
  USE headfile
  use precisi
  implicit none
!!$
  integer, intent(inout) :: n
  real(doubp), dimension(n+NFBins), intent(inout) :: Y
  integer j, i
!!$
  so_mass_temp=0.0_doubp
  call freez_prob(1)
!!$
  if (sum(prob(1:NABins)) == 0.0_doubp) return !nothing freezes
!!$
  freezing:  do  j=1, NABins
     q_help = 0
     i=NABins-(j-1)
!!$     Let's check the probability for freezing in size class i
!!$     If probability is o, we skip the ice main program; elseif probability is 1,
!!$     we change it to a little bit under 1
     if (prob(i) == 0.0_doubp) cycle freezing
     if (prob(i) == 1.0_doubp) prob(i)=(1-ower_limit)
!!$         if(prob(i) == 1) prob(i)=0.8
!!$     At first, we must check if number consentration in liquid phase for size class i
!!$     is too low. This option makes the model more stable.
     if (C(i) < C_min_limit) cycle freezing
     call calc_variables(i,1)
!!$     What about the frozen particle consentration, it should't be too small
     if (C(i)*prob(i) < dr_c) then 
        cycle freezing
     else
        n=NYConcs + NDyna+NY1+NY2+NY3+NY4
        call num_f_bins
        call Y_to_Y1(i,Y,n)

     end if
  end do freezing



  IF(cloud_mod ==1) THEN


     freezingY2:  do  j=1,NY2
       
        i=NaBins+j
!!$     Let's check the probability for freezing in size class j in Y2
!!$     If probability is o, we skip the ice main program; elseif probability is 1,
!!$     we change it to a little bit under 1
        if (prob(i) == 0.0_doubp) cycle freezingY2
        if (prob(i) == 1.0_doubp) prob(i)=(1-ower_limit)
!!$         if(prob(i) == 1) prob(i)=0.8
!!$     At first, we must check if number consentration in cloud droplets for size class j
!!$     is too low. This option makes the model more stable.
        if (C_Y2(j) < C_min_limit) cycle freezingY2
!!$     What about the frozen particle consentration, it should't be too small
!!     write(*,*) 'here',prob(i)
        if (C_Y2(j)*prob(i) < dr_c) then 
           cycle freezingY2
        else
           n=NYConcs + NDyna+NY1+NY2+NY3+NY4
           call Y2_to_Y1(j,Y,n)
           
        end if
     end do freezingY2
  ENDIF



end subroutine ice_bin_divider

!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine calc_variables(i,i_mass)
!
!  I'M NOT SURE IF THIS IS NEEDED OR NOT
!
  USE headfile
  use precisi
  implicit none

  integer, intent(in) :: i_mass, i
  integer ::  ISpec
!!$     Meaning of this program is to calculate some needed variables for the subroutine 
!!$     ic_bin_divider. Some of these cariables are for help and some go straight to 
!!$     solid phase


  so_mass_temp=0.0_doubp
  solidmass: do ISpec=1,NALiquids + NASolids
     if (ISpec == IHNO3l) then
        so_mass_temp = so_mass_temp+C(InC(ISpec,i))*WtMol(ISpec)

!!$           
!!$        if (method_HNO3 == 2 .or. method_HNO3 == 4) then
!!$           cycle solidmass
!!$        else
!!$           so_mass_temp = so_mass_temp+C(InC(ISpec,i))*WtMol(ISpec)
!!$        end if
     else
        so_mass_temp = so_mass_temp+C(InC(ISpec,i))*WtMol(ISpec)
     end if
  end do solidmass
!!$
  so_mass_temp = so_mass_temp/C(i)
!!$
  if (NFABins == 0) then
     ind_sma = (NFBins)
  else
     ind_sma = (NFBins+1-NFABins)
  end if
!!$
  ice2y = NYAero/NABins !how many substances are in vector Y
  istart_y = NYConcs  !starting point for solid phase H2O (in vector Y)
!!$
!  if (i_mass == 2) then
!     call test_mass(so_mass_temp)
!  end if
end subroutine calc_variables

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine  Y_to_Y1(i,Y,n)

!
!  This is a new subroutine to calculate the greezing and to move 
!  particles/droplets from liquid phase to ice phase
!
!
  USE headfile
  use precisi
  implicit none
  integer, intent(in) :: i
  integer :: ISpec,j_inx,j_step,j
  integer, intent(inout) :: n
  real(doubp), dimension(n+1), intent(inout) :: Y

  j_inx=NY1/2
  
!
! First look where freezing drop is placed
!
  finder: do j_step=1,NY1/2+1   
     if (R_Y1(j_inx) < Rp(i) .and. R_Y1(j_inx-1) > Rp(i) ) then
        j_inx = j_inx-1
        exit finder
      endif

     if (R_Y1(j_inx) > Rp(i)) then
        j_inx=max(2,j_inx-1)
     else
        j_inx=min(j_inx+1,NY1-1)
     endif
  end do finder

  j=max(1,j_inx)

  F_num(j) = 1

!
! Move numbers from liquid distribution to solir  
!
  C_Y1(j) = C_Y1(j) + C(i)*prob(i)
!
  C(i) = C(i)*(1.0_doubp-prob(i))
  T_ice(j) = T
!
! and then water and other compounds
!
  do  ISpec=1,NALiquids + NASolids
     C_Y1(InCi(ISpec,j)) = C_Y1(InCi(ISpec,j)) + C(InC(ISpec,i)) * prob(i)
     C(InC(ISpec,i)) = C(InC(ISpec,i)) * (1.0_doubp-prob(i))
  end do
!  write(*,*) j,prob(i),C_Y1(j),C(i)*prob(i)
!
! Update also Y
!
  Y(istart_y+j)=C_Y1(InCi(1,j))
  Y(i)=C(InC(1,i))


  CALL RADIUSICE2(j)


end subroutine Y_to_Y1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine  Y2_to_Y1(i,Y,n)

!
!  This is a new subroutine to calculate the greezing and to move 
!  particles/droplets from liquid phase to ice phase
!
!
  USE headfile
  use precisi
  implicit none
  integer, intent(in) :: i
  integer :: ISpec,j_inx,j_step,j
  integer, intent(inout) :: n
  real(doubp), dimension(n+1), intent(inout) :: Y

  j_inx=NY1/2
  
!
! First look where freezing drop is placed
!
  finder: do j_step=1,NY2/2+1   
     if (R_Y1(j_inx) < R_cl(i) .and. R_Y1(j_inx-1) > R_cl(i) ) then
        j_inx = j_inx-1
        exit finder
      endif

     if (R_Y1(j_inx) > R_cl(i)) then
        j_inx=max(2,j_inx-1)
     else
        j_inx=min(j_inx+1,NY2-1)
     endif
  end do finder

  j=max(1,j_inx)

  F_num(j) = 1

!
! Move numbers from liquid distribution to solir  
!
  C_Y1(j) = C_Y1(j) + C_Y2(i)*prob(i+NABins)
!
!  write(*,*)   C_Y2(i),prob(i+NABins),i
  C_Y2(i) = C_Y2(i)*(1.0_doubp-prob(i+NABins))
  T_ice(j) = T
!
! and then water and other compounds
!
  do  ISpec=1,NALiquids + NASolids
     C_Y1(InCi(ISpec,j)) = C_Y1(InCi(ISpec,j)) + C_Y2(InCl(ISpec,i)) * prob(i+NaBins)
     C_Y2(InCl(ISpec,i)) = C_Y2(InCl(ISpec,i)) * (1.0_doubp-prob(i+NaBins))
  end do
!  write(*,*) j,prob(i),C_Y1(j),C(i)*prob(i),'her2'
!
! Update also Y
!
  Y(istart_y+j)=C_Y1(InCi(1,j))
  Y(NYCONCS+NY1+i)=C_Y2(InCl(1,i))

  CALL RADIUSICE2(j)


end subroutine Y2_to_Y1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine num_f_bins
  USE headfile
  implicit none
  intrinsic count
  integer ii, count
!!$
!!$     Program calculates the amount of frozen bins
!!$     
  NFABins = 0
  do ii=1, NFBins
     if (F_num(ii) == 1) then
        NFABins = NFABins + 1
     end if
  end do
  IF(NFABINS.NE.0) NFABins=NY1
  
end subroutine num_f_bins

!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!!$      subroutine increase_Y_length(Y,n)
!!$      USE headfile
!!$      double precision Y(n+1)
!!$ c      n=n+1
!!$      Y((NYAero+NAGases+2):n)=
!!$     >     Y((NYAero+NAGases+1):(n-1))
!!$      Y(NYAero+NAGases+1)=0.C
!!$
!!$      end subroutine increase_Y_length
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 


subroutine freez_prob(k)

!!$ C     Calculates the freezing probability.

!!$ C     References:
!!$ C     MacKenzie A.R., A. Laaksonen, E. Batris and M. Kulmala (1998), 
!!$ C     "The Turnbull correlation and the freezing of stratospheric aerosol
!!$ C     droplets," J. Geophys. Res., Vol. 25, pp. 10875-10884.
!!$ C
!!$ C     Koop et al. Water activity as the determinant for homogeneous ice
!!$ C     nucleation in aqueous solutions (2000), Nature, Vol. 406, pp. 611-614.

  USE headfile
  use precisi
  implicit none
  integer i
  integer, intent(in) :: k
  real(doubp), dimension(NBins) :: v_d
  real(doubp), dimension(NY2) :: v_cl
!!$ C      DOUBLE PRECISION 
!!$ C     First we must choose the nucleation rate calculation method.
!!$ C     1 = Laaksonen et al. or 2 = Koop et al.
!!$  if (rate_method == 1) then
!!$     call ice_nuc_rate_class 
!!$  elseif (rate_method == 2) then
!!$     call ice_nuc_rate_thermo 
!!$  else
!!$     write(*,*) 'Unknown input in nucleation rate calculation, ', &
!!$          &     'rate_method =', rate_method, ', it should be 1 or 2', &
!!$          &     ', program stopped in subroutine freez_prob (can be ', &
!!$          &     'found in liq_ice.f)'
!!$     stop
!!$  end if
  select case(rate_method)
  case (1) !!$ classical ice nucleation
     call ice_nuc_rate_class(k)
  case (2) !!$ Koop's parameterization for ice nuc.
     call ice_nuc_rate_thermo(k)
  end select
!!$ C     Let's calculate the volume of liquid particles

  call ice_het_kiehl
 
  v_d(1:NABins) = (4./3.)*pi*((Rp(1:NABins)**3))
!  pause
!!$ C     It's showtime... and the probability is...
  prob(1:NABins)=0.0_doubp
  z_rate(FirstBinMode(I_MODE_CONTACT) : FirstBinMode(I_MODE_CONTACT)+NBinsPerMode(I_MODE_CONTACT)-1) = 0
  prob(1:NABins)=(1-exp(-DtOut_used*v_d(1:NABins)*(z_rate(1:NABins)+J_immersion(1:NABins))))!!!!!!!!!!
  sum_prob=sum(prob(1:NABins))

  IF(cloud_mod == 1)  THEN
     prob(NaBins+1:NaBins+NY2)=0.0_doubp
     v_cl(1:NY2) = (4./3.)*pi*((R_cl(1:NY2)**3))
     prob(NaBins+1:NaBins+NY2)=(1-exp(-DtOut_used*v_cl(1:NY2)*z_rate(NaBins+1:NaBins+NY2)))
     sum_prob=sum(prob(1:NABins+NY2))
     do i=1, NY2
         if (prob(i+NABins) < ower_limit) then 
        prob(i+NABins) = 0.0_doubp
     elseif (prob(i+NABins) > upper_limit) then
        prob(i+NABins) = upper_limit
     end if
     IF (C_Y2(i)*prob(i+NABINS)<1e-8) prob(i+NABINS)=0.0
     end do
  ENDIF
  
  do i=1, NABins
!!$ C     Here are the lower and upper limit of probability output
     if (prob(i) < ower_limit) then 
        prob(i) = 0.0_doubp
     elseif (prob(i) > upper_limit) then
        prob(i) = upper_limit
     end if
     IF (C(i)*prob(i)<1e-8) prob(i)=0.0
  end do
!!$ C
end subroutine freez_prob


!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine ice_nuc_rate_class(k)
!!$ C     Calculates the ice nucleation rate inside a liquid phase droplet
!!$ C     by classical homogenous nucleation theory for ice (CNT).
!!$ C
!!$ C     References:
!!$ C
!!$ C     [1] MacKenzie A.R., A. Laaksonen, E. Batris and M. Kulmala (1998), 
!!$ C     "The Turnbull correlation and the freezing of stratospheric aerosol
!!$ C     droplets," J. Geophys. Res., Vol. 25, pp. 10875-10884.
!!$ C
!!$ C
!!$ C     Output:
!!$ C     z_rate = nucleation rate [#/(s*m3)] 

  USE headfile
  use precisi
  implicit none
  
  real(doubp), dimension(NABins) :: delta_g_star
  real(doubp) :: delta_h_nucl, v_ice, sig_sl, delta_g
  real(doubp) :: cons_mol, rho_ice
  integer :: i, j,jj
  integer, intent(in) :: k
!!$ C     
!!$ C     First, let's cheack if the number of size classes is big enough
  if (NABins < 1) then
     print *, 'Liquid phase size classes are below 1 in ice_nuc_rate_class'
     write(*,*) 'Value of NABins =', NABins
     stop
  end if
!!$ C
  z_rate(1:NBins)=0.0_doubp
  do i=1, NABins
!!$ C     Next we calculate the nucleation ratio. If the temperature is over
!!$ C     T_ntp, the ice nucleation ratio is set to zero.
     if (T > T_ntp) then
        z_rate(i) = 0.0_doubp
     else
!!$ C
        call water_activity_for_icemodel(i,k) 
        call mass_frac(J)
!!!!        rho_w=rho_water(wtHNO3(i),wtSO4(i),T_drop(i))
!!$ C
!!!        v_h2o = 18.016e-3_doubp/rho_w
!!$ C
        cons_mol=0.0_doubp
        do J=1,(NYLiquids/NABins)
           cons_mol=cons_mol+C(Inc(IH2Ol,i)+(J-1)*NABins)
        end do
        cons_mol=cons_mol*Avog*1.e6_doubp
!!$ C
        delta_h_nucl = 6008.0_doubp + 41.7_doubp*(T_drop(i)-T_ntp)
!!$ C
        v_ice = 18.016e-3_doubp/(rho_ice(T))
!!$ C
        sig_sl = 0.32_doubp*delta_h_nucl/((Avog*v_ice**2)**(1.0/3.0))
!!$ C     Minus is because the integration limits were the opposite way :)
        delta_g=-(6008.0_doubp*(T_drop(i)/T_ntp-1)+T*41.7_doubp*log(T_drop(i)/ &
             &  T_ntp)+41.7_doubp*(T_ntp-T_drop(i)))
        if ((delta_g+R*T_drop(i)*log(a_w(i))) <= 0) then
           delta_g_star(i) = (16.0_doubp*pi*(sig_sl**3)*(v_ice**2))/1.e-36_doubp
        else
           delta_g_star(i) = (16.0_doubp*pi*(sig_sl**3)*(v_ice**2))/ &
                & (3.0*(delta_g+R*T_drop(i)*log(a_w(i)))**2)
        end if
        z_rate(i) = (0.1_doubp*cons_mol*Avog*Boltz*T_drop(i)/ &
             & h_planck)*exp(-delta_g_star(i)/(Boltz*T_drop(i)))  
     end if
  end do

end subroutine ice_nuc_rate_class


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine ice_nuc_rate_thermo(k)
!!$ C     This program calculates the nucleation ratio inside a liquid
!!$ C     phase droplet by thermodynamically consistent parameterization.
!!$ C
!!$ C     References:
!!$ C     Koop et al. Water activity as the determinant for homogeneous ice
!!$ C     nucleation in aqueous solutions (2000), Nature, Vol. 406, pp. 611-614.

  USE headfile
  use precisi
  implicit none
!!$ C
  real(doubp), dimension(NABins+NY2) :: pmax, dmu, aiw, vw0, vi0, grl, daw
  real(doubp) ::  p_G, c_kt0, dkt0dp, c_kti, dktidp, g_iw, g_ii, temperature
  integer :: i,j
  integer, intent(in) :: k
!!$ C     Convert p (mbar) -> p (GPa)
  p_G = TotP*1.0e-7_doubp
!!$ C     Isothermal compressibility of pure water at ambient pressure [GPa^-1]
  c_kt0 = 1.6_doubp
!!$ C     Pressure depedence of KT0 [GPa^-2]
  dkt0dp = -8.8_doubp

!!$ C     Isothermal compressibility of hexagonal ice at ambient pressure [GPa^-1]
  c_kti = 0.22_doubp

!!$ C     Pressure depedence of KTI [GPa^-2]
  dktidp = -0.17_doubp

!!$ C     Pressure depedent part of integral I
  g_iw = (p_G-0.5_doubp*c_kt0*p_G**2-(1.0/6.0)*dkt0dp*p_G**3)
  g_ii = (p_G-0.5_doubp*c_kti*p_G**2-(1.0/6.0)*dktidp*p_G**3)
!!$ C
  z_rate(1:NBins)=0.0_doubp
!!$ C
  do i = 1,NABins
!!$ C     Maximum pressure (GPa)
     pmax(i) = -0.93_doubp+T_drop(i)*(1.37e-2_doubp)-(T_drop(i)**2)*4.12e-5_doubp
     if (k ==2) then !this decreases the temperature in order to check if ice formation
        ! starts within few next time steps (if starts, time step is increased) 
        temperature = T_drop(i)-2.0_doubp*w/10.0_doubp
        pmax(i) = -0.93_doubp+temperature*(1.37e-2_doubp)-(temperature**2)*4.12e-5_doubp
     end if
!     pmax(i) = -0.93+T_drop(i)*1.37d-2-(T_drop(i)**2)*4.12d-5
!!$ C     Check pressure.
     z_rate(i) = 0.0_doubp
!     write(*,*) p_G ,pmax(i)
     if(p_G >= pmax(i)) then
        z_rate(i) = 0.0_doubp
     else
        if ( T <= 240.0_doubp .and. T >= 170.0_doubp ) then
           call water_activity_for_icemodel(i,k) 

!!$ C     Difference between ice and water chemical potentials
           dmu(i) = 210368.0_doubp+131.438_doubp*T_drop(i)-(3.32373e6_doubp)/T_drop(i) &
                &   -41729.1_doubp*log(T_drop(i))
!!$ C     Water activity in a solution in equilibrium with ice
           aiw(i) = exp(dmu(i)/(R*T_drop(i)))
!!$ C     Molar volume of pure liquid water
           vw0(i) = -230.76_doubp-0.1478_doubp*T_drop(i)+4099.2_doubp &
                & /T_drop(i)+48.8341_doubp*log(T_drop(i))
!!$ C     Molar volume of hexagonal ice
           vi0(i) = 19.43_doubp-(2.2e-3_doubp)*T_drop(i)+(1.08e-5_doubp)*(T_drop(i)**2)
!!$ C     Integral
           grl(i) = vw0(i)*g_iw - vi0(i)*g_ii
!!$ C
           if (a_w(i) > 1) a_w(i)=1.0_doubp
!!$ C     Water activity change
           daw(i) = a_w(i)*exp(grl(i)/(R*T_drop(i))) - aiw(i)
           if ( daw(i) >= 0.26_doubp .and. daw(i) <= 0.34_doubp ) then
!!$ C     Nucleation rate of ice in a solution. Check validity ranges. In Koop et at (2000),
!!$ C     validity range is 170K < TA < 240K and P < PMAX. Upper temperature limit check.
              z_rate(i) = 10.0_doubp**(-906.7_doubp+8502.0_doubp*daw(i)- &
                   & 26924.0_doubp*daw(i)**2+29180.0_doubp*daw(i)**3) 
           else
              z_rate(i) = 0.0_doubp
           end if
        else
           z_rate(i) = 0.0_doubp
        end if
     end if
  end do

!!$ C     Convert the homogeneous ice nucleation rate coefficient from [#/(s*cm3)]
!!$ C     to [#/(s*m3)]
  z_rate(1:NBins) = z_rate(1:NBins)*1e6_doubp



  IF (cloud_mod == 1) THEN
     
     koop:  do j = 1,NY2
        IF(C_Y2(j)<1e-10) cycle koop
        i=j+NaBins
!!$ C     Maximum pressure (GPa)
        pmax(i) = -0.93_doubp+T_drop(i)*(1.37e-2_doubp)-(T_drop(i)**2)*4.12e-5_doubp
        if (k ==2) then !this decreases the temperature in order to check if ice formation
        ! starts within few next time steps (if starts, time step is increased) 
           temperature = T_drop(i)-2.0_doubp*w/10.0_doubp
           pmax(i) = -0.93_doubp+temperature*(1.37e-2_doubp)-(temperature**2)*4.12e-5_doubp
        end if
!     pmax(i) = -0.93+T_drop(i)*1.37d-2-(T_drop(i)**2)*4.12d-5
!!$ C     Check pressure.
        z_rate(i) = 0.0_doubp
        !     write(*,*) p_G ,pmax(i)
        if(p_G >= pmax(i)) then
           z_rate(i) = 0.0_doubp
        else
           if ( T <= 240.0_doubp .and. T >= 170.0_doubp ) then
              call water_activity_for_icemodel(i,k) 
!!$ C     Difference between ice and water chemical potentials
              dmu(i) = 210368.0_doubp+131.438_doubp*T_drop(i)-(3.32373e6_doubp)/T_drop(i) &
                   &   -41729.1_doubp*log(T_drop(i))
!!$ C     Water activity in a solution in equilibrium with ice
              aiw(i) = exp(dmu(i)/(R*T_drop(i)))
!!$ C     Molar volume of pure liquid water
              vw0(i) = -230.76_doubp-0.1478_doubp*T_drop(i)+4099.2_doubp &
                   & /T_drop(i)+48.8341_doubp*log(T_drop(i))
!!$ C     Molar volume of hexagonal ice
              vi0(i) = 19.43_doubp-(2.2e-3_doubp)*T_drop(i)+(1.08e-5_doubp)*(T_drop(i)**2)
!!$ C     Integral
              grl(i) = vw0(i)*g_iw - vi0(i)*g_ii
!!$ C
              if (a_w(i) > 1) a_w(i)=1.0_doubp
!!$ C     Water activity change
              daw(i) = a_w(i)*exp(grl(i)/(R*T_drop(i))) - aiw(i)
              if ( daw(i) >= 0.26_doubp .and. daw(i) <= 0.34_doubp ) then
!!$ C     Nucleation rate of ice in a solution. Check validity ranges. In Koop et at (2000),
!!$ C     validity range is 170K < TA < 240K and P < PMAX. Upper temperature limit check.
                 z_rate(i) = 10.0_doubp**(-906.7_doubp+8502.0_doubp*daw(i)- &
                      & 26924.0_doubp*daw(i)**2+29180.0_doubp*daw(i)**3) *1e6_doubp
              else
                 z_rate(i) = 0.0_doubp
              end if

           else
              z_rate(i) = 0.0_doubp
           end if
        end if
     end do koop
  END IF



!!$ C
end subroutine ice_nuc_rate_thermo

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

function rho_water(xm22,xm33,T) result(rho)

!!$ c       Based on article, E.Martin C. George and P. Mirabel: Densities
!!$ c       and Surface Tensions of H2SO4/HNO3/H2O Solutions. Geophysical
!!$ c       Research Letters 27, no.2 197-200, 2000
!!$ c       xm2 = HNO3 weight fraction, xm3 = H2SO4 weight fraction, T =
!!$ c       temperature
  use precisi
  implicit none
  real(doubp), intent(in) :: xm22, xm33, T
  real(doubp) :: rho, y, x, rho293, rho253, kk, x2, x3, y2, y3
  real(doubp) :: xm2, xm3
!1  real, dimension(2) :: hh

  if (xm22+xm33 > 0.8) then
     xm3 = .8/(1.+xm22/xm33)
     xm2 = 0.8-xm3
  else
     xm2 = xm22
     xm3 = xm33
  end if

  y = 100.0_doubp*(xm2)
  x = 100.0_doubp*(xm3)
  x2 = x**2
  x3 = x**3
  y2 = y**2
  y3 = y**3
  
!!$ c       CCCCCCCCCCCCC  T = 293;   CCCCCCCCCCCCCC
  
  rho293 =   0.9982_doubp+ 5.8487e-3_doubp*y -2.3873e-5_doubp*y2 + &
       &     7.9815e-7_doubp*y3 + 7.9119e-3_doubp*x + 2.9369e-5_doubp*x*y&
       &     + 1.8949e-6_doubp*x*y2-3.6905e-8_doubp*x*y3 -7.6431e-5_doubp*x2 &
       &     + 2.8093e-6_doubp*x2*y  -6.4247e-8_doubp*x2*y2  + 2.2885e-6_doubp*x3 &
       &     -2.1422e-8_doubp*x3*y -1.4651e-8_doubp*x**4
!!$ c
!!$ c       CCCCCCCCCCCCCCCC  T = 253;   CCCCCCCCCCCCCCCCC
!!$ c
  rho253 = 1.0015_doubp + 9.7509e-3_doubp*y -1.834e-4_doubp*y2+ &
       &   2.9113e-6_doubp*y3 + 9.6589e-3_doubp*x-3.9433e-5_doubp*x*y&
       &   + 3.8149e-6_doubp*x*y2 - 6.3144e-8_doubp*x*y3 &
       &     -1.1562e-4_doubp*x2 + 4.3442e-6_doubp*x2*y - 7.0749e-8_doubp*x2*y2 + &
       &     2.6848e-6_doubp*x3 - 3.6871e-8_doubp*x3*y - 1.6015e-8_doubp*x**4
!!$ c
!!$ cc      CCCCCCCCCCCCCCCC  T = T  CCCCCCCCCCCCCCCCCC
!!$ c
  kk = (rho293-rho253)/40.0_doubp
  rho = kk*(T-253.0_doubp)+rho253
  rho = 1.0e3_doubp*rho
!!$ c       rho=1000.
end function rho_water


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



function rho_ice(T) result(rho)
!!$ C
!!$ C     This program calculates the density of ice in the specified temperature
!!$ C
!!$ C     Input:
!!$ C     T = Temperature [K]
!!$ C
!!$ C     Output:
!!$ C     rho = density of ice [kg/m3]
!!$ C
!!$ C     References:
!!$ C     Pruppacher and Klett (1997), Microphysics of clouds and precipitation,
!!$ C     Kluwer Academic Publishers, The Netherlands, 2nd edition, pages 79-80
!!$ C     
  use precisi
  implicit none
  real(doubp), intent(in) :: T
  real(doubp) :: rho, a0, a1, a2, T_C, rho_apu

  a0 = 0.9167_doubp
  a1 = -1.75e-4_doubp
  a2 = -5.0e-7_doubp
  
  T_C = T - 273.15_doubp
  rho_apu = a0 + a1*T_c+a2*T_C**2
!!$ C     [g/cm3] -> [kg/m3]
  rho = rho_apu*1.e3_doubp

end function rho_ice

!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine sat_vap_press_ice
  USE headfile
  use precisi
  implicit none
!  integer ind_s

!  ind_s = NFBins - NFAbins +1
  select case(e_sat_method)
  case(1) ! Magnus equation
     e_sat_i=6.1064_doubp*exp((21.88_doubp*(T-T_ntp))/(T-7.65_doubp))
!!$ C     [hPa] -> [pa]
     e_sat_i = 100.0_doubp*e_sat_i
!!$ C         
  case(2) ! Goff Gratch
     e_sat_i=10.0_doubp**(-9.09718_doubp*(273.16_doubp/T-1)-3.56654_doubp* &
          & log10(273.16_doubp/T)+ 0.876793_doubp*(1-T/273.16_doubp) + &
          & log10(6.1071_doubp))
!!$ C     [hPa] -> [pa]
     e_sat_i = 100.0_doubp*e_sat_i
  case default
     write(*,*) 'unknow method for determining the saturation ', &
          &     'vapor pressure over ice, e_sat_method =',  &
          &     e_sat_method, ' (it should be 1 or 2), error in subroutine', &
          &     ' sat_vap_press_ice'
     stop
  end select
  
end subroutine sat_vap_press_ice

!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine init_Y2C(Y,n)
  use headfile
  use precisi
  implicit none
  integer, intent(inout) :: n
  real(doubp), intent(inout), dimension(n) :: Y
  integer :: J, IFam, JJ, I, L3, ISpec
  real(doubp) :: so_mass_sum
!!$ C      
!!$ C     This program moves the values from Y to C before freezing program
!!$ C
  DO J=1,NYLiquids
     C(IY2C(J))=y(J)
  END DO
!!$ C     
  if (ice_mod == 1) then
     do  J=1,NY1
        C_Y1(InCi(IH2Ol,J)) = Y(NYConcs + J)
     end do
  end if
!!$ C
  J = NYLiquids
  DO IFam = 1,NAFl
     IF (NAFMbr(IFam) >= 2) THEN
        DO L3 = 1,NABins
           J = J+1
           TotF(IFam) = 0.0_doubp
           DO JJ = 1, NFMbr(IFam)
              IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
                 TotF(IFam) = TotF(IFam) + C(InC(IFMbr(IFam,JJ),L3))
              END IF
           END DO
           DO JJ = 1, NFMbr(IFam)
              IF(State(IFMbr(IFam,JJ)) /= 'G') THEN
                 C(InC(IFMbr(IFam,JJ),L3)) = C(InC(IFMbr(IFam,JJ),L3)) / &
                      & TotF(IFam)*y(J) 
              END IF
           END DO
        END DO
     END IF
  END DO
!!$ C     
  DO J = 1, NAGases
     C(NAero+J)=y(J + NYAero)
  END DO
!!$ C     
!  if (cloud_mod == 1) then
!     c_cl_tot      = y(NODEs-6)
!  endif
!  if (ice_mod == 1) then
!     C_ice_tot     = y(NODEs-5)
!  end if
!  Adia_fex         = y(NODEs-4)
!  z                = y(NODEs-3)
!  P                = y(NODEs-2)
!  Ctot             = y(NODEs-1)
!  T                = y(NODEs)
!!$ C
  DO I = 1, NABins
     CALL RADIUS(I)
  END DO
  DO I = 1,NY1
     CALL RADIUSICE2(I)
  END DO
  DO I = 1,NY2
     CALL RADIUS_Y2(I)
  END DO
!!$ C
  if (ice_mod == 1) then
     DO  I = 1,NY1
        so_mass_sum = 0.0_doubp
        IF(C_Y1(I)>0.0) THEN
           DO ISpec = 1, NALiquids + NASolids
              so_mass_sum = so_mass_sum + C_Y1(InCi(ISpec,I))*WtMol(ISpec)/C_Y1(I)
           END DO
!           CALL RADIUSICE2(I)
        ENDIF
        so_mass(I) = so_mass_sum
       
     END DO
  end if
!!$ C
end subroutine init_Y2C

!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine water_activity_for_icemodel(J,k)
  USE headfile
  use precisi
  implicit none
!!$ C
!  real(doubp), intent(in) :: ASOLV
  real(doubp) :: pwe, phno3, phcl, surface_tension, temperature
!  real(doubp) :: wat_act_old
  integer,intent(in) :: J, k
!!$ C
  !  write(*,*) 'kekelol'
  if (ice_mod ==  0) return
  select case(Ithermo_mod)
  case(1) ! Luo's thermodynamics
     !C     Water activity
     pwe   = 0.0_doubp 
     phno3 = 0.0_doubp 
     phcl  = 0.0_doubp 
     CALL s_luo(T_drop(J),0.0_doubp,0.0_doubp,0.0_doubp,0.0_doubp,pwe,phno3,phcl)
     call mass_frac(J)
     temperature= T_drop(J)
     a_w(J) = PSurf(IH2Ol,J)/exp(2.0_doubp*surface_tension(temperature, &
          & wtSO4(J)*100.0_doubp,wtHNO3(J)*100.0_doubp)*WtMol(IH2Ol)/(Rp(J)*R*T_drop(J)* &
          & 1.e6_doubp))/(pwe*100.0_doubp)
  case(2) !Jacobson thermo
     a_w(J)=watakti(J)
  case(3) !tomi
     a_w(J)=Water_activity(J)
  case(4) !AIM
     a_w(J)=Water_activity(J)
  end select
  if (wat_act_old(J) == 0) wat_act_old(J) = 1.0_doubp
  !this increases the water activity in order to check if ice formation
  !starts within few next time steps (if starts, decreases the time step)
  !if (k == 2) wat_act_old(J) = Water_activity(J)*Water_activity(J)/ &
  !     & wat_act_old(J)*(1.0_doubp+w/100.0_doubp)
!  write(*,*) k,Water_activity(J),J
  
  if (k == 2 .and. Water_activity(J)/wat_act_old(J) >= 1.0_doubp) then
     Water_activity(J) = Water_activity(J)*Water_activity(J)/ &
          & wat_act_old(J)*(1.0_doubp+w/1000.0_doubp)
  else if (k == 2 .and. Water_activity(J)/wat_act_old(J) < 1.0_doubp) then
     Water_activity(J) = Water_activity(J)*wat_act_old(J)/ &
          & Water_activity(J)*(1.0_doubp+w/1000.0_doubp)
end if
! write(*,*) k,Water_activity(J),J
!  write(*,*) a_w(J) , J
end subroutine water_activity_for_icemodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine mass_frac(J)
  USE headfile
  use precisi
  implicit none
!!$ C
  integer, intent(in) :: J
  real(doubp) :: CSVITEMP, CHNO3EMP, CH2OTEMP, WtTOT
  integer :: Itot, IFam

  CSVITEMP = 0.0_doubp
  CHNO3EMP = 0.0_doubp
  CH2OTEMP = 0.0_doubp
  CH2OTEMP = C(InC(IH2Ol,J))*WtMol(IH2Ol)
!     CNH3TEMP = 0.0_doubp
  WtTOT    = 0.0_doubp
  do Itot=1,NALiquids
     WtTot=WtTot+C(InC(Itot,J))*WtMol(Itot)
  end do

  LOOP1: do IFam=1,nafmbr(1)
     if (State(IFMbr(1,IFam)) == 'G') cycle LOOP1 
     CSVITEMP = CSVITEMP + C(InC(ifmbr(1,IFam),J))*WtMol(ifmbr(1,IFam))
  end do LOOP1

  LOOP2: do IFam=1,nafmbr(2)
     if (State(IFMbr(2,IFam)) == 'G') cycle LOOP2
     CHNO3EMP = CHNO3EMP + C(InC(ifmbr(2,IFam),J))*WtMol(ifmbr(2,IFam))
  end do LOOP2

  wtHNO3(J) = CHNO3EMP/WtTot
  WTSO4(J) = CSVITEMP/WtTot
  wtH2O(J) = CH2OTEMP/WtTot

end subroutine mass_frac

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


SUBROUTINE initialize
!!
  USE headfile
  use precisi
    implicit none
    integer I
!!$     This program initializes the frozen mode variables to zero
!!$     
  T_ice(1:NFBins) = 0.0_doubp
  R_ice(1:NFBins) = 0.0_doubp
  DO I=1,NFBins+NFBins*NLiquids+NFBins*NSolids
     C_i(I)= 0.0_doubp
     C_Y1(I)= 0.0_doubp
  ENDDO
!  I=NFBins+NFBins*NLiquids+NFBins*NSolids
!  C_i(1:I)  = 0.0
!  C_Y1(1:(NFBins+NFBins*NLiquids+NFBins*NSolids))  = 0.0_doubp

  C_ice_tot = 0.0_doubp
  F_num(1:NFBins) = 0.0_doubp
  liq_mass(1:NFBins) = 0.0_doubp
  so_mass(1:NFBins) = 0.0_doubp
  thetaHNO3(1:NFBins)=0.0_doubp
  theta_old=0.0_doubp
  NFABins=0
  ads_ice=0.0_doubp
  wat_act_old(1:NABins) = 0.0_doubp
  DtOut_used = DtOut

end subroutine initialize

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


SUBROUTINE halletmossop(dCdx,dC_Y1dx,dC_Y2dx,Ibin,I_ice,B_rate)

!
! SUBROUTINE TO CALCULATE ICE SPLINTERING WITH HALLET-MOSSOP (1974) 
! AND MOSSOP (1976) METHODS.
!
  USE headfile
  use precisi
  implicit none
  REAL(doubp), intent(inout), DIMENSION(NABins + NABins * NAerT)  :: dCdx
  REAL(doubp), intent(inout), DIMENSION(NY1 + NY1 * NAerT)  :: dC_Y1dx
  REAL(doubp), intent(inout), DIMENSION(NY2 + NY2 * NAerT)  :: dC_Y2dx
  REAL(doubp), intent(in) :: B_rate
  REAL(doubp) :: rho_ice
  integer, intent(in) :: Ibin,I_ice
  I_HM=I_HM_NEW
  IF(I_HM_METHOD == 1) THEN
     IF(T_drop(Ibin) < 270.16 .AND. T_drop(Ibin) > 268.16) THEN
        dC_Y1dx(I_HM)= dC_Y1dx(I_HM) + 3.5d5*dmdt(I_ice)*(T_drop(Ibin)-268.16)/2.0
        dC_Y1dx(I_HM+NY1)= dC_Y1dx(I_HM+NY1) + 3.5d5*dmdt(I_ice)*(T_drop(Ibin)-268.16)/2.0 & 
             &  *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
        dC_Y1dx(I_ice+NY1)=  dC_Y1dx(I_ice+NY1) - 3.5d5*dmdt(I_ice)*(T_drop(Ibin)-268.16)/2.0 & 
             &  *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
     ELSEIF(T_drop(Ibin) < 268.16 .AND. T_drop(Ibin) > 265.16) THEN
        dC_Y1dx(I_HM)= dC_Y1dx(I_HM)+3.5d5*dmdt(I_ice)*(T_drop(Ibin)-265.16)/3.0
        dC_Y1dx(I_HM+NY1)= dC_Y1dx(I_HM+NY1) + 3.5d5*dmdt(I_ice)*(T_drop(Ibin)-265.16)/3.0 & 
             &  *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
        dC_Y1dx(I_ice+NY1)=  dC_Y1dx(I_ice+NY1) - 3.5d5*dmdt(I_ice)*(T_drop(Ibin)-265.16)/3.0 & 
             &  *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
!        write(*,*)  dC_Y1dx(I_HM), dC_Y1dx(I_HM+NY1),dC_Y1dx(I_ice+NY1),'2'
     ENDIF
!     write(*,*)  dC_Y1dx(55)
     

  ELSE
     IF(IBin>NABins) THEN
        IF(T_drop(Ibin) < 270.16 .AND. T_drop(Ibin) > 268.16 .AND. R_cl(Ibin-NaBins)>12.0e-6) THEN
           dC_Y1dx(I_HM)=dC_Y1dx(I_HM) + 1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0
           dC_Y1dx(I_HM+NY1)= dC_Y1dx(I_HM+NY1)+1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
           dC_Y1dx(I_ice+NY1)=  dC_Y1dx(I_ice+NY1)-1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
        ELSEIF(T_drop(Ibin) < 268.16 .AND. T_drop(Ibin) > 265.16 .AND. R_cl(Ibin-NaBins)>12.0e-6) THEN
           dC_Y1dx(I_HM) = dC_Y1dx(I_HM) + 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0
           dC_Y1dx(I_HM+NY1) = dC_Y1dx(I_HM+NY1) + 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
           dC_Y1dx(I_ice+NY1) =  dC_Y1dx(I_ice+NY1) - 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL) 
        ENDIF
     ELSE
        IF(T_drop(Ibin) < 270.16 .AND. T_drop(Ibin) > 268.16 .AND. Rp(Ibin)>12.0e-6) THEN
           dC_Y1dx(I_HM)=dC_Y1dx(I_HM) + 1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0
           dC_Y1dx(I_HM+NY1)= dC_Y1dx(I_HM+NY1)+1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
           dC_Y1dx(I_ice+NY1)=  dC_Y1dx(I_ice+NY1)-1.0/250.0*B_rate*(T_drop(Ibin)-268.16)/2.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
        ELSEIF(T_drop(Ibin) < 268.16 .AND. T_drop(Ibin) > 265.16 .AND. Rp(Ibin)>12.0e-6) THEN
           dC_Y1dx(I_HM) = dC_Y1dx(I_HM) + 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0
           dC_Y1dx(I_HM+NY1) = dC_Y1dx(I_HM+NY1) + 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL)
           dC_Y1dx(I_ice+NY1) =  dC_Y1dx(I_ice+NY1) - 1.0/250.0*B_rate*(T_drop(Ibin)-265.16)/3.0 &
                &   *4./3.*pi*R_ice(I_HM)**3*rho_ice(T)*1.d3/WtMol(IH2OL) 
        ENDIF
        
     ENDIF
  ENDIF
END SUBROUTINE halletmossop


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


SUBROUTINE CONTACT_FREEZING(Ibin,Icontact,frac_ice,I_ice)
!
! This is a subroutine to calculate fraction of collision between 
! particles in mode I_MODE_CONTACT with cloud droplet leads to freezing.
!
! Parameterization from:
! MEYERS et al.: New primary ice-nucleation parameterizations in an explicit 
! cloud model. J. Appllied meteorol., 31, 708--, 1992
!
  USE headfile
  use precisi
  implicit none
  REAL(doubp) :: rho_ice
  REAL(doubp) :: N_cont,N_cont_min
  integer, intent(in) :: Ibin,Icontact
  real(doubp), intent(out) :: frac_ice
  integer ::III
  integer, intent(out) :: I_ice

  N_cont_min=1.0e-1
  N_cont = 0.0
  IF(T_drop(Ibin) <270.16) THEN
     N_cont=min(N_cont_min*(270.16-T_drop(Ibin))**1.3,  &
          &  sum(C(FirstBinMode(I_MODE_CONTACT) : FirstBinMode(I_MODE_CONTACT)+NBinsPerMode(I_MODE_CONTACT)-1)))
  ENDIF
  
  frac_ice = N_cont/sum(C(FirstBinMode(I_MODE_CONTACT) : FirstBinMode(I_MODE_CONTACT)+NBinsPerMode(I_MODE_CONTACT)-1))

  ice:DO III = 1, NY1
     IF (Rp(Ibin) > R_Y1(III) .AND. &
          &  Rp(Ibin) < R_Y1(III+1)) THEN
        I_ice=III
        cycle ice
     ENDIF
  ENDDO ice

END SUBROUTINE CONTACT_FREEZING


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


subroutine ice_het_kiehl
!C   
!C  Immersion freezing based on article:  
!C    
!C  K. DIEHL AND S. WURZLERS: Heterogeneous Drop Freezing in the Immersion Mode: 
!C  Model Calculations Considering Soluble and Insoluble Particles in the Drops
!C  J. Atmos. Sci 
!C

  USE headfile
  use precisi
  implicit none
!!$ C
  real(doubp) :: find_freez_dt,p_G
  integer :: i,j,k
  J_immersion(1:NBins)=0.0_doubp
!!$ C
  k=1
  B_kiehl=1.0e-1*1.0e6
  a_kiehl=1.0
  p_G = TotP*1.0e-7_doubp
  do i = 1,NABins

     call water_activity_for_icemodel(i,k) 

     
     IF(IModeBin(I) == I_MODE_IMMERSION.AND.IF_FREEZ_IMMERSION==1) J_immersion(i)=-a_kiehl*B_kiehl*  &
          &   exp(a_kiehl*(273.15-(T_drop(i)+find_freez_dt(p_G,a_w(i)))))*dtemp
  end do



!!$ C
end subroutine ice_het_kiehl


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


function find_freez_dt(p_G,aw1)
!
!
! This function calculates the freezing point depression due to solutes.
! Koop form for homogeneous freezing used
!

  use precisi
  implicit none

  real(doubp),intent(in) :: p_G,aw1
  real(doubp) :: TTMIN,TTMAX,TT, pmax, z_rate,find_freez_dt,TT2
  real(doubp) :: dmu, vw0, aiw, vi0, grl, daw,R,aw
  real(doubp) :: c_kt0, dkt0dp, c_kti,dktidp, g_iw, g_ii
  integer :: I
  R = 8.3145
  TTMIN=200.0
  TTMAX=240.0
 
  aw=1.0

  c_kt0 = 1.6_doubp
!!$ C     Pressure depedence of KT0 [GPa^-2]
  dkt0dp = -8.8_doubp

!!$ C     Isothermal compressibility of hexagonal ice at ambient pressure [GPa^-1]
  c_kti = 0.22_doubp

!!$ C     Pressure depedence of KTI [GPa^-2]
  dktidp = -0.17_doubp

!!$ C     Pressure depedent part of integral I
  g_iw = (p_G-0.5_doubp*c_kt0*p_G**2-(1.0/6.0)*dkt0dp*p_G**3)
  g_ii = (p_G-0.5_doubp*c_kti*p_G**2-(1.0/6.0)*dktidp*p_G**3)


  DO I = 1,15
     TT=(TTMIN+TTMAX)/2.0
     pmax = -0.93_doubp+TT*(1.37e-2_doubp)-(TT**2)*4.12e-5_doubp

     z_rate = 0.0_doubp

     if(p_G >= pmax) then
        z_rate = 0.0_doubp
     else
        if ( TT <= 240.0_doubp .and. TT >= 170.0_doubp ) then
!!$ C     Difference between ice and water chemical potentials
           dmu = 210368.0_doubp+131.438_doubp*TT-(3.32373e6_doubp)/TT &
                &   -41729.1_doubp*log(TT)
!!$ C     Water activity in a solution in equilibrium with ice
           aiw = exp(dmu/(R*TT))
!!$ C     Molar volume of pure liquid water
           vw0 = -230.76_doubp-0.1478_doubp*TT+4099.2_doubp &
                & /TT+48.8341_doubp*log(TT)
!!$ C     Molar volume of hexagonal ice
           vi0 = 19.43_doubp-(2.2e-3_doubp)*TT+(1.08e-5_doubp)*(TT**2)
!!$ C     Integral
           grl = vw0*g_iw - vi0 *g_ii
!!$ C
           if (aw > 1) aw=1.0_doubp
!!$ C     Water activity change
           daw = aw*exp(grl/(R*TT)) - aiw
           if ( daw >= 0.26_doubp .and. daw <= 0.34_doubp ) then
!!$ C     Nucleation rate of ice in a solution. Check validity ranges. In Koop et at (2000),
!!$ C     validity range is 170K < TA < 240K and P < PMAX. Upper temperature limit check.
              z_rate = 10.0_doubp**(-906.7_doubp+8502.0_doubp*daw- &
                   & 26924.0_doubp*daw**2+29180.0_doubp*daw**3)* 1e6_doubp
           else
              z_rate = 0.0_doubp
           end if
        else
           z_rate = 0.0_doubp
        end if
     end if
     IF((1-exp(-1.0*1.0e-15*z_rate))>0.5.OR.daw>0.34_doubp) THEN
        TTMIN=TT
     ELSE
        TTMAX=TT
     ENDIF
  END DO
  

  TTMIN=180.0
  TTMAX=240.0
  aw=aw1
  
 DO I = 1,15
     TT2=(TTMIN+TTMAX)/2.0
     pmax = -0.93_doubp+TT2*(1.37e-2_doubp)-(TT2**2)*4.12e-5_doubp

     z_rate = 0.0_doubp

     if(p_G >= pmax) then
        z_rate = 0.0_doubp
     else
        if ( TT2 <= 240.0_doubp .and. TT2 >= 170.0_doubp ) then
!!$ C     Difference between ice and water chemical potentials
           dmu = 210368.0_doubp+131.438_doubp*TT2-(3.32373e6_doubp)/TT2 &
                &   -41729.1_doubp*log(TT2)
!!$ C     Water activity in a solution in equilibrium with ice
           aiw = exp(dmu/(R*TT2))
!!$ C     Molar volume of pure liquid water
           vw0 = -230.76_doubp-0.1478_doubp*TT2+4099.2_doubp &
                & /TT2+48.8341_doubp*log(TT2)
!!$ C     Molar volume of hexagonal ice
           vi0 = 19.43_doubp-(2.2e-3_doubp)*TT2+(1.08e-5_doubp)*(TT2**2)
!!$ C     Integral
           grl = vw0*g_iw - vi0 *g_ii
!!$ C
           if (aw > 1) aw=1.0_doubp
!!$ C     Water activity change
           daw = aw*exp(grl/(R*TT2)) - aiw
           if ( daw >= 0.26_doubp .and. daw <= 0.34_doubp ) then
!!$ C     Nucleation rate of ice in a solution. Check validity ranges. In Koop et at (2000),
!!$ C     validity range is 170K < TA < 240K and P < PMAX. Upper temperature limit check.
              z_rate = 10.0_doubp**(-906.7_doubp+8502.0_doubp*daw- &
                   & 26924.0_doubp*daw**2+29180.0_doubp*daw**3)* 1e6_doubp
           else
              z_rate = 0.0_doubp
           end if
        else
           z_rate = 0.0_doubp
        end if
     end if
     IF((1-exp(-1.0*1.0e-15*z_rate))>0.5.OR.daw>0.34_doubp) THEN
        TTMIN=TT2
     ELSE
        TTMAX=TT2
     ENDIF
  END DO

  find_freez_dt=TT-TT
end function find_freez_dt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine ice_deposition(Y,n)
!C   
!C  This subroutine calculates the (ice)activated fraction 
!C
!C   
!C   
!C   

  USE headfile
  use precisi
  implicit none
!!$ C
  integer, intent(inout) :: n
  real(doubp), dimension(n+NFBins), intent(inout) :: Y
  real(doubp) :: f_act,S0,a_f
  integer I, J,ISpec,j_inx,j_step


!  Current equation type from:
!  Mohler et al.: Efficiency of the deposition mode ice nucleation on mineral dust particles
!  Atmos. Chem. Phys., 6, 3007-3021, 2006
!
  a_f=6.03
  S0 = 1.07

  CALL sat_vap_press_ice

  f_act=exp(a_f*(C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp/e_sat_i - S0))-1

  freezdep:  do j=FirstBinMode(I_MODE_DEPOS),FirstBinMode(I_MODE_DEPOS)+NBinsPerMode(I_MODE_DEPOS)-1
     
     if (C(InC(IH2Ogi,NABins)) * R * T * 1.0e6_doubp/e_sat_i < 1.0_doubp) cycle freezdep
     if (f_act<f_act_max) cycle freezdep
     

! Find the bin where ice is placed     
     j_inx=NY1/2
     
     finder: do j_step=1,NY1/2+1   
        if (R_Y1(j_inx) < Rp(j) .and. R_Y1(j_inx-1) > Rp(j) ) then
           j_inx = j_inx-1
           exit finder
        endif
        
        if (R_Y1(j_inx) > Rp(j)) then
           j_inx=max(2,j_inx-1)
        else
           j_inx=min(j_inx+1,NY1-1)
        endif
     end do finder
     
     I=max(1,j_inx)

! Change in number concetration due to freezing     
     C_Y1(I) = C_Y1(I) + C(j)*(f_act-f_act_dep)
     C(j) = C(j)- C(j)*(f_act-f_act_dep)
    
     do  ISpec=1,NALiquids + NASolids
        C_Y1(InCi(ISpec,I)) = C_Y1(InCi(ISpec,I)) + C(InC(ISpec,J)) * (f_act-f_act_dep)
        C(InC(ISpec,J)) = C(InC(ISpec,J))- C(InC(ISpec,J)) *(f_act-f_act_dep)
     end do
     
     Y(NYCONCS+j_inx)=C_Y1(InCi(1,j_inx))
     Y(j)=C(InC(1,j))
     
  end do freezdep

  f_act_dep=f_act
  if(f_act>f_act_max) f_act_max=f_act
  

end subroutine ice_deposition
