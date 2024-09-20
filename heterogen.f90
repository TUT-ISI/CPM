subroutine ice_nuc_rate_het(k)
!!$ C     Calculates the ice nucleation rate inside a liquid phase droplet
!!$ C     by classical heterogeneous nucleation rate
!!$ C
!!$ C     References:
!!$ C
!!$ C     KHVOROSTYANOV, V.I.,  CURRY, J.A.: The Theory of Ice Nucleation 
!!$ C     by Heterogeneous Freezing of Deliquescent Mixed CCN.
!!$ C     Part I: Critical Radius, Energy, and Nucleation Rate'
!!$ C     J. Atmos. Sci 61, 2676 -->
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
  integer :: i, j
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
 
        
        m_is=(sigma_ns-sigma_ni)/sigma_is
        
!!! EQ (2.11)
        
        y1 = (1.0_doubp-2.0_doubp*m_is*x1+x1**2)**(.5)
        ksii = (x1-m_is)/y1

        f_m = .5_doubp*(1.0_doubp+((1.0_doubp-m_is*x1)/y1)**3 &
             &  +x1**3*(2.0_doubp-3.0_doubp*ksii+ksii**3) &
             &  3.0_doubp*m_is*x1**2*(ksii-1.0_doubp))



!!! EQ (2.10)        
        
        delta_Fcr=(16.0_doubp*Pi/3.0_doubp)*sigma_is**3*f_m/ &
             &  (rho_i*L_m*log(T0/T1*S_w)-C1*epsilon**2 &
             &  -r_sc/r_d)**2 -alpha*r_n**2*(1.0_doubp-m_is)
















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

end subroutine ice_nuc_rate_het
