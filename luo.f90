subroutine s_luo(T,wt1,wt2,wt3,wt4,ph,phno3,phcl)
  use precisi
  implicit none
  
  REAL(doubp), intent(in) :: T,wt1,wt2,wt3,wt4
  REAL(doubp), intent(out) :: ph,phno3,phcl
  REAL(doubp) :: a_hno3, b_hno3,a_hcl,b_hcl,xm2,xm3,xm4,xm5
!!!  character(len=10) kk
!!$C Pressures from Luo et al. GRL 1995 or something.
!!$C Be care for units. 
  xm2 = wt3
  xm3 = wt1      
  xm5 = 0.0_doubp
  if ( xm2 > 0.0_doubp ) then
     a_hno3 = 22.74_doubp + 29.0_doubp*xm2 - 0.457_doubp*xm3 - 5.03_doubp*sqrt(xm2) &
          & + 1.447_doubp*sqrt(xm3) - 29.72_doubp*xm2*xm2  &
          & - 13.9_doubp*xm2*xm3 + 6.1_doubp*xm3*xm3
     b_hno3 = -7689.8_doubp - 2896.1_doubp*xm2 + 2859.8_doubp*xm3 - 274.2_doubp*sqrt(xm2) &
          &   - 389.5_doubp*sqrt(xm3) + 7281.0_doubp*xm2*xm2 + 6475.0_doubp*xm2*xm3 &
          &   + 801.0_doubp*xm3*xm3
     phno3 =exp(a_hno3+ b_hno3/T + log(xm2*(xm2 + 0.09_doubp*xm3)))
  else
     phno3 = 0.0_doubp
  end if
	
  if ( xm5 > 0.0_doubp ) then
     a_hcl = 21.0_doubp + 46.61_doubp*xm2 + 4.069_doubp*xm3 - 4.837_doubp*sqrt(xm2) &
          & + 2.186_doubp*sqrt(xm3) - 63.0_doubp*xm2*xm2 - 40.17_doubp*xm2*xm3 &
          & - 1.571_doubp*xm3*xm3
     b_hcl = -7437.0_doubp - 8327.8_doubp*xm2 + 1300.9_doubp*xm3 + 1087.2_doubp*sqrt(xm2) &
          &  - 242.71_doubp*sqrt(xm3) + 18749.0_doubp*xm2*xm2 + &
          &  18500.0_doubp*xm2*xm3 + 5632.0_doubp*xm3*xm3
     phcl = exp(a_hcl +b_hcl+log(xm5*(xm5+ xm2 +0.61_doubp*xm3)))
  else
     phcl =  0.0_doubp
  end if

!!$lc JusT for iniTial values (Tämä on Tässä, kun LUO:n sysTeemissä on se väli,jolla se
!!$lc Toimii. ITse asiassa siinä paperissa oli  0.0 < xm3 + xm2 < 0.7. Koska mulla
!!$c on alussa xm3 = 0.8 ja xm2 = 0.0 -> seuraava raja:
	
  if ( xm3 > 0.8_doubp ) xm3 = 0.8_doubp
  ph = exp(23.306_doubp-4.5261_doubp*xm2-5.3465_doubp*xm3+(xm2+1.4408_doubp*xm3) &
       &  * (7.451_doubp*xm2+12.0_doubp*xm3)-(xm2+1.4408_doubp*xm3)*(xm2+1.4408_doubp*xm3) &
       &  * (4.0_doubp*xm2+8.19_doubp*xm3)+(1.0_doubp/T)*(-5814.0_doubp+1033.0_doubp*xm2+928.9_doubp*xm3 &
       &  - (xm2 + 1.4408_doubp*xm3)*(2309.0_doubp*xm2 + 1876.7_doubp*xm3)))
!!$c       write(*,*) ph/exp(23.306+(1.0/T)*(-5814.0))
!!$c        ph = min(ph,exp(23.306+(1.0/T)*(-5814.0)))
end subroutine s_luo
!!$C
!!$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$C
subroutine rho_Mir(xm22,xm33,T,rho)
	        
!!$c       Based on article, E.Martin C. George and P. Mirabel: Densities
!!$c       and Surface Tensions of H2SO4/HNO3/H2O Solutions. Geophysical
!!$c       Research Letters 27, no.2 197-200, 2000
!!$c       xm2 = HNO3 weight fraction, xm3 = H2SO4 weight fraction, T =
!!$c       temperature 
	
  real :: xm2,xm3,T,rho,y,x,rho293,rho253,kk
  real,dimension(2) :: hh
  
  IF (xm22+xm33 > 0.8) THEN
     xm3 = .8/(1.+xm22/xm33)
     xm2 = 0.8-xm3
  else
     xm2 = xm22
     xm3 = xm33
  END IF
  y = 100.0*(xm2)
  x = 100.0*(xm3)
  x2 = x**2
  x3 = x**3
  y2 = y**2
  y3 = y**3

!!$c       CCCCCCCCCCCCC  T = 293;   CCCCCCCCCCCCCC

  rho293 =   0.9982+ 5.8487e-3*y -2.3873e-5*y2 + &
       &  7.9815e-7*y3 + 7.9119e-3*x + 2.9369e-5*x*y + 1.8949e-6*x*y2 &
       &  -3.6905e-8*x*y3 -7.6431e-5*x2 &
       &	+ 2.8093e-6*x2*y  -6.4247e-8*x2*y2  + 2.2885e-6*x3 &
       &	-2.1422e-8*x3*y -1.4651e-8*x**4
!!$c  
!!$c       CCCCCCCCCCCCCCCC  T = 253;   CCCCCCCCCCCCCCCCC
!!$c    
  rho253 = 1.0015 + 9.7509e-3*y -1.834e-4*y2+ 2.9113e-6*y3 + &
       &  9.6589e-3*x-3.9433e-5*x*y +  3.8149e-6*x*y2 - 6.3144e-8*x*y3 &
       &	-1.1562e-4*x2 + 4.3442e-6*x2*y - 7.0749e-8*x2*y2 + &
       &  2.6848e-6*x3 - 3.6871e-8*x3*y - 1.6015e-8*x**4
  
!!$c  
!!$cc      CCCCCCCCCCCCCCCC  T = T  CCCCCCCCCCCCCCCCCC
!!$c  
  kk = (rho293-rho253)/40.0
  rho = kk*(T-253.0)+rho253
  rho = 1.0e3*rho
END subroutine rho_Mir
