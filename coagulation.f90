SUBROUTINE COAGULATION(Y,n,x_v,x)

!!C     
!!C *****************************************************************
!!C *   SUBROUTINE COAGULATION CALCULATES COAGULATION RATE:         *
!!C *                                                              *
!!C *                                                               *
!!C *****************************************************************
!!C     
  USE headfile
  use precisi
  use aikaa_mittaavat_parametrit
  implicit none
  INTEGER, intent(inout)  :: n
  REAL(doubp), intent(inout), DIMENSION(n) :: Y
  REAL(doubp), intent(inout)  :: x_v,x
  INTEGER :: L3,II,I_P, ISpec,JJ,JJ2,I,Imode1,Imode2,Idrop,Icontact,I_ice,J_coag,I_COAG
  REAL(doubp), DIMENSION(NABins + NABins * NAerT)  :: dCdx,dCdx2
  REAL(doubp) :: CI_P, frac1,frac2, frac_ice,B_KoagKern, B_Koag,dt,dt2,dt_old
  REAL(doubp), DIMENSION(NY1 + NY1 * NAerT)  :: dC_Y1dx,dC_Y1dx2
  REAL(doubp), DIMENSION(NY2 + NY2 * NAerT)  :: dC_Y2dx,dC_Y2dx2
  REAL(doubp), DIMENSION(NY3 + NY3 * NAerT)  :: dC_Y3dx,dC_Y3dx2
  REAL(doubp), DIMENSION(NY4 + NY4 * NAerT)  :: dC_Y4dx,dC_Y4dx2
  REAL(doubp), DIMENSION(NaBins)  :: Rp_coag
!  INTEGER,  DIMENSION(NaLiquids+NaSolids) :: indexL3,indexII,indexIP
  !      IF(x_v.NE.x.OR.X.LT.1.1e-10) call find_sop(x)
!
!
! At the moment coagulation is solved with "middle step" method, e.g. slightly 
! better than simple Euler. To estimate derivate at the middle of internal time step
! the first step taken is very short.
!
! To make routine faster the coagulation Kernel as well as target bin is 
! calculated only at the beginning of the time step, so accuracy decreases 
! with longer time steps.
!
! Because coagulation is outside solver, the mean size during previous time step
! taken to solve condensation is used to estimate the kernel. Or should be....
!
! In the presence of large particles/drops internal time step is decreased
!
! NOTE: At the moment value of coagulation kernel is limited to maximum 1.0.
!
! NOTE: There is no coagulation for C_Y3 and C_Y4. Sorry about that :)
!       It should work with copy-paste, but I didn't have time to do that
! 
!
  
!  write(*,*) C_Y2(1:NY2),'here'
  eetta = 1.8325d-5*(416.16_doubp/(T+120.0_doubp))*(T/296.16_doubp)**1.5_doubp ![kg m-1 s-1]
  GASPEED2       = SQRT(8._doubp * R * T/(Pi *0.029_doubp))         ![m s-1]
  dt2=x-x_v
  J_coag=2
  J_coag=max(2,int(dt2/1.0)+1)
!  write(*,*) C_Y2(1:NY2)
!  IF(NFABins2 > 0 .OR. N_CL2 > 0 ) THEN
  DO II=1,NaBins
     IFCOAGL3(II)=1
     IF(C(II)<1.0e-5) IFCOAGL3(II)=0
  END DO

  DO II=1,NY1
     IF(R_ice(II)>4.0d-5 .AND. C_Y1(II)>1.0e-12) J_coag=int(dt2/.5)+1
     IF(R_ice(II)>8.0d-5 .AND. C_Y1(II)>1.0e-12) J_coag=int(dt2/.2)+1
     JJ=II+NABINS
     IFCOAGL3(JJ)=1
     IF(C_Y1(II)<1.0e-6) IFCOAGL3(JJ)=0
  END DO
  DO II=1,NY2
     IF(R_cl(II)>4.0d-5 .AND. C_Y2(II)>1.0e-12) J_coag=max(int(dt2/.5)+1, J_coag)
     IF(R_cl(II)>1.0d-4 .AND. C_Y2(II)>1.0e-12) J_coag=max(int(dt2/.1)+1,J_coag)
!     write(*,*) R_cl(II),C_Y2(II),int(dt2/.01)
     JJ=II+NABINS+NY1
     IFCOAGL3(JJ)=1
     IF(C_Y2(II)<1.0e-6) IFCOAGL3(JJ)=0
     IF(II>65) IFCOAGL3(JJ)=0
  END DO
!  ENDIF

  dt=1.0e-4_doubp
  CALL FIND_BIN

  DO I_COAG=1,J_COAG
!     CALL FIND_BIN
     dmdt(1:NY1)=0.0_doubp
     dCdx(1:NABins + NABins * NAerT)=0.0_doubp
     II=NY1 + NY1 * NAerT
     dC_Y1dx(1:II)=0.0_doubp
     II=NY2 + NY2 * NAerT
     dC_Y2dx(1:II)=0.0_doubp
!     II=NY3 + NY3 * NAerT
!     dC_Y3dx(1:II)=0.0_doubp
!     II=NY4 + NY4 * NAerT
!     dC_Y4dx(1:II)=0.0_doubp

     Rp_coag=(Rp(1:NaBins)+Rp_old(1:NaBins))/2.0_doubp

!  write(*,*) ' !!  COAGULATION SHOULD BE FIXXED SO THAT AEROSOL SIZE IS NOT THE ONE AT THE END OT TIME STEP BUT ONE IN THE MIDDLE !!'
!  write(*,*) ' !!  COAGULATION KERNELS FOR ICE PARTICLES NEED TO BE WRITTEN  !!'
!
!
! Find bins in which collision product is placed
!
 
     
     MIN_COAG=1
     MAX_COAG=NABins-1


! First calculate coagulation between aerosols


     DO  L3 = MIN_COAG,MAX_COAG
        
        
!        IFCOAGL3(L3)=1
!        IF(C(L3)*C(L3+NABINS)<1.0e-26_doubp) IFCOAGL3(L3)=0
!        IF(C(L3)<1.0e-6_doubp) IFCOAGL3(L3)=0
        DO  II = L3,MAX_COAG
!           IFCOAGL3(II)=1
           
!           IF(C(II)*C(II+NABINS)<1.0e-26_doubp) IFCOAGL3(II)=0    ! If concentration too small then no coagulation
!           IF(C(II)<1.0e-6_doubp) IFCOAGL3(II)=0
           
           IF(IFCOAGL3(L3)*IFCOAGL3(II).EQ.1) THEN  

              IF (Rp(L3) < rmin_cloud .OR. Rp(II) < rmin_cloud) THEN         ! This is to decide if drops should be moved to cloud droplets

                 frac_ice = 0.0_doubp
            
                 IF(IF_FREEZ_CONTACT == 1 .AND. & 
                      & (IModeBin(L3) == I_MODE_CONTACT .OR. IModeBin(II) == I_MODE_CONTACT) &
                      & .AND. (Rp(L3) > rmin_cloud .OR. Rp(II) > rmin_cloud)) THEN
                    Idrop=II
                    Icontact=L3
                    IF(Rp(L3)>Rp(II)) THEN
                       Idrop=L3
                       Icontact=II
                    ENDIF
                    CALL CONTACT_FREEZING(Idrop,Icontact,frac_ice,I_ice)
                 ENDIF
!C
!C--- Result of coagulation is divided to two bins I_P and I_P+1
!C--- frac1 goes to I_P and frac2 to I_P+1
!C
                 CI_P=min(float(MAX_COAG),(G_ind(L3,II)))
                 I_P = int(CI_P)
                 frac1 = (1.-G_ind(L3,II)+float(I_P))*(1.0-frac_ice)
                 frac2 = (G_ind(L3,II)-float(I_P))*(1.0-frac_ice)
!C
!C--- Calculate the new coagulation kernel
!C            
           
                 IF(I_COAG==1) B_coagul(L3,II) =  B_KoagKern(min(Rp_coag(L3), &
                      & Rp_coag(II)),max(Rp_coag(L3),Rp_coag(II)))
                 B_Koag =   B_coagul(L3,II)
!C
!C--- Check that kernel is reasonable
!C           
                 IF(B_Koag.LT.0) then
                    write(*,*) ' COAGULATION KERNEL IS WRONG1'
                    pause
                 ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bins I_P and I_P+1
!C--- If I_P=NRAD, THEN all goes to I_P
!C           
!c            IF(L3.NE.II) THEN
                 dCdx(L3) = dCdx(L3)-C(L3)*B_Koag*C(II)
                 dCdx(II) = dCdx(II)-C(L3)*B_Koag*C(II)
!C     
!c
                 IF(I_P<1) write(*,*) L3,II,I_P
                 dCdx(I_P) = dCdx(I_P)+C(L3)*B_Koag*C(II)*frac1
                 dCdx(I_P+1) = dCdx(I_P+1)+C(L3)*B_Koag*C(II)*frac2
                 IF(frac_ice > 0.0)  dC_Y1dx(I_ice) = dC_Y1dx(I_ice)+C(L3)*B_Koag*C(II)*frac_ice
!              IF(frac_ice>0) write(*,*) frac1,frac2,frac_ice
!C
!C--- Move also all substances just like the number
!c
                 coag1:   DO ISpec=1,NALiquids + NASolids

!                    IF(C(InC(ISpec,L3))<1.0d-50 .AND. C(InC(ISpec,II))<1.0d-50) cycle coag1

                    dCdx(InC(ISpec,L3)) = dCdx(InC(ISpec,L3))-B_Koag*C(II) &
                         &   *C(InC(ISpec,L3))
!                    write(*,*) C(InC(ISpec,L3))
                    
                    dCdx(InC(ISpec,II)) = dCdx(InC(ISpec,II))-C(L3)*B_Koag &
                         &   *C(InC(ISpec,II))
                    
                    dCdx(InC(ISpec,I_P)) = dCdx(InC(ISpec,I_P))+(B_Koag*C(II) &
                         &   *C(InC(ISpec,L3))+C(L3)*B_Koag &
                         &   *C(InC(ISpec,II)))*frac1
                    dCdx(InC(ISpec,I_P+1)) = dCdx(InC(ISpec,I_P+1)) &
                         &   +(B_Koag*C(II) &
                         &   *C(InC(ISpec,L3))+C(L3)*B_Koag &
                         &   *C(InC(ISpec,II)))*frac2
                    IF(frac_ice > 0.0) dC_Y1dx(InCi(ISpec,I_ice)) = dC_Y1dx(InCi(ISpec,I_ice)) &
                         &  +(C(L3)*B_Koag*C(II)/ &
                         &   C(L3)*C(InC(ISpec,L3))+C(L3)*B_Koag*C(II)/ &
                         &   C(II)*C(InC(ISpec,II)))*frac_ice
                    
                    
                    
                 END DO coag1  ! ISpec
                 
!                 indexL3=InC(1:(NaLiquids+NaSolids),L3)
!                 dCdx(indexL3)= dCdx(indexL3)-B_Koag*C(II) *C(indexL3)
!                 indexII=InC(1:(NaLiquids+NaSolids),II)
!                 dCdx(indexII)= dCdx(indexII)- B_Koag*C(L3) *C(indexII)
!                 indexIP=InC(1:(NaLiquids+NaSolids),I_P)
!                 dCdx(indexIP)= dCdx(indexIP)+(B_Koag*C(II)*C(indexL3)+C(L3)*B_Koag*C(indexII))*frac1
!                 dCdx(indexIP+1)= dCdx(indexIP+1)+(B_Koag*C(II)*C(indexL3)+C(L3)*B_Koag*C(indexII))*frac2
!                 write(*,*) C(InC(1,L3):InC(1,L3)+((NALiquids + NASolids-1)*NaBins):NaBins)
!                 pause
              ELSE

             
!C
!C--- Result of coagulation is placed to cloud droplet bin I_P
!C
!C
                 CI_P=min(float(NY2),(G_ind(L3,II)))
                 I_P = int(CI_P)

!C
!C--- Calculate the new coagulation kernel
!C            
           
                 IF(I_COAG==1) B_coagul(L3,II) =  B_KoagKern(min(Rp_coag(L3), &
                      & Rp_coag(II)),max(Rp_coag(L3),Rp_coag(II)))
                 B_Koag =   B_coagul(L3,II)
!C
!C--- Check that kernel is reasonable
!C           
                 IF(B_Koag.LT.0) then
                    write(*,*) ' COAGULATION KERNEL IS WRONG21'
                    pause
                 ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P 
!C
!C           
                 dCdx(L3) = dCdx(L3)-C(L3)*B_Koag*C(II)
                 dCdx(II) = dCdx(II)-C(L3)*B_Koag*C(II)
!C     
!c
                 IF(I_P<1) write(*,*) L3,II,I_P,'5+'
                 dC_Y2dx(I_P) = dC_Y2dx(I_P)+C(L3)*B_Koag*C(II)
!              IF(I_P==40) write(*,*) dC_Y2dx(I_P),C(L3)*B_Koag*C(II),'2'
!C
!C--- Move also all substances just like the number
!c
                 coag2: DO ISpec=1,NALiquids + NASolids
!                 IF(C(InC(ISpec,L3))<1.e-50.AND.C(InC(ISpec,II))<1.0e-50) cycle coag2
 
                    dCdx(InC(ISpec,L3)) = dCdx(InC(ISpec,L3))-B_Koag*C(II) &
                         &   *C(InC(ISpec,L3))
                    
                    dCdx(InC(ISpec,II)) = dCdx(InC(ISpec,II))-C(L3)*B_Koag &
                         &   *C(InC(ISpec,II))
                    
                    dC_Y2dx(InCl(ISpec,I_P)) = dC_Y2dx(InCl(ISpec,I_P))+(B_Koag*C(II) &
                         &   *C(InC(ISpec,L3))+C(L3)*B_Koag &
                         &   *C(InC(ISpec,II)))

                 END DO coag2  ! ISpec
              
              ENDIF
           ENDIF
        END DO  ! II
     END DO  ! L3
!
! Now do the same for particle - ice nuclei collisions
!
!  GOTO 100
     DO  L3 = MIN_COAG,MAX_COAG
     
        ice: DO  II = 1,NY1
           
          
           
!           IF(C_Y1(II)<1e-8_doubp.OR.C(L3)*C(L3+NABINS)<1.0e-26_doubp) cycle ice        ! Coagulation limited to concentration higher than
          
           IF(IFCOAGL3(L3)*IFCOAGL3(NABins+II)==1) THEN
              JJ=NABins+II
!C
!C--- Result of coagulation is put to bin I_P
!C--
!C
              CI_P=min(float(NY1),(G_ind(L3,JJ)))
              I_P = int(CI_P)
!C
!C--- Calculate the new coagulation kernel
!C
           
              IF(I_COAG==1) B_coagul(L3,JJ) =  B_KoagKern(min(Rp_coag(L3), &
                   & R_ice(II)),max(Rp_coag(L3),R_ice(II)))
              B_Koag =   B_coagul(L3,JJ)
!C
!C--- Check that kernel is reasonable
!C           
              
              IF(B_Koag.LT.0) then
                 write(*,*) II, R_ice(II)
                 write(*,*) ' COAGULATION KERNEL IS WRONG22'
                 write(*,*) C_Y1(II),C_Y1(II+NY1),C_Y1(II+2*NY1),&
                      & C_Y1(II+3*NY1),C_Y1(II+4*NY1),C_Y1(II+5*NY1)
                 
                 pause
              ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P 
!C           
              dCdx(L3) = dCdx(L3)-C(L3)*B_Koag*C_Y1(II)
              dC_Y1dx(II) = dC_Y1dx(II)-C(L3)*B_Koag*C_Y1(II)
!C     
!c
              dC_Y1dx(I_P) = dC_Y1dx(I_P)+C(L3)*B_Koag*C_Y1(II)
              dmdt(I_P) = C(L3)*B_Koag*C_Y1(II)/ &
                   &   C(L3)*C(InC(1,L3))*wtmol(IH2OL)
!C
!C--- Move also all substances just like the number
!c
              DO ISpec=1,NALiquids + NASolids
!         
                 dCdx(InC(ISpec,L3)) = dCdx(InC(ISpec,L3))-C(L3)*B_Koag*C_Y1(II)/ &
                      &   C(L3)*C(InC(ISpec,L3))
                 
                 dC_Y1dx(InCi(ISpec,II)) = dC_Y1dx(InCi(ISpec,II))-C(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II))
                 
                 dC_Y1dx(InCi(ISpec,I_P)) = dC_Y1dx(InCi(ISpec,I_P))+(C(L3)*B_Koag*C_Y1(II)/ &
                      &   C(L3)*C(InC(ISpec,L3))+C(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II)))

              END DO  ! ISpec
!              IF(Rp(L3) > r_min_hm) CALL halletmossop(dCdx,dC_Y1dx,dC_Y2dx,L3,I_P,C(L3)*B_Koag*C_Y1(II))
           ENDIF
        END DO ice
     END DO  ! L3
!C
!C
!C Now do the same for ICE - ICE nuclei collisions
!C 
!C
     ice2:DO  L3 = 1,NY1

        IF(IFCOAGL3(L3+NaBins)==0) cycle ice2

        DO  II = L3,NY1
           
           JJ=II+NABINS
           JJ2=L3+NABINS

!           IF(C_Y1(II)<1e-8_doubp .OR.C_Y1(L3)<1e-8_doubp ) cycle ice2   ! Coagulation limited to concentration higher than

           IF(IFCOAGL3(JJ2)*IFCOAGL3(JJ)==1) THEN
!C
!C--- Result of coagulation is but to bin I_P
!C
              CI_P=min(float(NY1),(G_ind(JJ2,JJ)))
              I_P = int(CI_P)
!C
!C--- Calculate the new coagulation kernel
!C            
           
               IF(I_COAG==1) B_coagul(JJ2,JJ) =  B_KoagKern(min(R_ice(L3), &
                   & R_ice(II)),max(R_ice(L3),R_ice(II)))
              B_Koag =   B_coagul(JJ2,JJ)
!C
!C--- Check that kernel is reasonable
!C           
              IF(B_Koag.LT.0) then
                 write(*,*) ' COAGULATION KERNEL IS WRONG3'
                 pause
              ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P
!C           
!c            IF(L3.NE.II) THEN
              dC_Y1dx(L3) = dC_Y1dx(L3)-C_Y1(L3)*B_Koag*C_Y1(II)
              dC_Y1dx(II) = dC_Y1dx(II)-C_Y1(L3)*B_Koag*C_Y1(II)
!C     
!c
              dC_Y1dx(I_P) = dC_Y1dx(I_P)+C_Y1(L3)*B_Koag*C_Y1(II)

!C
!C--- Move also all substances just like the number
!c
              DO ISpec=1,NALiquids + NASolids
                 
                 dC_Y1dx(InCi(ISpec,L3)) = dC_Y1dx(InCi(ISpec,L3))-C_Y1(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(L3)*C_Y1(InCi(ISpec,L3))
                 
                 dC_Y1dx(InCi(ISpec,II)) = dC_Y1dx(InCi(ISpec,II))-C_Y1(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II))
                 
                 dC_Y1dx(InCi(ISpec,I_P)) = dC_Y1dx(InCi(ISpec,I_P))+(C_Y1(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(L3)*C_Y1(InCi(ISpec,L3))+C_Y1(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II)))
              
              END DO   ! ISpec
!           IF(L3==41.OR.II==41.OR.I_P==41) write(*,*) L3,II,I_P,dC_Y1dx(41),B_Koag,R_ice(L3),R_ice(II)
           ENDIF
        END DO 
     END DO ice2 ! L3
!C
!C
!C coagulation for aerosol particle - cloud droplet collisions
!C 
!C
100  DO  L3 = MIN_COAG,MAX_COAG
        
        liq: DO  II = 1,NY2

           JJ=II+NABINS+NY1
           
!           IF(C_Y2(II)<1e-8_doubp) cycle liq        ! Coagulation limited to concentration higher than
           
           IF(IFCOAGL3(L3).NE.0.AND.IFCOAGL3(JJ).NE.0) THEN
              
              frac_ice = 0.0_doubp
              
              IF(IF_FREEZ_CONTACT == 1 .AND. & 
                   & (IModeBin(L3) == I_MODE_CONTACT) &
                   & .AND. Rp(L3) > rmin_cloud) THEN
                 Idrop=II
                 Icontact=L3
                 IF(Rp(L3)>Rp(II)) THEN
                    Idrop=L3
                    Icontact=II
                 ENDIF
                 CALL CONTACT_FREEZING(Idrop,Icontact,frac_ice,I_ice)
              ENDIF
           


!C
!C--- Result of coagulation is put to bin I_P
!C--
!C
              CI_P=min(float(NY2),(G_ind(L3,JJ)))
              I_P = int(CI_P)
              
              frac1 = 1.0-frac_ice
!C
!C--- Calculate the new coagulation kernel
!C
              
               IF(I_COAG==1) THEN
                  B_coagul(L3,JJ) =  B_KoagKern(min(Rp_coag(L3), &
                       & R_cl(II)),max(Rp_coag(L3),R_cl(II)))
!                  IF(B_coagul(L3,JJ)>1.d-4) write(*,*) B_coagul(L3,JJ)*C(L3), L3,JJ,I_P
                  B_coagul(L3,JJ) = min(B_coagul(L3,JJ),1./(sum(C(MIN_COAG:MAX_COAG))*C_Y2(II)))   !For numerical stability
               ENDIF
              B_Koag =   B_coagul(L3,JJ)
            
             !C
!C--- Check that kernel is reasonable
!C           
              
              IF(B_Koag.LT.0) then
                 write(*,*) II, R_cl(II)
                 write(*,*) ' COAGULATION KERNEL IS WRONG23'
                 write(*,*) C_Y2(II),C_Y2(II+NY2),C_Y2(II+2*NY2),&
                      & C_Y2(II+3*NY2),C_Y2(II+4*NY2),C_Y2(II+5*NY2)
                 
                 pause
              ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P 
!C           

              dCdx(L3) = dCdx(L3)-C(L3)*B_Koag*C_Y2(II)
              dC_Y2dx(II) = dC_Y2dx(II)-C(L3)*B_Koag*C_Y2(II)
!              if(B_Koag>1e-1.AND.I_COAG==1)  write(*,*) B_coagul(L3,NABINS+NY1+1:JJ)
!              if(B_Koag>2e-1.AND.I_COAG==1.AND.L3.EQ.10)  write(*,*) sum(B_coagul(L3,NABINS+NY1+1:JJ)*C_Y2(1:II)),'***'
!              if(B_Koag>1e-1.AND.I_COAG==1) pause
!C     
!c
              dC_Y2dx(I_P) = dC_Y2dx(I_P)+C(L3)*B_Koag*C_Y2(II)*frac1
              IF(frac_ice > 0.0)  dC_Y1dx(I_ice) = dC_Y1dx(I_ice)+C(L3)*B_Koag*C(II)*frac_ice
!           IF(II==40) write(*,*) dC_Y2dx(II),B_Koag,C_Y2(II),I_P,II,L3
!           IF(I_P==40) write(*,*) dC_Y2dx(I_P),C(L3)*B_Koag*C_Y2(II),'1',II,L3
!C
!C--- Move also all substances just like the number
!c
              DO ISpec=1,NALiquids + NASolids
                 
                 dCdx(InC(ISpec,L3)) = dCdx(InC(ISpec,L3))-C(L3)*B_Koag*C_Y2(II)/ &
                      &   C(L3)*C(InC(ISpec,L3))
                 
                 dC_Y2dx(InCl(ISpec,II)) = dC_Y2dx(InCl(ISpec,II))-C(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(II)*C_Y2(InCl(ISpec,II))
                 
                 dC_Y2dx(InCl(ISpec,I_P)) = dC_Y2dx(InCl(ISpec,I_P))+(C(L3)*B_Koag*C_Y2(II)/ &
                      &   C(L3)*C(InC(ISpec,L3))+C(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(II)*C_Y2(InCl(ISpec,II)))
              END DO   ! ISpec
           ENDIF
        END DO liq
     END DO  ! L3
!C
!C
!C Now do the same for DROP - DROP collisions
!C
!C

     liq2:DO  L3 = 1,NY2

        IF(IFCOAGL3(L3+NaBins+NY1)==0) cycle liq2

        DO  II = L3,NY2
           
           JJ=II+NABINS+NY1
           JJ2=L3+NABINS+NY1
!           IFCOAGL3(JJ)=1
!           IFCOAGL3(JJ2)=1
           
!           IF(C_Y2(II)<1e-8_doubp .OR.C_Y2(L3)<1e-8_doubp ) cycle liq2   ! Coagulation limited to concentration higher than

           IF(IFCOAGL3(JJ2)*IFCOAGL3(JJ)==1) THEN
              
!C
!C--- Result of coagulation is but to bin I_P
!C
              CI_P=min(float(NY2),(G_ind(JJ2,JJ)))
              I_P = nint(CI_P)
!              frac1 = (1.-G_ind(JJ2,JJ)+float(I_P))
!              frac2 = (G_ind(JJ2,JJ)-float(I_P))
!              write(*,*) I_P,G_ind(JJ2,JJ),frac1,frac2
!C
!C--- Calculate the new coagulation kernel
!C            
           
               IF(I_COAG==1) THEN 
                  B_coagul(JJ2,JJ) =  B_KoagKern(min(R_cl(L3), &
                       & R_cl(II)),max(R_cl(L3),R_cl(II)))
!                   B_coagul(JJ2,JJ) = min( B_coagul(JJ2,JJ)
!                  IF(B_coagul(JJ2,JJ)>1.d-4) write(*,*) B_coagul(JJ2,JJ)*C_Y2(L3),JJ2,JJ,I_P
                  B_coagul(JJ2,JJ) = min(B_coagul(JJ2,JJ),1.0/(sum(C_Y2(II:NY2))**2))   !For numerical stability
               ENDIF
              B_Koag =   B_coagul(JJ2,JJ)
!              write(*,*) L3,II,I_P,B_Koag,R_cl(L3),R_cl(II)
!C
!C--- Check that kernel is reasonable
!C           
              IF(B_Koag.LT.0) then
                 write(*,*) ' COAGULATION KERNEL IS WRONG3'
                 pause
              ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P
!C           
!c            IF(L3.NE.II) THEN
              dC_Y2dx(L3) = dC_Y2dx(L3)-C_Y2(L3)*B_Koag*C_Y2(II)
              dC_Y2dx(II) = dC_Y2dx(II)-C_Y2(L3)*B_Koag*C_Y2(II)
!C     
!c
              dC_Y2dx(I_P) = dC_Y2dx(I_P)+C_Y2(L3)*B_Koag*C_Y2(II)
          
!           IF(I_P==40) write(*,*) dC_Y2dx(I_P),C_Y2(L3)*B_Koag*C_Y2(II)
!C
!C--- Move also all substances just like the number
!c
              DO ISpec=1,NALiquids + NASolids
                 
                 dC_Y2dx(InCl(ISpec,L3)) = dC_Y2dx(InCl(ISpec,L3))-C_Y2(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(L3)*C_Y2(InCl(ISpec,L3))
                 
                 dC_Y2dx(InCl(ISpec,II)) = dC_Y2dx(InCl(ISpec,II))-C_Y2(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(II)*C_Y2(InCl(ISpec,II))
                 
                 dC_Y2dx(InCl(ISpec,I_P)) = dC_Y2dx(InCl(ISpec,I_P))+(C_Y2(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(L3)*C_Y2(InCl(ISpec,L3))+C_Y2(L3)*B_Koag*C_Y2(II)/ &
                      &   C_Y2(II)*C_Y2(InCl(ISpec,II)))
              
              END DO   ! ISpec
!              write(*,*) L3,II,I_P,C_Y2(I_P)
!           IF(I_BOX>10.aND.I_P==30) write(*,*) L3,II,dC_Y2dx(I_P)/C_Y2(I_P),dC_Y2dx(I_P+NY2)/C_Y2(I_P+NY2)
           ENDIF
        END DO
     END DO  liq2 ! L3
!C
!C
!C Finally DROP - ICE collisions
!C
!C

     liqice:DO  L3 = 1,NY2

        IF(IFCOAGL3(L3+NaBins+NY1)==0) cycle liqice

        DO  II = 1,NY1

           JJ=II+NABINS
           JJ2=L3+NABINS+NY1
!           IFCOAGL3(JJ)=1
!           IFCOAGL3(JJ2)=1
           
!           IF(C_Y1(II)<1e-8 .OR.C_Y2(L3)<1e-8 ) cycle liqice   ! Coagulation limited to concentration higher than
           
           IF(IFCOAGL3(JJ2)*IFCOAGL3(JJ)==1) THEN
!C
!C--- Result of coagulation is but to bin I_P
!C
              CI_P=min(float(NY1),(G_ind(JJ2,JJ)))
              I_P = int(CI_P)
!C
!C--- Calculate the new coagulation kernel
!C            
           
               IF(I_COAG==1) B_coagul(JJ2,JJ) =  B_KoagKern(min(R_cl(L3), &
                   & R_ice(II)),max(R_cl(L3),R_ice(II)))
              B_Koag =   B_coagul(JJ2,JJ)
!C
!C--- Check that kernel is reasonable
!C           
              IF(B_Koag.LT.0) then
                 write(*,*) ' COAGULATION KERNEL IS WRONG3'
                 pause
              ENDIF
!C
!C--- Number decreases in bins L3 and II
!C--- Number increases in bin I_P
!C           
!c            IF(L3.NE.II) THEN
              dC_Y2dx(L3) = dC_Y2dx(L3)-C_Y2(L3)*B_Koag*C_Y1(II)
              dC_Y1dx(II) = dC_Y1dx(II)-C_Y2(L3)*B_Koag*C_Y1(II)
!C     
!c
              dC_Y1dx(I_P) = dC_Y1dx(I_P)+C_Y2(L3)*B_Koag*C_Y1(II)
              dmdt(I_P) = C_Y2(L3)*B_Koag*C_Y1(II)/ &
                   &   C_Y2(L3)*C_Y2(InCl(1,L3))*wtmol(IH2OL)
!C
!C--- Move also all substances just like the number
!c
              DO ISpec=1,NALiquids + NASolids
              
                 dC_Y2dx(InCl(ISpec,L3)) = dC_Y2dx(InCl(ISpec,L3))-C_Y2(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y2(L3)*C_Y2(InCl(ISpec,L3))
                 
                 dC_Y1dx(InCi(ISpec,II)) = dC_Y1dx(InCi(ISpec,II))-C_Y2(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II))
                 
                 dC_Y1dx(InCi(ISpec,I_P)) = dC_Y1dx(InCi(ISpec,I_P))+(C_Y2(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y2(L3)*C_Y2(InCl(ISpec,L3))+C_Y2(L3)*B_Koag*C_Y1(II)/ &
                      &   C_Y1(II)*C_Y1(InCi(ISpec,II)))
              
              END DO   ! ISpec
!              IF(R_CL(L3) > r_min_hm) CALL halletmossop(dCdx,dC_Y1dx,dC_Y2dx,L3+NABins,I_P,C_Y2(L3)*B_Koag*C_Y1(II))
           ENDIF
        END DO
     END DO  liqice ! L3

!C
!C Make changes due to cagulation to concentrations
!
!  DO I=1,NY1*2
!     IF(abs(dC_Y1dx(I)/C_Y1(I)).GT.1e-10) write(*,*) dC_Y1dx(I)/C_Y1(I),C_Y1(I),I
!  END DO
!  
!   DO I=1,NaBins*2
!     IF(abs(dCdx(I)/C(I)).GT.1e-10) write(*,*) dCdx(I)/C(I),C(I),I,'2'
!  END DO
     
     DO I=1,NY2-1
        IF  (C_Y2(I) < -dC_Y2dx(I)*dt.OR.C_Y2(I+NY2) < -dC_Y2dx(I+NY2)*dt) THEN

!           write(*,*) C_Y2(I),dC_Y2dx(I),dC_Y2dx2(I),(dC_Y2dx(I)-dC_Y2dx2(I))/dt_old
!           write(*,*) C_Y2(I+NY2),dC_Y2dx(I+NY2)*dt,dC_Y2dx2(I+NY2)*dt,dt_old
!           write(*,*) I,I_coag,R_cl(I),R_Y2(I+1)
           
!           pause
!           dC_Y2dx(I+1)=dC_Y2dx(I+1)+(dC_Y2dx(I)+C_Y2(I)*0.99/dt)
!           dC_Y2dx(I) = -C_Y2(I)*0.99/dt
!           DO ISpec=1,NALiquids + NASolids
!              dC_Y2dx(InCl(ISpec,I+1)) = dC_Y2dx(InCl(ISpec,I+1))+(dC_Y2dx(InCl(ISpec,I))+C_Y2(InCl(ISpec,I))*0.99/dt)
!              dC_Y2dx(InCl(ISpec,I)) = -C_Y2(InCl(ISpec,I))*0.99/dt
!           END DO
!        write(*,*) I,dC_Y2dx(I+1),C_Y2(I+1),dC_Y2dx(I+1+NY2),C_Y2(I+1+NY2)
!        pause
        ENDIF
     END DO

     DO I=1,NY1-1
        IF  (C_Y1(I) < -dC_Y1dx(I)*dt) THEN

!           write(*,*) C_Y1(I),dC_Y1dx(I)*dt,C_Y1(I+NY1),dC_Y1dx(I+NY1)*dt,I,dt,R_ice(I),R_Y1(I+1)
!           pause
!           dC_Y1dx(I+1)=dC_Y1dx(I+1)+(dC_Y1dx(I)+C_Y1(I)*0.99/dt)
!           dC_Y1dx(I) = -C_Y1(I)*0.99/dt
!           DO ISpec=1,NALiquids + NASolids
!              dC_Y1dx(InCi(ISpec,I+1)) = dC_Y1dx(InCi(ISpec,I+1))+(dC_Y1dx(InCi(ISpec,I))+C_Y1(InCi(ISpec,I))*0.99/dt)
!              dC_Y1dx(InCi(ISpec,I)) = -C_Y1(InCi(ISpec,I))*0.99/dt
!           END DO
        ENDIF
     END DO

!  
!  Chance in concentration due to coagulation dC/dt is simply multiplied by dt to get new concentrations after time step
!
!     write(*,*)  C_Y2(38),dC_Y2dx(38),I_coag

!     IF(x>6000) THEN
!        write(*,*) dC_Y2dx(20:50),'here',I_COAG
!        write(*,*) 
!        write(*,*) dC_Y2dx(20:50) + .5*dt*(dC_Y2dx(20:50)-dC_Y2dx2(20:50))/dt_old
!        pause
!     ENDIF
     IF(I_COAG>0) THEN         !!!!!!!!!!CHECK!!!!!!!!!!
        C(1:NABINS*(NALiquids + NASolids+1))=C(1:NABINS*(NALiquids + NASolids+1))+ &
             &  dCdx(1:NABINS*(NALiquids + NASolids+1))*dt
        C_Y1(1:NY1*(NALiquids + NASolids+1))=C_Y1(1:NY1*(NALiquids + NASolids+1))+ & 
             &  dC_Y1dx(1:NY1*(NALiquids + NASolids+1))*dt
        C_Y2(1:NY2*(NALiquids + NASolids+1))=C_Y2(1:NY2*(NALiquids + NASolids+1))+ &
             &  dC_Y2dx(1:NY2*(NALiquids + NASolids+1))*dt
        
        dCdx2(1:NABINS*(NALiquids + NASolids+1))=dCdx(1:NABINS*(NALiquids + NASolids+1))
        dC_Y1dx2(1:NY1*(NALiquids + NASolids+1))= dC_Y1dx(1:NY1*(NALiquids + NASolids+1))
        dC_Y2dx2(1:NY2*(NALiquids + NASolids+1))= dC_Y2dx(1:NY2*(NALiquids + NASolids+1))
        dt_old=dt
     ELSE
        C(1:NABINS*(NALiquids + NASolids+1))=C(1:NABINS*(NALiquids + NASolids+1))+ & 
             &  dCdx(1:NABINS*(NALiquids + NASolids+1))*dt &
             &  + .5*dt**2*(dCdx(1:NABINS*(NALiquids + NASolids+1))- &
             &  dCdx2(1:NABINS*(NALiquids + NASolids+1)))/dt_old

        C_Y1(1:NY1*(NALiquids + NASolids+1))=C_Y1(1:NY1*(NALiquids + NASolids+1))+ & 
             &  dC_Y1dx(1:NY1*(NALiquids + NASolids+1))*dt &
             &  + .5*dt**2*(dC_Y1dx(1:NY1*(NALiquids + NASolids+1))-  &
             &  dC_Y1dx2(1:NY1*(NALiquids + NASolids+1)))/dt_old


        C_Y2(1:NY2*(NALiquids + NASolids+1))=C_Y2(1:NY2*(NALiquids + NASolids+1))+ & 
             &  dC_Y2dx(1:NY2*(NALiquids + NASolids+1))*dt &
             &  + .5*dt**2*(dC_Y2dx(1:NY2*(NALiquids + NASolids+1))- & 
             &  dC_Y2dx2(1:NY2*(NALiquids + NASolids+1)))/dt_old

!        write(*,*) dC_Y2dx(10:NY2)/ dC_Y2dx2(10:NY2)
!        write(*,*) C_Y2(10:NY2)
!        write(*,*)
       
        dCdx2(1:NABINS*(NALiquids + NASolids+1))=dCdx(1:NABINS*(NALiquids + NASolids+1))
        dC_Y1dx2(1:NY1*(NALiquids + NASolids+1))= dC_Y1dx(1:NY1*(NALiquids + NASolids+1))
       
        dC_Y2dx2(1:NY2*(NALiquids + NASolids+1))= dC_Y2dx(1:NY2*(NALiquids + NASolids+1))
        dt_old=dt
     ENDIF
    
!  write(*,*) sum(dCdx(NABINS+1:NABINS*2))+sum(dC_Y2dx(NY2+1:NY2*2)), sum(dCdx(22:NABINS-5)),sum(dC_Y2dx(1:NY2))


!     DO II = 1,NY1
!        CALL RADIUSICE2(II)
!     END DO
!     DO II = 1,NY2
!        CALL RADIUS_Y2(II)
!     END DO
     dt=dt2/(J_coag-1)
  END DO

!  pause
  !        IF(I_BOX>10) write(*,*) R_cl(27:35),'2'
!  DO II = 1,NaBIns
!     CALL RADIUS(II)
!  END DO

  Ctot=sum(C(1:NaBins))
  Y(NODES-1) = Ctot
  
  IF (ice_mod == 1) THEN
     C_ice_tot=sum(C_Y1(1:NY1))
     y(NODEs-5)=C_ice_tot
  ENDIF
  IF (cloud_mod == 1) THEN
     c_cl_tot=sum(C_Y2(1:NY2))
     Y(NODES-6) = C_cl_tot
  ENDIF
!     
! Change also Y.  
! CHECK THAT Y FOR LIQUIDS IS CORRECTED SOMEWHERE ELSE !!!!!!!!!!!!!!!!!!! 
!
!
  DO I=1,NY1
     IF (C_Y1(I) > 0.0) THEN        
        Y(NYConcs + I)= C_Y1(NY1+I)           
     ENDIF
  END DO
   DO I=1,NY2
     IF (C_Y2(I) > 0.0) THEN        
        Y(NYConcs + NY1 + I)= C_Y2(NY2+I) 
     ENDIF
  END DO
  DO I=1,NaBins
     Y(I)= C(NaBins+I) 
  END DO

!  write(*,*) dmdt,wtmol(IH2OL)

END SUBROUTINE COAGULATION





Function B_KoagKern(r_1,r_2)

!!C *****************************************************************
!!C         FUNCTION CALCULATES COAGULATION KERNEL                  *
!!C                                                                 *
!!C *****************************************************************
!!C
!!C r_1 and r_2 should be in meters so that r_1 < r_2
!!C
!!C
!!C
  USE headfile
  use precisi
  implicit none
  REAL(doubp), intent(in) :: r_1,r_2
  REAL(doubp) :: BC_shear,  BC_inertia, BC_grav, BBB,B_KoagKern
  REAL(doubp) :: C_Kn,D1,D2,v1,v2,sigma1,sigma2
  REAL(doubp) :: velo1,velo2,Re_2,Sc_1, Rho1,eps_k,set_vel,col_ef

!  eetta = 1.8325d-5*(416.16/(T+120.0))*(T/296.16)**1.5 ![kg m-1 s-1]
!  GASPEED2       = SQRT(8. * R * T/(Pi *0.029))         ![m s-1]
  B_KoagKern=0.0
!  Rho = Rho3(1)*1000.
  TotP             = P + C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  Rho1=TotP / (10000.0_doubp  * R / WtAir * T)*1000.0_doubp        ![kg/m3] CHECK THIS   
  C_Kn = 2.*eetta/(Rho1*GASPEED2*r_1)
  D1 = Boltz*T/(6*Pi*r_1*eetta) &
       & *(1+C_Kn*(1.249+0.42*exp(-0.87/C_Kn)))  ![m s-1]
  
  C_Kn = 2._doubp*eetta/(Rho1*GASPEED2*r_2)
  D2 = Boltz*T/(6.0_doubp*Pi*r_2*eetta) &
       & *(1.0_doubp+C_Kn*(1.249_doubp+0.42_doubp*exp(-0.87_doubp/C_Kn)))

  v1 = SQRT(8.0_doubp * Boltz * T/(Pi*4.0_doubp/3.0_doubp*Pi*r_1**3*1000._doubp))
  v2 = SQRT(8.0_doubp * Boltz * T/(Pi*4.0_doubp/3.0_doubp*Pi*r_2**3*1000._doubp))  ! CHECK THIS

  sigma1 = ((2.0_doubp*r_1+8.0_doubp*D1/(Pi*v1))**3 &
       &     -(4.0_doubp*r_1**2+(8.0_doubp*D1/(Pi*v1))**2)**(3.0_doubp/2.0_doubp))/ &
       &     (6.0_doubp*r_1*8.0_doubp*D1/(Pi*v1))-2.0_doubp*r_1
    
  sigma2 = ((2.0_doubp*r_2+8.0_doubp*D2/(Pi*v2))**3 &
       &     -(4.0_doubp*r_2**2+(8.0_doubp*D2/(Pi*v2))**2)**(3.0_doubp/2.0_doubp))/ &
       &     (6.0_doubp*r_2*8.0_doubp*D2/(Pi*v2))-2.0_doubp*r_2
!C
!C     BROWNIAN COAGULATION
!C
  B_KoagKern = 4.0_doubp*Pi*(r_1+r_2)*(D1+D2)/( &
       &     (r_1+r_2)/(r_1+r_2+(sigma1**2+sigma2**2)**.5_doubp) + &
       &     4.0_doubp*(D1+D2)/((v1**2+v2**2)**.5_doubp*(r_1+r_2)))*1.0d6 

!C
!C     CONVECTIVE BROWNIAN DIFFUSION ENHANCEMENT
!C      
  velo1 = set_vel(r_1)
  velo2 = set_vel(r_2)
  Re_2 = 2.0_doubp*r_2*velo2/(eetta/Rho1)
  Sc_1 = eetta/rho1/D1
  IF(Re_2.LT.1.0_doubp) THEN
     BBB=B_KoagKern*0.45_doubp*Re_2**(1.0_doubp/3.0_doubp)*Sc_1**(1.0_doubp/3.0_doubp)
  ELSE
     BBB=B_KoagKern*0.45_doubp*Re_2**(1.0_doubp/2.0_doubp)*Sc_1**(1.0_doubp/3.0_doubp)
  ENDIF
!C

  B_KoagKern = B_KoagKern + BBB

!C
!C Gravitational settling    
!C     
  
  BC_grav = col_ef(r_1,r_2,velo1,velo2)* &
       &     Pi*(r_1+r_2)**2*(velo2-velo1)*1.0d6
  

!C Turbulent inertial motion

  eps_k = 10.0d-4  ! Turbulent energy dissibation rate
  
  BC_inertia = Pi*eps_k**(3.0_doubp/4.0_doubp)/(g*(eetta/Rho1)**(1.0_doubp/4.0_doubp)) &
       &  *(r_1+r_2)**2*(velo2-velo1)*1.0d6
  
!C
!C Turbulent shear
!C

  BC_shear =  SQRT(8.0_doubp*Pi*eps_k/(15.0_doubp*(eetta/Rho1)))* &
       &          (r_1+r_2)**3*1.0d6
  
  B_KoagKern =  B_KoagKern + SQRT(BC_grav**2+BC_inertia**2 &
       &  + BC_shear**2)


!  if(B_KoagKern>0.1) write(*,*) r_1,r_2,B_KoagKern,BC_grav
 
  B_KoagKern=min(1.0_doubp,B_KoagKern)
  

END Function B_KoagKern


!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


Function set_vel(r_1)

!C ********************************************************************
!C          Function CALCULATES GRAVITATIONAL SETTLING VELOCITY       *
!C                                                                    *
!C ********************************************************************
  USE headfile
  use precisi
  implicit none
  REAL(doubp), intent(in) :: r_1
  REAL(doubp) :: C_Kn,set_vel,Re,C_D,set_vel_old, Rho1,C_c
  INTEGER :: III

!  eetta = 1.8325d-5*(416.16/(T+120.0))*(T/296.16)**1.5 ![kg m-1 s-1]
!  GASPEED2 = SQRT(8. * R * T/(Pi *0.029))
  
  TotP             = P + C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  Rho1=TotP / (10000.0_doubp  * R / WtAir * T)*1000.0_doubp        
  
  
  C_Kn = 2.*eetta/(Rho1*GASPEED2*r_1)
  
  set_vel = 1._doubp/18._doubp*(2._doubp*r_1)**2*1000._doubp* &
       &  g*(1._doubp+C_Kn*(1.249_doubp+0.42_doubp*exp(-0.87_doubp/C_Kn)))/eetta
  
  C_c=1.0_doubp+C_Kn*(1.249_doubp+0.42_doubp*exp(-0.87_doubp/C_Kn))
  III=0
20 Re = 2.0_doubp*r_1*set_vel/(eetta/(Rho1))
  III = III+1
  IF(Re.LT.0.1) THEN
     C_D = 24.0_doubp/Re
  ELSEIF(Re.LT.1.0.AND.Re.GT.0.1) THEN
     C_D = 24.0_doubp/Re*(1.0_doubp+ min(1.0_doubp,100.0_doubp*(Re-.1_doubp))* &
          & (3.0_doubp/16.0_doubp*Re+9.0_doubp/160.0_doubp*Re**2*log(2.0_doubp*Re)))
  ELSEIF(Re.GT.1.0.AND.Re.LT.3.0) THEN
     C_D = 24.0_doubp/Re*((3.0_doubp-Re)/2.0_doubp*(1._doubp+ min(1.0_doubp,100.0_doubp*(Re-.1_doubp))* &
          & (3.0_doubp/16.0_doubp*Re+9.0_doubp/160.0_doubp*Re**2*log(2.0_doubp*Re)))+ &
          & (Re-1.0_doubp)/2.0_doubp*(1.0_doubp+0.15_doubp*Re**0.687_doubp))
  ELSEIF(Re.LT.500..AND.Re.GT.3.0) THEN
     C_D = 24.0_doubp/Re*(1.0_doubp+0.15*Re**0.687_doubp)
  ELSE
     C_D = 0.44_doubp
  ENDIF
   
  set_vel_old = set_vel
  set_vel = (8.0_doubp*9.81_doubp*r_1*C_c*1000.0_doubp/(3.0_doubp*C_D*rho1))**.5_doubp
  IF(abs(set_vel/set_vel_old-1.0_doubp).GT.1d-4.AND.III.LT.25) GOTO 20
  
END Function set_vel


Function col_ef(r_1,r_2,velo1,velo2)

!C ********************************************************************
!C    Function calculates raindrop-aerosol collision efficiency       *
!C                                                                    *
!C ********************************************************************
     
  USE headfile
  use precisi
  REAL(doubp), intent(in) ::r_1,r_2,velo1,velo2
  REAL(doubp) :: K_12, Rho1,Re_2, E_A12,E_v12,col_ef


!  eetta = 1.8325d-5*(416.16/(T+120.0))*(T/296.16)**1.5 ![kg m-1 s-1]
  
  TotP             = P + C(InC(IH2Ogi,NABins)) * R * T * 1.e4_doubp
  Rho1=TotP / (10000.0_doubp  * R / WtAir * T)*1000.0_doubp        

  Re_2 = 2._doubp*r_2*velo2/(eetta/(Rho*1000))                         
  
  K_12 = velo1*abs(velo2 - velo1)/(r_2*9.81_doubp)

  IF (K_12.GT.1.214) THEN
     E_v12 = (1._doubp + 0.75_doubp*log(2._doubp*K_12)/(K_12-1.214_doubp))**(-2)
  ELSE
     E_v12 = 0.0_doubp
  ENDIF
  
  E_A12 = K_12**2/(K_12 +0.5_doubp)**2
  
  col_ef = (60.0_doubp*E_v12 + E_A12*Re_2)/(60.0_doubp+Re_2)
  
END Function col_ef

!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


SUBROUTINE FIND_BIN
  
!C ********************************************************************
!C    Subroutine calculates target bins for collision products        *
!C                                                                    *
!C ********************************************************************
!
!  At the end something added for C_Y3 and C_Y4, but it is not tested.
!
!

  USE headfile
  use precisi
  implicit none
  
  
  REAL(doubp),DIMENSION(NABins) :: Vlo, Vhi, d_p, Vtemp
  INTEGER :: I_best,L3,II,III,JJ,I_mode,JJ2,I_mode2,IND_MIN
  REAL(doubp),DIMENSION(NY1+1) :: Vtemp_ice
  REAL(doubp),DIMENSION(NY1) :: Vtemp_ice2
  REAL(doubp),DIMENSION(NY2) :: Vtemp_CL2
  REAL(doubp),DIMENSION(NY2+1) :: Vtemp_CL
  REAL(doubp),DIMENSION(NY3+1) :: Vtemp_C3
  REAL(doubp),DIMENSION(NY3) :: Vtemp_C3_2
  REAL(doubp),DIMENSION(NY4+1) :: Vtemp_C4
  REAL(doubp),DIMENSION(NY4) :: Vtemp_C4_2

  Vtemp=4.0_doubp/3.0_doubp*Pi*Rp(1:NaBins)**3
  Vtemp_ice=4.0_doubp/3.0_doubp*Pi*R_Y1**3
  Vtemp_ice2=4.0_doubp/3.0_doubp*Pi*R_ice**3
  Vtemp_CL=4.0_doubp/3.0_doubp*Pi*R_Y2**3
  Vtemp_C3=4.0_doubp/3.0_doubp*Pi*R_Y3**3
  Vtemp_C3_2=4.0_doubp/3.0_doubp*Pi*R_C3**3
  Vtemp_C4=4.0_doubp/3.0_doubp*Pi*R_Y4**3
  Vtemp_C4_2=4.0_doubp/3.0_doubp*Pi*R_C4**3

!
! FIND TARGET BIN IN AEROSOL PARTICLE -AEROSOL PARTICLE COLLISIONS
!

!!
!! System to decide if particles should be moved to cloud drops is missing
!!

  IF(NaModes.EQ.1) THEN
     DO L3 = 1,NABins
        small:   DO II = L3,NABins
           !        IFC_COAG(L3,II) = 0
           DO III = L3, NABins-1
              IF (Vtemp(L3)+Vtemp(II) > Vtemp(III) .AND. &
                   &  Vtemp(L3)+Vtemp(II) < Vtemp(III+1)) THEN
                 G_ind(L3,II) = float(III)+ &
                      & (((Vtemp(L3)+Vtemp(II))-Vtemp(III))/(Vtemp(III+1)-Vtemp(III)))
                 cycle small
              ENDIF
           END DO
           G_ind(L3,II) = NABins-1.0e-8
           
        END DO small
     END DO
  ELSE
     JJ=0
     DO I_mode = 1,NaModes
        DO L3 = 1, NBinsPerMode(I_mode)
           
           small_mod:   DO II = L3, NBinsPerMode(I_mode)
              !        IFC_COAG(L3,II) = 0
                 
                 DO III = L3, NBinsPerMode(I_mode)-1
!                    write(*,*) NBinsPerMode(I_mode), III,JJ,L3,II
                    IF (Vtemp(L3+JJ)+Vtemp(II+JJ) > Vtemp(III+JJ) .AND. &
                         &  Vtemp(L3+JJ)+Vtemp(II+JJ) < Vtemp(III+1+JJ)) THEN
                       G_ind(L3+JJ,II+JJ) = float(III+JJ)+ &
                            & (((Vtemp(L3+JJ)+Vtemp(II+JJ))-Vtemp(III+JJ))/(Vtemp(III+1+JJ)-Vtemp(III+JJ)))
                       cycle small_mod
                    ENDIF
                 END DO
                 G_ind(L3+JJ,II+JJ) =  NBinsPerMode(I_mode)+JJ-1.0e-8

           END DO small_mod
         
        END DO
        JJ=JJ+NBinsPerMode(I_mode)
     END DO
     
!!!     write(*,*) '!!  FIX THE TARGET MODE TO BE CALCULATED BASED ON DRY VOLUME !!'

     JJ=0
     JJ2=0
     DO I_mode = 1,NaModes-1
        JJ2=sum(NBinsPerMode(1:I_mode))
        DO I_mode2 =  I_mode+1,NaModes
           DO L3 = 1, NBinsPerMode(I_mode)
             

              small_mixmod:   DO II = 1, NBinsPerMode(I_mode2)
                 !        IFC_COAG(L3,II) = 0

                    
                    IF (Rp(L3+JJ)>Rp(II+JJ2)) THEN
                       
                       DO III = L3, NBinsPerMode(I_mode)-1
                          
                          IF (Vtemp(L3+JJ)+Vtemp(II+JJ2) > Vtemp(III+JJ) .AND. &
                               &  Vtemp(L3+JJ)+Vtemp(II+JJ2) < Vtemp(III+1+JJ)) THEN
                             
                             G_ind(L3+JJ,II+JJ2) = float(III+JJ)+ &
                                  & (((Vtemp(L3+JJ)+Vtemp(II+JJ2))-Vtemp(III+JJ))/(Vtemp(III+1+JJ)-Vtemp(III+JJ)))
                             cycle small_mixmod
                          ENDIF
              
                       
                          
                       END DO
                       G_ind(L3+JJ,II+JJ2) =  NBinsPerMode(I_mode)+JJ-1.0e-8
                       
                    ELSE
                       
                       DO III = II, NBinsPerMode(I_mode2)-1
                          
                          IF (Vtemp(L3+JJ)+Vtemp(II+JJ2) > Vtemp(III+JJ2) .AND. &
                               &  Vtemp(L3+JJ)+Vtemp(II+JJ2) < Vtemp(III+1+JJ2)) THEN
                             
                             G_ind(L3+JJ,II+JJ2) = float(III+JJ2)+ &
                                  & (((Vtemp(L3+JJ)+Vtemp(II+JJ2))-Vtemp(III+JJ2))/(Vtemp(III+1+JJ2)-Vtemp(III+JJ2)))
                             
                             cycle small_mixmod
                          ENDIF
             
                       END DO
                 
                       G_ind(L3+JJ,II+JJ2) =  NBinsPerMode(I_mode2)+JJ2-1.0e-8
                    ENDIF

              END DO small_mixmod
           
           END DO
          
           JJ2=JJ2+NBinsPerMode(I_mode2)
        END DO
        JJ=JJ+NBinsPerMode(I_mode)
     END DO
     


  ENDIF
 
!
! Make correction to case there cloud droplets are present in different distribution C_Y2
!
  IFCLOUD=1
  IF(IFCLOUD.EQ.1) THEN 
     DO L3 = 1,NABins
        small_cl:   DO II = L3,NABins
           !        IFC_COAG(L3,II) = 0
           IF (Rp(II) > rmin_cloud .AND. Rp(L3) > rmin_cloud) THEN
              DO III = 1, NY2
                 IF (Vtemp(L3)+Vtemp(II) > Vtemp_CL(III) .AND. &
                      &  Vtemp(L3)+Vtemp(II) < Vtemp_CL(III+1)) THEN
                    
                    G_ind(L3,II) = float(III)
                    G_ind(II,L3) = float(III)
                    cycle small_cl
                 ENDIF
                 G_ind(L3,II) = NY2-1.0e-8
                 G_ind(II,L3) = NY2-1.0e-8
              END DO
           END IF
           
        END DO small_cl
     END DO
  ENDIF


!
! FIND TARGET BIN IN AEROSOL PARTICLE - ICE PARTICLE COLLISIONS.
! TWO POSSIBILITIES, FIRST THAT ICE BIN CONTAINS PARTICLES AND SECOND ITS EMPTY
!
  DO L3 = 1,NABins
     small2:   DO II = 1,NY1
        IF(C(L3)*C_Y1(II)<1.0e-50) cycle small2
        DO III = II, NY1
           IF (Vtemp(L3)+Vtemp_ice2(II) > Vtemp_ice(III) .AND. &
                &  Vtemp(L3)+Vtemp_ice2(II) < Vtemp_ice(III+1)) THEN
              JJ=II+NABins
              G_ind(L3,JJ) = float(III)
              G_ind(JJ,L3) = float(III)
!
! NOTE: This is added to avoid negative concentrations in coagulation routine. 
!
              IF(Vtemp(L3)<Vtemp_ice2(II)*0.02) THEN
                 G_ind(L3,JJ) = float(II)
                 G_ind(JJ,L3) = float(II)
              ENDIF
              cycle small2
           ENDIF
        END DO
        G_ind(L3,II+NABins) = float(NY1)
        G_ind(II+NABins,L3) = float(NY1)
     END DO small2
  END DO

!
! FIND TARGET BIN IN ICE PARTICLE - ICE PARTICLE COLLISIONS
!

  DO L3 = 1,NY1
     small3:   DO II = L3,NY1
        IF(C_Y1(L3)*C_Y1(II)<1.0e-50) cycle small3
        DO III = 1, NY1
           IF (Vtemp_ice2(L3)+Vtemp_ice2(II) > Vtemp_ice(III) .AND. &
                &  Vtemp_ice2(L3)+Vtemp_ice2(II)< Vtemp_ice(III+1)) THEN
              JJ=II+NABins
              G_ind(L3+NABINS,JJ) = float(III)
              cycle small3
           ENDIF
        END DO
        G_ind(L3+NABINS,II+NABINS) = float(NY1)
        
     END DO small3
  END DO


!
! FIND TARGET BIN IN aerosol particle - cloud droplet coalescence
!

  Vtemp_CL2=4.0_doubp/3.0_doubp*Pi*R_CL**3
  
  DO L3 =  1,NABINS
     small4:   DO II = 1,NY2
        IF(C(L3)*C_Y2(II)<1.0e-50) cycle small4
        DO III = II, NY2
!           IF (C_Y2(II)>1.e-10) THEN
              IF (Vtemp(L3)+Vtemp_CL2(II) > 1.0*Vtemp_CL(III) .AND. &
                   &  Vtemp(L3)+Vtemp_CL2(II) <1.0*Vtemp_CL(III+1)) THEN
                 JJ=II+NABins+NY1
                 G_ind(L3,JJ) = float(III)
                 G_ind(JJ,L3) = float(III)
                 IF(Vtemp(L3)<Vtemp_CL2(II)*0.02) THEN
                    G_ind(L3,JJ) = float(II)
                    G_ind(JJ,L3) = float(II)
                 ENDIF
                 cycle small4
              ENDIF
!           ELSE
!              IF (Vtemp(L3)+.5*Vtemp_CL(II)+.5*Vtemp_CL(II+1) > Vtemp_CL(III) .AND. &
!                   &  Vtemp(L3)+.5*Vtemp_CL(II)+.5*Vtemp_CL(II+1) < Vtemp_CL(III+1)) THEN
!                 JJ=II+NABins+NY1
!                 G_ind(L3,JJ) = float(III)
!                 G_ind(JJ,L3) = float(III)
!                 cycle small4
!              ENDIF
!           ENDIF

        END DO
        G_ind(L3,II+NABINS+NY1) = float(NY2)
        G_ind(II+NABINS+NY1,L3) = float(NY2)
        IF(Vtemp(L3)+Vtemp_CL2(II) < Vtemp_CL(1))  G_ind(L3,II+NABINS+NY1) = 1.0
        IF(Vtemp(L3)+Vtemp_CL2(II) < Vtemp_CL(1))  G_ind(II+NABINS+NY1,L3) = 1.0
     END DO small4
  END DO

!
! FIND TARGET BIN IN ice nuclei - cloud droplet coalescence
!

  DO L3 =  1,NY1
     small5:   DO II = 1,NY2
        IF(C_Y1(L3)*C_Y2(II)<1.0e-50) cycle small5
        DO III = 1, NY1
           
           IF (Vtemp_ice2(L3)+Vtemp_CL2(II) > Vtemp_ice(III) .AND. &
                &  Vtemp_ice2(L3)+Vtemp_CL2(II) < Vtemp_ice(III+1)) THEN
              JJ=II+NABins+NY1
              G_ind(JJ,L3+NABins) = float(III)
              G_ind(L3+NABins,JJ) = float(III)
              IF(Vtemp_CL2(II)<Vtemp_ice2(L3)*0.02) THEN
                 G_ind(JJ,L3+NABins) = float(L3)
                 G_ind(L3+NABins,JJ) = float(L3)
              ENDIF
           
              cycle small5
           ENDIF
           


        END DO
        G_ind(L3+NABINS,II+NABINS+NY1) = float(NY1) 
        G_ind(II+NABINS+NY1,L3+NABINS) = float(NY1) 
     END DO small5
  END DO

!
! FIND TARGET BIN IN cloud droplet - cloud droplet coalescence
!

  DO L3 =  1,NY2
     small6:   DO II = 1,NY2
        !        IFC_COAG(L3,II) = 0
        IF(C_Y2(L3)*C_Y2(II)<1.0e-50) cycle small6
        IND_MIN=min(L3,II)
        DO III = IND_MIN, NY2
           
           IF (Vtemp_CL2(L3)+Vtemp_CL2(II) > Vtemp_CL(III) .AND. &
                &  Vtemp_CL2(L3)+Vtemp_CL2(II) < Vtemp_CL(III+1)) THEN
              JJ=II+NABins+NY1
              G_ind(L3+NABins+NY1,JJ) = float(III)
              G_ind(JJ,L3+NABins+NY1) = float(III)
              cycle small6
           ENDIF
           
        END DO
        G_ind(L3+NABINS+NY1,II+NABINS+NY1) = float(NY2)
        G_ind(II+NABINS+NY1,L3+NABINS+NY1) = float(NY2)
        IF(Vtemp_CL2(L3)+Vtemp_CL2(II) < Vtemp_CL(1))  G_ind(L3+NABINS+NY1,II+NABINS+NY1) = 1.0
        IF(Vtemp_CL2(L3)+Vtemp_CL2(II) < Vtemp_CL(1))  G_ind(II+NABINS+NY1,L3+NABINS+NY1) = 1.0
     END DO small6
  END DO


!
! FIND TARGET BIN IN aerosol particle - ice (graupel) particle coalescence
!
  
  DO L3 =  1,NABINS
     small7:   DO II = 1,NY3
        IF(C(L3)*C_Y2(II)<1.0e-50) cycle small7
        DO III = 1, NY3
           
           IF (Vtemp(L3)+Vtemp_C3_2(II) > Vtemp_C3(III) .AND. &
                &  Vtemp(L3)+Vtemp_C3_2(II) < Vtemp_C3(III+1)) THEN
              JJ=II+NABins+NY1+NY2
              G_ind(L3,JJ) = float(III)
              G_ind(JJ,L3) = float(III)
              cycle small7
           ENDIF

           G_ind(L3,II+NABINS+NY1+NY2) = float(NY3)
           G_ind(II+NABINS+NY1+NY2,L3) = float(NY3)
           
           IF(Vtemp(L3)+Vtemp_C3_2(II) < Vtemp_C3(1))  THEN
              G_ind(L3,II+NABINS+NY1+NY2) = 1.0
              G_ind(II+NABINS+NY1+NY2,L3) = 1.0
           ENDIF
        END DO
     END DO small7
  END DO

  !
  ! FIND TARGET BIN IN ice (graubel) nuclei - cloud droplet coalescence
  !

  DO L3 =  1,NY3
     small8:   DO II = 1,NY2
        IF(C_Y3(L3)*C_Y2(II)<1.0e-50) cycle small8
        DO III = 1, NY3
           
           IF (Vtemp_C3_2(L3)+Vtemp_CL2(II) > Vtemp_C3(III) .AND. &
                &  Vtemp_C3_2(L3)+Vtemp_CL2(II) < Vtemp_C3(III+1)) THEN
              JJ=II+NABins+NY1
              G_ind(JJ,NABins+NY1+NY2+L3) = float(III)
              G_ind(NABins+NY1+NY2+L3,JJ) = float(III)
              cycle small8
           ENDIF
           
        END DO
        G_ind(II+NABins+NY1,NABins+NY1+NY2+L3) = float(NY3) 
        G_ind(NABins+NY1+NY2+L3,II+NABins+NY1) = float(NY3) 
     END DO small8
  END DO

  !
  ! FIND TARGET BIN IN ice (graubel) nuclei - ice particle coalescence
  !

  DO L3 =  1,NY3
     small9:   DO II = 1,NY1
        IF(C_Y3(L3)*C_Y1(II)<1.0e-50) cycle small9
        DO III = 1, NY3
           
           IF (Vtemp_C3_2(L3)+Vtemp_ice2(II) > Vtemp_C3(III) .AND. &
                &  Vtemp_C3_2(L3)+Vtemp_ice2(II) < Vtemp_C3(III+1)) THEN
              JJ=II+NABins
              G_ind(JJ,NABins+NY1+NY2+L3) = float(III)
              G_ind(NABins+NY1+NY2+L3,JJ) = float(III)
              cycle small9
           ENDIF
           
        END DO
        G_ind(II+NABins,NABins+NY1+NY2+L3) = float(NY3) 
        G_ind(NABins+NY1+NY2+L3,II+NABins) = float(NY3) 
     END DO small9
  END DO



  !
  ! FIND TARGET BIN IN ice (graubel) nuclei - ice (graubel particle coalescence
  !

  DO L3 =  1,NY3
     small10:   DO II = L3,NY3
        IF(C_Y3(L3)*C_Y1(II)<1.0e-50) cycle small10
        JJ=NABins+NY1+NY2
        DO III = L3, NY3
           
           IF (Vtemp_C3_2(L3)+Vtemp_C3_2(II) > Vtemp_C3(III) .AND. &
                &  Vtemp_C3_2(L3)+Vtemp_C3_2(II) < Vtemp_C3(III+1)) THEN
              G_ind(JJ+L3,JJ+II) = float(III)
              G_ind(JJ+II,JJ+L3) = float(III)
              cycle small10
           ENDIF
           
        END DO
        G_ind(II+JJ,L3+JJ) = float(NY3) 
        G_ind(JJ+L3,II+JJ) = float(NY3) 
     END DO small10
  END DO

  !
! FIND TARGET BIN IN aerosol particle - ice (snow) particle coalescence
!

  DO L3 =  1,NABINS
     small11:   DO II = 1,NY4
        IF(C(L3)*C_Y2(II)<1.0e-50) cycle small11
        JJ=II+NABins+NY1+NY2+NY3
        DO III = II, NY4
           
           IF (Vtemp(L3)+Vtemp_C4_2(II) > Vtemp_C4(III) .AND. &
                &  Vtemp(L3)+Vtemp_C4_2(II) < Vtemp_C4(III+1)) THEN
              JJ=II+NABins+NY1+NY2+NY3
              G_ind(L3,JJ) = float(III)
              G_ind(JJ,L3) = float(III)
              cycle small11
           ENDIF

           G_ind(L3,II+NABINS+NY1+NY2) = float(NY4)
           G_ind(II+NABINS+NY1+NY2,L3) = float(NY4)
           
           IF(Vtemp(L3)+Vtemp_C4_2(II) < Vtemp_C4(1))  THEN
              G_ind(L3,JJ) = 1.0
              G_ind(JJ,L3) = 1.0
           ENDIF
        END DO
     END DO small11
  END DO

  !
  ! FIND TARGET BIN IN ice (snow) nuclei - cloud droplet coalescence
  !  NOT DONE FROM THIS ON

!  DO L3 =  1,NY3
!     small8:   DO II = 1,NY2
!        IF(C_Y3(L3)*C_Y2(II)<1.0e-50) cycle small8
!        DO III = 1, NY3
!           
!           IF (Vtemp_C3_2(L3)+Vtemp_CL2(II) > Vtemp_C3(III) .AND. &
!                &  Vtemp_C3_2(L3)+Vtemp_CL2(II) < Vtemp_C3(III+1)) THEN
!              JJ=II+NABins+NY1
!              G_ind(JJ,NABins+NY1+NY2+L3) = float(III)
!              G_ind(NABins+NY1+NY2+L3,JJ) = float(III)
!              cycle small8
!           ENDIF
!           
!        END DO
!        G_ind(II+NABins+NY1,NABins+NY1+NY2+L3) = float(NY3) 
!        G_ind(NABins+NY1+NY2+L3,II+NABins+NY1) = float(NY3) 
!     END DO small8
!  END DO

!  !
!  ! FIND TARGET BIN IN ice (graubel) nuclei - ice particle coalescence
!  !
!
!  DO L3 =  1,NY3
!     small9:   DO II = 1,NY1
!        IF(C_Y3(L3)*C_Y1(II)<1.0e-50) cycle small9
!        DO III = 1, NY3
!           
!           IF (Vtemp_C3_2(L3)+Vtemp_ice2(II) > Vtemp_C3(III) .AND. &
!                &  Vtemp_C3_2(L3)+Vtemp_ice2(II) < Vtemp_C3(III+1)) THEN
!              JJ=II+NABins
!              G_ind(JJ,NABins+NY1+NY2+L3) = float(III)
!              G_ind(NABins+NY1+NY2+L3,JJ) = float(III)
!              cycle small9
!           ENDIF
!           
!        END DO
!        G_ind(II+NABins,NABins+NY1+NY2+L3) = float(NY3) 
!        G_ind(NABins+NY1+NY2+L3,II+NABins) = float(NY3) 
!     END DO small9
!  END DO



!  !
!  ! FIND TARGET BIN IN ice (graubel) nuclei - ice (graubel particle coalescence
!  !

!  DO L3 =  1,NY3
!     small10:   DO II = L3,NY3
!        IF(C_Y3(L3)*C_Y1(II)<1.0e-50) cycle small10
!        JJ=NABins+NY1+NY2
!        DO III = L3, NY3
!           
!           IF (Vtemp_C3_2(L3)+Vtemp_C3_2(II) > Vtemp_C3(III) .AND. &
!                &  Vtemp_C3_2(L3)+Vtemp_C3_2(II) < Vtemp_C3(III+1)) THEN
!              G_ind(JJ+L3,JJ+II) = float(III)
!              G_ind(JJ+II,JJ+L3) = float(III)
!              cycle small10
!           ENDIF
!           
!        END DO
!        G_ind(II+JJ,L3+JJ) = float(NY3) 
!        G_ind(JJ+L3,II+JJ) = float(NY3) 
!     END DO small10
!  END DO


END SUBROUTINE FIND_BIN
