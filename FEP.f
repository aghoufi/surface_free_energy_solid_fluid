         program   FEP
   
         implicit none 
         
         include 'ener.inc'

!--------------------------------------------
	 
         Myunit = 2
         
	 open(Myunit,file='vdw.dat',status='unknown')
	 read(Myunit,*) Nconf 
	 read(Myunit,*)nmoltype
	 do i=1,nmoltype
	    read(Myunit,*)Nbmolec(i),Lengt(i)
	 enddo
	 do i=1,nmoltype
	    read(Myunit,*) Ntypevdw(i)
	    do j=1,Ntypevdw(i)
	       read(Myunit,*)atnamevdw(i,j),epsilont(i,j),sigmat(i,j)
	    enddo
	 enddo
	 read(Myunit,*) kmaxx,kmaxy,kmaxz,alpha_ew
	 read(Myunit,*)cutoff	
	 read(Myunit,*)Nfile_Error
         read(Myunit,*)typePertu 
         read(Myunit,*)Nlambda,Nwindows
         do i=1,Lengt(typePertu)
         read(Myunit,*)epsilon_ini(i)
         enddo
	 close(Myunit)   
	 epsilon_TA=0.1
!--------------------------------------------
	 calcul=.true.
  	 avtest=0.
         avterm2=0.
         avterm1=0. 
         Myunit = 3
	 Nfile = 0       	 
	 Etot_ww = 0.0
         Etot_wp = 0.0
         Etot_wm = 0.0
         Etot_pp = 0.0
         Etot_mm = 0.0
         Etot_pm = 0.0
         
         hist_rdf(:,:)= 0.
         av_volume=0.
         Reci_elec =0.
         avEtot_vdw=0.
	 avEtot_vdw2=0.
         avEtot_elec=0.
         Eelectot_self=0.
         Eelectot_self1=0.
         Eelectot_self2=0.
         Eelectot_self3=0.
         
         Reci_elec1=0.
         Reci_elec2=0.
         Reci_elec3=0.
         
         ppnn=0.
         pptt=0.
         
         hist_z(:,:)=0.
	 pn_z(:)=0.
	 pt_z(:)=0.
	 pt2_z(:)=0.
	 gamma_z(:)=0.
         
	 pn_z2(:)=0.
	 pt_z2(:)=0.
	 pt2_z2(:)=0.
	 gamma_z2(:)=0.
	 
	 avPxx_IK(:)=0.
	 avPyy_IK(:)=0.
	 avpzz_IK(:)=0.
	 
	 avPxx_Gh(:)=0.
	 avPyy_Gh(:)=0.
	 avpzz_Gh(:)=0.
	 
	 av_diffU=0.
	 av_diff_E(:)=0.
	 av_diff_EGh(:)=0.
	     
         open(Myunit,file='HISTORY', status='unknown')                         
        	

!------------------------------------------------------------------
         do i=1,Nfile_Error
            read(Myunit,*)
         enddo 
	 read(Myunit,*)
	 read(Myunit,*) imcon, ikey, Natms 
         do i=1, Nconf   
            read(Myunit,*)          
            read(Myunit,*) boxl(1), boxl(2),boxl(3)
            read(Myunit,*) boxl(4), boxl(5),boxl(6)
            read(Myunit,*) boxl(7), boxl(8),boxl(9) 
	    Nfile=Nfile + 1                                      
            do j=1, Natms
               read(Myunit,*) atname(j), indice, masse(j)!, charge(j)
               read(Myunit,*) xxx(j), yyy(j), zzz(j)
               if(atname(j).eq.'Cg') masse(j)=12.
            enddo
         !   diff =100.-boxl(9)
	    boxl9h=boxl(9)!+diff

              	     	                 

	     
	    
	   	   
	  
	     
	     
1000   FORMAT('ATOM  ',i5,2x,a1,8x,i4,4x,3f8.3,16x,a2) 
2000   FORMAT('REMARK 2000',i5)
	     
!------------------------------------------------------------             
	     
             volume=abs(boxl(1)*(boxl(5)*boxl(9)-boxl(6)*boxl(8))+
     &boxl(2)*(boxl(6)*boxl(7)-boxl(4)*boxl(9))+
     &boxl(3)*(boxl(4)*boxl(8)-boxl(5)*boxl(7)))
             
             av_volume=av_volume+volume

!----------------------------------------------------------------- Transformation
!-----------------------------------------------------------------

	  if(Nfile.eq.1)then	     
             compt=0
             compt_molec=0
	     do ii=1,nmoltype
	        do j=1,Nbmolec(ii)
                   compt_molec=compt_molec+1
		   do k=1,Lengt(ii)
		      compt=compt+1
		       
		      indice_molec(compt)= compt_molec
		      indice_Lengt(compt)=k
		      indice_type (compt)=ii
		      indice_atm  (ii,j,k)=compt
		      
		      do l=1,Ntypevdw(ii)
		          if(atname(compt).eq.atnamevdw(ii,l))then
		             sigma(ii,k)   = sigmat(ii,l)
			     epsilon(ii,k) = epsilont(ii,l)
			   endif   
		      enddo
                   enddo
		 enddo
             enddo
	  endif  		
	  
!!----------------------------------------------
             boxl9h=boxl(9)!+diff
	     call reconstruction()	    	     		 

            open(643,file='visualisation.pdb',status='unknown')
            do j=1, Natms
            write(643,1000) j,atname(j),j, xxx (j), yyy (j), zzz(j),atname(j)                          
            enddo
            close(643)

             call init_fact_ew()  
	     call Inter_Realspace()
             call Inter_Reci(xxx,yyy,zzz,0)
            
	     U_LJ_0=Etot_tot_0
	     Reci1=Reci_elec 
             sigma_old(:,:)=sigma(:,:) 
             epsilon_old(:,:)=epsilon(:,:)
             charge_old(:)=charge(:)
             
	     call trans_Area( X_new,Y_new,Z_new,Xmol_New,Ymol_New,Zmol_New,1)
             
             call init_fact_ew()
             call energy_TA(1,X_new,Y_new,Z_new,Xmol_New,Ymol_New,Zmol_New,U_LJ_1,Energy_bin1)
             call Inter_Reci(X_new,Y_new,Z_new,1)
             Reci2=Reci_elec   
             sigma(:,:)=sigma_old(:,:)
             epsilon(:,:)=epsilon_old(:,:)
             charge(:)=charge_old(:) 
!!----------------------------------------------
!!----------------------------------------------
!!---------------------------------------------- Macroscopique
             av_diffU=av_diffU + exp(-(U_LJ_1+Reci2 - U_LJ_0-Reci1)/(298.15*8.31/1000.0))!/epsilon_TA
            
             gamma=- ((8.31*298.15)/1000.00)*log(av_diffU/real(Nfile))! (av_diffU/real(Nfile))

!--------------------------------	     
	     
	         
             write(*,'(15g20.8)') real(i),U_LJ_0,U_LJ_1,(gamma/(2*boxl(1)*boxl(5)))*166.279072
!              write(*,'(15g20.8)') real(i),gamma,xmolcm2(250),ymolcm2(250),zmolcm2(250)   
     

10              format(a5,5x,i10)
20              format(3g20.10)
     
!!----------------------------------------------
     
           
    
!!----------------------------------------------
         
	 

!!----------------------------------------------
         
	
          
   
         enddo
         close(Myunit)
	 
	
        
        end      

!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
!!----------------------------------------------------------------------------------
               
      subroutine conditions_perio (rx,ry,rz,rx1,ry1,rz1)         
      
      implicit none
      
      include 'ener.inc'
      
      real(kind=8)                       :: b(10), d, r, ssx, ssy, ssz,Boxlength_Ma(10),bbb(10)
      real(kind=8),intent(in)            :: rx,ry,rz
      real(kind=8),intent(out)           :: rx1,ry1,rz1

    
     
                  
	     Boxlength_Ma(:) = boxl(:)	                  	         	     
             b(1) = Boxlength_Ma(5) * Boxlength_Ma(9) -  Boxlength_Ma(6) * Boxlength_Ma(8)
             b(2) = Boxlength_Ma(3) * Boxlength_Ma(8) -  Boxlength_Ma(2) * Boxlength_Ma(9)
             b(3) = Boxlength_Ma(2) * Boxlength_Ma(6) -  Boxlength_Ma(3) * Boxlength_Ma(5)
             b(4) = Boxlength_Ma(6) * Boxlength_Ma(7) -  Boxlength_Ma(4) * Boxlength_Ma(9)
             b(5) = Boxlength_Ma(1) * Boxlength_Ma(9) -  Boxlength_Ma(3) * Boxlength_Ma(7)
             b(6) = Boxlength_Ma(3) * Boxlength_Ma(4) -  Boxlength_Ma(1) * Boxlength_Ma(6)
             b(7) = Boxlength_Ma(4) * Boxlength_Ma(8) -  Boxlength_Ma(5) * Boxlength_Ma(7)
             b(8) = Boxlength_Ma(2) * Boxlength_Ma(7) -  Boxlength_Ma(1) * Boxlength_Ma(8)
             b(9) = Boxlength_Ma(1) * Boxlength_Ma(5) -  Boxlength_Ma(2) * Boxlength_Ma(4)
             
             d = Boxlength_Ma(1)*b(1) + Boxlength_Ma(4)*b(2) + Boxlength_Ma(7)*b(3)
             r=0.
             if(abs(d).gt.0.) r=1./d
             
             b(1) =  b(1) * r
             b(2) =  b(2) * r
             b(3) =  b(3) * r
             b(4) =  b(4) * r
             b(5) =  b(5) * r
             b(6) =  b(6) * r
             b(7) =  b(7) * r
             b(8) =  b(8) * r
             b(9) =  b(9) * r
                             
             ssx =  b(1)*rx + b(4)*ry + b(7)*rz
             ssy =  b(2)*rx + b(5)*ry + b(8)*rz 
             ssz =  b(3)*rx + b(6)*ry + b(9)*rz                 
             ssx = ssx - ANINT(ssx)
             ssy = ssy - ANINT(ssy)
             ssz = ssz - ANINT(ssz)                   
             rx1 = Boxlength_Ma(1)*ssx + Boxlength_Ma(4)*ssy + Boxlength_Ma(7)*ssz
             ry1 = Boxlength_Ma(2)*ssx + Boxlength_Ma(5)*ssy + Boxlength_Ma(8)*ssz
             rz1 = Boxlength_Ma(3)*ssx + Boxlength_Ma(6)*ssy + Boxlength_Ma(9)*ssz
             
!----------------------
c     calculate lengths of cell vectors

               bbb(1)=
     &sqrt(Boxlength_Ma(1)*Boxlength_Ma(1)+Boxlength_Ma(2)*Boxlength_Ma(2)+Boxlength_Ma(3)*Boxlength_Ma(3))
               bbb(2)=
     &sqrt(Boxlength_Ma(4)*Boxlength_Ma(4)+Boxlength_Ma(5)*Boxlength_Ma(5)+Boxlength_Ma(6)*Boxlength_Ma(6))
               bbb(3)=
     &sqrt(Boxlength_Ma(7)*Boxlength_Ma(7)+Boxlength_Ma(8)*Boxlength_Ma(8)+Boxlength_Ma(9)*Boxlength_Ma(9))

c     calculate cosines of cell angles
           
               bbb(4)=
     &(Boxlength_Ma(1)*Boxlength_Ma(4)+Boxlength_Ma(2)*Boxlength_Ma(5)+Boxlength_Ma(3)*Boxlength_Ma(6))
     &/(bbb(1)*bbb(2))
               bbb(5)=
     &(Boxlength_Ma(1)*Boxlength_Ma(7)+Boxlength_Ma(2)*Boxlength_Ma(8)+Boxlength_Ma(3)*Boxlength_Ma(9))
     &/(bbb(1)*bbb(3))
               bbb(6)=
     &(Boxlength_Ma(4)*Boxlength_Ma(7)+Boxlength_Ma(5)*Boxlength_Ma(8)+Boxlength_Ma(6)*Boxlength_Ma(9))
     &/(bbb(2)*bbb(3))
              
	      volume=abs(boxl(1)*(boxl(5)*boxl(9)-boxl(6)*boxl(8))+
     &boxl(2)*(boxl(6)*boxl(7)-boxl(4)*boxl(9))+
     &boxl(3)*(boxl(4)*boxl(8)-boxl(5)*boxl(7)))
      
	      Unit_a=sin(Acos(bbb(4)))*bbb(1)*bbb(2)/volume
	      Unit_b=sin(Acos(bbb(5)))*bbb(1)*bbb(3)/volume
	      Unit_c=sin(Acos(bbb(6)))*bbb(2)*bbb(3)/volume
              
	      
           
                    
!--------------------------------------------------------------
             
          end subroutine conditions_perio

      
********************************************************************************
      SUBROUTINE reconstruction()
********************************************************************************
      
      IMPLICIT NONE
      INCLUDE 'ener.inc'
      
      real(8)  :: dx, dy,dz, dxx, dyy, dzz, rijsq
      real(8)  ::xmolcmtot,ymolcmtot,zmolcmtot,masstotot
      integer  ::  np, mp,ind_atm(maxmol,maxlength)
      	
	
      
      
      xmolcm(:)  = 0.0
      ymolcm(:)  = 0.0
      zmolcm(:)  = 0.0
      masstot(:) = 0.0
      
      xmolcm2(:)  = 0.0
      ymolcm2(:)  = 0.0
      zmolcm2(:)  = 0.0
      masstot2(:) = 0.0
      
      xmolcm3(:,:)  = 0.0
      ymolcm3(:,:)  = 0.0
      zmolcm3(:,:)  = 0.0
      masstot3(:,:) = 0.0
      
      
      
!-------reconstruction
    	
      do k=1,nmoltype
        do j=1,Nbmolec(k)
	       do i =2,Lengt(k)
	    
	       dx =xxx(indice_atm(k,j,i)) - xxx(indice_atm(k,j,1))!xxx(ind_atm(j,np))
	       dy =yyy(indice_atm(k,j,i)) - yyy(indice_atm(k,j,1))!yyy(ind_atm(j,np)) 
	       dz =zzz(indice_atm(k,j,i)) - zzz(indice_atm(k,j,1))!zzz(ind_atm(j,np)) 
	       call conditions_perio (dx,dy,dz,dx,dy,dz)
	       xxx(indice_atm(k,j,i))=dx+xxx(indice_atm(k,j,1))!xxx(ind_atm(j,np))
	       yyy(indice_atm(k,j,i))=dy+yyy(indice_atm(k,j,1))!yyy(ind_atm(j,np)) 
	       zzz(indice_atm(k,j,i))=dz+zzz(indice_atm(k,j,1))!zzz(ind_atm(j,np)) 
	    enddo
	 enddo
	enddo
	
       xmolcmtot=0.
       ymolcmtot=0.
       zmolcmtot=0.
       masstotot=0.
       do j=1,Natms
          xmolcmtot=xmolcmtot+xxx(j)*masse(j)
	  ymolcmtot=ymolcmtot+yyy(j)*masse(j)
	  zmolcmtot=zmolcmtot+zzz(j)*masse(j)
	  masstotot=masstotot+ masse(j)
       enddo
       xmolcmtot=xmolcmtot/masstotot
       ymolcmtot=ymolcmtot/masstotot
       zmolcmtot=zmolcmtot/masstotot
      
	
!-------centre of masse
  		      
       
      
      compt=0
      do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
            do i =1,Lengt(k)
	           xmolcm2(compt)=xmolcm2(compt) + xxx(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
	           ymolcm2(compt)=ymolcm2(compt) + yyy(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
	           zmolcm2(compt)=zmolcm2(compt) + zzz(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
               masstot2(compt) = masstot2(compt) + masse(indice_atm(k,j,i))
	      
	          xmolcm3(k,j)=xmolcm3(k,j) + xxx(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
	          ymolcm3(k,j)=ymolcm3(k,j) + yyy(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
	          zmolcm3(k,j)=zmolcm3(k,j) + zzz(indice_atm(k,j,i))*masse(indice_atm(k,j,i))
              masstot3(k,j) = masstot3(k,j) + masse(indice_atm(k,j,i))
	      
	        enddo
	        xmolcm2(compt)  = xmolcm2(compt)/masstot2(compt)
            ymolcm2(compt)  = ymolcm2(compt)/masstot2(compt)
            zmolcm2(compt)  = zmolcm2(compt)/masstot2(compt)
	    
	        xmolcm3(k,j)=xmolcm3(k,j)/masstot3(k,j)
            ymolcm3(k,j)=ymolcm3(k,j)/masstot3(k,j)
            zmolcm3(k,j)=zmolcm3(k,j)/masstot3(k,j)
	    


          enddo 
       enddo     

!------------------------------- 

       compt=0
       do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
	        xmolcm2(compt)=xmolcm2(compt) -xmolcmtot
	        ymolcm2(compt)=ymolcm2(compt) -ymolcmtot
	        zmolcm2(compt)=zmolcm2(compt) -zmolcmtot
	       enddo
       enddo	    
       
       
       do j=1,Natms
          xxx(j)=xxx(j) -xmolcmtot
	      yyy(j)=yyy(j) -ymolcmtot
	      zzz(j)=zzz(j) -zmolcmtot	  
       enddo
       
        
       
    
       
  

!-------------------------------
              
      end   

!!----------------------------------------------------------------------------------
               
      subroutine Inter_Realspace         
      
      implicit none
      
      include 'ener.inc'
           
      real(kind=8) :: dx,dy,dz,rr,rr_rad,etmp
      real(kind=8) :: U_LJ,sum_LJ,U_LJ2,sum_LJ2,Pi
      integer      :: tpi,Lei,tpj,Lej,comptmol1,comptmol2,kk,iii,kkk,jjj
      real(kind=8) :: epsilon_Mix,sigma_Mix,DPDErfc
      real(kind=8) :: sum_LJ_ww,sum_ew_ww,sum_LJ_wp,sum_ew_wp,sum_LJ_wm
      real(kind=8) :: sum_ew_wm,sum_LJ_mm,sum_ew_mm,sum_LJ_pp,sum_ew_pp,sum_LJ_pm,sum_ew_pm
      real(kind=8) :: Eelectot,Evdwtot,Npaire_wp,Npaire_wm
      real(kind=8) :: ppnn_tmp_ew,ppnn_tmp_LJ,pptt_tmp_ew,pptt_tmp_LJ,dxm,dym,dzm,rrm
      real(kind=8) :: stress_LJ2,ppnn_tmp_LJ2,pptt_tmp_LJ2
      
      real(kind=8) :: Ash,Bsh,Csh,Ash2,Bsh2,Csh2,repulsif,attractif
      real(kind=8) :: epsilon_Mix0 ,sigma_Mix0,epsilon_Mix1,sigma_Mix1 
      
      Pi=Acos(-1.0)
      
!-------------------------------

      
        Eelectot    =0.
        Evdwtot     =0.

        comptmol1=0
         Do ii=1,1
            Do jj=1,Nbmolec(ii)
            comptmol1=comptmol1+1
            comptmol2=0
                Do iii=2,2
                   Do jjj=1,Nbmolec(iii)
                       comptmol2=comptmol2+1
                      
                          Do kk=1,Lengt(ii)
                             Do kkk=1,Lengt(iii)
                                 i=indice_atm(ii,jj,kk)
                                 j=indice_atm(iii,jjj,kkk)
	                         dx =xxx(i) - xxx(j)
	                         dy =yyy(i) - yyy(j) 
	                         dz =zzz(i) - zzz(j) 
	                         call conditions_perio (dx,dy,dz,dx,dy,dz)
	                         rr=sqrt(dx**2+dy**2+dz**2)
	                         tpi=indice_type(i)
	                         Lei=indice_Lengt(i)
	                         tpj=indice_type(j)
	                         Lej=indice_Lengt(j)
	       		
	                         if(rr.lt.cutoff) then
                                    epstmp= (1.0/(real(Nlambda)-1))*(Nwindows-1) !+ epsilon_TA
                                 
                                    epsilon_Mix  =epsilon_ini(kk)*(1-epstmp)+epsilon_fin(kk)*epstmp
                                    sigma_Mix   =(sigma(tpi,Lei)+sigma(tpj,Lej))/2.0
                                  
	                                repulsif=  4.d0*epsilon_Mix*(sigma_Mix**12)
	                                attractif= 4.d0*epsilon_Mix*(sigma_Mix**6)	    
                                    Evdwtot =Evdwtot+ repulsif*((1./(rr**12))) -attractif*((1./(rr**6)))
		                  
		                    
                                    dxm =xmolcm2(indice_molec(i)) - xmolcm2(indice_molec(j))
                                    dym =ymolcm2(indice_molec(i)) - ymolcm2(indice_molec(j))
                                    dzm =zmolcm2(indice_molec(i)) - zmolcm2(indice_molec(j))
                                    call conditions_perio (dxm,dym,dzm,dxm,dym,dzm)
                                    rrm=sqrt(dxm**2+dym**2+dzm**2)
	                             
                                    if(rrm.lt.cutoff) then   
                                       Eelectot=Eelectot+charge(i)*charge(j)* DPDErfc(rr*alpha_ew)/rr
		                    
		                            endif  
	                        	 endif
                               enddo
                             enddo
                           
                         enddo
                       enddo
                     enddo
                   enddo
      
      
                  avEtot_elec=avEtot_elec+Eelectot*(1.667563183d5*8.31/1000.0)
                  avEtot_vdw =avEtot_vdw +Evdwtot
                  Etot_tot_0=(Evdwtot+Eelectot*(1.667563183d5*8.31/1000.0))
!-------------------------------      	     
      end   




********************************************************************************

      REAL(KIND=8) FUNCTION DPDErfc(xx)

      REAL(KIND=8) :: xx
      REAL(KIND=8) :: zz, tt

      zz = abs(xx)
      tt = 1.d0 / (1.d0+0.5d0*zz)
      
      DPDErfc = tt*EXP(-zz*zz-1.26551223+tt*(1.00002368+tt*(0.37409196+              
     &tt*(0.09678418+tt*(-0.18628806+tt*(0.27886807+tt*(-1.13520398+              
     &tt*(1.48851587+tt*(-0.82215223+tt*0.17087277)))))))))

      IF (xx < 0.d0) DPDErfc = 2.0 - DPDErfc
           
      RETURN
      END

********************************************************************************
********************************************************************************
      SUBROUTINE init_fact_ew ()

!--------- Initialisation pour la sommation d'Ewald    
      
********************************************************************************
      
      IMPLICIT NONE
      INCLUDE 'ener.inc'
  
      
      INTEGER :: kx,ky,kz,symmetry_x,x,y,z
      REAL(KIND=8) :: sqrtPi,k2t,kmax,Pi
      REAL(KIND=8),DIMENSION(3) :: tpl
      
      
********************************************************************************
       
      Pi = Acos(-1.)
      sqrtPi = sqrt(Pi) 
      
     
      tpl(1) = 4*Pi*Pi / (boxl(1)**2)
      tpl(2) = 4*Pi*Pi / (boxl(5)**2)
      tpl(3) = 4*Pi*Pi / (boxl(9)**2)
      
      kmax = min(tpl(1)*kmaxx*kmaxx, tpl(2)*kmaxy*kmaxy,tpl(3)*kmaxz*kmaxz)     
      
      DO kx = 0, kmaxx
         IF (kx == 0) THEN
            symmetry_x = 1
         ELSE
            symmetry_x = 2
         ENDIF
         DO ky = -kmaxy, kmaxy
            DO kz = -kmaxz, kmaxz
               k2t = tpl(1)*kx*kx + tpl(2)*ky*ky + tpl(3)*kz*kz             
               IF( k2t > 0 .and. k2t <= kmax) THEN
                   fact_ew(kx,ky,kz) = exp( -k2t /(4*alpha_ew**2)) / k2t 
     &* symmetry_x                  
               ELSE
               fact_ew(kx,ky,kz) = 0.
               ENDIF
            ENDDO
         ENDDO   
      ENDDO          
!      fact_ew(:, -1:-kmaxy:-1, :           )             = fact_ew(:, 1: kmaxy, :       )
!      fact_ew(:, :           , -1:-kmaxz:-1)             = fact_ew(:,        :, 1: kmaxz)
!      fact_ew(:, -1:-kmaxy:-1, -1:-kmaxz:-1)             = fact_ew(:, 1: kmaxy, 1: kmaxz)
      
      
      RETURN
      END


********************************************************************************
      SUBROUTINE Inter_Reci(X_Ntmp,Y_Ntmp,Z_Ntmp,pert)

!--------- Calculs des contributions dans l'espace reciproque      
      
********************************************************************************
      
      IMPLICIT NONE
      INCLUDE 'ener.inc'

      
      INTEGER :: kx,ky,kz,count,pert
      REAL(KIND=8) , DIMENSION( 3) :: fack 
      REAL(KIND=8) :: fact   ,Pi,Energy_Elec_reci,Energy_Elec_reci1,Energy_Elec_reci2,Energy_Elec_reci3
      COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:) :: eikr    
      real(kind=8), dimension(maxatm)         :: X_Ntmp,Y_Ntmp,Z_Ntmp
      COMPLEX(8) :: sum_zeikr,sum_zeikr1,sum_zeikr2,sum_zeikr3
      COMPLEX(8) :: zeikr  
      
      ALLOCATE(eikr (maxatm, 3, -k_max:k_max))
      
********************************************************************************        Fiare une version moleculaire       
      
      Pi = Acos(-1.)         
      fack(1)    = 2.*Pi/boxl(1) 
      fack(2)    = 2.*Pi/boxl(5) 
      fack(3)    = 2.*Pi/boxl(9) 
     
      eikr(:,:,0) = 1               
!      if(pert.eq.0) then                            
!      eikr(:,1,1) = exp( cmplx(0.d0,fack(1))*xxx(:) )         
!      eikr(:,2,1) = exp( cmplx(0.d0,fack(2))*yyy(:) )         
!      eikr(:,3,1) = exp( cmplx(0.d0,fack(3))*zzz(:) )
!      elseif(pert.eq.1) then
      eikr(:,1,1) = exp( cmplx(0.d0,fack(1))*X_Ntmp(:) )
      eikr(:,2,1) = exp( cmplx(0.d0,fack(2))*Y_Ntmp(:) )
      eikr(:,3,1) = exp( cmplx(0.d0,fack(3))*Z_Ntmp(:) )
!      endif      
   
      
      DO i=2, k_max
          eikr(:,:,i)    = eikr   (:,:,i-1) * eikr(:,:,1)
      ENDDO
      eikr(:,:,-1:-k_max:-1) = conjg(eikr(:,:,1:k_max))
      Energy_Elec_reci =0.0
      DO kz=-kmaxz,kmaxz  
         DO ky=-kmaxy,kmaxy 
            DO kx=0,kmaxx
               fact = fact_ew(kx,ky,kz)	        
               IF(fact/=0.) THEN
                  sum_zeikr              = 0.0 
                  DO i=1,Natms
                        sum_zeikr = sum_zeikr +  charge(i)*eikr(i,1,kx)*eikr(i,2,ky)*eikr(i,3,kz)
                  ENDDO
		  Energy_Elec_reci = Energy_Elec_reci +sum_zeikr*conjg(sum_zeikr)*fact  
                                 
	        ENDIF
	     ENDDO
	 ENDDO
       ENDDO
       volume=boxl(1)*boxl(5)*boxl(9)
       Energy_Elec_reci =(Energy_Elec_reci *2.0*Pi /volume)*(1.667563183d5*8.31/1000.0)!/(4.0*Acos(-1.0))
 
       Reci_elec = Energy_Elec_reci
       
       END	     		  
       



!!----------------------------------------------------------------------------------
               
      subroutine Energy_TA (typepert,X_Ntmp,Y_Ntmp,Z_Ntmp,Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp,U_LJ_tot,Energy_bin)         
      
      implicit none
      
      include 'ener.inc'
           
      real(kind=8) :: dx,dy,dz,rr,rr_rad
      real(kind=8) :: U_LJ,sum_LJ,U_LJ2,sum_LJ2,Pi
      integer      :: tpi,Lei,tpj,Lej,comptmol1,comptmol2,kk,iii,kkk,jjj
      real(kind=8) :: epsilon_Mix,sigma_Mix,DPDErfc
      real(kind=8) :: sum_LJ_ww,sum_ew_ww,sum_LJ_wp,sum_ew_wp,sum_LJ_wm
      real(kind=8) :: sum_ew_wm,sum_LJ_mm,sum_ew_mm,sum_LJ_pp,sum_ew_pp,sum_LJ_pm,sum_ew_pm
      real(kind=8) :: Eelectot,Evdwtot,Npaire_wp,Npaire_wm
      real(kind=8) :: ppnn_tmp_ew,ppnn_tmp_LJ,pptt_tmp_ew,pptt_tmp_LJ,dxm,dym,dzm,rrm
      real(kind=8) :: stress_LJ2,ppnn_tmp_LJ2,pptt_tmp_LJ2
      
      real(kind=8) :: Ash,Bsh,Csh,Ash2,Bsh2,Csh2,repulsif,attractif
      
      integer                                 :: typepert
      real(kind=8), dimension(maxatm)         :: X_Ntmp,Y_Ntmp,Z_Ntmp
      real(kind=8)                            :: U_LJ_tot,etmp
      real(kind=8), dimension(maxmol)         :: Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp
      real(kind=8) :: epsilon_Mix0 ,sigma_Mix0,epsilon_Mix1,sigma_Mix1 
      Pi=Acos(-1.0)
      
!-------------------------------
      
      Evdwtot=0.
      Eelectot=0.
     
      
       comptmol1=0
         Do ii=1,1
            Do jj=1,Nbmolec(ii)
               comptmol1=comptmol1+1
               comptmol2=0
                Do iii=2,2
                   Do jjj=1,Nbmolec(iii)
                       comptmol2=comptmol2+1
                      
                          Do kk=1,Lengt(ii)
                             Do kkk=1,Lengt(iii)
                                 i=indice_atm(ii,jj,kk)
                                 j=indice_atm(iii,jjj,kkk)
	                         dx =X_Ntmp(i) - X_Ntmp(j)
	                         dy =Y_Ntmp(i) - Y_Ntmp(j) 
	                         dz =Z_Ntmp(i) - Z_Ntmp(j) 
	                         call conditions_perio (dx,dy,dz,dx,dy,dz)
	                         rr=sqrt(dx**2+dy**2+dz**2)
	                         tpi=indice_type(i)
	                         Lei=indice_Lengt(i)
	                         tpj=indice_type(j)
	                         Lej=indice_Lengt(j)
	       		
                                 if(rr.lt.cutoff) then 
		                  
!                                  epstmp=   1-(1.0/(real(Nlambda)-1))*(Nwindows-1) + epsilon_TA       
                                   epstmp= (1.0/(real(Nlambda)-1))*(Nwindows-1) + epsilon_TA
                                   epsilon_Mix  =epsilon_ini(kk)*(1-epstmp)+epsilon_fin(kk)*epstmp
                                   sigma_Mix   =(sigma(tpi,Lei)+sigma(tpj,Lej))/2.0
              
                                   repulsif=  4.d0*epsilon_Mix*(sigma_Mix**12)
                                   attractif= 4.d0*epsilon_Mix*(sigma_Mix**6)
                                   Evdwtot =Evdwtot+ repulsif*((1./(rr**12))) -attractif*((1./(rr**6)))
                                   etmp=             repulsif*((1./(rr**12))) -attractif*((1./(rr**6)))
		  
                  
                                   dxm =Xmol_Ntmp(indice_molec(i)) - Xmol_Ntmp(indice_molec(j))
	                               dym =Ymol_Ntmp(indice_molec(i)) - Ymol_Ntmp(indice_molec(j))
	                               dzm =Zmol_Ntmp(indice_molec(i)) - Zmol_Ntmp(indice_molec(j)) 
                                   call conditions_perio (dxm,dym,dzm,dxm,dym,dzm)
	                               rrm=sqrt(dxm**2+dym**2+dzm**2) 
		   
		                
		  
                                   if(rrm.lt.cutoff) then      
	                                Eelectot=Eelectot+charge(i)*charge(j)* DPDErfc(rr*alpha_ew)/rr	         
			                       
			 
			 
                                   endif  
                                endif
                               enddo
                             enddo
                          
                         enddo
                       enddo
                     enddo
                   enddo
      
      
                   U_LJ_tot=(Eelectot*(1.667563183d5*8.31/1000.0)+Evdwtot)
                  
      
   
!-------------------------------      	     
      end   


!---------------------------------------------------------------

      subroutine  trans_Area(X_Ntmp,Y_Ntmp,Z_Ntmp,Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp,delta)
      
      IMPLICIT NONE
            
      include 'ener.inc'
      
      integer :: delta
      real(kind=8), dimension(maxatm)         :: X_Ntmp,Y_Ntmp,Z_Ntmp 
      real(kind=8), dimension(maxmol)         :: Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp
      REAL(KIND=8) :: rappx,rappy,rappz,dx,dy,dz
      REAL(KIND=8) :: incre1,incre2,boxl_pl(10)
      
    
      
     
      
        boxl_pl(1) =boxl(1)
        boxl_pl(5) =boxl(5)
        boxl_pl(9) =boxl(9)
        
       compt=0
       do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
            Xmol_Ntmp(compt)=xmolcm2(compt)
	        Ymol_Ntmp(compt)=ymolcm2(compt)
	        Zmol_Ntmp(compt)=zmolcm2(compt)	       
	     enddo
       enddo       
         	          
      compt=0
      do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
            do i =1,Lengt(k)
	           X_Ntmp(indice_atm(k,j,i))=xxx(indice_atm(k,j,i))
               Y_Ntmp(indice_atm(k,j,i))=yyy(indice_atm(k,j,i))
	           Z_Ntmp(indice_atm(k,j,i))=zzz(indice_atm(k,j,i))
               
               if(k.eq.typePertu) then
                   epstmp=   (1.0/(real(Nlambda)-1))*(Nwindows-1) + epsilon_TA
                             
              endif
	       enddo
	     enddo
       enddo	    


      end subroutine trans_Area



