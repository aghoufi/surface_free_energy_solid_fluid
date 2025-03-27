!!!!------ Fortran code for a system of charged particules with the interface along xy plane and z is the normel of the interface


     do i=1,N-1
        do j=i+1,N 
           dx =xxx(i) - xxx(j)
           dy =yyy(i) - yyy(j)
           dz =zzz(i) - zzz(j)
           call conditions_perio (dx,dy,dz,dx,dy,dz)
           rr=sqrt(dx**2+dy**2+dz**2)
 
!---Stress - Lennard-Jonnes  
           stress_LJ=-24.d0*epsilon_Mix*((sigma_Mix/rr)**6)*(2.d0*((sigma_Mix/rr)**6)-1.d0)/(rr**2)
!---Stress - Ewald-Real
           stress_ew3 =charge(i)*charge(j)* (Erfc(rr*alpha_ew)
     &            + ((2.*alpha_ew*rr)/(Pi**0.5))* exp(-(alpha_ew*rr)**2)) /(rr**3)
!-----Pressure     
     
           pressz= (-((dz *dz))*stress_LJ)+(((dzm *dz))*stress_ew3*(1.667563183d5*8.31/1000.0) )
           pressx= (-((dx *dx))*stress_LJ)+(((dxm *dx))*stress_ew3*(1.667563183d5*8.31/1000.0) )
           pressy= (-((dy *dy))*stress_LJ)+(((dym *dy))*stress_ew3*(1.667563183d5*8.31/1000.0) )
       
           call press_ik(zzz(i),zzz(j))
  
        enddo
      enddo   
  


!--------------------------------------------------------------------------------
         SUBROUTINE press_ik (pos1,pos2)
!--------------------------------------------------------------------------------
!-                           --- Calcule les composantes du tenseur de pression -
!-                           --- (méthode Irving-Kirkwood).                     -
!--------------------------------------------------------------------------------           
!    boxl(9) is the dimension the bow length according to z direction (normal)
!    maxbin - number of bin allowing the slicing along z

!--------------------------------------------------------------------------------
            deltaz =  boxl(9)/(2.0*real(maxbin)+1.)         
            xfx=pressx
            yfy=pressy
            zfz=pressz  

            rzi = pos1-boxl(9)*anint(pos1/boxl(9))
            rzj = pos2-boxl(9)*anint(pos2/boxl(9))
            IF (rzi.gt.rzj) THEN
              rzi = pos2-boxl(9)*anint(pos2/boxl(9))
              rzj = pos1-boxl(9)*anint(pos1/boxl(9))
            ENDIF
            zij = rzj-rzi
            zlim = abs(zij)
            zij = zij-boxl(9)*anint(zij/boxl(9))
            zij = abs(zij)

            iz  = anint(rzi/(deltaz))
            jz  = anint(rzj/deltaz)
            jzlim = anint(0.5*boxl(9)/deltaz)
            izlim = -jzlim
            IF (iz.eq.jz) THEN
              Pxx_IK(iz) = Pxx_IK(iz)+(xfx+yfy)/2.
              Pyy_IK(iz) = Pyy_IK(iz)+(yfy+xfx)/2.
              Pzz_IK(iz) = Pzz_IK(iz)+zfz       
            ELSE
              IF (zlim.lt.(0.5*boxl(9))) THEN
                a = (REAL(iz)*deltaz+0.5*deltaz)-rzi
                b = rzj-(REAL(jz)*deltaz-0.5*deltaz)
              ELSE
                a = rzi-(REAL(iz)*deltaz-0.5*deltaz)
                b = (REAL(jz)*deltaz+0.5*deltaz)-rzj
              ENDIF

              Pxx_IK(iz) = Pxx_IK(iz)+a*(xfx+yfy)/(2*zij)
              Pyy_IK(iz) = Pyy_IK(iz)+a*(yfy+xfx)/(2*zij)
              Pzz_IK(iz) = Pzz_IK(iz)+a*zfz/zij
              Pxx_IK(jz) = Pxx_IK(jz)+b*(xfx+yfy)/(2*zij)
              Pyy_IK(jz) = Pyy_IK(jz)+b*(yfy+xfx)/(2*zij)
              Pzz_IK(jz) = Pzz_IK(jz)+b*zfz/zij 
              IF (zlim.lt.(0.5*boxl(9))) THEN
                DO k=iz+1,jz-1
                  Pxx_IK(k) = Pxx_IK(k)+deltaz*(xfx+yfy)/(2*zij)
                  Pyy_IK(k) = Pyy_IK(k)+deltaz*(yfy+xfx)/(2*zij)
                  Pzz_IK(k) = Pzz_IK(k)+deltaz*zfz/zij          
                ENDDO
              ELSE
                DO k=izlim,iz-1
                  Pxx_IK(k) = Pxx_IK(k)+deltaz*(xfx+yfy)/(2*zij)
                  Pyy_IK(k) = Pyy_IK(k)+deltaz*(yfy+xfx)/(2*zij)
                  Pzz_IK(k) = Pzz_IK(k)+deltaz*zfz/zij          
                ENDDO
                DO k=jz+1,jzlim
                  Pxx_IK(k) = Pxx_IK(k)+deltaz*(xfx+yfy)/(2*zij)
                  Pyy_IK(k) = Pyy_IK(k)+deltaz*(yfy+xfx)/(2*zij)
                  Pzz_IK(k) = Pzz_IK(k)+deltaz*zfz/zij          
                ENDDO
              ENDIF
            ENDIF

      end subroutine  press_ik
      
      
                 
      