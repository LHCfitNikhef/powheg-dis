!     DIS phase space taken from nlojet++
!     dPhi = dx dOmega/32 pi^2
!          = dx dcth/16pi * dphi/2pi       [cth = 1-2y]
!          = dx 2dy/16 pi * dphi/2pi  
!     with Q^2 = x y S      
      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'

      real * 8 xborn(ndiminteg-4)
      real * 8 vec(3), beta
! DIS variables
      real * 8 x,Q2,s,y
      real * 8 jac, Eh, El
!     Limits on Q2 and x
      real * 8 Q2min, Q2max, xmin, xmax, ymin, ymax,ymn,ymx
      integer mu,k,pow
      logical ini
      data ini/.true./
      save ini, s, El, Eh, Q2min, Q2max, xmin, xmax, ymin, ymax

      real*8 powheginput,sqrts,theta,phi,cost,sint
      real*8 yaxis(3),zaxis(3),q(0:3)
      parameter (yaxis = (/0d0,1d0,0d0/))
      parameter (zaxis = (/0d0,0d0,1d0/))


C -   N.B: To check the phase space volume is correct,
C -   replace res(j) in sigborn.f by res(j)=kn_jacborn
C -   then check that the cross section is equal to,
C -   hc^2/(8*pi) = 1.54929*10^7.
c~       print*, 'xborn', xborn
      ! Initialise stuff
      if(ini) then
C -   Set initial- and final-state masses for Born and real
         do k=1,nlegborn
            kn_masses(k)=0
         enddo
         kn_masses(nlegreal)=0

!     Read in bounds on x,y,Q2
         Q2min = (powheginput('Qmin'))**2
         Q2max = (powheginput('Qmax'))**2
         xmin = powheginput('xmin')
         xmax = powheginput('xmax')
         ymin = powheginput('ymin')
         ymax = powheginput('ymax')
         
!     Set centre of mass energy
         El = kn_beams(0,1)
         Eh = kn_beams(0,2)
         s = 4d0 * Eh * El

!     Sanity checks
         if(Q2min.lt.0d0) stop 'Q2min negative'
         if(xmin.lt.0d0) stop 'xmin negative'
         if(ymin.lt.0d0) stop 'ymin negative'
         if(Q2min.gt.Q2max) stop 'Q2min > Q2max'
         if(xmin.gt.xmax) stop 'xmin > xmax'
         if(ymin.gt.ymax) stop 'ymin > ymax'
         if(Q2max.gt.s) Q2max = s

         if(Q2min.eq.Q2max.and.xmin.eq.xmax.and.ymin.eq.ymax) stop
     $        'Only two variables can be constrained'

         if(s*xmax*ymax.lt.Q2min .or. s*xmin*ymin.gt.Q2max) stop
     $        'No phase space available'

         if(Q2min.eq.Q2max.and.xmin.eq.xmax) then
            ymin = Q2min/(s*xmin)
            ymax = ymin
         else
            if(xmin*(ymax*s) .lt. Q2min) xmin = Q2min/(ymax*s)
            if(xmax*(ymin*s) .gt. Q2max) xmax = Q2max/(ymin*s)
         endif
!         if(s*xmin*ymin .gt. Q2min) Q2min = s*xmin*ymin
!         if(s*xmax*ymax .gt. Q2max) Q2max = s*xmax*ymax

         print*, '****************************************'
         print*, 'Doing DIS with'
         print*, 'xmin, xmax', xmin, xmax
         print*, 'ymin, ymax', ymin, ymax
         print*, 'Q2min, Q2max', Q2min, Q2max
         print*, '****************************************'

         !read(*,*)
         ini=.false.
      endif
c~       xborn = (/ 0.1d0, 0.2d0, 0.3d0/)
c~       print*, 'xborn', xborn
!     Initial jacobian pre factor
c~       jac = 1d0/(16d0*pi)
      jac = 1d0/(32d0*pi**2)                     !!!! Phi = dX dOmega /(32pi^2)
!     Generate first momentum fraction x
      if(xmin.eq.xmax) then
         x = xmin
         jac = 1d0 * jac
!      elseif(xmin.lt.ymax) then
      elseif(xmin.lt.xmax) then
         x = xmin * exp(xborn(1)*log(xmax/xmin))
         jac = jac * log(xmax/xmin) * x         !!!! Include the jacobian for dX /dxborn(1)
      endif
      kn_sborn = x * s
      sqrts = sqrt(kn_sborn)
      kn_xb1 = 1d0              ! Lepton
      kn_xb2 = x                ! Quark
      
      kn_cmpborn = 0d0
!     Initial state lepton in CM frame
      kn_cmpborn(0,1) = 1d0
      kn_cmpborn(3,1) = 1d0
!     Initial state quark in CM frame
      kn_cmpborn(0,2) =  1d0
      kn_cmpborn(3,2) = -1d0
      
!     Now generate the lepton-quark system back to back in the CM frame
!     Through generation of y
!     Check available phase space for y/x
      ymn = ymin
      ymx = ymax

      if(ymin.lt. Q2min/(s*x)) ymn = Q2min/(s*x)
      if(ymax.gt. Q2max/(s*x)) ymx = Q2max/(s*x)


      if(ymn.eq.ymx) then
         y = ymn
         jac = 1d0 * jac
         if(Q2min .eq. Q2max) jac = jac/(x*s)   !! If I fix Q2, I need the jacobian dy/dQ2
      elseif(ymn.lt.ymx) then
         y = ymn * exp(xborn(2)*log(ymx/ymn))
         jac = jac * log(ymx/ymn) * y         !!!! Include the jacobian for dY /dxborn(2)     
      else
         print*, ymn,ymx,x,s,Q2min,Q2max
         stop 'No phase space available for y'
      endif
! For numerical stability
      if(Q2min.eq.Q2max) then
         Q2 = Q2min
      else
         Q2 = y * x * s
      endif
c~       print*, 'Rolled Q2 in Born phase space', Q2
      if(Q2.lt.Q2min) then
         print*, 'Q2 too low', Q2, Q2min
      elseif(Q2.gt.Q2max) then
         print*, 'Q2 too high', Q2, Q2max
      endif
!      print*, Q2
!     theta = acos(1-2d0*y)
      cost = 1-2d0*y
      sint = sqrt(1d0 - cost**2)
      jac = 2d0 * jac ! y = (1-cos(theta))/2   !!!! Include the jacobian for dcth /dY

      phi = xborn(3) * 2d0 * pi
      jac = jac * 2d0 * pi

!     Outgoing lepton in CM frame
      kn_cmpborn(0,3) = 1d0
      kn_cmpborn(1,3) = Cos(phi) * sint
      kn_cmpborn(2,3) = Sin(phi) * sint
      kn_cmpborn(3,3) = cost
!     Rotate around y and z axis
c~       call mrotate(yaxis,sint,cost,kn_cmpborn(1:3,3))
c~       call mrotate(zaxis,sin(phi),cos(phi),kn_cmpborn(1:3,3))
!      call mrotate(yaxis,sin(theta),cos(theta),kn_cmpborn(1:3,3))

!     Quark back to back
      kn_cmpborn(0,4) = 1d0
      kn_cmpborn(1:3,4) = -kn_cmpborn(1:3,3)
!     Now rescale to total CM energy
      kn_cmpborn = sqrts * 0.5d0 * kn_cmpborn

!     Boost into lab frame
      beta=(El-x*Eh)/(x*Eh+El)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegborn,vec,beta,kn_cmpborn(0,1),kn_pborn(0,1))

      q = kn_pborn(:,1) - kn_pborn(:,3)
      call mlab2breit2(4,q,kn_pborn,kn_bfpborn,.False.)

      kn_jacborn = jac

c$$$      call checkmomzero(nlegborn,kn_pborn)
c$$$      if(x.lt.xmin.or.x.gt.xmax) print*, 'x out of bounds', x, xmax, xmin
c$$$      if(y.lt.ymin.or.y.gt.ymax) print*, 'y out of bounds', y, ymax, ymin
c$$$      if(Q2.lt.Q2min.or.Q2.gt.Q2max) print*, 'Q2 out of bounds', Q2, Q2max, Q2min
c$$$      print*, 'DIS variables x, y, Q', x, y, sqrt(Q2)
c$$$      print*, 'Incoming particles in lab frame'
c$$$      print*, kn_pborn(:,1)
c$$$      print*, kn_pborn(:,2)
c$$$      print*, kn_pborn(:,3)
c$$$      print*, kn_pborn(:,4)
c$$$      print*, ''
c$$$      print*, 'Incoming particles in CM frame'
c$$$      print*, kn_cmpborn(:,1)
c$$$      print*, kn_cmpborn(:,2)
c$$$      print*, kn_cmpborn(:,3)
c$$$      print*, kn_cmpborn(:,4)
c$$$      print*, ''
c$$$      print*, 'kn_jacoborb=', kn_jacborn
c$$$      print*, 'Event generated successfully'
c$$$  read(*,*)
      end



      subroutine born_suppression(fact)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real *8 fact
      real *8, save :: q2suppr
      real *8 q2
      real *8, external :: dotp, powheginput
      if(ini) then
         ini = .false.
         q2suppr = max(powheginput("#bornsuppfact"), 0d0)
         q2suppr = q2suppr**2
      endif
      q2 = 2d0*dotp(kn_cmpborn(:,1), kn_cmpborn(:,3))
      fact = q2/(q2 + q2suppr)
      end

      subroutine regular_suppression(fact)
      implicit none
      real * 8 fact
      call rmn_suppression(fact)
      end


      subroutine global_suppression(c,fact)
      implicit none
      character * 1 c
      real * 8 fact
      fact=1d0
      end

      subroutine rmn_suppression(fact)
      implicit none    
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      logical ini
      data ini/.true./
      real *8 fact
      real *8, save :: q2suppr
      real *8 q2
      real *8, external :: dotp, powheginput
      if(ini) then
         ini = .false.
         q2suppr = max(powheginput("#bornsuppfact"), 0d0)
         q2suppr = q2suppr**2
      endif
      q2 = 2d0*dotp(kn_cmpborn(:,1), kn_cmpborn(:,3))
      fact = q2/(q2 + q2suppr)
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'PhysPars.h'
      include 'pwhg_st.h'
      include 'pwhg_flg.h'
      real * 8 muf,mur
      real * 8 q(0:3), Q2
      logical, save :: runningscales, ini
      real * 8, external:: powheginput
      data ini/.true./
      if (ini) then
      runningscales = (powheginput('#runningscales').ne.0)
      endif
      if(runningscales) then
         if ((flg_btildepart.eq.'r').or.(flg_btildepart.eq.'R')) then
            q = kn_preal(:,1) - kn_preal(:,3)
         else
            q = kn_pborn(:,1) - kn_pborn(:,3) 
         endif
         Q2 = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2
         Q2 = -Q2
         muf=max(sqrt(Q2),1d0)
         mur=muf
      else
         muf = ph_Zmass
         mur = muf
      endif
      end


      subroutine born_phsp_silvia(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'

      real * 8 xborn(ndiminteg-4)
      real * 8 vec(3), beta
! DIS variables
      real * 8 x,Q2,s,y
      real * 8 jac, Eh, El
!     Limits on Q2 and x
      real * 8 Q2min, Q2max, xmin, xmax, ymin, ymax,ymn,ymx, cth, sth
      real * 8 Q2maxloc, xminloc
      integer mu,k,pow
      logical ini
      data ini/.true./
      save ini, s, El, Eh, Q2min, Q2max, xmin, xmax, ymin, ymax
      real*8, external:: powheginput, dotp

      jac = 1d0/(16d0*pi)

!     Generate Q2
      Q2min = 20d0
      Q2max = kn_sbeams
      Q2 = Q2min + (Q2max-Q2min) * xborn(1)
      jac = jac *(Q2max-Q2min)
    
!     Generate x      
      kn_xb1 = 1d0     
      xmin   = Q2/kn_sbeams
      xmax   = 1d0
      kn_xb2 = xmin +(xmax-xmin)*xborn(2)
      if(xmax > xmin) jac = jac *(xmax-xmin)

      
!     Partonic s
      kn_sborn = kn_sbeams * kn_xb2

      kn_pborn(:,1) = kn_beams(:,1)
      kn_pborn(:,2) = kn_xb2 * kn_beams(:,2)
      
      kn_cmpborn(:, 1) = sqrt(kn_sborn)/2d0 * (/ 1d0, 0d0, 0d0,  1d0 /)
      kn_cmpborn(:, 2) = sqrt(kn_sborn)/2d0 * (/ 1d0, 0d0, 0d0, -1d0 /)


!     ! Q^2 = s/2(1-cth) ---> (1-cth) = 2 Q^2/s--> cth = 1 - 2Q^2/s
     

      cth = 1d0 -2d0/kn_sborn * Q2
      sth = sqrt(max(1d0-cth**2,0d0))
      jac = jac * 2d0/kn_sborn !* Q2
      
      kn_cmpborn(:, 3) = sqrt(kn_sborn)/2d0 * (/ 1d0,  sth, 0d0,  cth /)
      kn_cmpborn(:, 4) = sqrt(kn_sborn)/2d0 * (/ 1d0, -sth, 0d0, -cth /)
      
      beta = (kn_pborn(0,1) - kn_pborn(0,2))/(kn_pborn(0,1) + kn_pborn(0,2))
      
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(1,vec,beta,kn_cmpborn(0,3),kn_pborn(0,3))
      call mboost(1,vec,beta,kn_cmpborn(0,4),kn_pborn(0,4))

      kn_jacborn = jac

      !print *, "jac =", kn_jacborn
c$$$
c$$$      print *, "check "
c$$$      print *, "cmp1 = ", kn_cmpborn(:,1), "m2 = ", dotp(kn_cmpborn(:,1),kn_cmpborn(:,1))
c$$$      print *, "cmp2 = ", kn_cmpborn(:,2), "m2 = ", dotp(kn_cmpborn(:,2),kn_cmpborn(:,2))
c$$$      print *, "cmp3 = ", kn_cmpborn(:,3), "m2 = ", dotp(kn_cmpborn(:,3),kn_cmpborn(:,3))
c$$$      print *, "cmp4 = ", kn_cmpborn(:,4), "m2 = ", dotp(kn_cmpborn(:,4),kn_cmpborn(:,4))
c$$$      print *, "cmpin  = ", kn_cmpborn(:,1)+kn_cmpborn(:,2)
c$$$      print *, "cmpout = ", kn_cmpborn(:,3)+kn_cmpborn(:,4)
c$$$      print *, "s = ", 2d0 * dotp(kn_cmpborn(:,1),kn_cmpborn(:,2))
c$$$      print *, "t = ", -2d0 * dotp(kn_cmpborn(:,1),kn_cmpborn(:,3))
c$$$      print *, "u = ", -2d0 * dotp(kn_cmpborn(:,3),kn_cmpborn(:,2))
c$$$
c$$$      print *, "p1 = ", kn_pborn(:,1)
c$$$      print *, "p2 = ", kn_pborn(:,2)
c$$$      print *, "p3 = ", kn_pborn(:,3)
c$$$      print *, "p4 = ", kn_pborn(:,4)
c$$$      print *, "pin  = ", kn_pborn(:,1)+kn_pborn(:,2)
c$$$      print *, "pout = ", kn_pborn(:,3)+kn_pborn(:,4)
c$$$      print *, "s = ", 2d0 * dotp(kn_pborn(:,1),kn_pborn(:,2))
c$$$      print *, "t = ", -2d0 * dotp(kn_pborn(:,1),kn_pborn(:,3))
c$$$      print *, "u = ", -2d0 * dotp(kn_pborn(:,3),kn_pborn(:,2))
      end
      
