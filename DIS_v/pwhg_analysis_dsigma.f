!     This analysis computes differential cross section for DIS in bins
!     of  x, log x , Q^2
!     This is a modification of the inclusive analysis

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_rwl.h'
      integer xnbins, Q2nbins
      
      real * 8 Q2binsize, Q2min, Q2max, logQ2min, logQ2max
      real * 8 xbinsize, logxmin, logxmax, sbeams
      call inihists
      call setupmulti(rwl_num_weights+1)

      ! bins in x
      xnbins = 10
      xbinsize = 1d0 / dble(xnbins)
      call bookupeqbins('dsigma_dx',xbinsize,0d0,1d0)
      
      ! bins in log(x)
      logxmax = log(1d0)
      logxmin = -10d0
      xbinsize = (logxmax - logxmin) / dble(xnbins)
      call bookupeqbins('dsigma_dlx',xbinsize,logxmin,logxmax)

      ! bins in Q^2
      Q2nbins = 12
      Q2binsize = 5000d0
      call bookupeqbins('dsigma_dQ2',Q2binsize,0d0,5d4)

      ! bins in log Q^2
      sbeams = 4d0 * ebmup(1) * ebmup(2)
      logQ2max = log(sbeams)
      logQ2min = log(1d0)
      Q2binsize = (logQ2max - logQ2min) / dble(Q2nbins)
      call bookupeqbins('dsigma_dlQ2',Q2binsize,logQ2min,logQ2max)
 
      end
     
      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'nlegborn.h'
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include 'LesHouches.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_flst.h'
      include 'pwhg_st.h'
      include 'pwhg_weights.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwl.h'
      include 'PhysPars.h'
      real * 8 dsig(weights_max)

      logical ini
      data ini/.true./
      save ini
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=128)
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      data WHCPRG/'NLO   '/
      integer nlep, npartons                         !lepton, quark, lepton initialstate, parton initialstate 
      real * 8 plephard(4),plep(4,maxtrack),pquark(4,maxtrack),plis(4),ppis(4),q(4)
!     momenta of lepton, quark, lepton initialstate, parton is, momentum transfere
      real * 8 plep03(0:3,maxtrack),pquark03(0:3,maxtrack),plis03(0:3),ppis03(0:3)
     $     ,q03(0:3) 
      real * 8 Q2,sbeams,x,y, delx !momentum transfer variable Q2=-q^2
      real * 8, parameter :: alphaem = 1d0/137d0 ! The lhef_analysis does not fill the PhysPars.h apparently

      integer i

      real * 8, external :: phepdot,getpseudorapidity0

      real *8 sigma_r

      dsig=0d0
      call multi_plot_setup(dsig0,dsig, weights_max)
!      print*, 'weights', dsig
      
      if(ini) then         
         if(whcprg.eq.'NLO') then
            write(*,*) '       NLO analysis'
         elseif(WHCPRG.eq.'LHE   ') then
            write(*,*) '       LHE analysis'
         elseif(WHCPRG.eq.'HERWIG') then
            write (*,*) '           HERWIG ANALYSIS            '
         elseif(WHCPRG.eq.'PYTHIA') then
            write (*,*) '           PYTHIA ANALYSIS            '
         endif
         ini=.false.
      endif

!      if (nhep == 4) return
      nlep = 0
      npartons = 0
      pquark = 0d0
      plep = 0d0
      ppis = 0d0
      plis = 0d0
      if(whcprg.eq.'NLO   '.or.whcprg.eq.'LHE   '.or.whcprg.eq.'PYTHIA')
     $     then
         do i=1,nhep
!     Initial states
            if(isthep(i).eq.-1.or.isthep(i).eq.21) then
               ! if(abs(idhep(i)).eq.11) then
               if(abs(idhep(i)).eq.12) then
                  plis(1:4) = phep(1:4,i)
               else if(abs(idhep(i)).le.5.or.idhep(i).eq.21) then
                  ppis(1:4) = phep(1:4,i)
               else
                  stop 'Wrong initial state'
               endif
!     Final states
            else if(isthep(i).eq.1) then
               if(abs(idhep(i)).eq.11) then
                  nlep = nlep + 1
                  plep(1:4,nlep) = phep(1:4,i)
                  ! For now everything not an electron is assumed a
                  ! parton. If we run with QED shower we have to be
                  ! careful
               else 
                  npartons = npartons + 1
                  pquark(1:4,npartons) = phep(1:4,i)
               endif
            endif
         enddo
      else
         print*, 'Wrong program', whcprg
         stop
      endif


      if(nlep.ne.1) then
         print*, nlep
         stop 'Wrong number of leptons'
      endif

      sbeams = 4d0 * ebmup(1) * ebmup(2)
      
      q(:) = plis(:) - plep(:,nlep)
      
      Q2 = phepdot(q,q)
      Q2 = -Q2

      if(Q2.lt.0d0 .or. Q2.gt.sbeams) then
         print*, 'Q2', Q2
         stop 'Q2 out of bounds'
      endif
      
      y = phepdot(ppis,q) / phepdot(ppis,plis)

      if(y.gt.1d0 .or. y.lt.0d0) stop 'y out of bounds'

      x = Q2 / (sbeams * y)

      if(x.gt.1d0 .or. x.lt.0d0) stop 'x out of bounds'


      ! dsigma / dx
      call filld('dsigma_dx', x, dsig)
      ! dsigma / dlx
      call filld('dsigma_dlx', log(x), dsig)      
      ! dsigma / dQ2
      call filld('dsigma_dQ2', Q2, dsig)
      ! dsigma / dlQ2
      call filld('dsigma_dlQ2', log(Q2), dsig)

      end



      function phepDot(p_A,p_B)
      implicit none
      real * 8  phepDot
      real * 8  p_A(4),p_B(4)
      phepDot=p_A(4)*p_B(4)-p_A(1)*p_B(1)
     1       -p_A(2)*p_B(2)-p_A(3)*p_B(3)
      end

      function kt2(p)
      implicit none
      real * 8 kt2, p(0:3)

      kt2 = p(1)**2 + p(2)**2
      end

      function eta(p)
      implicit none
      real * 8 eta, p(0:3), normp, norm

      normp = norm(p)
      if(normp.gt.p(3)) then
         eta = 0.5d0 * log((normp + p(3)) / (normp - p(3)))
      else
         eta = sign(1d100,p(3)) 
      endif
      end

      function norm(p)
      implicit none
      real * 8 norm, p(0:3)

      norm = p(1)**2 + p(2)**2 + p(3)**2
      norm = sqrt(max(0d0,norm))
      end

      function getrapidity0(p)
      implicit none
      real * 8 p(0:3),getrapidity0
      getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end
      
      function getpseudorapidity0(p)
      implicit none
      real * 8 p(0:3),getpseudorapidity0,norm      
      getpseudorapidity0=0.5d0*log((norm(p)+p(3))/(norm(p)-p(3)))
      end

      subroutine getrapidity(p,y)
      implicit none
      real * 8 p(4),y
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function azi(p)
      implicit none
      real * 8 pi
      parameter(pi = 3.141592653589793D0)
      real * 8 azi,p(0:3)
      azi = atan(p(2)/p(1))
      if (p(1).lt.0d0) then
         if (azi.gt.0d0) then               
            azi = azi - pi
         else
            azi = azi + pi
         endif
      endif    
      end
      

