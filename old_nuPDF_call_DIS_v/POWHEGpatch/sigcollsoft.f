c subroutines soft, collfsr, softcollfsr, collisrp, collisrm,
c             softcollisrp, softcollisrm 
c are called by btildereal; they return (1-y) csi^2 R_alr (if kn_emitter
c is >2, i.e. final state) or (1-y^2) csi^2 R (if kn_emitter <=2) in the
c collinear final state                      (collfsr)
c soft-collinear final state                 (softcollfsr)
c collinear initial state in the + direction (collisrp)
c collinear initial state in the - direction (collisrm)
c soft-collinear initial state +             (softcollisrp)
c soft-collinear initial state -             (softcollisrm)
c The soft limit is taken keeping y (i.e. kn_y) fixed
c The collinear limit is taken keeping csitilde (kn_csitilde) fixed
c
c The Born and real phase space must have been set before calling them.
c The subroutine fills the array argument rc(maxalr) with the alr
c contributions that have as emitter the current emitter kn_emitter.
c All the others are set to zero.
c
c
c Functions: collbtl and softbtl are required to compute the damping
c            factor in bornzerodamp. They are called by sigreal, and
c            sigreal is used by btilde (hence the btl suffix)

c softbtl: computes the same quantity as soft, but does not include 
c          the luminosity
c
c collbtl: if kn_emitter is in the final state, computes the same quantity as
c          collfsr, but does not include the luminosity
c          if kn_emitter is 1 or 2  computes the same quantity as
c          collisrp or collisrm, without luminosity
c          if kn_emitter is 0 computes the average of the return values
c          of collisrp and collisrm, weighted with the factors 
c          (1+y)/2 and (1-y)/2.
c
c 
c
c Functions: collrad, softrad  are required to compute the damping
c            factor in bornzerodamp. They are called by sigreal_rad,
c            i.e. they are used in the generation of radiation.
c            They compute the same quantities as collbtl and softbtl,
c            except that they only fill alr values corresponding to
c            the current radiation underlying Born (i.e. the alr list
c            rad_alr_list(1...rad_alr_nlist)) and that share the current
c            radiation kinematics (i.e. the rad_kinreg value). Furthermore
c            they are divided by kn_csi^2(1-kn_y^2) (for initial state
c            kinematic regions) or  kn_csi^2(1-kn_y) (for final state
c            kinematic regions).

      subroutine collfsr(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      real * 8 rc(maxalr)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac
      integer alr
      call collfsrnopdf(rc)
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            rc(alr)=rc(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end

      subroutine collfsrnopdf(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      real * 8 rc(maxalr)
      integer alr,em,kres
      real * 8 kperp(0:3),kperp2,q0,xocsi,csi,x,phi
      em=kn_emitter
      csi=kn_csi
      call buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            call collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,rc(alr))
         else
            rc(alr)=0
         endif
      enddo
      end

      subroutine buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
      implicit none
      integer em
      real * 8 q0,xocsi,x,kperp(0:3),kperp2
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      integer kres
      real * 8 pem(0:3),pres(0:3),phi
      real * 8 dotp
      external dotp
      kres=flst_bornrescurr(em)
      if(kres.eq.0) then
         pres=kn_cmpborn(:,1)+kn_cmpborn(:,2)
         pem=kn_cmpborn(:,em)
      else
         call boost2reson(kn_cmpborn(:,kres),1,kn_cmpborn(:,em),pem)
         pres(1:3)=0
         pres(0)=sqrt(dotp(kn_cmpborn(:,kres),kn_cmpborn(:,kres)))
      endif
      q0=pres(0)
      xocsi=q0/2/pem(0)
c for fsr csi is y independent
      phi=kn_azi
      x=kn_csi*xocsi
      call buildkperp(phi,pem,kperp,kperp2)
c The correctness of this has not been tested yet in any specific process
c (typically a process with a gluon splitting in resonance decay).
c The question is: should the kperp vector be boosted back from
c the resonance frame to the CM frame?
c      if(kres.ne.0) then
c         call boost2reson(kn_cmpborn(:,kres),1,kperp,kperp)
c      endif
      end

      subroutine buildkperp(phi,pem,kperp,kperp2)
      implicit none
      real * 8 phi,pem(0:3),kperp(0:3),kperp2,dir(1:3)
c     Construct kperp; First construct a vector in the plane of p_em and
c     the third axis, orthogonal to p_em.
c     Then rotate it counterclockwise around the em direction
      kperp(1)=pem(1)
      kperp(2)=pem(2)
      kperp(3)=-(pem(2)**2+pem(1)**2)
     #/pem(3)
      kperp2=(kperp(1)**2+kperp(2)**2+kperp(3)**2)
      dir(1)=pem(1)/pem(0)
      dir(2)=pem(2)/pem(0)
      dir(3)=pem(3)/pem(0)
      call mrotate(dir,sin(phi),cos(phi),kperp(1))
      kperp(0)=0
      end

      subroutine softcollfsr(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      real * 8 rc(maxalr)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      integer alr
      real * 8 tmp,rescfac
      tmp=kn_csi
      kn_csi=0
      call collfsrnopdf(rc)
      kn_csi=tmp
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            rc(alr)=rc(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end



      subroutine collisrp(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 rc(maxalr)
      call collisr(1,rc)
      end

      subroutine collisrm(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 rc(maxalr)
      call collisr(2,rc)
      end

      subroutine softcollisrp(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      real * 8 rc(maxalr)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      integer alr
      real * 8 tmp,rescfac
      tmp=kn_csip
      kn_csip=0
      call collisrnopdf(1,rc)
      kn_csip=tmp
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            rc(alr)=rc(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end

      subroutine softcollisrm(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      real * 8 rc(maxalr)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      integer alr
      real * 8 tmp,rescfac
      tmp=kn_csim
      kn_csim=0
      call collisrnopdf(2,rc)
      kn_csim=tmp
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            rc(alr)=rc(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end

      subroutine collisr(i,rc)
      implicit none
      integer i
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      include 'pwhg_kn.h'
      real * 8 rc(maxalr)
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton),
     1     x1,x2
      integer alr
      real * 8 rescfac
      call collisrnopdf(i,rc)
      if(i.eq.1) then
         x1=kn_xb1/(1-kn_csip)
         x2=kn_xb2
      elseif(i.eq.2) then
         x1=kn_xb1
         x2=kn_xb2/(1-kn_csim)
      endif
      if(.not.flg_minlo) then
         call pdfcall(1,x1,pdf1)
         call pdfcall(2,x2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,x1,pdf1)
               call pdfcall(2,x2,pdf2)
            endif
            rc(alr)=rc(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end

      subroutine collisrnopdf(i,rc)
      implicit none
      integer i
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      real * 8 rc(maxalr)
      integer alr
      real * 8 kperp(0:3),csi,phi
      if(i.eq.1) then
         csi=kn_csip
      else
         csi=kn_csim
      endif
      phi=kn_azi
c Construct kperp
      kperp(1)=sin(phi)
      kperp(2)=cos(phi)
      kperp(3)=0
      kperp(0)=0
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            call collisralr(alr,i,csi,kperp,rc(alr))
         else
            rc(alr)=0
         endif
      enddo
      end



      subroutine soft(r0)
c blegs:        integer, number of legs of born
c bflav(nlegs): integer, flavours of the incoming partons, according to PDG conventions,
c               MODIFIED TO HAVE 0 FOR GLUONS (instead of 21)
c p(0:3,nleg):  real * 8, momenta, 0 is time component
c softvec(0:3): real * 8, 4-vector of soft gluon normalized to softvec(0)=1
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_br.h'
      include 'pwhg_pdf.h'
      real * 8 r0(maxalr)
      integer alr,em
      real * 8 y
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac
c Boost Born momenta to their rest frame
c find boost velocity
      em=kn_emitter
      y=kn_y
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1        flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            call softalr(alr,em,y,r0(alr))
         else
            r0(alr) = 0
         endif
      enddo
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1        flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            if(flg_minlo) then
               call setlocalscales(flst_alr2born(alr),2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            r0(alr)=r0(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1           *rescfac
         endif
      enddo         
      end


      subroutine collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_br.h'
      include 'pwhg_par.h'
      include 'pwhg_em.h'
      integer alr
      real * 8 csi,xocsi,x,q0,kperp(0:3),kperp2,res
      integer lalr,iub,em,emflav,raflav,mu,nu,kres
      real * 8 gtens(0:3,0:3),ap,q2
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      real * 8 chargeofparticle
      external chargeofparticle
      logical is_em
      is_em = .false.
      lalr = flst_alrlength(alr)
      iub=flst_alr2born(alr)
      em=flst_emitter(alr)
c no collinear subtraction for massive particle
      if(kn_masses(em).gt.0) then
         res=0
         return
      endif
      kres=flst_bornres(em,iub)
      if(kres.eq.0) then
         q2=kn_cmpborn(0,em)**2
      else
         q2= ( kn_cmpborn(0,kres)*kn_cmpborn(0,em)
     1        -kn_cmpborn(1,kres)*kn_cmpborn(1,em)
     2        -kn_cmpborn(2,kres)*kn_cmpborn(2,em)
     3        -kn_cmpborn(3,kres)*kn_cmpborn(3,em) )**2/
     4            (   kn_cmpborn(0,kres)**2
     5               -kn_cmpborn(1,kres)**2
     6               -kn_cmpborn(2,kres)**2
     7               -kn_cmpborn(3,kres)**2 )
      endif
      emflav=flst_alr(em,alr)
      raflav=flst_alr(lalr,alr)
      if(emflav.eq.0.and.raflav.eq.0) then
         ap=0
         do mu=0,3
            do nu=0,3
               ap=ap+(-gtens(mu,nu)*(csi*x/(1-x)+(1-x)/xocsi)
     #+2*x*(1-x)*csi*
     #(kperp(mu)*kperp(nu))/kperp2)*br_bmunu(mu,nu,em,iub)
            enddo
         enddo
         ap=ap*2*ca
c     In case of two equal gluon we also supply e E_em^p/(E_em^p+E_rad^p)
c     factor, and divide by 2 for the identical particles
c         ap=ap*(1-x)**par_2gsupp/((1-x)**par_2gsupp+x**par_2gsupp)
c Commented out: now this is done at the end of the if block
      elseif(emflav+raflav.eq.0) then
         ap=0d0
         do mu=0,3
            do nu=0,3
               ap=ap+(-gtens(mu,nu)
     #-4*x*(1-x)*(kperp(mu)*kperp(nu))/kperp2)*br_bmunu(mu,nu,em,iub)
            enddo
         enddo
         ap=ap*tf*csi
      elseif(raflav.eq.0.and.emflav.ne.0) then
         ap=cf*(1+(1-x)**2)/xocsi*br_born(iub)
      elseif(raflav.ne.0.and.emflav.eq.0) then
         ap=cf*(1+x**2)/(1-x)*csi*br_born(iub)
      elseif(raflav.eq.22.and.emflav.ne.0) then
         is_em = .true.
         ap=(1+(1-x)**2)/xocsi*br_born(iub)
     1        *chargeofparticle(emflav)**2
      elseif(raflav.ne.0.and.emflav.eq.22) then
         is_em = .true.
         ap=(1+x**2)/(1-x)*csi*br_born(iub)
     1        *chargeofparticle(raflav)**2
      else
         write(*,*) 'coll (fsr): unammissible flavour structure'
         call pwhg_exit(-1)
      endif
      if(flg_doublefsr.or.(emflav.eq.0.and.raflav.eq.0)) then
         ap=ap*(1-x)**par_2gsupp/((1-x)**par_2gsupp+x**par_2gsupp)
      endif
c     1/(p_em . p_ra) = 1/(p_bar_em(0,em)**2* x * (1-x) * (1-y);
c     we multiply everything by (1-y) csi^2; one csi is included
c     above; the other here.
      res=ap/q2/(xocsi*(1-x))
      if(is_em) then
         res = res * (4*pi*em_alpha)
      else
         res = res * (4*pi*st_alpha)
      endif
c provide multiplicity of emitter in underlyng Born
      res=res*flst_ubmult(alr)
      end


      subroutine collisralr(alr,i,csi,kperp,res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_br.h'
      integer alr,i,em
      real * 8 csi,kperp(0:3),res
      integer lalr,iub,emflav,raflav,mu,nu
      real * 8 gtens(0:3,0:3),ap,x
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      real * 8 chargeofparticle
      logical is_em
      external chargeofparticle
      lalr = flst_alrlength(alr)
      is_em=.false.
      x=1-csi
      iub=flst_alr2born(alr)
      em=flst_emitter(alr)
      emflav=flst_alr(i,alr)
      raflav=flst_alr(lalr,alr)
      if(emflav.eq.0.and.raflav.eq.0) then
         ap=0
         do mu=0,3
            do nu=0,3
               ap=ap+(-gtens(mu,nu)*(x+x*(1-x)**2)
     #+2*(1-x)**2/x*
     #(kperp(mu)*kperp(nu)))*br_bmunu(mu,nu,i,iub)
            enddo
         enddo
         ap=ap*2*ca
      elseif(raflav-emflav.eq.0) then
         ap=0d0
         do mu=0,3
            do nu=0,3
               ap=ap+(-gtens(mu,nu)*x
     #+4*(1-x)/x*(kperp(mu)*kperp(nu)))*br_bmunu(mu,nu,i,iub)
            enddo
         enddo
c gluon radiation case
         if(flst_born(i,iub).eq.0) then
            ap=ap*cf*(1-x)
c photon radiation case
         elseif(flst_born(i,iub).eq.22) then
            is_em=.true.
            ap=ap*(1-x)*chargeofparticle(emflav)**2
         else
            write(*,*) ' collisralr: unammissible ub flavour'
            call pwhg_exit(-1)
         endif
      elseif(raflav.eq.0.and.emflav.ne.0) then
         ap=cf*(1+x**2)*br_born(iub)
      elseif(emflav.eq.0.and.raflav.ne.0) then
         ap=tf*(x**2+(1-x)**2)*br_born(iub)*(1-x)
      elseif(emflav.ne.0.and.raflav.eq.22) then
c QED, consider photon emission
         is_em=.true.
         ap=(1+x**2)*br_born(iub)*chargeofparticle(emflav)**2
      elseif(emflav.eq.22.and.raflav.ne.0) then
c QED, consider emission by photon 
         is_em=.true.
         ap=(x**2+(1-x)**2)*br_born(iub)*(1-x)
     +        *3d0*chargeofparticle(raflav)**2
         if(abs(raflav)>6) then
            ap = ap/3
         endif
      else
         write(*,*) 'coll (isr): unammissible flavour structure'
      endif
c     In the real CM frame:
c     1/(p_ra . p_i)=1/(p^0_ra*p^0_i*(1-y))=1/(p^0_i)^2 /[(1-x)*(1-y)]
c     where 1/(p^0_i)^2=1/(p^0_1 * p^0_2)=1/(pborn^0_1 * pborn^0_2/x),
c     the last expression being boost invariant.
c     Supplying che csi^2 (1-y^2) factor in the collinear limit,
c     using csi=1-x we get
      if(is_em) then
c QED, electromagnetic coupling
         res=ap/(kn_pborn(0,1)*kn_pborn(0,2)/x) * 2
c     ew coupling:
     # *(4*pi*em_alpha)
      else
         res=ap/(kn_pborn(0,1)*kn_pborn(0,2)/x) * 2
c     strong coupling:
     # *(4*pi*st_alpha)
c     The remaining csi=1-x factor has been applied earlier
      endif
      end


      function colcorr(j,iub,res)
c     If the flag flg_zerowidth is .true.,
c     it returns true if parton j, in the underlying Born flavour
c     structure iub, belongs to a group of colour correlated particles
c     arising from the decay of the resonance res. This group is formed
c     by all coloured particles that are sons of the resonance res
c     (including eventually other coloured resonances) and by the resonance
c     itself if coloured.
c     The case res=0 corresponds to the partons produced promptly in the
c     hard process, included coloured resonances produced promptly.
c     If flg_zerowidth is .true., only initial and final state partons can
c     be colour correlated.
      implicit none
      logical colcorr
      integer iub,res,j
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      logical is_coloured,isfinalstate
      colcorr= is_coloured(flst_born(j,iub))
      if(flg_zerowidth) then
c Look out: if res = 0 the condition flst_bornres(j,iub).eq.res holds
c also for initial state partons.
         colcorr=colcorr .and.
     1        (j.eq.res .or. flst_bornres(j,iub).eq.res)
      else
c Should be coloured, and initial or final state
         colcorr=colcorr .and. (j.le.2 .or.
     1        isfinalstate(j,flst_bornlength(iub),flst_bornres(:,iub)) )
      endif
      end



      function colcorrem(j,iub,em)
c Returns true if parton j, in the underlying Born flavour structure
c iub, belongs to a group of colour correlated partons relevant for the
c emitter em. This is essentially as the colcorr function, except that
c it deals with the special case em=0 (that in POWHEG means radiation
c from either initial state partons)
      implicit none
      logical colcorrem
      integer iub,em,res,j
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      logical colcorr
      if(em.eq.0) then
         res=0
      else
         res=flst_bornres(em,iub)
      endif
      colcorrem=colcorr(j,iub,res)
      end


      subroutine softalr(alr,em,y,res)
      implicit none
      integer alr,em
      real * 8 y,res
c last argument=1: traditional POWHEG BOX behaviour;
c kn_softvec is assume to have 0th component = 1 in the
c resonance rest frame
      call softalr0(alr,em,y,res,1)
      end

      subroutine softalr0(alr,em,y,res,imode)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_br.h'
      integer alr,em,imode
      real * 8 y,res,chj,chk
      integer iub,j,k,kres,alrlen,bornlen
      real * 8 pjsq,sumdijinv,result,q2,rresweight,bresweight
      real * 8 dotp,chargeofparticle
      integer sg_j,sg_k
      logical colcorrem
      external dotp,chargeofparticle,colcorrem
      alrlen = flst_alrlength(alr)
      bornlen = alrlen - 1
      if(flst_emitter(alr).eq.em) then
         if(flst_alr(alrlen,alr).eq.0) then
            iub=flst_alr2born(alr)
            if(imode.eq.1) then
               if(em.eq.0) then
                  kres=0
               else
                  kres=flst_bornres(em,iub)
               endif
               if(kres.eq.0) then
                  q2=4*kn_cmpborn(0,1)*kn_cmpborn(0,2)
               else
                  q2=abs(kn_cmpborn(0,kres)**2
     1                 -kn_cmpborn(1,kres)**2
     2                 -kn_cmpborn(2,kres)**2
     3                 -kn_cmpborn(3,kres)**2)
               endif
            else
c if imode=0, kn_softvec is assumed to be the dimensionful soft momentum,
c and the return value (result) is not multiplied by csi^2.
c No 4/q2 factor needed, just set it to 1 (look for 4/q0^2 comment below)
               q2 = 4
            endif
c     loop over pairs of coloured particles
            result=0
            do j=1,bornlen
               if(colcorrem(j,iub,em)) then
                  do k=j+1,bornlen
                     if(colcorrem(k,iub,em)) then
                        result=result+
     1                       dotp(kn_cmpborn(0,j),kn_cmpborn(0,k))
     2                       /(dotp(kn_cmpborn(0,j),kn_softvec)*
     3                       dotp(kn_cmpborn(0,k),kn_softvec)) * br_bornjk(j,k,iub)
                     endif
                  enddo
               endif
            enddo
c     the previous sum should run over all indexes. Since br_bornjk is
c     symmetric, multiply by 2
            result = result*2
            do j=1,bornlen
               if(colcorrem(j,iub,em) .and.
     1              flst_born(j,iub).ne.0) then
                  pjsq=kn_cmpborn(0,j)**2-kn_cmpborn(1,j)**2-
     1                 kn_cmpborn(2,j)**2-kn_cmpborn(3,j)**2
                  result=result-pjsq/dotp(kn_cmpborn(0,j),
     1                 kn_softvec)**2*br_born(iub)*cf
               endif
            enddo
c     having chosen a soft four momentum of energy 1, supply
c     1/esoft^2. Multiply by csi^2, csi=2*esoft/q0, so net
c     factor 4/q0^2
            if(em.gt.2) then
               res=result*4/q2*(1-y)
            else
               res=result*4/q2*(1-y**2)
            endif
c     Coupling:
            res=res*(4*pi*st_alpha)
         elseif(flst_alr(alrlen,alr).eq.22) then
c QED: soft photon emission
            iub=flst_alr2born(alr)
            if(imode.eq.1) then
               if(em.eq.0) then
                  kres=0
               else
                  kres=flst_bornres(em,iub)
               endif
               if(kres.eq.0) then
                  q2=4*kn_cmpborn(0,1)*kn_cmpborn(0,2)
               else
                  q2=abs(kn_cmpborn(0,kres)**2
     1                 -kn_cmpborn(1,kres)**2
     2                 -kn_cmpborn(2,kres)**2
     3                 -kn_cmpborn(3,kres)**2)
               endif
            else
c if imode=0, kn_softvec is assumed to be the dimensionful soft momentum,
c and the return value (result) is not multiplied by csi^2.
c No 4/q2 factor needed, just set it to 1 (look for 4/q0^2 comment below)
               q2 = 4
            endif
c     loop over pairs of charged particles
            result=0
            do j=1,bornlen
               chj=chargeofparticle(flst_born(j,iub))
               if( j > 2 ) then
                  sg_j=-1
               else
                  sg_j=1
               endif
               
               if(chj.ne.0) then
                  do k=j+1,bornlen
                     chk=chargeofparticle(flst_born(k,iub))
                     if( k > 2 ) then
                        sg_k=-1
                     else
                        sg_k=1
                     endif
                     if(chk.ne.0) then
                        result= result+
     1                          dotp(kn_cmpborn(0,j),kn_cmpborn(0,k))
     2                         /( dotp(kn_cmpborn(0,j),kn_softvec)*
     3                            dotp(kn_cmpborn(0,k),kn_softvec) )
     4                           * chj * chk * (-1) * sg_j * sg_k
                     endif
                  enddo
               endif
            enddo
c     symmetric, multiply by 2
            result = result*2
            do j=1,bornlen
               chj=chargeofparticle(flst_born(j,iub))
               if(chj.ne.0.and.kn_masses(j).ne.0) then
                  pjsq   = dotp(kn_cmpborn(0,j),kn_cmpborn(0,j)) 
                  result = result
     +                    -pjsq/dotp(kn_cmpborn(0,j),kn_softvec)**2
     +                     *chj**2
               endif
            enddo
            result=result * br_born(iub)
c     having chosen a soft four momentum of energy 1, supply
c     1/esoft^2. Multiply by csi^2, csi=2*esoft/q0, so net
c     factor 4/q0^2
            if(em.gt.2) then
               res=result*4/q2*(1-y)
            else
               res=result*4/q2*(1-y**2)
            endif
c     Coupling:
            res=res*(4*pi*em_alpha)
c     The case of the emitter being a gluon requires no special treatment here!
c     the extra (1-x) factor is simply 1!
         else
            res=0
         endif
      else
         res=0
      endif
c     Multiply soft result by (soft limit Sij appropriate factors)
      if(res.ne.0) then
         call realresweight_soft(alr,rresweight)
         call bornresweight(flst_alr2born(alr),bresweight)
         res = res * rresweight/bresweight
      endif
      res=res*flst_ubmult(alr)
      end



      subroutine softbtl(r0)
c blegs:        integer, number of legs of born
c bflav(nlegs): integer, flavours of the incoming partons, according to PDG conventions,
c               MODIFIED TO HAVE 0 FOR GLUONS (instead of 21)
c p(0:3,nleg):  real * 8, momenta, 0 is time component
c softvec(0:3): real * 8, 4-vector of soft gluon normalized to softvec(0)=1
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_br.h'
      real * 8 r0(maxalr)
      integer alr,em
      real * 8 y
c Boost Born momenta to their rest frame
c find boost velocity
      em=kn_emitter
      y=kn_y
      do alr=1,flst_nalr
         if(flst_emitter(alr) == kn_emitter .and.
     1      flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
            call softalr(alr,em,y,r0(alr))
         else
            r0(alr) = 0
         endif
      enddo
      end


c This returns in rc the collinear approximation to
c the real cross section (multiplied by csi^2(1-y^2) or (1-y))
c to be used to construct the damping factor in the real
c cross section used in Btilde
      subroutine collbtl(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 rc(maxalr)
      integer alr,em,i
      real * 8 kperp(0:3),kperp2,q0,xocsi,csi,x,phi,r1,r2
      em=kn_emitter
      if(em.gt.2) then
c     for fsr csi is y independent
         csi=kn_csi
         call buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
         do alr=1,flst_nalr
            if(flst_emitter(alr) == kn_emitter .and.
     1           flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
               call collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,rc(alr))
            else
               rc(alr)=0
            endif
         enddo
      else
         phi=kn_azi
c Construct kperp
         kperp(1)=sin(phi)
         kperp(2)=cos(phi)
         kperp(3)=0
         kperp(0)=0
         do alr=1,flst_nalr
            if(flst_emitter(alr) == kn_emitter .and.
     1           flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
               if(em.ne.2) then
                  print*, 'em =', em
                  stop 'em should only be 2 in collrad'
                  i=1
                  csi=kn_csi*kn_csimaxp/kn_csimax
                  call collisralr(alr,i,csi,kperp,r1)
               endif
               if(em.ne.1) then
                  i=2
                  csi = kn_csitilde * kn_csimaxm
                  ! csi=kn_csi*kn_csimaxm/kn_csimax
                  call collisralr(alr,i,csi,kperp,r2)
               endif
               if(em.eq.0) then
                  rc(alr)=(r1*(1+kn_y)+r2*(1-kn_y))/2
               elseif(em.eq.1) then
                  rc(alr)=r1
               elseif(em.eq.2) then
                  rc(alr)=r2
               endif
            else
               rc(alr)=0
            endif
         enddo
      endif
      end

c This returns in rcs the soft-collinear approximation to
c the real cross section (multiplied by csi^2(1-y^2) or (1-y))
c to be used to construct the damping factor in the real
c cross section used in Btilde
      subroutine collsoftbtl(rcs)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      real * 8 rcs(maxalr)
      integer alr,em,i
      real * 8 kperp(0:3),kperp2,q0,xocsi,csi,x,phi,r1,r2,kncsisave

      kncsisave = kn_csi
      kn_csi = 0
      em=kn_emitter
      if(em.gt.2) then
c     for fsr csi is y independent
         csi=kn_csi
         call buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
         do alr=1,flst_nalr
            if(flst_emitter(alr) == kn_emitter .and.
     1           flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
               call collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,rcs(alr))
            else
               rcs(alr)=0
            endif
         enddo
      else
         phi=kn_azi
c Construct kperp
         kperp(1)=sin(phi)
         kperp(2)=cos(phi)
         kperp(3)=0
         kperp(0)=0
         do alr=1,flst_nalr
            if(flst_emitter(alr) == kn_emitter .and.
     1           flst_bornresgroup(flst_alr2born(alr)) == flst_ibornresgroup  ) then
               if(em.ne.2) then
                  i=1
!     kn_csitilde is zero: csi = 0   (ignore the branch cut // assume it is always -)
                  csi=0d0
                  call collisralr(alr,i,csi,kperp,r1)
               endif
               if(em.ne.1) then
                  i=2
!     kn_csitilde is zero: csi = 0   (ignore the branch cut // assume it is always - )
                  csi=0d0
                  call collisralr(alr,i,csi,kperp,r2)
               endif
               if(em.eq.0) then
                  rcs(alr)=(r1*(1+kn_y)+r2*(1-kn_y))/2
               elseif(em.eq.1) then
                  rcs(alr)=r1
               elseif(em.eq.2) then
                  rcs(alr)=r2
               endif
            else
               rcs(alr)=0
            endif
         enddo
      endif
      kn_csi = kncsisave
      end




c This returns in rc the collinear approximation to
c the real cross section to be used to construct the damping factor
c in the real radiation cross section
      subroutine collrad(rc)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      real * 8 rc(flst_nalr)
      integer alr,em
      real * 8 q0,xocsi,csi,phi,x,kperp(0:3),kperp2,r1,r2
      integer j
      rc=0
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         em=flst_emitter(alr)

         kn_csitilde = kn_csi / kn_csimax
         
c     check if emitter corresponds to current radiation region (i.e. rad_kinreg):
         if(rad_kinreg.eq.1.and.em.le.2) then
            if (kn_branchfactor > 0 .and. flg_dis_isr ) then
               kn_csitilde = (2d0*kn_csimax_minus-kn_csi)/kn_csimax
            else
               kn_csitilde = kn_csi / kn_csimax
            endif
            phi=kn_azi
c     Construct kperp
            kperp(1)=sin(phi)
            kperp(2)=cos(phi)
            kperp(3)=0
            kperp(0)=0
            if(em.ne.2) then
               csi=kn_csitilde*kn_csimaxp               
               call collisralr(alr,1,csi,kperp,r1)
            endif

            if(em.ne.1) then
               ! csi=kn_csi*kn_csimaxm/kn_csimax
               csi=kn_csitilde * kn_csimaxm
               call collisralr(alr,2,csi,kperp,r2)
            endif
    
            if(em.eq.0) then
               rc(alr)=(r1*(1+kn_y)+r2*(1-kn_y))/2
            elseif(em.eq.1) then
               rc(alr)=r1
            elseif(em.eq.2) then
               rc(alr)=r2
            endif
            
            rc(alr)=rc(alr)/((kn_csitilde*kn_csimax)**2*(1-kn_y**2))
           
         elseif(flst_lightpart+rad_kinreg-2.eq.em) then
            csi=kn_csi
            call buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
            call collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,rc(alr))
            rc(alr)=rc(alr)/csi**2/(1-kn_y)
         else
            rc(alr)=0
         endif
      enddo
      end


      subroutine softrad(r0)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      real * 8 r0(flst_nalr)
      integer alr,em,j
      real * 8 y,csi
      y=kn_y
      r0=0
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         em=flst_emitter(alr)
         csi = kn_csi      
         if(rad_kinreg.eq.1.and.em.le.2) then
            if (kn_branchfactor > 0 .and. flg_dis_isr) csi = 2d0*kn_csimax_minus-csi 
            call softalr(alr,em,y,r0(alr))
            r0(alr)=r0(alr)/csi**2/(1-y**2)
         elseif(flst_lightpart+rad_kinreg-2.eq.em) then
            call softalr(alr,em,y,r0(alr))
            r0(alr)=r0(alr)/csi**2/(1-y)
         else
            r0(alr)=0
         endif
      enddo
      end



c This returns in rcs the soft-collinear approximation to
c the real cross section to be used to construct the damping factor
c in the real radiation cross section
      subroutine collsoftrad(rcs)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      real * 8 rcs(flst_nalr)
      integer alr,em
      real * 8 q0,xocsi,csi,phi,x,kperp(0:3),kperp2,r1,r2,kncsisave
      integer j
      kncsisave = kn_csi
      rcs = 0
      
      do j=1,rad_alr_nlist
         csi=0d0
         kn_csi=0d0

         alr=rad_alr_list(j)
         em=flst_emitter(alr)

         
!     kn_csi has been generated, and stored in kncsisave,
!     but for ISR we need kn_csitilde, whose expression
!     depends on the branch cut
         if (kn_branchfactor > 0 .and. em==2 .and. flg_dis_isr) then
            kn_csitilde = (2d0*kn_csimax_minus-kncsisave)/kn_csimax
         else
            kn_csitilde = kncsisave/kn_csimax
         endif
         
c     check if emitter corresponds to current radiation region (i.e. rad_kinreg):
         if(rad_kinreg.eq.1.and.em.le.2) then
            phi=kn_azi
c     Construct kperp
            kperp(1)=sin(phi)
            kperp(2)=cos(phi)
            kperp(3)=0
            kperp(0)=0
            if(em.ne.2) then
               call collisralr(alr,1,csi,kperp,r1)
            endif
            if(em.ne.1) then
               call collisralr(alr,2,csi,kperp,r2)
            endif
            if(em.eq.0) then
               rcs(alr)=(r1*(1+kn_y)+r2*(1-kn_y))/2
            elseif(em.eq.1) then
               rcs(alr)=r1
            elseif(em.eq.2) then
               rcs(alr)=r2
            endif
            !use kn_csitilde * kn_csimax instead of kncsisave
            rcs(alr)=rcs(alr)/((kn_csitilde*kn_csimax)**2*(1-kn_y**2))
         elseif(flst_lightpart+rad_kinreg-2.eq.em) then
            call buildfsrvars(em,q0,xocsi,x,kperp,kperp2)
            call collfsralr(alr,csi,xocsi,x,q0,kperp,kperp2,rcs(alr))
            rcs(alr)=rcs(alr)/kncsisave**2/(1-kn_y)
         else
            rcs(alr)=0
         endif
      enddo

      
      kn_csi = kncsisave
      end

      subroutine softmismatchres(xrad,res,www)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      include 'pwhg_br.h'
      include 'pwhg_math.h'
      real * 8 xrad(3),res(flst_nborn),www
      integer alr,em,iborn,ires
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1     pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac
      real * 8 xjac,csi,csires,y,phi,rrr,rrrct,kemdkres,kem0,ksoft0,
     1     m2res,omcthse,cosphi,sinphi,costh,sinth,softvec(0:3),ctfac
      real * 8 ccf(-5:5)
      real * 8 dotp,momsq03
      external dotp,momsq03
      real * 8 tinyscale,tinycsi,mintheta2
      parameter (tinyscale=0.01d0)
      real * 8 chargeofparticle
      procedure() :: chargeofparticle
      ccf = (/ cf,cf,cf,cf,cf,ca,cf,cf,cf,cf,cf /)

      res = 0

      xjac = kn_sborn/(4*pi)**3 * kn_jacborn

c set up tiny protection for small csi and y->1,
c corresponding to tinyscale
      tinycsi = tinyscale/sqrt(kn_sborn)

      csi = tinycsi + xrad(1)/(1-xrad(1))
      xjac = xjac /(1-xrad(1))**2
      xjac = xjac * csi
c mintheta * csi*sqrt(kn_sborn)=tinyscale
      mintheta2 = (tinyscale/csi)**2/kn_sborn
c mintheta2 is at most 1-> y~0
      y = 1 - mintheta2 - (2-2*mintheta2)*xrad(2)**2
      xjac = xjac * 4 * xrad(2)
      costh = y
      sinth = sqrt(abs(1-y**2))
      phi = 2*pi*xrad(3)
      xjac = xjac * 2*pi
      sinphi = sin(phi)
      cosphi = cos(phi)

      softvec(0) = sqrt(kn_sborn) * csi/2
      softvec(1) = softvec(0)*sinth*cosphi
      softvec(2) = softvec(0)*sinth*sinphi
      softvec(3) = softvec(0)*costh

      do alr=1,flst_nalr
         if(.not.  (flst_alr(flst_alrlength(alr),alr) == 0 .or.
     1              flst_alr(flst_alrlength(alr),alr) == 22 ) ) cycle
         iborn = flst_alr2born(alr)
         if(flst_bornresgroup(iborn).eq.flst_ibornresgroup) then
            em = flst_emitter(alr)
            if(em.gt.2) then
               ires = flst_alrres(em,alr)
            else
               ires = 0
            endif
            if(ires.eq.0) cycle
            kn_softvec = softvec
            call rotate3tovec(kn_cmpborn(1:3,em),kn_softvec(1:3))
c if imode=0, kn_softvec is assumed to be the dimensionful soft momentum,
c and the return value rrr is not multiplied by csi^2.
            call softalr0(alr,em,y,rrr,0)
            m2res = momsq03(kn_cmpborn(:,ires))
            csires = 2*dotp(kn_softvec,kn_cmpborn(:,ires))/m2res
            rrr = rrr*(exp(-csires) - exp(-csi))/(1-y)
            res(iborn) = res(iborn) + rrr
c 1/(1-y)_+ subtraction; only for massless emitter
            if(kn_masses(em) == 0) then
               kemdkres = dotp(kn_cmpborn(:,em),kn_cmpborn(:,ires))
               kem0 = kn_cmpborn(0,em)
               ksoft0 = kn_softvec(0)
               if(flst_alr(flst_alrlength(alr),alr) == 0) then
                  ctfac = 32*pi*st_alpha * ccf(flst_alr(em,alr))
               elseif(flst_alr(flst_alrlength(alr),alr) == 22) then
                  ctfac = 32*pi*em_alpha * chargeofparticle(flst_alr(em,alr))**2
               endif
               rrrct = ctfac/kn_sborn/csi**2
     1              * br_born(iborn) *
     2              ( exp(-2*kemdkres/m2res * ksoft0/kem0)-exp(-csi) )
     3              /(1-y)
               res(iborn) = res(iborn) - rrrct
            endif
c            if(omcthse.lt.0.001) then
c               write(*,*) ' softmismatchres check:', rrrct/rrr,(rrrct-rrr)/rrr,omcthse/(1-y)
c               write(*,*) (exp(-csires) - exp(-csi))/ ( exp(-2*kemdkres/m2res * ksoft0/kem0)-exp(-csi) )
c               write(*,*) csires/(2*kemdkres/m2res * ksoft0/kem0),br_born(iborn)
c            endif
         endif
      enddo
      
      if(.not.flg_minlo) then
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
         rescfac=1
      endif
      do iborn=1,flst_nborn
         if(res(iborn).ne.0) then
            if(flg_minlo) then
               call setlocalscales(iborn,2,rescfac)
               call pdfcall(1,kn_xb1,pdf1)
               call pdfcall(2,kn_xb2,pdf2)
            endif
            res(iborn)=res(iborn)*pdf1(flst_born(1,iborn))*pdf2(flst_born(2,iborn))
     1           * rescfac * xjac
         endif
      enddo         
      
      end