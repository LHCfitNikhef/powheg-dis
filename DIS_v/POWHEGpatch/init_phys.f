      subroutine init_phys
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      include 'pwhg_rad.h'
      include 'pwhg_dbg.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_rnd.h'
      include 'pwhg_physpar.h'
      character * 5 scheme
      character * 3 whichpdfpk
      real * 8 powheginput
      integer iorder,iret,iun,j
      external whichpdfpk,powheginput
      real * 8 Ebeam1, Ebeam2
c Initialization of default values for common block
c variables. These may be overridden by the user program
c init_processes.

      par_diexp=powheginput("#par_diexp")
      par_dijexp=powheginput("#par_dijexp")
      par_2gsupp=powheginput("#par_2gsupp")
      if(par_diexp.lt.0) par_diexp=1
      if(par_dijexp.lt.0) par_dijexp=1
      if(par_2gsupp.lt.0) par_2gsupp=1
      if(par_diexp.ne.par_dijexp) then
         write(*,*) 'par_dijexp not equal to par_diexp;'
         write(*,*) 'not possible at present!!! Program exits'
         call exit(-1)
      endif
c
      par_isrtinycsi = 1d-6
      par_isrtinyy = 1d-6
      par_fsrtinycsi = 1d-5
      par_fsrtinyy = 1d-6
c
      rad_branching=1
c this is set to true in processes where the FSR jacobian
c can become singular (massless recoil particle)
      flg_jacsing=.false.
c flag to use importance sampling in the x variable in
c collinear remnant generation. Needed for charm at LHC
      flg_collremnsamp=.false.
c End initialization of common block defaults.
c     pdf_cutoff_fact^2, multiplied by pdf_q2min, modulates
c     the minimum q2 at which PDF's are called, namely
c     PDF's are frozen (lhapdf6if.f, for example) at scales < pdf_q2min*pdf_cutoff_fact**2.
c     During generation of radiation, st_mufact2 is never below
c     q2min, so setting this to 1 does not alter the traditional behaviour
c     During the Bbar computation differences  may arise if very small
c     scales are reached, as in MiNLO applications.
      pdf_cutoff_fact=1
      pdf_ih1=powheginput('ih1')
      pdf_ih2=powheginput('ih2')
      if(whichpdfpk().eq.'lha') then
         pdf_ndns1=powheginput('lhans1')
         pdf_ndns2=powheginput('lhans2')
      elseif(whichpdfpk().eq.'mlm') then
         pdf_ndns1=powheginput('ndns1')
         pdf_ndns2=powheginput('ndns2')
      else
         write(*,*) ' unimplemented pdf package',whichpdfpk()
         stop
      endif
c     the following are needed only if fullrwgt is on, and are
c     set in due time. They are better initialized as follows
      pdf_ndns1lhe = pdf_ndns1
      pdf_ndns2lhe = pdf_ndns2
      if(pdf_ndns1.ne.pdf_ndns2) then
         st_lambda5MSB=powheginput('QCDLambda5')
      else
         call genericpdfpar(pdf_ndns1,pdf_ih1,st_lambda5MSB,
     1                      scheme,iorder,iret)
         if(iret.ne.0) then
            write(*,*) ' faulty pdf number ',pdf_ndns1
            stop
         endif
      endif
      Ebeam1=powheginput('ebeam1')
      Ebeam2=powheginput('ebeam2')
      if(powheginput("#fixed_target") .ne. 1d0) then
!     Default behaviour: two massless beams colliding beams
         kn_beams(0,1)=Ebeam1
         kn_beams(0,2)=Ebeam2
      else
         write(*,*) "Colliding a beam with energy ", Ebeam1,
     $        " against a fixed target with mass ",  Ebeam2
!     Find two massless four-vectons kn_beams(:,1:2) such that
!     kn_beams(:,1) is parallel to the incoming lepton momentum = (Ebeam1,0,0,Ebeam1)
!     kn_beams(:,2) is anti-parallel to kn_beams(:,1) and
!     kn_beams(:,1)+kn_beams(:,2) = (Ebeam1 + Ebeam2, 0, 0, Ebeam1) (Ebeam2 is the mass of nucleon)
         kn_beams(0,1)=Ebeam1+0.5*Ebeam2
         kn_beams(0,2)=0.5*Ebeam2
      endif
      kn_beams(1,1)=0
      kn_beams(1,2)=0
      kn_beams(2,1)=0
      kn_beams(2,2)=0
      kn_beams(3,1)=kn_beams(0,1)
      kn_beams(3,2)=-kn_beams(0,2)
      kn_sbeams=4*kn_beams(0,1)*kn_beams(0,2)

c generation cut: see Gen_born_phsp.f
      kn_ktmin=powheginput("#bornktmin")
      if(kn_ktmin.lt.0) kn_ktmin=0

c masses for light fermions, used in momentum reshuffling
      do j=1,6
         physpar_mq(j)=0
      enddo
      do j=1,3
         physpar_ml(j)=0
      enddo

c thresholds 
      rad_ptsqmin=powheginput('#ptsqmin')
      if(rad_ptsqmin.lt.0) rad_ptsqmin=0.8d0
      rad_charmthr2=powheginput('#charmthr')
      if(rad_charmthr2.lt.0) rad_charmthr2=1.5d0
      rad_charmthr2=rad_charmthr2**2
      rad_bottomthr2=powheginput('#bottomthr')
      if(rad_bottomthr2.lt.0) rad_bottomthr2=5d0
      rad_bottomthr2=rad_bottomthr2**2
c scale factors
      st_renfact=powheginput('#renscfact')
      st_facfact=powheginput('#facscfact')
      if(st_facfact.lt.0) st_facfact=1
      if(st_renfact.lt.0) st_renfact=1

c     if true, perform the check that there are no coloured light
c     partons before flst_lightpart
      flg_lightpart_check=.true.
c

      flg_evenmaxrat = .false.
      if(powheginput("#evenmaxrat").eq.1) flg_evenmaxrat = .true.

c initialize Lambda values for radiation
      call init_rad_lambda
c

      call init_couplings
      call init_processes
      call setup_reson

      call genflavreglist
c initialize arrays containing all possible resonance
c histories for each alr
      call fill_res_histories

      dbg_softtest=.true.
      dbg_colltest=.true.
      if(powheginput("#softtest").eq.0) dbg_softtest=.false.
      if(powheginput("#colltest").eq.0) dbg_colltest=.false.
      if(flg_bornonly) then
         write(*,*)
     $        ' POWHEG: no soft and coll. tests if bornonly is set'
         dbg_softtest=.false.
         dbg_colltest=.false.
      endif
      if (dbg_softtest.or.dbg_colltest) then         
         call newunit(iun)
         if(rnd_cwhichseed == 'none') then
            open(unit=iun,file='pwhg_checklimits')
         else
            open(unit=iun,file='pwhg_checklimits-'//trim(rnd_cwhichseed))
         endif
         call checklims(iun)
         call flush(iun)
         write(*,*) ' POWHEG:  '
         write(*,*) ' Check of soft/collinear limits performed'
         write(*,*) ' Results in file pwhg_checklimits'
      endif
      end


      subroutine setup_reson
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      integer j,k,res
      flst_nreson=0
      do j=1,flst_nreal
         res=flst_realres(flst_reallength(j),j)
         do k=1,flst_nreson
            if(flst_reslist(k).eq.res) exit
         enddo
         if(k.eq.flst_nreson+1) then
c it didn't find the resonance on the list; add it up
            flst_nreson=flst_nreson+1
            flst_reslist(flst_nreson)=res
         endif
      enddo
      if(flst_nreson.eq.1.and.flst_reslist(flst_nreson).eq.0) then
         flg_withresrad=.false.
      else
         flg_withresrad=.true.
      endif
      end
