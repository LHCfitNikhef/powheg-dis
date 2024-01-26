MODULE ddxsec

   IMPLICIT NONE

   PRIVATE :: ddxs,is_present,write_to,dp
   PUBLIC  :: dists,bookupeqdd,init_dists,fill
   
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   INTERFACE ddxs
      MODULE PROCEDURE :: new_ddxs
   END INTERFACE ddxs

   TYPE ddxs
      CHARACTER(LEN=80) :: str
      REAL(KIND=dp)     :: lx,hx,ly,hy,sx,sy,lz,hz,sz
      INTEGER(KIND=8)   :: nbinx,nbiny,nbinz
      INTEGER(KIND=8)   :: nds,nevents
      REAL(KIND=dp),DIMENSION(:,:,:,:),ALLOCATABLE :: arr
      CONTAINS
         PROCEDURE :: fill
         PROCEDURE :: write_to
         PROCEDURE :: is_present
         PROCEDURE :: get_element
   END TYPE
   
   TYPE(ddxs),DIMENSION(:),POINTER,SAVE :: dists

   CONTAINS
      

      SUBROUTINE init_dists
         INTEGER(KIND=4) :: i
         ALLOCATE(dists(1))
         DO i=LBOUND(dists,DIM=1),UBOUND(dists,DIM=1)
            dists(i)%str=""
         END DO
      END SUBROUTINE init_dists

      SUBROUTINE filldd(str,sig,x,y,z)
         CHARACTER(LEN=80),INTENT(IN) :: str
         REAL(KIND=dp),INTENT(IN) :: x,y
         REAL(KIND=dp),INTENT(IN),OPTIONAL :: z
         REAL(KIND=dp),DIMENSION(*),INTENT(IN) :: sig
         INTEGER(KIND=4) :: i
         i=get_dist_index(str)
         WRITE(*,*) "Recieved dsig: ",sig(1)
         IF(PRESENT(z)) THEN
            CALL dists(i)%fill(sig,x,y,z)
         ELSE
            CALL dists(i)%fill(sig,x,y)
         END IF
      END SUBROUTINE filldd

      INTEGER FUNCTION get_dist_index(str) RESULT(ind)
         CHARACTER(LEN=80),INTENT(IN) :: str
         INTEGER(KIND=4) :: i
         DO i=LBOUND(dists,DIM=1),UBOUND(dists,DIM=1)
            IF (str.EQ.dists(i)%str) ind=i
         END DO
      END FUNCTION get_dist_index


      PURE FUNCTION get_element(self,i,j,k,l) RESULT(elem)
         CLASS(ddxs),INTENT(IN) :: self
         INTEGER(KIND=4),INTENT(IN) :: i,j,k,l
         REAL(KIND=8) :: elem
         elem=self%arr(i,j,k,l)
      END FUNCTION get_element

! Constructor for the ddxs derived type.
      TYPE(ddxs) FUNCTION new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt,lz,hz,sz)
         CHARACTER(LEN=80),INTENT(IN) :: str
         REAL(KIND=dp),INTENT(IN) :: lx,hx,ly,hy,sx,sy
         REAL(KIND=dp),INTENT(IN),OPTIONAL :: lz,hz,sz
         INTEGER(KIND=4),INTENT(IN) :: nwgt
         REAL(KIND=dp) :: nlx,nly,nhx,nhy,nhz,nlz
         INTEGER(KIND=8) :: nbinx,nbiny,nbinz,nds

         IF(PRESENT(lz).AND.PRESENT(hz).AND.PRESENT(sz)) THEN
            nds=3
         ELSE
            nds=2
         END IF

         CALL check_bounds(lx,hx,nlx,nhx)
         nbinx=get_nbin(lx,hx,sx)
         nhx=lx+nbinx*sx
         CALL check_bounds(ly,hy,nly,nhy)
         nbiny=get_nbin(ly,hy,sy)
         nhy=ly+nbiny*sy
         IF (nds.EQ.3) THEN
            CALL check_bounds(lz,hz,nlz,nhz)
            nbinz=get_nbin(lz,hz,sz)
            nhz=lz+nbinz*sz
            new_ddxs%lz=nlz
            new_ddxs%hz=nhz
            new_ddxs%sz=sz
            new_ddxs%nbinz=nbinz
         ELSE
            new_ddxs%lz=0.0_dp
            new_ddxs%hz=0.0_dp
            new_ddxs%sz=1.0_dp
            nbinz=1
            new_ddxs%nbinz=nbinz
         END IF

         new_ddxs%str=str 
         new_ddxs%lx=nlx
         new_ddxs%hx=nhx
         new_ddxs%sx=sx
         new_ddxs%ly=nly
         new_ddxs%hy=nhy
         new_ddxs%sy=sy
         new_ddxs%nbinx=nbinx
         new_ddxs%nbiny=nbiny
         new_ddxs%nds=nds
         new_ddxs%nevents=0
         ALLOCATE(new_ddxs%arr(nbinx,nbiny,nbinz,nwgt))
         new_ddxs%arr=0.0_dp
         CONTAINS
            SUBROUTINE check_bounds(lb,hb,nlb,nhb)
               REAL(KIND=dp),INTENT(IN) :: lb,hb
               REAL(KIND=dp),INTENT(OUT) :: nlb,nhb
               IF (lb.GT.hb) THEN
                  nlb=hb
                  nhb=lb
               ELSE
                  nlb=lb
                  nhb=hb
               END IF
            END SUBROUTINE check_bounds

            INTEGER FUNCTION get_nbin(lb,hb,step) RESULT(nbin)
               REAL(KIND=dp),INTENT(IN) :: lb,hb,step
               REAL(KIND=dp) :: tmp
               tmp=(hb-lb)/step
               nbin=SIGN(1.0_dp,tmp)*INT(ABS(tmp)+0.5_dp)
            END FUNCTION get_nbin
      END FUNCTION new_ddxs

      PURE INTEGER FUNCTION is_present(self,str) 
         CLASS(ddxs),INTENT(IN) :: self
         CHARACTER(LEN=80),INTENT(IN) :: str
         IF(TRIM(self%str).EQ."") THEN
            is_present=0
         ELSE IF (self%str.EQ.str) THEN
            is_present=1
         ELSE
            is_present=2
         END IF
      END FUNCTION is_present
      
      SUBROUTINE fill(self,xsec,x,y,z)
         CLASS(ddxs),INTENT(INOUT) :: self
         REAL(KIND=dp),INTENT(IN)  :: x,y
         REAL(KIND=dp),INTENT(IN),OPTIONAL  :: z
         REAL(KIND=dp),DIMENSION(LBOUND(self%arr,DIM=4):UBOUND(self%arr,DIM=4)-1),INTENT(IN) :: xsec
         INTEGER(KIND=4) :: i,j,k
         REAL(KIND=dp) :: dxdydz
         k=1
         self%nevents=self%nevents+1
         IF(PRESENT(z)) k=findbin(z,self%lz,self%sz,self%hz)
         i=findbin(x,self%lx,self%sx,self%hx)
         j=findbin(y,self%ly,self%sy,self%hy)
         IF (i.EQ.0.OR.j.EQ.0) return
         IF (PRESENT(z).AND.k.EQ.0) return
         dxdydz=self%sx*self%sy*self%sz
         dxdydz=self%sy
         self%arr(i,j,k,LBOUND(self%arr,DIM=4):UBOUND(self%arr,DIM=4)-1)=self%arr(i,j,k,:)+xsec(:)/dxdydz
         self%arr(i,j,k,UBOUND(self%arr,DIM=4))=self%arr(i,j,k,UBOUND(self%arr,DIM=4))+1.0_dp
      END SUBROUTINE fill

      INTEGER FUNCTION findbin(val,low,step,high) 
         REAL(KIND=dp),INTENT(IN) :: val,low,step,high
         INTEGER(KIND=8) :: i,nbins
         nbins=INT((high-low)/step,KIND=8)
         findbin=0
         DO i=1,nbins
            IF(val.GE.low.AND.val.LT.low+REAL(i,KIND=dp)*step) THEN
               findbin=i
               EXIT
            END IF
         END DO
      END FUNCTION findbin
         

      SUBROUTINE write_to(self,filename)
         CLASS(ddxs),INTENT(IN)       :: self
         CHARACTER(LEN=80),INTENT(IN) :: filename
         INTEGER(KIND=4) :: un
         INTEGER(KIND=4) :: i,j,k,l
         REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: x,y,z
         REAL(KIND=dp) :: nhits

         IF (ALLOCATED(y)) DEALLOCATE(y)
         ALLOCATE(y(self%nbiny+1))
         IF (ALLOCATED(x)) DEALLOCATE(x)
         ALLOCATE(x(self%nbinx+1))
         IF (ALLOCATED(z)) DEALLOCATE(z)
         ALLOCATE(z(self%nbinz+1))

         DO i=1,self%nbinx
            x(i)=self%lx+(i-1)*self%sx
         END DO
         x(self%nbinx+1)=self%hx
         DO i=1,self%nbiny
            y(i)=self%ly+(i-1)*self%sy
         END DO
         y(self%nbiny+1)=self%hy
         DO i=1,self%nbinz
            z(i)=self%lz+(i-1)*self%sz
         END DO
         z(self%nbinz+1)=self%hz

         nhits=0.0_dp
         DO i=LBOUND(self%arr,DIM=1),UBOUND(self%arr,DIM=1)
            DO j=LBOUND(self%arr,DIM=2),UBOUND(self%arr,DIM=2)
               nhits=nhits+SUM(self%arr(i,j,:,UBOUND(self%arr,DIM=4)))
            END DO
         END DO
         IF(nhits.EQ.0) nhits=0.0_dp
         CALL nunit(un)
         OPEN(UNIT=un,FILE=TRIM(filename),STATUS="UNKNOWN")
         WRITE(un,"(A)") "Distribution "//TRIM(self%str)
         IF (self%nds.EQ.2) THEN
            DO i=LBOUND(self%arr,DIM=1),UBOUND(self%arr,DIM=1)
               DO j=LBOUND(self%arr,DIM=2),UBOUND(self%arr,DIM=2)
                  WRITE(un,*) x(i),x(i+1),y(j),y(j+1),&
                &(self%arr(i,j,1,k)/self%nevents,k=LBOUND(self%arr,DIM=4),UBOUND(self%arr,DIM=4))
               END DO 
            END DO
         ELSE
            DO i=LBOUND(self%arr,DIM=1),UBOUND(self%arr,DIM=1)
               DO j=LBOUND(self%arr,DIM=2),UBOUND(self%arr,DIM=2)
                  DO k=LBOUND(self%arr,DIM=3),UBOUND(self%arr,DIM=3)
                  WRITE(un,*) x(i),x(i+1),y(j),y(j+1),z(k),z(k+1),&
                            & (self%arr(i,j,k,l), &
                            & l=LBOUND(self%arr,DIM=4),UBOUND(self%arr,DIM=4))
                  END DO
                END DO 
            END DO
         END IF
         CLOSE(un)

         IF (ALLOCATED(x)) DEALLOCATE(x)
         IF (ALLOCATED(y)) DEALLOCATE(y)
         IF (ALLOCATED(z)) DEALLOCATE(z)
         ! write the content of the array to a file
      END SUBROUTINE write_to

      SUBROUTINE nunit(u)
         INTEGER(KIND=4),INTENT(INOUT) :: u
         INTEGER(KIND=4) :: i
         LOGICAL :: isopen
         DO i=1,1000
            INQUIRE(UNIT=i,OPENED=isopen)
            IF(.NOT.isopen) THEN
               u=i
               EXIT
            END IF
         END DO
      END SUBROUTINE nunit


      SUBROUTINE bookupeqdd(str,lx,hx,sx,ly,hy,sy,lz,hz,sz,nweights)
         CHARACTER(LEN=80),INTENT(IN) :: str
         REAL(KIND=dp),INTENT(IN)     :: lx,hx,sx,ly,hy,sy
         REAL(KIND=dp),INTENT(IN),OPTIONAL :: lz,hz,sz
         INTEGER(KIND=4),INTENT(IN),OPTIONAL   :: nweights
         INTEGER(KIND=4) :: nwgt
         INTEGER(KIND=8)              :: i
         LOGICAL :: booked
         TYPE(ddxs),DIMENSION(:),POINTER :: tmp
         booked=.FALSE.
         
         IF (PRESENT(nweights)) THEN 
            nwgt=nweights+1
         ELSE
            nwgt=2
         END IF
         DO i=1,UBOUND(dists,DIM=1)
            IF (dists(i)%is_present(str).EQ.0) THEN 
               IF(PRESENT(lz).AND.PRESENT(hz).AND.PRESENT(sz)) THEN
                  dists(i)=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt,lz,hz,sz)
               ELSE
                  dists(i)=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt)
               END IF
               booked=.TRUE.
               EXIT
            ELSE IF (dists(i)%is_present(str).EQ.1) THEN
               WRITE(*,*) "Warning: Overwriting existing double differential&
                           & distribution"
               IF(PRESENT(lz).AND.PRESENT(hz).AND.PRESENT(sz)) THEN
                  dists(i)=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt,lz,hz,sz)
               ELSE
                  dists(i)=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt)
               END IF
               booked=.TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.booked) THEN
            WRITE(*,*) "Allocating new distribution"
            ALLOCATE(tmp(UBOUND(dists,DIM=1)))
            tmp=dists
            DEALLOCATE(dists)
            ALLOCATE(dists(UBOUND(tmp,DIM=1)+1))
            dists(1:UBOUND(tmp,DIM=1))=tmp(:)
            DEALLOCATE(tmp)
            IF(PRESENT(lz).AND.PRESENT(hz).AND.PRESENT(sz)) THEN
               dists(UBOUND(dists,DIM=1))=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt,lz,hz,sz)
            ELSE
               dists(UBOUND(dists,DIM=1))=new_ddxs(str,lx,hx,sx,ly,hy,sy,nwgt)
            END IF
         END IF
      END SUBROUTINE bookupeqdd
   
END MODULE ddxsec
