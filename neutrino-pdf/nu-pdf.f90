PROGRAM nupdf
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: dat1,dat2,dat3,datm,datp
   REAL(KIND=dp) :: s
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: E,r
   CHARACTER(LEN=80) :: filename1,filename2,filename3,tmp,outname
   CHARACTER(LEN=80) :: fmt,fname
   INTEGER :: argc,u,status,nline,i,j
   LOGICAL :: exists,replica
   INTEGER(KIND=4) :: nrep
   REAL(KIND=dp) :: xsec,evrate
   INTERFACE 
      SUBROUTINE readflux(filename,dat)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT) :: dat
      END SUBROUTINE readflux
      SUBROUTINE readuncert(filename,datm,datp)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT) :: datm,datp
      END SUBROUTINE readuncert
      SUBROUTINE box_muller(r)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         REAL(KIND=dp),DIMENSION(:),INTENT(INOUT) :: r
      END SUBROUTINE box_muller
   END INTERFACE
   argc=COMMAND_ARGUMENT_COUNT()
   IF(argc.LT.5)THEN
      WRITE(*,*) "Provide a file and a centre of mass energy and lhapdf &
        &file name as input."
      STOP
   END IF
   replica=.FALSE.
   IF(argc.EQ.8)THEN
      replica=.TRUE.
      CALL GETARG(8,tmp)
      READ(tmp,*,IOSTAT=status) nrep
      CALL checkread(status)
   END IF

   CALL GETARG(1,outname)
   CALL GETARG(2,tmp)
   CALL GETARG(3,filename1)
   CALL GETARG(4,filename2)
   CALL GETARG(5,filename3)

   status=0
   READ(tmp,*,IOSTAT=status) s
   IF(status.NE.0)THEN
      WRITE(*,*) "Could not parse value given for energy"
      STOP
   END IF
   CALL GETARG(6,tmp)
   READ(tmp,*,IOSTAT=status) xsec
   CALL checkread(status)
   CALL GETARG(7,tmp)
   READ(tmp,*,IOSTAT=status) evrate
   CALL checkread(status)

   CALL checkfilename(filename1)
   CALL checkfilename(filename2)
   CALL checkfilename(filename3)
   CALL getnlines(filename1,nline)
   ALLOCATE(dat1(4,nline))
   CALL getnlines(filename2,nline)
   ALLOCATE(dat2(4,nline))
   CALL getnlines(filename3,nline)
   ALLOCATE(dat3(4,nline))
   ALLOCATE(datp(4,nline))
   ALLOCATE(datm(4,nline))
   ALLOCATE(r(nline))
   CALL readflux(filename1,dat1)
   CALL readflux(filename2,dat2)
   CALL readflux(filename3,dat3)

   ! should check that the binning is the same.

   WRITE(tmp,"(I4.4)") nrep
   outname=ADJUSTL(TRIM(outname))//"_"// &
      ADJUSTL(TRIM(tmp))//".dat"
   OPEN(NEWUNIT=u,FILE=outname,STATUS="NEW",ACTION="WRITE")
   ALLOCATE(E(LBOUND(dat1,1):UBOUND(dat1,2)))
   E(:)=dat1(2,:)+dat1(1,:)
   E=E/2.0_dp/s
   dat1(3,:)=dat1(3,:)/2.0_dp*E(:)/xsec*evrate
   dat2(3,:)=dat2(3,:)/2.0_dp*E(:)/xsec*evrate
   dat3(3,:)=dat3(3,:)/2.0_dp*E(:)/xsec*evrate
   CALL readuncert(filename1,datm,datp)
   datm(3,:)=datm(3,:)/2.0_dp*E(:)/xsec*evrate
   datp(3,:)=datp(3,:)/2.0_dp*E(:)/xsec*evrate
   IF(replica)THEN
      CALL box_muller(r)
      IF(nrep.EQ.0) r=0.0_dp
      dat1(3,:)=dat1(3,:)+r(:)*(datp(3,:)-datm(3,:))/2.0_dp
   END IF
   CALL readuncert(filename2,datm,datp)
   datm(3,:)=datm(3,:)/2.0_dp*E(:)/xsec*evrate
   datp(3,:)=datp(3,:)/2.0_dp*E(:)/xsec*evrate
   IF(replica)THEN
      CALL box_muller(r)
      IF(nrep.EQ.0) r=0.0_dp
      dat2(3,:)=dat2(3,:)+r(:)*(datp(3,:)-datm(3,:))/2.0_dp
   END IF
   CALL readuncert(filename3,datm,datp)
   datm(3,:)=datm(3,:)/2.0_dp*E(:)/xsec*evrate
   datp(3,:)=datp(3,:)/2.0_dp*E(:)/xsec*evrate
   IF(replica)THEN
      CALL box_muller(r)
      IF(nrep.EQ.0) r=0.0_dp
      dat3(3,:)=dat3(3,:)+r(:)*(datp(3,:)-datm(3,:))/2.0_dp
   END IF
   IF(.NOT.replica)THEN
      WRITE(u,"(A)") "PdfType: central"
   ELSE
      WRITE(u,"(A)") "PdfType: replica"
   END IF
   WRITE(u,"(A)") "Format: lhagrid1"
   WRITE(u,"(A)") "---"
   WRITE(tmp,*) INT(SIZE(E),KIND=4)
   fmt="("//ADJUSTL(TRIM(tmp))
   fmt=ADJUSTL(TRIM(fmt))//"E14.7)"
   WRITE(u,fmt) (E(i),i=LBOUND(E,1),UBOUND(E,1))
   WRITE(u,"(2E14.7)") 0.1_dp,1.0_dp
   WRITE(u,"(6I4)") -16,-14,-12,12,14,16
   DO i=LBOUND(dat1,2),UBOUND(dat1,2)
      WRITE(u,"(6E16.7)") dat3(3,i),dat2(3,i),dat1(3,i),dat1(3,i), &
         dat2(3,i),dat3(3,i)
      WRITE(u,"(6E16.7)") dat3(3,i),dat2(3,i),dat1(3,i),dat1(3,i), &
         dat2(3,i),dat3(3,i)
   END DO
   WRITE(u,"(A)") "---"
   DEALLOCATE(E)
   DEALLOCATE(dat1)
   DEALLOCATE(dat2)
   DEALLOCATE(dat3)
   DEALLOCATE(datp)
   DEALLOCATE(datm)
   CLOSE(u)
END PROGRAM nupdf

SUBROUTINE readflux(filename,dat)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT) :: dat
   INTEGER(KIND=4) :: u,status,i,j
   OPEN(NEWUNIT=u,FILE=TRIM(filename),STATUS="OLD",ACTION="READ")
   DO i=1,UBOUND(dat,2)
      READ(u,*,IOSTAT=status) (dat(j,i),j=1,UBOUND(dat,1))
      IF(status.NE.0)EXIT
   END DO
   CLOSE(u)
END SUBROUTINE readflux

SUBROUTINE readuncert(filename,datm,datp)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT) :: datm,datp
   CHARACTER(LEN=80) :: fname
   INTEGER(KIND=4) :: u,status,i,j,ind
   ind=INDEX(filename,"11")
   fname=filename(1:ind-1)//"min"//filename(ind+2:)
   OPEN(NEWUNIT=u,FILE=TRIM(fname),STATUS="OLD",ACTION="READ")
   DO i=1,UBOUND(datm,2)
      READ(u,*,IOSTAT=status) (datm(j,i),j=1,UBOUND(datm,1))
      IF(status.NE.0)EXIT
   END DO
   CLOSE(u)
   fname=filename(1:ind-1)//"max"//filename(ind+2:)
   OPEN(NEWUNIT=u,FILE=TRIM(fname),STATUS="OLD",ACTION="READ")
   DO i=1,UBOUND(datp,2)
      READ(u,*,IOSTAT=status) (datp(j,i),j=1,UBOUND(datp,1))
      IF(status.NE.0)EXIT
   END DO
   CLOSE(u)
END SUBROUTINE readuncert

SUBROUTINE getnlines(filename,nline)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   INTEGER(KIND=4),INTENT(INOUT) :: nline
   CHARACTER(LEN=80) :: tmp
   INTEGER(KIND=4) :: u,status
   nline=0
   OPEN(NEWUNIT=u,FILE=TRIM(filename),STATUS="OLD",ACTION="READ")
   DO
      READ(u,*,IOSTAT=status) tmp
      IF(status.NE.0)EXIT
      nline=nline+1
   END DO
   CLOSE(u)
END SUBROUTINE getnlines

SUBROUTINE checkfilename(filename)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   LOGICAL :: exists
   INQUIRE(FILE=TRIM(filename),EXIST=exists)
   IF(.NOT.exists)THEN
      WRITE(*,*) "File "//ADJUSTL(TRIM(filename))//" not found."
      STOP
   END IF
END SUBROUTINE checkfilename

SUBROUTINE checkread(status)
   INTEGER(KIND=4),INTENT(IN) :: status
   IF(status.NE.0)THEN
      WRITE(*,*) "Could not parse value given for energy"
      STOP
   END IF
END SUBROUTINE checkread


SUBROUTINE box_muller(r)
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   REAL(KIND=dp),DIMENSION(:),INTENT(INOUT) :: r
   INTEGER(KIND=4) :: i
   REAL(KIND=dp) :: u1,u2
   DO i=1,UBOUND(r,1)
      CALL RANDOM_NUMBER(u1)
      CALL RANDOM_NUMBER(u2)
      r(i)=DSQRT(-2*DLOG(u1))*DCOS(2*PI*u2)
   END DO
END SUBROUTINE box_muller 
