PROGRAM nupdf
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: dat1,dat2,dat3
   REAL(KIND=dp) :: s
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: E
   CHARACTER(LEN=80) :: filename1,filename2,filename3,tmp,outname
   CHARACTER(LEN=80) :: fmt
   INTEGER :: argc,u,status,nline,i,j
   LOGICAL :: exists
   INTERFACE 
      SUBROUTINE readflux(filename,dat)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT) :: dat
      END SUBROUTINE readflux
   END INTERFACE
   argc=COMMAND_ARGUMENT_COUNT()
   IF(argc.LT.3)THEN
      WRITE(*,*) "Provide a file and a centre of mass energy and lhapdf &
        &file name as input."
      STOP
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
   CALL checkfilename(filename1)
   CALL checkfilename(filename2)
   CALL checkfilename(filename3)
   CALL getnlines(filename1,nline)
   ALLOCATE(dat1(4,nline))
   CALL getnlines(filename2,nline)
   ALLOCATE(dat2(4,nline))
   CALL getnlines(filename3,nline)
   ALLOCATE(dat3(4,nline))
   CALL readflux(filename1,dat1)
   CALL readflux(filename2,dat2)
   CALL readflux(filename3,dat3)
   ! should check that the binning is the same.
   OPEN(NEWUNIT=u,FILE=outname,STATUS="NEW",ACTION="WRITE")
   WRITE(u,"(A)") "PdfType: central"
   WRITE(u,"(A)") "Format: lhagrid1"
   WRITE(u,"(A)") "---"
   ALLOCATE(E(LBOUND(dat1,1):UBOUND(dat1,2)))
   E(:)=dat1(2,:)+dat1(1,:)
   E=E/2.0_dp/s
   dat1(3,:)=dat1(3,:)/2.0_dp*E(:)
   dat2(3,:)=dat2(3,:)/2.0_dp*E(:)
   dat3(3,:)=dat3(3,:)/2.0_dp*E(:)
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
