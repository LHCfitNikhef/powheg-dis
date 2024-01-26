PROGRAM mergetable
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   INTEGER(KIND=4) :: un
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: lx,ly,lz,hx,hy,hz
   REAL(KIND=dp),DIMENSION(:,:,:),ALLOCATABLE :: tables
   REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: res
   CHARACTER(LEN=80) :: filename
   INTEGER(KIND=4) :: i,j,nfiles,nbin,ncol,tmpbin,tmpcol
   LOGICAL :: exists
   CHARACTER(LEN=30) FMT

   INTERFACE
      SUBROUTINE filltable(filename,arr,nrow,ncol,header)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         INTEGER(KIND=4),INTENT(IN) :: nrow,ncol
         REAL(KIND=dp),DIMENSION(nrow,ncol),INTENT(INOUT) :: arr
         INTEGER(KIND=4),INTENT(IN),OPTIONAL :: header
      END SUBROUTINE filltable
      SUBROUTINE get_nbin(filename,nbin,header)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         INTEGER(KIND=4),INTENT(INOUT) :: nbin
         INTEGER(KIND=4),OPTIONAL :: header
      END SUBROUTINE get_nbin
      SUBROUTINE get_ncol(filename,ncol,header)
         INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
         CHARACTER(LEN=80),INTENT(IN) :: filename
         INTEGER(KIND=4),INTENT(INOUT) :: ncol
         INTEGER(KIND=4),INTENT(IN),OPTIONAL :: header
      END SUBROUTINE get_ncol
   END INTERFACE

   nfiles=iargc()
   DO i=1,nfiles
      CALL getarg(i,filename)
      INQUIRE(FILE=TRIM(filename),EXIST=exists)
      IF(.NOT.exists) THEN
         WRITE(*,*) "File ",TRIM(filename)," not found."
         STOP
      END IF
      CALL get_nbin(filename,tmpbin,1)
      IF(i.GT.1) THEN
         IF (nbin.NE.tmpbin) THEN
            WRITE(*,*) "Files have different length."
            STOP
         END IF
      END IF
      nbin=tmpbin
      CALL get_ncol(filename,tmpcol,1)
      IF(i.GT.1) THEN
         IF (ncol.NE.tmpcol) THEN
            WRITE(*,*) "Files have different length."
            STOP
         END IF
      END IF
      ncol=tmpcol
   END DO

   ALLOCATE(tables(nfiles,nbin,ncol))
   ALLOCATE(res(nbin,ncol))
   res=0
   DO i=1,nfiles
      CALL getarg(i,filename)
      CALL filltable(filename,tables(i,:,:),nbin,ncol,1) 
      res(:,:)=res(:,:)+tables(i,:,:)
   END DO
   DEALLOCATE(tables)
   res=res/REAL(nfiles,KIND=dp)
   WRITE(FMT,100) ncol
 100 FORMAT(I3)
   DO i=1,nbin
      WRITE(*,"("//ADJUSTL(FMT)//"E16.9"//")") (res(i,j),j=1,ncol)
   END DO
   DEALLOCATE(res)
END PROGRAM mergetable

SUBROUTINE filltable(filename,arr,nrow,ncol,header)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   INTEGER(KIND=4),INTENT(IN) :: nrow,ncol
   REAL(KIND=dp),DIMENSION(nrow,ncol),INTENT(INOUT) :: arr
   INTEGER(KIND=4),INTENT(IN),OPTIONAL :: header
   INTEGER(KIND=4) :: i,j,u,io
   
   OPEN(NEWUNIT=u,FILE=TRIM(filename),STATUS="old",IOSTAT=io)
   IF (io.NE.0) THEN
      WRITE(*,*) "File not found. File:",TRIM(filename)
      STOP
   END IF

   IF(PRESENT(header)) THEN
      DO i=1,header
         READ(u,*) 
      END DO
   END IF

   DO i=1,nrow
      READ(u,*) (arr(i,j),j=1,ncol)
   END DO
END SUBROUTINE filltable


SUBROUTINE get_nbin(filename,nbin,header)
   IMPLICIT NONE
   CHARACTER(LEN=80),INTENT(IN) :: filename
   INTEGER(KIND=4),INTENT(INOUT) :: nbin
   INTEGER(KIND=4),OPTIONAL :: header
   INTEGER(KIND=4) :: i,io,u
   CHARACTER(LEN=200) :: string
   
   OPEN(NEWUNIT=u,FILE=TRIM(filename),STATUS="old",IOSTAT=io)
   IF (io.NE.0) THEN
      WRITE(*,*) "File not found. File:",TRIM(filename)
      STOP
   END IF

   IF(PRESENT(header)) THEN
      DO i=1,header
         READ(u,"(A)") string
      END DO
   END IF

   nbin=0
   DO
      READ(u,FMT="(A)",IOSTAT=io) string
      IF (io.LT.0) EXIT
      nbin=nbin+1
   END DO
   REWIND(u)
   CLOSE(u)
END SUBROUTINE get_nbin

SUBROUTINE get_ncol(filename,ncol,header)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   CHARACTER(LEN=80),INTENT(IN) :: filename
   INTEGER(KIND=4),INTENT(INOUT) :: ncol
   INTEGER(KIND=4),INTENT(IN),OPTIONAL :: header
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: arr
   INTEGER(KIND=4) :: i,j,io,u
   CHARACTER(LEN=400) :: string
   LOGICAL :: readingnum


   OPEN(NEWUNIT=u,FILE=TRIM(filename),STATUS="old",IOSTAT=io)
   IF (io.NE.0) THEN
      WRITE(*,*) "File not found. File:",TRIM(filename)
      STOP
   END IF

   ncol=0
   IF(PRESENT(header)) THEN
      DO i=1,header
         READ(u,"(A)") string
      END DO
   END IF
   READ(u,"(A)") string
   readingnum=.FALSE.
   DO i=1,LEN_TRIM(string)
      IF(string(i:i).NE." ") THEN
         IF (.NOT.readingnum) ncol=ncol+1
         readingnum=.TRUE.
      ELSE
         readingnum=.FALSE.
      END IF
   END DO
   REWIND(u)
   CLOSE(u)
END SUBROUTINE get_ncol


