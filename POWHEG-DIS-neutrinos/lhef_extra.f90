SUBROUTINE lhefwriteevrw_new(nlf)
   IMPLICIT NONE
   INCLUDE "nlegborn.h"
   INCLUDE "pwhg_flst.h"
   INCLUDE "pwhg_rad.h"
   INCLUDE "pwhg_rwl.h"
   CHARACTER(LEN=132) :: string
   INTEGER :: id1,id2,iret
   REAL(KIND=8) :: x1,x2,xf1,xf2,xmufact
   CHARACTER(LEN=100) :: buffer

   CALL pdfreweightinfo(id1,id2,x1,x2,xmufact,xf1,xf2)
   WRITE(buffer,111)'#pdf ',id1,id2,x1,x2,xmufact,xf1,xf2

   CALL pwhg_io_write(nlf,trim(buffer))
   WRITE(string,*)"rwgt",rwl_type,rwl_index,rwl_weight,rwl_seed,&
                   & rwl_n1,rwl_n2,rad_kinreg,rad_ubfailrad
   CALL pwhg_io_write(nlf,TRIM(ADJUSTL(string)))
   111  FORMAT(a,2(1x,i2),5(1x,e14.8))
END SUBROUTINE lhefwriteevrw

SUBROUTINE lhefreadev_new(nlf)
   IMPLICIT NONE
   INTEGER(KIND=8),INTENT(IN) :: nlf
   CHARACTER(LEN=200) :: string
   INCLUDE "LesHouches.h"
   INCLUDE "pwhg_kn.h"
   INTEGER(KIND=8) :: i,j,iret
   INTEGER(KIND=8) :: flav1,flav2,lhanum
   REAL(KIND=8) ::  x1,x2,scale1,scale2
   string=" "
   iret=-9999
   DO WHILE(iret.NE.0)
      CALL pwhg_io_read(nlf,string,iret)
      IF(iret.NE.0.OR.string.EQ."</LesHouchesEvents>") THEN
         nup=0
         RETURN
      END IF 
      IF(string(1:6).EQ."<event") THEN
         CALL pwhg_io_read(nlf,string,iret)
         IF(iret.NE.0) THEN
            nup=0
            RETURN
         END IF
         READ(string,*,ERR=1) nup,idprup,xwgtup,scalup,aqedup,aqcdup
         DO i=1,nup
            CALL pwhg_io_read(nlf,string,iret)
            IF(iret.NE.0) THEN
               nup=0
               RETURN
            END IF
            READ(string,*,ERR=1) idup(i),istup(i),mothup(1,i), &
                               & mothup(2,i),icolup(1,i),icolup(2,i),&
                               & (pup(j,i),j=1,5),vtimup(i),spinup(i)
         END DO
         CALL lhefreadextra(nlf,iret)
         ! Read the line with the pdf information.
         CALL pwhg_io_read(nlf,string,iret)
         IF(string(1:4.EQ."#pdf") THEN
            READ(string,*) flav1,flav2,x1,x2,scale1,scale2,lhanum
         END IF
         END IF
      END IF
   END DO


   END IF
END SUBROUTINE lhefreadev_new
