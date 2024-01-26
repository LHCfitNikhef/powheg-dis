PROGRAM test
   USE ddxsec
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   CALL init_dists
   CALL bookupeqdd("mydd1",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd2",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd3",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd4",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd1",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd2",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd3",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   CALL bookupeqdd("mydd4",1.0_dp,10.0_dp,2.0_dp,20.0_dp,40.0_dp,3.0_dp,2)
   !CALL dists(1)%fill(1.0_dp,24.0_dp,(/10.0_dp,2.0_dp,44.0_dp,22.0_dp,22.0_dp/))
   !CALL dists(1)%fill(2.0_dp,22.0_dp,(/10.0_dp,2.0_dp,44.0_dp,22.0_dp,22.0_dp/))
   !CALL dists(1)%fill(3.0_dp,35.0_dp,(/3.0_dp,3.0_dp,44.0_dp,22.0_dp,22.0_dp/))
   !CALL dists(1)%fill(4.0_dp,30.0_dp,(/3.0_dp,2.0_dp,44.0_dp,22.0_dp,22.0_dp/))
   !CALL dists(1)%fill(1.0_dp,20.0_dp,(/2.0_dp,2.0_dp,44.0_dp,22.0_dp,22.0_dp/))
   CALL filldd("mydd1",2.0_dp,10.0_dp,(/1.0_dp,2.9_dp/))
   CALL filldd("mydd2",2.0_dp,10.0_dp,(/1.0_dp,2.9_dp/))
   CALL filldd("mydd3",2.0_dp,10.0_dp,(/1.0_dp,2.9_dp/))
   CALL filldd("mydd1",2.0_dp,10.0_dp,(/1.0_dp,2.9_dp/))
   CALL filldd("mydd4",2.0_dp,10.0_dp,(/1.0_dp,2.9_dp/))
   CALL dists(2)%write_to("mynewfile")
      
END PROGRAM test
