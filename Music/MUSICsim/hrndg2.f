*********************************************************************
* Return a random point (x1,x2) distributed according to the contents 
*  of the 2D Matrix ff(,) generated by the function "gaisser"
*********************************************************************
      Subroutine hrndg2(ff,ii,a1,b1,jj,a2,b2,sum,x1,x2,pro)
      implicit none
      real*8 gaisser
      EXTERNAL gaisser
      integer i,j,k,ii,jj,flag_p(2),flag_bin(2), m,n
      real*8 an1,bn1,an2,bn2
      real*8 abin1,bbin1,abin2,bbin2,bsum    ! for integer a bin
      real*8 a1,b1,a2,b2, x,y
      real*8 sum,sumbox1,binx,biny,binz,x1,x2,pp
      real*8 sumbox2,binxn,binyn
      real*8 vec(3)
      real pro
      real*8  sta1                          ! sta1: flag to determine->if need to redo integral
      Real*8 ff(ii,jj)                      ! when "lowmu"(a1) changes, the flux changes too.
      Real*8 bb(100,20)
      
      logical first
      data first/.TRUE./
      data m,n/100,3/
      data sta1/0.0/

**********************************
      IF(FIRST.or.sta1.ne.a1) THEN
         first = .FALSE.
      binx=(b1-a1)/ii
      biny=(b2-a2)/jj

*--> book the 2d matrix & do integral
**********************************
       sumbox1 = 0.d0
***      print *,'*--> do integral !  ---   Start!!!   --- '
       do i=1,ii
***        if(mod(i,500).eq.0) print*,'^^^',i,'   ','sumbox1: ',sumbox1
        do j=1,jj
*----------------------
*--> calculate the bin integral
* Which is only need to calculate for one time
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          bsum=0.d0
          abin1=a1 + binx*(i-1)
          bbin1=abin1 + binx
          abin2=a2 + biny*(j-1)
          bbin2=abin2 + biny
           if(abin1.gt.80)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,5,4)
           else if(abin1.gt.40)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,8,4)
           else if(abin1.gt.10)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,40,4)
           else if(abin1.gt.5)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,200,4)
           else if(abin1.gt.1)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,800,3)
           else if(abin1.gt.0.6)then
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,6000,3)
           else
            call d2_integral(abin1,bbin1,abin2,bbin2,bsum,40000,3)
           endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*----------------------
         sumbox1=bsum+sumbox1
         ff(i,j)=sumbox1
        enddo
       enddo
      sum = sumbox1
      print *,'*--> do integral !  ---  Result!!!  --- ',sum,'hz'

       If (pro.gt.10.)then
        print *,'double precision !'
       Else
        print *,'single precision !'
       Endif

      sta1=a1
      ENDIF                        !   only once called  for the same "lowmu"

**********************************


*--> generate three variables according to the muon_flux(,) 
**********************************
       call RM48(vec,1)
       pp=vec(1)*sumbox1
      do i=1,ii
       do j=1,jj
         if(ff(i,j).ge.pp) then
          flag_p(1)=i
          flag_p(2)=j
          goto 17
         endif
       enddo
      enddo 

*-> do above again in the selected bin
*    to get better precision
***************
17    IF (pro.gt.10.)then
       an1=a1 + binx*(flag_p(1)-1)
       bn1=an1 + binx
       an2=a2 + biny*(flag_p(2)-1)
       bn2=an2 + biny
       binxn=binx/m
       binyn=biny/n

        sumbox2=0.d0
        do i=1,m
         do j=1,n
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          bsum=0.d0
          abin1=an1 + binxn*(i-1)
          bbin1=abin1 + binxn
          abin2=an2 + binyn*(j-1)
          bbin2=abin2 + binyn
           call d2_integral(abin1,bbin1,abin2,bbin2,bsum,2,2)
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          sumbox2=bsum+sumbox2
          bb(i,j)=sumbox2
         enddo
        enddo

        call RM48(vec,1)
        pp=vec(1)*sumbox2
       do i=1,m
        do j=1,n
          if(bb(i,j).ge.pp) then
           flag_bin(1)=i
           flag_bin(2)=j
           goto 18
          endif
        enddo
       enddo 

* Give two random numbers to get the (x1,x2).
**********************************
18      call RM48(vec,2)
        x1=an1+binxn*(flag_bin(1)-1) + binxn*vec(1)
        x2=an2+binyn*(flag_bin(2)-1) + binyn*vec(2)
 
      ELSE 
       call RM48(vec,2)
        x1=a1+binx*(flag_p(1)-1) + binx*vec(1)
        x2=a2+biny*(flag_p(2)-1) + biny*vec(2)
      ENDIF
***************
 
*      print'("vec(4) : ",4f12.5)',(vec(i),i=1,4)
*       if (flag_p(1).gt.0)then
*        print'("flag_p : ",3i12)',(flag_p(i),i=1,3)
*        print'("coordi : ",3f12.3)',x1,x2,x3
*       endif
**********************************
      END


