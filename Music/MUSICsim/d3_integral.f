*********************************************************************
* do integral of the Funcation || 
*  method Divide to sub mini bins & add content ||
*********************************************************************
      Subroutine d3_integral(a1,a2,b1,b2,c1,c2,sum,m,n,l)
      implicit none
      real*8 gaisser
      EXTERNAL gaisser
      integer i,j,k, m,n,l
      real*8 a1,a2,b1,b2,c1,c2
      real*8 sum,sumbox,binx,biny,binz
      real*8 x,y,z
      Real*8 fff(m,n,l)
      
**********************************
      binx=(a2-a1)/m
      biny=(b2-b1)/n
      binz=(c2-c1)/l

*--> book 3d matrix, Get the sum
**********************************
      sumbox=0.d0
      do i=1,m
       do j=1,n
        do k=1,l
         x=a1+binx*(i-0.5)
         y=b1+biny*(j-0.5)
         z=c1+binz*(k-0.5)
         fff(i,j,k)=gaisser(x,y,z)
         sumbox=fff(i,j,k)+sumbox
         enddo
       enddo
      enddo

**********************************

*--> do integral
**********************************
      SUM = sumbox*binx*biny*binz
*
*      print*,'binx,y,z; sum: ',binx,biny,binz,sum
**********************************

      end


