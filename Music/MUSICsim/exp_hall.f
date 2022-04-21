************************************
!Consider experiment hall, Because it is empty, no rock, So it may influence the
!muon flux, average energy. 
************************************
      Subroutine exp_hall(hall_len,costhet,phi,wh_hall)
      implicit none
      Real*8 costhet, aaa
      Real  phi, thetar, phir, hall_len
      Real  xpro, ypro, zpro, yz, xy,xz,yx,zx,zy
      Real x0,y0,rot,pi,degpi,radpi
      Real x1, y1, z1, x2, y2, z2                     !the root of muon line equation
      Real b,delta,a                                  !the answer of the equation
      Logical flag
      integer wh_hall
      parameter (pi=3.1415926536,degpi=180./pi,radpi=pi/180.)
      thetar=dacos(costhet)
      phir=phi
      IF (wh_hall.eq.3) then      !Far
        rot = 29.5*radpi
        phir = -rot + phir
      else if (wh_hall.eq.2) then      !LA
        rot = 10.4*radpi
        phir = -rot + phir
      else if (wh_hall.eq.1) then     !DYB
        rot = 33.3*radpi
        phir = -rot+phir
      else 
        print*,'ERROR! You did not choise one site!'
        print*,'3=Far, 2=LA, 1=Near'
      ENDIF
!      print*,'costhet, phi ',costhet, phi
!      print*,'thetar, phir ',thetar*degpi, phir*degpi
!      phir=2*pi
!      aaa=12.235/sqrt(8*8+12.235*12.235)
!      thetar=dacos(aaa)
!      print*,'thetar=',thetar*degpi
*--->Unit vector projection
      xpro=sin(thetar)*cos(phir)
      ypro=sin(thetar)*sin(phir)
      zpro=cos(thetar)
      yz=ypro/zpro
      xy=xpro/ypro
      xz=xpro/zpro
      yx=ypro/xpro
      zx=zpro/xpro
      zy=zpro/ypro
!      print*,'xpro, ypro, zpro=',xpro,ypro,zpro
!      print*,'yz, xy=',yz,xy
*--->muon line equation
*       x/xpro=y/ypro=z/zpro
        
*--->Far site
*     x/xpro=y/ypro=z/zpro
*     (y+0.5)^2+(z)^2=16.55^2    
*     two equation together get the cross point.
*     condition 12.793<z<16.551, -11<y<10, -22.3<x<11 
      IF (wh_hall.eq.3) then
        flag=.false.
!       muon line cross with half column.
        delta=(yz)**2-4*(yz**2+1.)*(0.25-16.55**2)
        if (delta.lt.0) then
          print*,'ERROR! delta <0'
        else
          b=yz
          a=yz**2+1
          z1=(-b+sqrt(delta))/(2.*a)
          y1=yz*z1
          x1=xy*y1
!          print*,'x1, y1, z1',x1,y1,z1
*         limit the range
          if((12.793.le.z1)
     &      .and.(-22.3.le.x1.and.x1.le.11.)) then
            flag=.true.
            if (z1.gt.16.551) then
              print*,'x1, y1, z1',x1,y1,z1
              print*,'!!!!!!ERROR, cross point z1 is lager than 21.55m '
            endif
!            print*,'Through the column above'
            go to 101
          else
            y1=-11.   !y=-11 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-22.3.and.x1.le.11.
     &         .and.z1.ge.0.) then
              flag=.true.
!              print*,'x1, y1, z1',x1,y1,z1
!              print*,'Through y=-11 plane'
              go to 101
            endif
            y1=10.    !y=10 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-22.3.and.x1.le.11.
     &         .and.z1.ge.0.) then
              flag=.true.
!              print*,'x1, y1, z1',x1,y1,z1
!              print*,'Through y=10 plane'
              go to 101
            endif
            x1=-22.3    !x=-22.3 plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-11.and.y1.le.10.
     &          .and.z1.ge.0.) then
              flag=.true.
!              print*,'Through x=-22.3 plane'
              go to 101
            endif
            x1=11.      !x=11. plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-11.and.y1.le.10.
     &          .and.z1.ge.0.) then
              flag=.true.
!              print*,'Through x=11 plane'
              go to 101
            endif
          endif 
        endif
      
101     if (flag.and.wh_hall.eq.3) then
          z2=0.
          x2=xpro/zpro*z2
          y2=ypro/zpro*z2
!          print*,'x2,y2,z2',x2,y2,z2
!          if(x2.ge.-22.3.and.x2.le.11..and.y2.ge.-11..and.y2.le.10.)
          hall_len=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        else
          print*,'theta, phi ', thetar*degpi, phir*degpi
          hall_len=0.0
        endif
      ENDIF
*--->LA site
*     x/xpro=y/ypro=z/zpro
*     (x+0.5)^2+(z)^2=14.35^2    
*     two equation together get the cross point.
*     condition 17.234<z<19.351, -22.3<y<11, -8<x<7 
      IF (wh_hall.eq.2) then
        flag=.false.
!       muon line cross with half column.
        delta=(xz)**2-4*(xz**2+1.)*(0.25-14.35**2)
        if (delta.lt.0) then
          print*,'ERROR! delta <0'
        else
          b=xz
          a=xz**2+1
          z1=(-b+sqrt(delta))/(2.*a)
          y1=yz*z1
          x1=xz*z1
!          print*,'x1, y1, z1',x1,y1,z1
*         limit the range
          if((12.234.le.z1)
     &      .and.(-22.3.le.y1.and.y1.le.11.)) then
            flag=.true.
            if (z1.gt.14.351) then
              print*,'x1, y1, z1',x1,y1,z1
              print*,'!!!!!!ERROR, cross point z1 is lager than 19.35m '
            endif
!            print*,'Through the column above'
            go to 201
          else
            y1=11.   !y=11 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-8.and.x1.le.7.
     &         .and.z1.ge.0) then
              flag=.true.
!              print*,'Through y=11 plane'
              go to 201
            endif
            y1=-22.3    !y=-22.3 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-8.and.x1.le.7.
     &         .and.z1.ge.0) then
              flag=.true.
!              print*,'Through y=-22.3 plane'
!              print*,'*******'
              go to 201
            endif
            x1=-8    !x=-8 plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-22.3.and.y1.le.11.
     &          .and.z1.ge.0) then
              flag=.true.
!              print*,'Through x=-8 plane'
              go to 201
            endif
            x1=7.      !x=7. plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-22.3.and.y1.le.11.
     &          .and.z1.ge.0) then
              flag=.true.
!              print*,'Through x=7. plane'
              go to 201
            endif
          endif 
        endif
      
201     if (flag.and.wh_hall.eq.2) then
          z2=0
          x2=xpro/zpro*z2
          y2=ypro/zpro*z2
!          if(y2.ge.-22.3.and.y2.le.11..and.x2.ge.-8..and.x2.le.7.)
          hall_len=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        else
          print*,'theta, phi ', thetar*degpi, phir*degpi
          hall_len=0.0
        endif
      ENDIF
*--->DYB site  ********************************************************
*     x/xpro=y/ypro=z/zpro
*     (x-0.5)^2+(z)^2=14.35^2    
*     two equation together get the cross point.
*     condition 12.234<z<14.35, -11<y<16, -7<x<8 
      IF (wh_hall.eq.1) then
        flag=.false.
!       muon line cross with half column.
        delta=(xz)**2-4*(xz**2+1.)*(0.25-14.35**2)
        if (delta.lt.0) then
          print*,'ERROR! delta <0'
        else
          b=-xz
          a=xz**2+1
          z1=(-b+sqrt(delta))/(2.*a)
          y1=yz*z1
          x1=xz*z1
!          print*,'x1, y1, z1',x1,y1,z1
*         limit the range
          if((12.234.le.z1)
     &      .and.(-11.le.y1.and.y1.le.16.)) then
            flag=.true.
            if (z1.gt.14.351) then
              print*,'x1, y1, z1',x1,y1,z1
              print*,'!!!!!!ERROR, cross point z1 is lager than 14.35m '
            endif
!            print*,'Through the column above'
            go to 301
          else
            y1=-11.   !y=-11 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-7.and.x1.le.8.
     &         .and.z1.ge.0) then
              flag=.true.
!              print*,'Through y=-11 plane'
              go to 301
            endif
            y1=16.    !y=16 plane
            x1=xy*y1
            z1=zy*y1
            if (x1.ge.-7.and.x1.le.8.
     &         .and.z1.ge.0) then
              flag=.true.
!              print*,'Through y=16 plane'
              go to 301
            endif
            x1=-7    !x=-7 plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-11.and.y1.le.16.
     &          .and.z1.ge.0) then
              flag=.true.
!              print*,'Through x=-7 plane'
              go to 301
            endif
            x1=8.      !x=8. plane
            y1=yx*x1
            z1=zx*x1
            if (y1.ge.-11.and.y1.le.16.
     &          .and.z1.ge.0) then
              flag=.true.
!              print*,'Through x=8 plane'
              go to 301
            endif
          endif 
        endif
      
301     if (flag.and.wh_hall.eq.1) then
          z2=0
          x2=xpro/zpro*z2
          y2=ypro/zpro*z2
!          if(x2.ge.-22.3.and.x2.le.11..and.y2.ge.-11..and.y2.le.10.)
          hall_len=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        else
          print*,'theta, phi ', thetar*degpi, phir*degpi
          hall_len=0.0
        endif
      ENDIF
!       hall_len=100.0005
!        print*,'hall_len=',hall_len
!        print*,'   '
!       print*,'wh_hall=',wh_hall
      END
