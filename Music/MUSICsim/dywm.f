
      program main
*******************
      implicit none
      integer nmuon, imu,iprint,lrec,istat,icycle,npass ,ne,nth, ii,jj
      real*8 lowmu,upmu,th1,th2,dim_sum, emu0, costh0
      real depth0, th0, phi0,  costh,phi_tr
      real x0, y0, z0, tmu0
      real shortest,  Trate, emean,  bine, binth
      real xxx,yyy, zzz                                 ! the site coordinate
      real*8 x,y,z,theta,phi,emu,depth,rho,tmu     ! changed to transfer visual parameters 

      integer evtId, nsub, listId, spId
      integer evtsubn(1037393),evtlistId(1037393),evtspId(1037393)
      real muE(2005812),muth(2005812),muphi(2005812)
      real muposx(2005812),muposy(2005812)
      integer sId
      integer isub    !loop sub event num
      integer subln   !sub event live number
      integer npsub   !npass sub all events
      real posx0,posy0

                                                   ! to suit MUSIC parameters
      real a,b,c,gaisser
      real hall_lenth
      real bsize, bloc, bwgt, logup
      EXTERNAL gaisser, hrndg2, muon_transport, exp_hall
      real  qqpi,TWO_PI, hpi, qqradpi, log5
      parameter(qqpi=3.1415926536,qqradpi = 180./qqpi, log5=0.69897)
      parameter(TWO_PI = 2.0*qqpi, hpi=qqpi/2.0)

      real*8 muon_flux(4000,100)
*
      integer i,j,k, iranlux, whi_hall
      real  vec(1)
*
      integer NWPAWC
      real H
      parameter (NWPAWC =5000000)
      common /pawc/H(nwpawc)
************************************************************************
************************************************************************
*      call hlimit(NWPAWC)
      call initialize_music
      iranlux=1
      call rluxgo(3,iranlux,0,0)

*--->  in function "depth_cal",  parameters:
* 1,2 : the angle 
* 3,4 : the site coordinate 
* 5,6 : sample Layer level   &    the const of Weighed average 
**********************************************************
      open(35, file='mountain.NearS', status='unknown')
**      open(36, file='outgen-.dat', status='unknown')
**********************************************************
*---> calculate the muon Flux in different location
*
**********************************************************
*        xxx = 1102.79004   ! n1
*        yyy = 647.22998
*        xxx = 1626.        ! n2
*        yyy = 1376.
*        xxx = 1147         ! m 
*        yyy = 1280
**        xxx = 769.0         ! f
**        yyy = 2364.0 

*******************
*the new scheme:
*******************
*        xxx = 1114.77      ! N1
*        yyy = 494.043
*       xxx = 1601.19      ! N2 change
*        yyy = 1364.88
*        xxx = 746.691-18.0      ! F
*        yyy = 2107.22
*********************************************************
* New Map
*         xxx=60894.42      !DYB
*         yyy=5525.81
*         zzz=-20
*         whi_hall=1

*         xxx=61382.73      !LA
*         yyy=6391.07
*         zzz=-16.62
*         whi_hall=2

*        xxx=60500.25      !Far
*        yyy=7141.09
*        zzz=-15.15
*        whi_hall = 3             !near=1, LA=2, Far=3 if equal 0,NO exp. hall
*****************************************************
* All Map --eliipse earth
*Near nuclear plant in new map:x=60753.2131 y=5186.6915 
*         xxx=60894.42-60753.2131      !DYB
*         yyy=5525.81-5186.6915
*         xxx=61382.73-60753.2131      !LA
*         yyy=6391.07-5186.6915
*         xxx=60500.25-60753.2131      !Far
*         yyy=7141.09-5186.6915

        xxx=141.21                 !DYB
        yyy=339.12
        zzz=-20
        whi_hall=1


*         xxx=629.52             !LA
*         yyy=1204.38
*         zzz=-16.62
*         whi_hall=2



*         xxx=-252.96          !Far origin
*         yyy=1954.40
*         zzz=-15.15
*         whi_hall = 3  

*         xxx=-252.96+41.3          !Far
*         yyy=1954.40-73.1
*         zzz=-15.15-0.25
*         whi_hall = 3
**********************************************************
        iprint = 1
        nmuon  = 100 !  changed the numbers to be tracked 
        rho  = 2.6d0
        tmu0 = 0.
        x0   = 0.
        y0   = 0.
        z0   = 0.

* --> book the muon Flux with energy less than 5000GeV.  
**********************************
        a=100
        b=0.0
        c= 0.37  ! 0.57   ! 0.8   !  0.37
        shortest = 1000;
        print*,'--- Shortest OverBurden: ',shortest,'  meter;  w_f: ',c 
        lowmu = 500.0*(exp(shortest*rho/2500.0)-1.0)       ! // see Gaisser's book.
        lowmu = lowmu-(lowmu/45.0)*4.0
*         lowmu=int(500.0*(exp(shortest*rho/2500.0)-1.0)/5.0)*5.0 - 3
        if(lowmu.lt.0) lowmu = 0
        lowmu = 40.0
        upmu  = lowmu+10000.d0                !  // change
        th1 = 0.00366518 !89.79 !0.25881906667075   ! 75deg;   --->  cos(70deg) = 0.3420202d0
        th2 = 1.d0
        ne = 400
        nth = 10
*         call hrndg2(muon_flux,ne,lowmu,upmu,nth,th1,th2
*     &            ,dim_sum, emu0, costh0, 1.)  
*          Trate = dim_sum*qqpi*2.0*1.0e4    ! intergral of phi give a 2pi, and convert to /m^2
         print*,'Energy range, Flux ', lowmu,upmu, Trate

        Trate = 184*0.0003256 !if lowmuE=0.105, 177Hz, if lowmu=0, 184Hz, 400GeV-0.0003256
         print*,'Energy range, Flux ', lowmu,upmu, Trate
************ read muon sample on surface from corsika ***************

        open(80, file='muon.txt',status='unknown')
        do i=1, 2005812
         read (80,*) evtId, nsub, listId, spId, 
     &       muE(i),muth(i),muphi(i),muposx(i),muposy(i)
         evtsubn(evtId) = nsub 
         evtlistId(evtId) = listId
         evtspId(evtId) = spId
        enddo
        close(80)
        print*, 'read muon.txt OK!'
*       do i=1, 1037393
*        evtId = evtlistId(i)-evtsubn(i)+1
*        print*,i,evtsubn(i),evtId,muE(evtId),muth(evtId)
*       enddo


**********************************************************************

        emean = 0.0
        npass = 0
        npsub = 0
        imu = 0
301     imu = imu + 1
        if (mod(imu,100000).eq.0)then
      print'("^^^^^^^^^^^^^^^^",2i10,"    ^^^^^^^^^^^^^^^^")',imu,npass
        endif

*--> get random numbers from the distribution
**********************************************
*        call hrndg2(muon_flux,ne,lowmu,upmu,nth,th1,th2
*     &            ,dim_sum, emu0, costh0, 111.)           ! the last par --> precision (>10?)
*        call RANLUX(vec,1)
******* get muon samples from the file ***************
401     call RANLUX(vec,1)
        sId = INT(vec(1)*1037393)
        if(sId.lt.1.or.sId.gt.1037393) then
          goto 401
        endif
        evtId = evtlistId(sId)-evtsubn(sId)+1
*        print*,sId,evtsubn(sId),evtId,muE(evtId),muth(evtId)
        subln = 0
        do isub=1 , evtsubn(sId)
*        print*,'isub: ',imu,isub,evtsubn(sId)
        if((evtId+isub-1).gt.2005812) then
         print*,'ERROR:',evtId+isub-1,' is larger than 2005812'
        endif
        emu0 = muE(evtId+isub-1)
        costh0 = cos(muth(evtId+isub-1)/qqradpi)
        phi0 = muphi(evtId+isub-1)/qqradpi
        posx0 = muposx(evtId+isub-1)
        posy0 = muposy(evtId+isub-1)
        
*        phi0 = vec(1)*TWO_PI
  
102     format(1x, 3f12.5)
**       write(36,102) emu0,costh0   !     ,phi0     

*-->  calculate the muon track length through the mountain
*   muon direction changed with different coordinate system
************************************
        costh = costh0
        phi_tr = phi0
        if(phi_tr.le.3.0*hpi) then 
          phi_tr = 3.0*hpi - phi_tr     ! phi_tr<=270 , phi_tr+= 270-phi_tr
        else
          phi_tr = 7.0*hpi - phi_tr     ! phi_tr >270 , phi_tr+= 640-phi_tr
        endif  
*--------> the last parameter stand for    the sample Layer level: "C" 
*------->! xxx,yyy  given the site coordinate
        depth0 = 1000;
        if(depth0.lt.0.0)then
          depth0 = 10
          print*,"------ Error! depth is little than Zero!" 
        endif
**       print*,'depth0,costh0, phi_tr ',depth0,costh0, phi_tr*qqradpi
 
        hall_lenth=0.0
!         print*,'costh0==========================================',costh0
        call exp_hall(hall_lenth,costh0,phi_tr,whi_hall)
!         print*,'costh==========================================',costh
!         print*,'hall_lenth=',hall_lenth
!         print*,'                    '
!         print*,'                    '
        depth0 = depth0 - hall_lenth 
*--> call MUSIC transport
***********************************

        x   = x0
        y   = y0
        z   = z0
        emu = emu0
        tmu = tmu0                    ! ns
        depth = depth0*100.           ! change to cm
        theta = dacos(costh0)
        phi   = phi0 
  
******************

        call muon_transport
     &    (x,y,z,theta,phi,emu,depth,rho,tmu)

  	     theta = theta*qqradpi
  	     phi   = phi*qqradpi
         th0   = dacos(costh0)*qqradpi
         phi0  = phi0*qqradpi
         depth = depth/100.

*--> save data
************************************
101     format(1x,4i10, 9f12.5)
        if(emu.gt.0.0) then !emu arrive to 0
          subln = subln + 1
        write(35,101) npass+1,evtsubn(sId),subln,evtspId(sId),emu,
     &        theta,phi,emu0,th0,phi0,posx0,posy0,depth
          emean = emean + emu
          npsub = npsub + 1
        endif
        enddo
        if(subln.gt.0) then
          npass = npass + 1
        endif

        if(npass.lt.nmuon) goto 301     ! changed to track 300k effective muon event 

        Trate = Trate/imu
        if (whi_hall.eq.3) then
          write(35,*) ' Which site        :  Far'
        else if (whi_hall.eq.2) then
          write(35,*) ' which site        :  LA'
        else if (whi_hall.eq.1) then
          write(35,*) ' which site        :  DYB'
        endif

        write(35,*) 'Total events       : ', nmuon,imu
        write(35,*) 'Emu   bins, range  : ', ne, lowmu, upmu
        write(35,*) 'costh bins, range  : ', nth, th1,th2
        write(35,*) 'Normalization      : ', Trate
        write(35,*) 'Total muon rate    : ', npass*trate
        write(35,*) 'Mean energy of muon: ', emean/npsub
        close(35)
*        close(36)
************************************************************************
        if (whi_hall.eq.3) then
          write(*,*) ' Which site        :  Far'
        else if (whi_hall.eq.2) then
          write(*,*) ' which site        :  LA'
        else if (whi_hall.eq.1) then
          write(*,*) ' which site        :  DYB'
        endif
        write(*,*) ' Total events      : ', nmuon,imu
        write(*,*) ' Emu   bins, range : ',  ne, lowmu, upmu
        write(*,*) ' costh bins, range : ',  nth, th1,th2
        write(*,*) ' Normalization     : ', Trate
        write(*,*) ' Total muon rate    : ', npass*trate
        write(*,*) ' Mean energy of muon: ', emean/npsub


      END

************************************************************************
************************************************************************
