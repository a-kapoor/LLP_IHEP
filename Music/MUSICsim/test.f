
      program main
      implicit none
      real*8 x,y,z,theta,phi,emu,depth,rho,tmu     ! changed to transfer visual parameters 
      rho = 2.6d0
      depth = 100.
      theta = 0 
      phi = 0
      emu = 2000

      print*,'--- start an empty test code : ' 

      print*,'--- initialize music : ' 
      call initialize_music
      print*,'--- call music : ' 
      call muon_transport
     &    (0,0,0,theta,phi,emu,depth,rho,0)
      print*,'--- finish music : ' 

      END
