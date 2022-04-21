      real*8 function gaisser(emu, costh0)
      implicit none
      real*8 costh,costh0, emu,en1, t115, t850
      real*8 A0, gamma
      real*8 p1,p2,p3,p4,p5;
      parameter( A0 = 0.14d0, gamma = 2.7d0)
C     parameter( A0 = 0.26d0, gamma = 2.78d0)

      p1 = 0.102573;
      p2 = -0.068287;
      p3 = 0.958633;
      p4 = 0.0407253;
      p5 = 0.817285;
      costh = sqrt((costh0*costh0+p1*p1+p2*costh0**p3+p4*costh0**p5)
     &                                /(1+p1*p1+p2+p4));

      en1= emu*(1+3.63698/(emu*costh**1.29685));

      t115 = 1.d0  /( 1.d0 + 1.1d0*emu*costh/115.d0 )
      t850 = 0.054d0/( 1.d0 + 1.1d0*emu*costh/850.0d0 )

      gaisser = A0 * en1**(-gamma) * ( t115 + t850 )
* --> cm^2 s sr GeV

      END
