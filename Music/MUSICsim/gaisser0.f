      real*8 function gaisser(emu, costh)
      implicit none
      real*8 costh, emu, t115, t850
      real*8 A0, gamma
      parameter( A0 = 0.14d0, gamma = 2.7d0)
C     parameter( A0 = 0.26d0, gamma = 2.78d0)

      t115 = 1.d0  /( 1.d0 + 1.1d0*emu*costh/115.d0 )
      t850 = 0.054d0/( 1.d0 + 1.1d0*emu*costh/850.0d0 )

      gaisser = A0 * emu**(-gamma) * ( t115 + t850 )
* --> cm^2 s sr GeV

      END
