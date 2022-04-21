g++ -c music.cpp

g++ -o music MUSICsim/gaisser.o MUSICsim/d2_integral.o MUSICsim/hrndg2.o MUSICsim/exp_hall.o MUSICsim/music-sr.o MUSICsim/ranlux.o MUSICsim/corset.o MUSICsim/corgen.o MUSICsim/rnorml.o MUSICsim/ranmar.o MUSICsim/rm48.o music.o -lgfortran

./music
