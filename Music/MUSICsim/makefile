
user_src = dywm gaisser d2_integral hrndg2 exp_hall music-sr ranlux corset corgen rnorml ranmar rm48

#user_src = test music-sr ranlux corset corgen rnorml ranmar

main = sim
CLEANUP = *.o *~

FFLAGS  = -fno-automatic  -fbounds-check

# create user list
user_object = $(foreach u_object,$(user_src),$(u_object).o)

# create executable
$(main): $(user_object) $(user_inc) makefile
	gfortran $(FFLAGS) -Wall -o $(main) $(user_object) $(user_lib) -lm

graf: $(user_object) $(user_inc) makefile
	gfortran -c $(FFLAGS) -DGRAF_TYPE -Wall 
	gfortran $(FFLAGS) -Wall -o $(main) $(user_object) $(user_lib) -lm

# resolve the dependencies
.f.o:
	gfortran -c $(FFLAGS) -Wall $*.f

# cleanup
clear:
	/bin/rm -f $(CLEANUP)


