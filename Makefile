#
# Solar system simulation using particle-particle method
# CIS 455/555 WQ'06
# Configured for IBM p690
#					

# Uncomment one of these lines to specify a compiler -- use g++
# for the sequential program, mpic++ for parallel.

CXX = g++
# CXX = mpic++

nbody:	nbody.o Vector.o Body.o
	$(LINK.cc) -o nbody nbody.o Vector.o Body.o

mpinbody:	mpinbody.o Vector.o Body.o
	$(LINK.cc) -o mpinbody mpinbody.o Vector.o Body.o

vdemo:	vdemo.o Vector.o
	$(LINK.cc) -o vdemo vdemo.o Vector.o

clean:
	/bin/rm -rf core *.o *~ *.bak ii_files


