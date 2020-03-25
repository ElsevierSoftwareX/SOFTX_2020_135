#Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %\n
#Junior University.  All rights reserved.                             %\n
#
#The programs in GSLIB are distributed in the hope that they will be  %\n
#useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %\n
#responsibility to anyone for the consequences of using them or for   %\n
#whether they serve any particular purpose or work at all, unless he  %\n
#says so in writing.  Everyone is granted permission to copy, modify  %\n
#and redistribute the programs in GSLIB, but only under the condition %\n
#that this notice and the above copyright notice remain intact.       %\n
#
# Makefile for generating "libgs.a"
#
LIBGS = libgs.a
SRCS = acorni.f backtr.f beyond.f blue.f chknam.f chktitle.f cova3.f \
	dlocate.f dpowint.f dsortem.f gauinv.f gcum.f getindx.f getz.f \
	green.f hexa.f ksol.f ktsol.f locate.f nscore.f numtext.f \
	ordrel.f picksupr.f powint.f psfill.f psline.f pstext.f rand.f \
	red.f resc.f scal.f setrot.f setsupr.f sortem.f sqdist.f \
	srchsupr.f strlen.f srchsuprsisim.f
OBJS = $(SRCS:.f=.o)

FFLAGS =-c -m64 -march=k8 -fno-second-underscore -ffixed-form -O3 -ffast-math
FFLAGS += -fbounds-check -Wline-truncation -Wuninitialized
FFLAGS += -fopenmp -fdefault-integer-8 -I./

all: $(LIBGS)

$(LIBGS): $(OBJS)
	/usr/bin/ar -q $@ $(OBJS)
	/usr/bin/ar -ts $@

.f.o:
	gfortran $(FFLAGS) $*.f

install:
	install -s -m 644 $(LIBGS) $(LIBDIR)

clean:
	/bin/rm -rf $(OBJS) $(LIBGS) core
