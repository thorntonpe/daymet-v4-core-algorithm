# makefile for tools.d directory
# revised to produce the minimal code set for US daymet runs

LDFLAGS += -L/usr/lib
LIBDIR = ../library.d
INCDIR = ../include.d
BINDIR = ../bin
CC = gcc
CFLAGS = -O2 -I/usr/include 

STDLIBS = -lm ${LIBDIR}/libfio.a ${LDFLAGS} -lnetcdf  
LIBS1 = ${LIBDIR}/libgeo.a 

EXECS = xylocate_nc

OBJS5 = xylocate_nc.o

INC1 = ${INCDIR}/metsrc2.h
INC2 = ${INCDIR}/cproj.h ${INCDIR}/proj.h
INC3 = ${INCDIR}/dailywx.h

all : ${EXECS}

xylocate_nc : ${OBJS5}
	${CC} -I${INCDIR} -O -o $@ ${OBJS5} ${LIBS1} ${STDLIBS} ${LDFLAGS}
	mv $@ ${BINDIR}

xylocate_nc.o : ${INC1} ${INC2} ${INC3}
xylocate_nc.o : xylocate_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -O -c xylocate_nc.c

clean:
	rm -f *.o
