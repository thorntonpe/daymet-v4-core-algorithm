# makefile for interpolate.d directory
# this version generates the minimal code set for US daymet runs

LDFLAGS += -L/usr/lib 

# - No modification necessary below this line ------------------------------

LIBDIR = ../library.d
INCDIR = ../include.d
BINDIR = ../bin
CC = gcc
CFLAGS = -O2 -I/usr/include 

STDLIBS = -lm ${LIBDIR}/libfio.a ${LDFLAGS} -lnetcdf

EXECS = interpolate2_nc

OBJS2 = interpolate2_nc.o interpolate_func.o

INC1 = ${INCDIR}/metsrc.h
INC2 = ${INCDIR}/cproj.h ${INCDIR}/proj.h

all : ${EXECS}

interpolate2_nc : ${OBJS2}
	$(CC) -I${INCDIR} -O -o $@ ${OBJS2} ${STDLIBS} ${LDFLAGS}
	mv $@ ${BINDIR}

interpolate2_nc.o : ${INC1}
interpolate2_nc.o : interpolate2_nc.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c interpolate2_nc.c

interpolate_func.o : ${INC1}
interpolate_func.o : interpolate_func.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c interpolate_func.c

clean:
	rm -f *.o
