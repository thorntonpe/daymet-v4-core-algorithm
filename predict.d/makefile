# makefile for predict directory
# revised to produce the minimal code set for US daymet runs
LDFLAGS += -L/usr/lib

# - No modification necessary below this line ------------------------------

LIBDIR = ../library.d
INCDIR = ../include.d
BINDIR = ../bin
CC = gcc
CFLAGS = -O2 -I/usr/include

STDLIBS = -lm ${LIBDIR}/libfio.a ${LDFLAGS} -lnetcdf
LIBS1 = ${LIBDIR}/libgeo.a 

EXECS = predict_daily_srad_vp_nc predict_daily_tair3_nc predict_daily_prcp3_nc

INC2 = ${INCDIR}/cproj.h ${INCDIR}/proj.h
INC3 = ${INCDIR}/metsrc2.h

OBJS4 = predict_daily_srad_vp_nc.o predict_daily_srad_vp_func.o yearday_tools.o

OBJS7 = predict_daily_tair3_nc.o predict_daily_tair_func3.o yearday_tools.o

OBJS8 = predict_daily_prcp3_nc.o predict_daily_prcp_func3.o yearday_tools.o

all : ${EXECS}

predict_daily_srad_vp_nc : ${OBJS4}
	$(CC) -I${INCDIR} -O3 -o $@ ${OBJS4} ${STDLIBS} ${LIBS1}
	mv $@ ${BINDIR}

predict_daily_tair3_nc : ${OBJS7}
	$(CC) -I${INCDIR} -O -o $@ ${OBJS7} ${STDLIBS} 
	mv $@ ${BINDIR}

predict_daily_prcp3_nc : ${OBJS8}
	$(CC) -I${INCDIR} -O -o $@ ${OBJS8} ${STDLIBS}
	mv $@ ${BINDIR}

predict_daily_srad_vp_nc.o : ${INC1} ${INC2}
predict_daily_srad_vp_nc.o : predict_daily_srad_vp_nc.c
	$(CC) -I${INCDIR} ${CFLAGS} -O3 -c predict_daily_srad_vp_nc.c

predict_daily_srad_vp_func.o : ${INC1}
predict_daily_srad_vp_func.o : predict_daily_srad_vp_func.c
	$(CC) -I${INCDIR} ${CFLAGS} -O3 -c predict_daily_srad_vp_func.c

predict_daily_tair3_nc.o : ${INC3}
predict_daily_tair3_nc.o : predict_daily_tair3_nc.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c predict_daily_tair3_nc.c

predict_daily_tair_func3.o : ${INC3}
predict_daily_tair_func3.o : predict_daily_tair_func3.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c predict_daily_tair_func3.c

predict_daily_prcp3_nc.o : ${INC3}
predict_daily_prcp3_nc.o : predict_daily_prcp3_nc.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c predict_daily_prcp3_nc.c

predict_daily_prcp_func3.o : ${INC3}
predict_daily_prcp_func3.o : predict_daily_prcp_func3.c
	$(CC) -I${INCDIR} ${CFLAGS} -O -c predict_daily_prcp_func3.c

yearday_tools.o : yearday_tools.c
	$(CC) -O3 -c yearday_tools.c

clean:
	rm -f *.o
