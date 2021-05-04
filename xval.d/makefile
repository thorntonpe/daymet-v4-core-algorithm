# makefile for cross validation (xval) directory
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

EXECS = cmd_predict_daily_tair_xval3_nc cmd_predict_daily_prcp_xval3_nc cmd_prcp_ndays_filter_nc \
fill_daily_tair3_nc fill_daily_prcp3_nc cmd_predict_daily_tair_xval_filter3_nc cmd_tair_ndays_filter_nc

INC0 = ${INCDIR}/metsrc2.h
INC1 = ${INCDIR}/metsrc.h
INC2 = ${INCDIR}/cproj.h ${INCDIR}/proj.h
INC3 = ${INCDIR}/dailywx.h

OBJS19 = ../interpolate.d/interpolate_func.o \
../tools.d/convolve.o xval_daily_tair_func3.o cmd_predict_daily_tair_xval3_nc.o

OBJS20 = ../interpolate.d/interpolate_func.o ../predict.d/predict_daily_prcp_func3.o \
../tools.d/convolve.o xval_daily_prcp_func3.o cmd_predict_daily_prcp_xval3_nc.o

OBJS21 = ../interpolate.d/interpolate_func.o ../predict.d/predict_daily_tair_func3.o \
../tools.d/convolve.o xval_daily_tair_func3.o fill_daily_tair3_nc.o

OBJS22 = ../interpolate.d/interpolate_func.o ../predict.d/predict_daily_prcp_func3.o \
../tools.d/convolve.o xval_daily_prcp_func3.o fill_daily_prcp3_nc.o

OBJS23 = ../interpolate.d/interpolate_func.o ../predict.d/predict_daily_tair_func3.o \
../tools.d/convolve.o xval_daily_tair_func3.o cmd_predict_daily_tair_xval_filter3_nc.o

all : ${EXECS}

cmd_predict_daily_tair_xval3_nc : interp_obj  convolve_obj \
xval_daily_tair_func3.o cmd_predict_daily_tair_xval3_nc.o
	${CC} -I${INCDIR} -o $@ ${OBJS19} ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

cmd_predict_daily_prcp_xval3_nc : interp_obj predict_day_prcp3_obj convolve_obj \
xval_daily_prcp_func3.o cmd_predict_daily_prcp_xval3_nc.o
	${CC} -I${INCDIR} -o $@ ${OBJS20} ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

cmd_prcp_ndays_filter_nc : cmd_prcp_ndays_filter_nc.o
	${CC} -I${INCDIR} -o $@ cmd_prcp_ndays_filter_nc.o ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

cmd_tair_ndays_filter_nc : cmd_tair_ndays_filter_nc.o
	${CC} -I${INCDIR} -o $@ cmd_tair_ndays_filter_nc.o ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

fill_daily_tair3_nc : interp_obj predict_day_tair3_obj convolve_obj \
xval_daily_tair_func3.o fill_daily_tair3_nc.o
	${CC} -I${INCDIR} -o $@ ${OBJS21} ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

fill_daily_prcp3_nc : interp_obj predict_day_prcp3_obj convolve_obj \
xval_daily_prcp_func3.o fill_daily_prcp3_nc.o
	${CC} -I${INCDIR} -o $@ ${OBJS22} ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

cmd_predict_daily_tair_xval_filter3_nc : interp_obj predict_day_tair3_obj convolve_obj \
xval_daily_tair_func3.o cmd_predict_daily_tair_xval_filter3_nc.o
	${CC} -I${INCDIR} -o $@ ${OBJS23} ${LIBS1} ${STDLIBS}
	mv $@ ${BINDIR}

interp_obj :
	cd ../interpolate.d ; make interpolate_func.o

predict_day_prcp3_obj :
	cd ../predict.d ; make predict_daily_prcp_func3.o

predict_day_tair3_obj :
	cd ../predict.d ; make predict_daily_tair_func3.o

convolve_obj :
	cd ../tools.d ; make convolve.o

gctp_tools_nc.o : ${INC1} ${INC2} ${INC3}
gctp_tools_nc.o : gctp_tools_nc.c
	$(CC) -I${INCDIR} ${CFLAGS} -c gctp_tools_nc.c

cmd_predict_daily_tair_xval3_nc.o : ${INC0}
cmd_predict_daily_tair_xval3_nc.o : cmd_predict_daily_tair_xval3_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c cmd_predict_daily_tair_xval3_nc.c

cmd_predict_daily_prcp_xval3_nc.o : ${INC0}
cmd_predict_daily_prcp_xval3_nc.o : cmd_predict_daily_prcp_xval3_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c cmd_predict_daily_prcp_xval3_nc.c

cmd_prcp_ndays_filter_nc.o : ${INC0}
cmd_prcp_ndays_filter_nc.o : cmd_prcp_ndays_filter_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c cmd_prcp_ndays_filter_nc.c

cmd_tair_ndays_filter_nc.o : ${INC0}
cmd_tair_ndays_filter_nc.o : cmd_tair_ndays_filter_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c cmd_tair_ndays_filter_nc.c

xval_daily_tair_func3.o : ${INC0}
xval_daily_tair_func3.o : xval_daily_tair_func3.c
	${CC} -I${INCDIR} ${CFLAGS} -c xval_daily_tair_func3.c
	
xval_daily_prcp_func3.o : ${INC0}
xval_daily_prcp_func3.o : xval_daily_prcp_func3.c
	${CC} -I${INCDIR} ${CFLAGS} -c xval_daily_prcp_func3.c

fill_daily_tair3_nc.o : ${INC0}
fill_daily_tair3_nc.o : fill_daily_tair3_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c fill_daily_tair3_nc.c

fill_daily_prcp3_nc.o : ${INC0}
fill_daily_prcp3_nc.o : fill_daily_prcp3_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c fill_daily_prcp3_nc.c

cmd_predict_daily_tair_xval_filter3_nc.o : ${INC0}
cmd_predict_daily_tair_xval_filter3_nc.o : cmd_predict_daily_tair_xval_filter3_nc.c
	${CC} -I${INCDIR} ${CFLAGS} -c cmd_predict_daily_tair_xval_filter3_nc.c

clean:
	rm -f *.o

