# top-level makefile for daymet_compute_src directory

all : gctp_lib libfio tools interpolate predict xval

libfio : 
	cd libfio.d ; make

gctp_lib : 
	cd libgeo.d ; make 

tools : gctp_lib libfio
	cd tools.d ; make all

interpolate : libfio
	cd interpolate.d ; make all

predict : gctp_lib libfio
	cd predict.d ; make all

xval : libfio
	cd xval.d ; make all

.PHONY : clean
clean :
	rm -f bin/*
	rm -f interpolate.d/*.o
#	rm -f libfio.d/*.o
#	rm -f libgeo.d/*.o
	rm -f predict.d/*.o
	rm -f tools.d/*.o
	rm -f xval.d/*.o

