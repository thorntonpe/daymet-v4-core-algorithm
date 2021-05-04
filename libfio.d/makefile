
ARCHIVE = libfio.a
INCS    = fileio_tools.h 
CC	= gcc
CFLAGS = -O2

OBJECTS=$(ARCHIVE)(fileio_tools.o)

SOURCES = fileio_tools.c 

$(ARCHIVE): $(OBJECTS)
	ar rc $(ARCHIVE) fileio_tools.o
	mv $@ ../library.d/

$(ARCHIVE)(fileio_tools.o): $(INCS)
