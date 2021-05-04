#ifndef __FILEIO_TOOLS_HEADER__
  #define __FILEIO_TOOLS_HEADER__

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations for fileio_tools.c */

/* structure definition for custom file handling */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>


typedef struct
{
  char name[128];
  FILE * ptr;
}
file;


extern int file_open( file * target, char mode );

extern int scan_value( file ini, void * var, char mode );

extern int scan_open( file ini, file * target, char mode );

extern int bin_out( file * fptr, char * outprefix, char * outsuffix );


#ifdef __cplusplus
}
#endif


#endif //__FILEIO_TOOLS_HEADER__
