#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "fileio_tools.h"

int file_open( file * target, char mode );
int file_open (file *target, char mode)
/*  'r' for read binary
	'i' for read ascii
	'w' for write binary
	'o' for write ascii */
{
	int flag = 0;

	switch (mode)
	{
        case 'r':
            if ((target->ptr = fopen(target->name,"rb")) == NULL)
            {
                printf("Can't open %s for binary read ... Exiting\n",target->name);
                flag = 1;
            }
            break;

        case 'i':
            if ((target->ptr = fopen(target->name,"r")) == NULL)
            {
                printf("Can't open %s for ascii read ... Exiting\n",target->name);
                flag = 1;
            }
            break;

        case 'w':
            if ((target->ptr = fopen(target->name,"wb")) == NULL)
            {
                printf("Can't open %s for binary write ... Exiting\n",target->name);
                flag = 1;
            }
            break;

        case 'o':
            if ((target->ptr = fopen(target->name,"w")) == NULL)
            {
                printf("Can't open %s for ascii write ... Exiting\n",target->name);
                flag = 1;
            }
            break;

        default:
            printf("Invalid mode specification for file_open ... Exiting\n");
            flag = 1;
    }
    return(flag);
}


int scan_value (
    file ini,       /* initialization file containing ascii data */
    void *var,      /* pointer for returning value */
    char type)      /* type to read from ini:
                        'i' for integer
                        'd' for double
                        's' for string */

{
    int ok_scan;
    int flag = 0;

    switch (type)
    {
        case 'i':
            ok_scan = fscanf(ini.ptr, "%d%*[^\n]",var);
            if (ok_scan == 0 || ok_scan == EOF)
			{
				printf("Error reading int value from %s ... exiting\n",ini.name);
				flag = 1;
			}
            break;

        case 'd':
            ok_scan = fscanf(ini.ptr, "%lf%*[^\n]",var);
            if (ok_scan == 0 || ok_scan == EOF)
			{
				printf("Error reading double value from %s... exiting\n",ini.name);
				flag = 1;
			}
            break;

        case 's':
            ok_scan = fscanf(ini.ptr, "%s%*[^\n]",var);
            if (ok_scan == 0 || ok_scan == EOF)
			{
				printf("Error reading string value from %s... exiting\n",ini.name);
				flag = 1;
			}
            break;

        default:
            printf("Invalid type specifier for scan_value ... Exiting\n");
            flag = 1;
    }
    return(flag);
}


int scan_open (
	file ini,		/* FILE struct pointer to initialization file*/
	file *target,	/* pointer to FILE struct: file to open */
	char mode)		/*  'r' for read binary
						'i' for read ascii
						'w' for write binary
						'o' for write ascii */
{
	int flag = 0;

	if (scan_value(ini,target->name,'s'))
	{
		printf("Error reading filename from %s... Exiting\n",ini.name);
		flag = 1;
		return(flag);
	}
	flag = file_open(target,mode);
	return(flag);
}

int bin_out(file* fptr, char* outprefix, char* outsuffix)
/* create filenames and open files for binary output */
{
	int ok=1;
	strcpy(fptr->name, outprefix);
	strcat(fptr->name, outsuffix);
	if (file_open(fptr, 'w')) ok=0;

	return (!ok);
}
