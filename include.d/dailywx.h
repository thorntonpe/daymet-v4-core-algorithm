/*
dailywx.h
Peter Thornton
12/30/96

Header file for NCDC and SNOTEL extract, append, etc. code
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>

#define PHDRTAG "PHDRTAG"

int mondays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

char *mon_name[12] = 
{"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"};

char *list_name[14] = 
{"TMAX","TMIN","TOBS","PRCP","EVAP","MXPN","MNPN","WDMV","SNOW","SNWD",
"WTEQ","SX","SN","SO"};

char *stname[13] = 
{"ak","az","ca","co","id","mt","nm","nv","or","sd","ut","wa","wy"};

char *list_desc[14] = 
{"Daily maximum temperature","Daily minimum temperature",
"Temperature at observation time","Daily precipitation",
"Daily evaporation",
"Daily maximum temperature of water in evaporation pan",
"Daily minimum temperature of water in evaporation pan",
"24-hour wind movement","Daily snowfall","Snow depth at observation time",
"Water equivalent of snow depth","Daily maximum soil temperature",
"Daily minimum soil temperature","Soil temperature at time of observation"};

/* the qual field indicates the degree of agreement between the
station history file and the station data...
1 = year and month match
2 = year matches
3 = year matches within +/- 5 years */
typedef struct
{
	char name[24];
	char id[8];
	short year;
	char month;
	char type;
	char quality;
	float stnlat;
	float stnlon;
	float stnelev;
} monmetastr;

typedef struct
{
	char isgood[32];
	short daydat[32];
} mondatastr;
	
typedef struct
{
	monmetastr meta;
	mondatastr data;
} monrecstr;

typedef struct
{
	char tag[8];
	char elem_name[5];
	int year;
	long int record_size;
	int nstns[12];
} avsmhdrstr;

struct state_link
{
	int stateid;
	struct station_list *stnlst;
	struct state_link *next;
};

struct state_list
{
	struct state_link *head;
	struct state_link *tail;
};

struct station_link
{
	int stnid;
	int month_offset[12];
	char textid[8];
	struct station_link *prev;
	struct station_link *next;
};

struct station_list
{
	struct station_link *head;
	struct station_link *tail;
	struct station_link *try;
};

struct permout_hdr
{
	char tag[8];
	char element[5];
	int year;
	int nstations;
	int nstations_state[100];
	int state_offset[100];
	int nstations_month[12];
	int month_offset[12];
};

struct state_hdr
{
	int stateid;
	int nstations;
};

struct index
{
	char stnid[8];
	int offset[12];
};
