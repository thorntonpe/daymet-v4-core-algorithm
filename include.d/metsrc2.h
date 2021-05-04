/*
metsrc.h
Peter Thornton
9/4/96

Last modified : 1/24/97

Structures and constants used by various daymet and avgclim programs
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#define IMGTAG "IMGTAG0"
#define LSTTAG "LSTTAG0"
#define STNTAG "STNTAG0"
#define DAYTAG "DAYTAG0"
#define AVGTAG "AVGTAG0"
#define TAIRTAG "TAIRTAG0"
#define PRCPTAG "PRCPTAG0"
#define SRADTAG "SRADTAG0"
#define VPTAG "VPTAG0"
#define DAYLTAG "DAYLTAG0"
#define SWETAG "SWETAG0"

#define DFGSP 1.5        /* density filter gaussian shape parameter */
#define DFNIT 3          /* number of iterations for density filter */
#define DFNSM 1.0        /* number of stations multiplier for density filter */
#define PI 	3.141592653589793238

#define DEFAULT_TVSZ -0.005   /* default air temperature lapse rate (C/m) */


/* structure definition for custom file handling */
typedef struct
{
	char name[128];
	FILE *ptr;
} file;

/* General Cartographic Transformation Package (GCTP) interface structure */
typedef struct
{
	long outsys;
	long outzone;
	long outdatum;
	double outparm[15];
	char fn27[64];
	char fn83[64];
} gctp_struct;

/* image header structure definition */
typedef struct
{
	char tag[8];      /* used to verify proper header */
	int ncols;        /* number of columns in the image */
	int nrows;        /* number of rows in the image */
	long typesize;    /* sizeof(type) used to store image */
	double ulx;       /* upper-left X projection coordinate (usu. meters) */
	double uly;       /* upper-left Y projection coordinate (usu. meters) */
	double cellsize;  /* gridcell width (meters) */
	gctp_struct gctp; /* GCTP projection parameter structure */
} imghdr_struct;

/* (x,y,z) list header structure definition */
typedef struct
{
	char tag[8];      /* used to verify proper header */
	int ncells;       /* number of cells stored in file */
	int nvals;        /* number of values stored per cell */
	long typesize;    /* sizeof(type) used to store list */
	gctp_struct gctp; /* GCTP projection parameter structure */
} lsthdr_struct;

/* station meta-data header structure */
typedef struct
{
	char tag[8];      /* used to verify proper header */
	int nstns;        /* number of stations in the station metadata list */
	int xyloc_flag;   /* indicates if xylocate has been run on stations */
	int zfilt_flag;   /* indicates if zfilter has been run on stations */
	int xvalfilt_flag; /* indicates if xval filtering has been run */
	double xy_width;  /* width for relocation filter (m) */
	double z_width;   /* width for elevation smoothing filter (m) */
	double critp;     /* proportion of days with good data to retain */
	gctp_struct gctp; /* GCTP projection parameter structure */
	imghdr_struct dem; /* structure describing DEM used for filtering */
} stnmetahdr_struct;

typedef struct
{
	/* original (official) location and elevation data */
	double lon;       /* recorded longitude, decimal degrees, +=E, -=W */
	double lat;       /* recorded latitude, decimal degrees, +=N, -=S */
	double x;         /* projection coordinate from recorded location */
	double y;         /* projection coordinate from recorded location */
	double z;         /* recorded elevation, meters */

	/* filtered location and elevation data */
	double fx;        /* relocated projection coordinate */
	double fy;        /* relocated projection coordinate */
	double fz;        /* filtered elevation, meters */
	
	/* station classification data */
	int type;         /* 0=NWS, 1=SNOTEL, others undefined */
	int uclass1;      /* first user-defined classification */
	int uclass2;      /* second user-defined classification */
	
	/* station identification data */
	char code[8];     /* official code or ID string, null terminated */
	char name[64];    /* name string, null terminated */
} stnmeta_struct;

typedef struct
{
	char tag[8];      /* used to verify proper header */
	int nstns;        /* number of stations in data file */
	int ndays;        /* number of days per station, if daily data */
	int start_day;    /* yearday of first day of data (base 0) */
	int start_year;   /* year of first day of data */
	int fill_flag;    /* 1=missing data filled in, 0=missing data unfilled */
} daydathdr_struct;

typedef struct
{
	char tag[8];      /* used to verify proper header */
	int nstns;        /* number of stations in data file */
} avgdathdr_struct;

typedef struct
{
	char tag[9];      /* used to verify proper header */
	char noz;         /* set to 1 for output with no elevation corrections */
	double valmax;    /* maximum data value for scaling algorithm */
	double valmin;    /* minimum data value for scaling algorithm */
	int bytmax;       /* maximum scaled byte value */
	int bytmin;       /* minimum scaled byte value */
	int start_yday;   /* yearday of first day of data (1-366) */
	int start_year;   /* year of first day of data (e.g. 1993) */
	int ndays;        /* days of data stored in the output file */
	int ioflag;       /* 1 == MAP io, 0 == LIST io */
	imghdr_struct maskhdr;  /* header information for the mask image used */
	lsthdr_struct listhdr;  /* header information for the list used */
} dayouthdr_struct;

typedef struct
{
	double value;     /* data element (either day or averaging period) */
	int code;         /* integer code for value (yday, month, etc.) */
	int mflag;        /* 1 = original data missing, 0 = original data */
} dat_struct;

/* interpolation data structure */
typedef struct
{
	double* stnx;     /* list of station x coordinates (meters) */
	double* stny;     /* list of station y coordinates */
	double* sqdist;   /* list of squared distances to each station */
	double* wt;       /* list of interpolation weights */
	short* id;        /* list of interpolation station identifiers */
	double x;         /* x-coordinate for interpolation point */
	double y;         /* y-coordinate for interpoaltion point */
	double gsp;       /* unitless gaussian shape parameter */
	double ans;       /* average number of stations to use in interpoaltion */
	double trunc;     /* filter truncation parameter */
	double sqrad;     /* squared radius of filter */
	double dfsqrad;   /* squared radius of density filter */
	double inv_dfsqrad;  /* inverse of dfsqrad */
	double dfarea;    /* area under density filter */
	double dftrunc;   /* density filter truncation parameter */
	double dfavgwt;   /* average weight under density filter */
	int nstns;        /* total number of stations to search */
	short count;      /* number of stations returned in interpolation lists */
} interpolate_struct;

typedef struct
{
	double* stnobs;   /* list of observed average (or total) precipitation */
	double* stnz;     /* list of precip station elevations */
	double* listobs;  /* interpolation list observed precip */
	double* listz;    /* interpolation list elevations */
	double* regx;     /* weighted regression x values */
	double* regy;     /* weighted regression y values */
	double* regwt;    /* weighted regression weights */
	double* wt;       /* interpolation list weights */
	short* id;        /* interpoaltion list identifiers */
	double pp;        /* predicted precipitation */
	double ppnoz;     /* predicted precipitation without elevation correction */
	double elev;      /* elevation of prediction point */
	double slope;     /* slope of the precip-to-elevation regression */
	double f_max;     /* parameter limiting elevation extrapolations */
	int nregr;        /* number of regression points */
	int nstns;        /* total number of stations */
	int switch_init;  /* initial value for dual sign-switching algorithm */
	short count;      /* number of station in interpolation list */
	char noz;         /* set to 1 for output with no elevation corrections */
} avgclim_prcp_struct;

typedef struct
{
	double* stnobs;   /* list of observed daily precipitation */
	double* stnsmobs; /* list of smoothed daily precip observations */
	double* stnx;     /* list of precip station x coordinates */
	double* stny;     /* list of precip station y coordinates */
	double* stnz;     /* list of precip station elevations */
	double* listobs;  /* interpolation list observed precip */
	double* listsmobs;/* interpolation list of smoothed observed precip */
	double* listx;    /* interpolation list x-coordinates */
	double* listy;    /* interpolation list y-coordinates */
	double* listz;    /* interpolation list elevations */
	double** regx;     /* weighted regression x values */
	double* regy;     /* weighted regression y values */
	double* regwtall; /* weighted regression weights (all points) */
	double* regwt;    /* weighted regression weights (points with precip) */
	double* wt;       /* interpolation list weights */
	short* id;        /* interpoaltion list identifiers */ 
	double pp;        /* predicted precipitation */
	double ppnoz;     /* predicted precipitation without elevation correction */
	double mapx;      /* projected x-coordinate of prediction point */
	double mapy;      /* projected y-coordinate of prediction point */
	double elev;      /* elevation of prediction point */
	double coef[4];   /* slopes and int of the t_air-to-(x,y,z) regression */
	double f_max;     /* parameter limiting elevation extrapolations */
	long int dam_off; /* offset into day-major data file */
	int nregr;        /* number of regression points */
	int ndays;        /* number of days of simulation */
	int nstns;        /* total number of stations */
	int switch_init;  /* initial value for dual sign-switching algorithm */
	short count;      /* number of station in interpolation list */
	char noz;         /* set to 1 for output with no elevation corrections */
} daymet_prcp_struct;

typedef struct
{
	double* stnobs;   /* list of observed daily air temperature */
	double* stnsmobs; /* list of smoothed daily t_air observations */
	double* stnx;     /* list of t_air station projected x coordinates */
	double* stny;     /* list of t_air station projected y coordinates */
	double* stnz;     /* list of t_air station elevations */
	double* listobs;  /* interpolation list observed t_air */
	double* listsmobs;/* interpolation list of smoothed observed t_air */
	double* listx;    /* interpolation list x-coordinates */
	double* listy;    /* interpolation list y-coordinates */
	double* listz;    /* interpolation list elevations */
	double** regx;     /* weighted regression x values */
	double** regdx;     /* weighted regression x values */
	double* regy;     /* weighted regression y values */
	double* regwtall; /* weighted regression weights (all points) */
	double* regwt;    /* weighted regression weights (points in regression) */
	double* wt;       /* interpolation list weights */
	double x1m;
	double x2m;
	double x3m;
	double wm1;
	double wm;
	double ssx1;
	double ssx2;
	double ssx3;
	double** inv;
	double sigx[3];
	short* id;        /* interpoaltion list identifiers */  
	double ta;        /* predicted daily air temperature */
	double tanoz;     /* predicted t_air, no elevation correction */
	double mapx;      /* projected x-coordinate of prediction point */
	double mapy;      /* projected y-coordinate of prediction point */
	double elev;      /* elevation of prediction point */
	double coef[4];   /* slopes and int of the t_air-to-(x,y,z) regression */
	long int dam_off; /* offset into day-major data file */
	int nregr;        /* number of regression points */
	int ndays;        /* number of days of simulation */
	int nstns;        /* total number of stations */
	int switch_init;  /* initial value for dual sign-switching algorithm */
	short count;      /* number of station in interpolation list */
	char noz;         /* set to 1 for output with no elevation corrections */
} daymet_tair_struct;
