/*******************************************************************************
NAME                           FOR_INIT 

PURPOSE:	Initializes forward projection transformation parameters

PROGRAMMER                DATE		REASON
----------                ----		------
T. Mittan		  3-09-93	Initial Development
S. Nelson		  11-94		Added Clarke spheroid default to UTM
Raj Gejjagaraguppe(ARC)   8-30-96       Landsat Ratio is removed as hard
                                        coded value.  Now this ratio can be
                                        an input from the user through the
                                        projection parameter array element
                                        number 9.
Raj Gejjagaraguppe(ARC)   1-07-97       Added a new projection type called
                                        Integerized Sinusoidal Grid to 
                                        support MODIS level 3 datasets.
D. Wynne(ARC)		  3-24-97	Added Support for Power Challenge
					(R10000 Processor Chip Revision: 2.5)
					Long is 8 bytes, on all other currently
					supported platforms Long is 4 bytes.

ALGORITHM REFERENCES

1.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
    Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
    State Government Printing Office, Washington D.C., 1987.

2.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
    U.S. Geological Survey Professional Paper 1453 , United State Government
    Printing Office, Washington D.C., 1989.
*******************************************************************************/
#include "cproj.h"
#include "proj.h"


for_init(outsys,outzone,outparm,outdatum,fn27,fn83,iflg,for_trans)

#if defined(SGI64)
int outsys;		/* output system code				*/
int outzone;		/* output zone number				*/
double *outparm;	/* output array of projection parameters	*/
int outdatum;		/* output datum					*/
char *fn27;		/* NAD 1927 parameter file			*/
char *fn83;		/* NAD 1983 parameter file			*/
int *iflg;		/* status flag					*/
int (*for_trans[])();	/* forward function pointer			*/
#else
long outsys;		/* output system code				*/
long outzone;		/* output zone number				*/
double *outparm;	/* output array of projection parameters	*/
long outdatum;		/* output datum					*/
char *fn27;		/* NAD 1927 parameter file			*/
char *fn83;		/* NAD 1983 parameter file			*/
long *iflg;		/* status flag					*/
long (*for_trans[])();	/* forward function pointer			*/
#endif
{
long zone;		/* zone number					*/
double azimuth;		/* azimuth					*/
double alf;		/* SOM angle					*/
double angle;		/* rotation anlge				*/
double lon1;		/* longitude point in utm scene			*/
double lon2;		/* 2nd longitude 				*/
double lat1;		/* 1st standard parallel			*/
double lat2;		/* 2nd standard parallel			*/
double center_long;	/* center longitude				*/
double center_lat;	/* center latitude				*/
double h;		/* height above sphere				*/
double lon_origin;	/* longitude at origin				*/
double lat_origin;	/* latitude at origin				*/
double r_major;		/* major axis in meters				*/
double r_minor;		/* minor axis in meters				*/
double scale_factor;	/* scale factor					*/
double false_easting;	/* false easting in meters			*/
double false_northing;	/* false northing in meters			*/
double shape_m;		/* constant used for Oblated Equal Area		*/
double shape_n;		/* constant used for Oblated Equal Area		*/
long   start;		/* where SOM starts beginning or end		*/
double time;		/* SOM time					*/
double radius;		/* radius of sphere				*/
double paksz();		/* function to convert DMS to DEG		*/
long tmpdatum;		/* temporary datum for UTM			*/
long path;		/* SOM path number				*/
long satnum;		/* SOM satellite number				*/
long mode;		/* which initialization method  to use A or B	*/
double sat_ratio;       /* satellite ratio which specify the start point*/
double dzone;           /* number of longitudinal zones in ISG          */
double djustify;        /* justify flag in ISG projection               */

long thing;		/* used to initialize 8 byte pointer, added	*/
			/* for Power Challenge		 		*/
long *iflg64;		/* 8 byte status flag, for Power Challenge	*/

long (*for_trans64[MAXPROJ + 1])();
			/* 8 byte forward function pointer for Power 	*/
			/* Challenge */



 
	/* Function declarations for function pointer use
	-----------------------------------------------*/
long utmfor();
long stplnfor();
long alberfor();
long lamccfor();
long merfor();
long psfor();
long polyfor();
long eqconfor();
long tmfor();
long sterfor();
long lamazfor();
long azimfor();
long gnomfor();
long orthfor();
long gvnspfor();
long sinfor();
long equifor();
long millfor();
long vandgfor();
long omerfor();
long somfor();
long hamfor();
long robfor();
long goodfor();
long molwfor();
long imolwfor();
long alconfor();
long wivfor();
long wviifor();
long obleqfor();
long isinusfor();

thing = 0;			/* These lines are to initialize the 	*/
iflg64 = &thing;		/* the 8-byte pointer address           */

/* Initialize forward transformations
-----------------------------------*/
  /* find the correct major and minor axis
  --------------------------------------*/
  sphdz(outdatum,outparm,&r_major,&r_minor,&radius);
  false_easting  = outparm[6];
  false_northing = outparm[7];

  if (outsys == UTM)
    {
    /* this is the call to initialize U T M
    -------------------------------------*/
     /* set Clarke 1866 spheroid if negative datum code  
        ----------------------------------------------*/
     if (outdatum < 0)
        {
        tmpdatum = 0;
	sphdz(tmpdatum,outparm,&r_major,&r_minor,&radius);
        }
    zone = outzone;
    if (zone == 0)
      {
      lon1 = paksz(outparm[0],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
      if (*iflg64 != 0)
        return ERROR;
      lat1 = paksz(outparm[1],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
      if (*iflg64 != 0)
        return ERROR;
      zone = calc_utm_zone(lon1 * R2D);
      if (lat1 < 0)
         zone = -zone;
      }
    scale_factor = .9996;
    *iflg64 = utmforint(r_major,r_minor,scale_factor,zone);
    *iflg = *iflg64;
    for_trans64[outsys] = utmfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == SPCS)
    {
    /* this is the call to initialize STATE PLANE 
    -------------------------------------------*/
    *iflg64 = stplnforint(outzone,outdatum,fn27,fn83);
    *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    for_trans64[outsys] = stplnfor;
        #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == ALBERS)
    {
    /* this is the call to initialize ALBERS CONICAL EQUAL AREA 
    ----------------------------------------------------------*/
    lat1 = paksz(outparm[2],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat2 = paksz(outparm[3],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat_origin = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = alberforint(r_major,r_minor,lat1,lat2,center_long,lat_origin,
		       false_easting, false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = alberfor;
        #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == LAMCC)
    {
    /* this is the call to initialize LAMBERT CONFORMAL CONIC 
    --------------------------------------------------------*/
    lat1 = paksz(outparm[2],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat2 = paksz(outparm[3],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat_origin  = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = lamccforint(r_major,r_minor,lat1,lat2,center_long,lat_origin,
	   	       false_easting, false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = lamccfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == MERCAT)
    {
    /* this is the call to initialize MERCATOR
    ----------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat1   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = merforint(r_major,r_minor,center_long,lat1,false_easting,
		     false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = merfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == PS)
    {
    /* this is the call to initialize POLAR STEREOGRAPHIC 
    ----------------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
       *iflg = (int) *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat1  = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;

    *iflg64 = psforint(r_major,r_minor,center_long,lat1,false_easting,
		    false_northing);
      *iflg = *iflg64;

    for_trans64[outsys] = psfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == POLYC)
    {
    /* this is the call to initialize POLYCONIC
    -----------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat_origin   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = polyforint(r_major,r_minor,center_long,lat_origin,false_easting,
		      false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = polyfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == EQUIDC)
    {
    /* this is the call to initialize EQUIDISTANT CONIC 
    -------------------------------------------------*/
    lat1 = paksz(outparm[2],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat2 = paksz(outparm[3],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat_origin   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    if (outparm[8] == 0)
       mode = 0;
    else 
       mode = 1;
    *iflg64 = eqconforint(r_major,r_minor,lat1,lat2,center_long,lat_origin,
		false_easting,false_northing,mode);
      *iflg = *iflg64;
    for_trans64[outsys] = eqconfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == TM)
    {
    /* this is the call to initialize TRANSVERSE MECTAR
    -------------------------------------------------*/
    scale_factor = outparm[2];
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat_origin   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = tmforint(r_major,r_minor,scale_factor,center_long,lat_origin,
		    false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = tmfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == STEREO)
    {
    /* this is the call to initialize STEREOGRAPHIC
    ---------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = sterforint(radius,center_long,center_lat,false_easting,
		      false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = sterfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == LAMAZ)
    {
    /* this is the call to initialize LAMBERT AZIMUTHAL
    -------------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat  = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = lamazforint(radius,center_long, center_lat,false_easting,
		       false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = lamazfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == AZMEQD)
    {
    /* this is the call to initialize AZIMUTHAL EQUIDISTANT
    -----------------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = azimforint(radius,center_long,center_lat,false_easting,
		      false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = azimfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == GNOMON)
    {
    /* this is the call to initialize GNOMONIC 
    ----------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = gnomforint(radius,center_long,center_lat,false_easting,
		      false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = gnomfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == ORTHO)
    {
    /* this is the call to initalize ORTHOGRAPHIC
    -------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = orthforint(radius,center_long,center_lat,false_easting,
		      false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = orthfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == GVNSP)
    {
    /* this is the call to initalize GENERAL VERTICAL NEAR-SIDE PERSPECTIVE
    ----------------------------------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    h = outparm[2];
    *iflg64 = gvnspforint(radius,h,center_long,center_lat,false_easting,
		       false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = gvnspfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == SNSOID)
    {
    /* this is the call to initialize SINUSOIDAL 
    -------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = sinforint(radius, center_long,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = sinfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == EQRECT)
    {
    /* this is the call to initialize EQUIRECTANGULAR
    -----------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    lat1   = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = equiforint(radius,center_long,lat1,false_easting,false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = equifor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == MILLER)
    {
    /* this is the call to initialize MILLER CYLINDRICAL 
    --------------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64) * 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = millforint(radius, center_long,false_easting,false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = millfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == VGRINT)
    {
    /* this is the call to initialize VAN DER GRINTEN 
    -----------------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = vandgforint(radius, center_long,false_easting,false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = vandgfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == HOM)
     {
     /* this is the call to initialize HOTLINE OBLIQUE MERCATOR
     ---------------------------------------------------------*/
     scale_factor = outparm[2];
     lat_origin = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
     if (*iflg64 != 0)
        return ERROR;
     if (outparm[12] != 0)
        {
        mode = 1;
        azimuth = paksz(outparm[3],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        lon_origin = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        }
     else
        {
        mode = 0;
        lon1 = paksz(outparm[8],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        lat1 = paksz(outparm[9],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        lon2 = paksz(outparm[10],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        lat2 = paksz(outparm[11],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
        if (*iflg64 != 0)
           return ERROR;
        }
     *iflg64 = omerforint(r_major,r_minor,scale_factor,azimuth,lon_origin,
                        lat_origin,false_easting, false_northing,lon1,lat1,
                        lon2,lat2,mode);
      *iflg = *iflg64;
     for_trans64[outsys] = omerfor;
     #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
     }
  else
  if (outsys == SOM)
    {
    /* this is the call to initialize SOM 
    -----------------------------------*/
    path = outparm[3];
    satnum = outparm[2];
    if (outparm[12] == 0)
       {
       mode = 1;
       alf = paksz(outparm[3],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
       if (*iflg64 != 0)
          return ERROR;
       lon1 = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
       if (*iflg64 != 0)
          return ERROR;
       time = outparm[8];
       sat_ratio = outparm[9];
       start = outparm[10];
       }
    else
       mode = 0;
/*
    *iflg64 = somforint(r_major,r_minor,satnum,path,false_easting,false_northing);
      *iflg = *iflg64;
*/
    *iflg64 = somforint(r_major,r_minor,satnum,path,alf,lon1,false_easting,
		      false_northing,time,start,mode,sat_ratio);
      *iflg = *iflg64;
    for_trans64[outsys] = somfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == HAMMER)
    {
    /* this is the call to initialize HAMMER 
    --------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = hamforint(radius,center_long,false_easting,false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = hamfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == ROBIN)
    {
    /* this is the call to initialize ROBINSON 
    ----------------------------------------*/
    center_long  = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = robforint(radius,center_long,false_easting,false_northing); 
      *iflg = *iflg64;
    for_trans64[outsys] = robfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == GOODE)
    {
    /* this is the call to initialize GOODE'S HOMOLOSINE
    ---------------------------------------------------*/
    *iflg64 = goodforint(radius);
      *iflg = *iflg64;
    for_trans64[outsys] = goodfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == MOLL)
    {
    /* this is the call to initialize MOLLWEIDE
    ------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = molwforint(radius, center_long,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = molwfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == IMOLL)
    {
    /* this is the call to initialize INTERRUPTED MOLLWEIDE
    -----------------------------------------------------*/
    *iflg64 = imolwforint(radius);
      *iflg = *iflg64;
    for_trans64[outsys] = imolwfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == ALASKA)
    {
    /* this is the call to initialize ALASKA CONFORMAL 
    ------------------------------------------------*/
    *iflg64 = alconforint(r_major,r_minor,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = alconfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == WAGIV)
    {
    /* this is the call to initialize WAGNER IV 
    -----------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = wivforint(radius, center_long,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = wivfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == WAGVII)
    {
    /* this is the call to initialize WAGNER VII 
    ------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = wviiforint(radius, center_long,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = wviifor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == OBEQA)
    {
    /* this is the call to initialize OBLATED EQUAL AREA 
    ---------------------------------------------------*/
    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    center_lat  = paksz(outparm[5],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    shape_m = outparm[2];
    shape_n = outparm[3];
    angle = paksz(outparm[8],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    *iflg64 = obleqforint(radius,center_long,center_lat,shape_m, shape_n, 
		angle,false_easting,false_northing);
      *iflg = *iflg64;
    for_trans64[outsys] = obleqfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
  else
  if (outsys == ISINUS)
    {
    /* this is the call to initialize INTEGERIZED SINUSOIDAL GRID
    ------------------------------------------------------------*/

    center_long = paksz(outparm[4],iflg64)* 3600 * S2R;
      *iflg = *iflg64;
    if (*iflg64 != 0)
       return ERROR;
    dzone = outparm[8];
    djustify = outparm[10];

    *iflg64 = isinusforinit(radius, center_long, false_easting, false_northing,
                    dzone, djustify);
      *iflg = *iflg64;

    for_trans64[outsys] = isinusfor;
    #if defined(SGI64)
    for_trans[outsys] = (int (*)()) for_trans64[outsys];
    #else
    for_trans[outsys] = for_trans64[outsys];
    #endif
    }
       
return OK;
}
