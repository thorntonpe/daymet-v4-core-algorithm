/*
predict_daily_srad_vp_nc.c
Peter Thornton
2/13/01

Based on The MTCLIM v4.3 code for simultaneous estimation of radiation
and humidity.  A part of the Daymet code, relying on inputs from the
predict_daily_tair and predict_daily_prcp routines.

Revision notes: 
2/12/01, PET: Adding support for east horizon, west horizon, and average
horizon input files. For now this means that the LIST inputs will not
be consistent with the MAP inputs, so remember not to use this version 
with the LIST option. To make it easier to update the LIST option later
if necessary, I'm leaving the LIST handling as is.
The three input horizon angle grids are assumed to be stored as IEEE
floats (32-bit), in units of radians, measuring the positive angle above
the horizon. For consistency with the rest of the Daymet code, these 
grids are required to be registered, using the register_image routine.
The projection information for these grids must match that for the DEM,
mask, slope and aspect files, as well as the input daily metdata files.

2/15/01, PET: Adding flag that makes it possible to avoid the snow correction,
since there is some question whether the parameters from the Austrian
study are applicable for the U.S. Also adding flag to turn off output of
daylength, since it doesn't vary between years.
Also modified wet-day correction, so that the daily total prcp has to 
be greater than 2 mm in order to initiate radiation reduction.

July 2005 : netcdf grids - CSH

1/12/10, PET: fixed error in typecast of slope and aspect grids,
  introduced during conversion to netcdf grids.

Usage:
For now the only supported option is "MAP" input, which takes a command
line with the following structure
arg# description
0		executable name
1		"MAP" (string)
2		good mask value (integer) - the mask value of the simulation region
3		mask filename (string)
4		dem filename (string)
5		slope filename (string)
6		aspect filename (string)
7		east horizon filename (string)
8		west horizon filename (string)
9		average horizon filename (string)
10		tmax filename (string)
11		tmin filename (string)
12		prcp filename (string)
13		outprefix (string)
14		dt (float) - sub-daily timestep, in seconds
15		daylength output flag (integer) - 0=no, 1=yes
16		snowpack correction (and output) flag (integer) - 0=no, 1=yes

*/

#include "metsrc2.h"
#include "cproj.h"
#include "proj.h"
#include "netcdf.h"

#define SECPERRAD 13750.9871            /* seconds per radian of hour angle */
#define RADPERDAY 0.017214              /* radians of orbit per julian day */
#define RADPERDEG 0.01745329            /* radians per degree */
#define TRADPERDEG 17.45329             /* thousandths of radians per degree */
#define MINDECL -0.4092797              /* minimum declination (radians) */
#define DAYSOFF 11.25                   /* julian day offset winter solstice */

/* parameters for the Tair algorithm */
#define TDAYCOEF     0.45  /* (dim) daylight air temperature coefficient (dim) */

/* parameters for the snowpack algorithm */
#define SNOW_TCRIT   -6.0  /* (deg C) critical temperature for snowmelt */   
#define SNOW_TRATE  0.042  /* (cm/degC/day) snowmelt rate */

/* parameters for the radiation algorithm */
#define TBASE       0.870  /* (dim) max inst. trans., 0m, nadir, dry atm */
#define ABASE     -6.1e-5  /* (1/Pa) vapor pressure effect on transmittance */
#define C             1.5  /* (dim) radiation parameter */
#define B0          0.013  /* (dim) radiation parameter */
#define B1          0.201  /* (dim) radiation parameter */
#define B2          0.185  /* (dim) radiation parameter */
#define RAIN_SCALAR  0.75  /* (dim) correction to trans. for rain day */
#define DIF_ALB       0.6  /* (dim) diffuse albedo for horizon correction */
#define SC_INT       1.32  /* (MJ/m2/day) snow correction intercept */
#define SC_SLOPE    0.096  /* (MJ/m2/day/cm) snow correction slope */

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define R        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
#define EPS      0.62196351      /* (MW/MA) unitless ratio of molec weights */
#define TWOPI    (2.0*PI)        /* 2*pi */
#define HALFPI   (PI/2.0)        /* pi/2 */

#define SRAD_VALMAX 800.0
#define SRAD_VALMIN 0.0
#define SRAD_BYTMAX 250
#define SRAD_BYTMIN 0
#define DAYL_VALMAX 86400.0
#define DAYL_VALMIN 0.0
#define DAYL_BYTMAX 250
#define DAYL_BYTMIN 0
#define VP_VALMAX 10000.0
#define VP_VALMIN 0.0
#define VP_BYTMAX 250
#define VP_BYTMIN 0
#define SWE_VALMAX 100.0
#define SWE_VALMIN 0.0
#define SWE_BYTMAX 250
#define SWE_BYTMIN 0

/* function prototypes */
double calc_pet(double rad, double ta, double pa, double dayl);

int main(int argc, char* argv[])
{
	file ini;
	file mask_f;
	file dem_f;
	file asp_f;
	file slp_f;
	file ehoriz_f, whoriz_f, avghoriz_f;
	file list_f;
	file tmax_f;
	file tmin_f;
	file prcp_f;
	file dayl_out_f, srad_out_f, vp_out_f, swe_out_f;
	file dayl_tavg_f, dayl_savg_f;
	file srad_tavg_f, srad_savg_f;
	file vp_tavg_f, vp_savg_f;
	file swe_tavg_f, swe_savg_f;
	
	file mask_nc_f;
	file dem_nc_f;
	file asp_nc_f;
	file slp_nc_f;
	file ehoriz_nc_f, whoriz_nc_f, avghoriz_nc_f;
	imghdr_struct mask_nc_hdr, dem_nc_hdr, asp_nc_hdr, slp_nc_hdr, ehoriz_nc_hdr, whoriz_nc_hdr,avghoriz_nc_hdr;
	long int *mask_array;
	size_t lenp;
	size_t index[2];
	int *dem_array;
	float *slp_array, *asp_array;
	float slp, asp;
	float *ehoriz_array, *whoriz_array, *avghoriz_array;
	int ncid_mask, colsid_mask, rowsid_mask, imageid_mask;
	int ncid_dem, colsid_dem, rowsid_dem, imageid_dem;
	int ncid_slp, colsid_slp, rowsid_slp, imageid_slp;
	int ncid_asp, colsid_asp, rowsid_asp, imageid_asp;
	int ncid_ehoriz, colsid_ehoriz, rowsid_ehoriz, imageid_ehoriz;
	int ncid_whoriz, colsid_whoriz, rowsid_whoriz, imageid_whoriz;
	int ncid_avghoriz, colsid_avghoriz, rowsid_avghoriz, imageid_avghoriz;
	int status;

	dayouthdr_struct tmaxhdr, tminhdr, prcphdr, sradhdr, daylhdr, vphdr, swehdr;
	imghdr_struct maskhdr, demhdr, asphdr, slphdr, ehorizhdr, whorizhdr,avghorizhdr;
	lsthdr_struct listhdr;
	gctp_struct gctp;
	
	gctp_struct gstr;
	long iflg;
	long (*inv_trans[MAXPROJ + 1])();
	
	double celldat[5];

	double sin_decl[365];                /* coded on yearday */
	double cos_decl[365];                /* coded on yearday */
	double sin_0to2pi[6284];             /* coded on thousandths of a radian */
	double cos_0to2pi[6284];             /* coded on thousandths of a radian */
	double sin_negpitopi[6284];          /* coded on thousandths of a radian */
	double cos_negpitopi[6284];          /* coded on thousandths of a radian */
	double trans;
	double dir_beam_topa,dir_beam_bota,dir_flat_topa,dir_flat_bota;
	double dir_slope,dif_slope,total_slope;
	double sum_dir_slope, sum_dif_slope;
	double max_rate;
	double t_dif,dif_flat_bota;
	double cbsa;
	double cosegeom,sinegeom;
	double cosdecl,sindecl;
	double coslat,sinlat;
	double cosh,sinh;
	double cosslope,sinslope,coshalfslope;
	double cosasp,sinasp;
	float ehoriz, whoriz, avghoriz;
	double coszeh, coszwh;
	double dt; /* timestep in seconds */
	double dh;
	double dayl;
	double save1,save2,save3;
	double coshss,hss;
	double h,starth,stoph;
	double sc,cza,oam,am; 
	double dtr_max, dtr_min, scale_ratio;
	double cst, pcst, b, c, t, dtr;
	double angle;
	double total,trate,maxr;
	double elev;
	double ttmax0[366];
	double flat_potrad[366];
	double slope_potrad[366];
	double daylength[366];
	double sum_trans, sum_flat_potrad, sum_slope_potrad;

	/* solar constant by month (W/m^2) */
	/* original values, average = 1403, too high by about 35 W/m^2
	double solcon[12] = {1445.0,1431.0,1410.0,1389.0,1368.0,1354.0,1354.0,
	1375.0,1403.0,1424.0,1438.0,1445.0};
	*/
	/* new values subtracting 35 from old values */
	double solcon[12] = {1410.0,1396.0,1375.0,1354.0,1333.0,1319.0,1319.0,
	1340.0,1368.0,1389.0,1403.0,1410.0};

	/* optical airmass by degrees */
	double optam[21] = {2.90,3.05,3.21,3.39,3.69,3.82,4.07,4.37,4.72,5.12,5.60,
	6.18,6.88,7.77,8.90,10.39,12.44,15.36,19.79,26.96,30.00};
	
	double *tmax_ts, *btmin_ts, *bprcp_ts;
	int *yday_ts;
	double *tmin_ts, *prcp_ts, *dtr_ts, *smooth_dtr_ts, *dayl_ts, *srad_ts, *vp_ts;
	double srad_tavg, dayl_tavg, vp_tavg, swe_tavg;
	double *srad_savg, *dayl_savg, *vp_savg, *swe_savg;
	double map_x, map_y, longitude, latitude;
	double slope, aspect;
	double map_ulx, map_uly, cellsize;
	double t1, t2, pratio, trans1,trans2;
	double avg_horizon, horizon_scalar;
	double slope_excess, slope_scalar, sky_prop;
	double tmax, tmin, tmean;
	double avgdtr;
	double pva, t_tmax, t_final, pdif, pdir, srad1, srad2;
	double snowpack, newsnow, snowmelt, sum;
	double *parray, *window, *t_fmax, *tdew, *swe_ts;
	
	int i,j,k,l,day;
	int bdtr_max, bdtr_min;
	int flat_flag,yday;
	int ami,month,hh;
	int ncols, nrows;
	int row,col;
	int ncells;
	int ndays;
	int start_yday, end_yday, prev_yday, count, start_year;
	int slp_ind, asp_ind, lat_ind;
	int coszeh_ind, coszwh_ind;
	int n_in;
	int test;
	int junk;
	int add;
	
	int z;
	
	char outprefix[80], round[16];
	char endkey[8];
	char iocase[16];
	char ioc,map;
	long int mask;
	double tmax_brange, tmin_brange, tmax_bmin, tmin_bmin;
	double prcp_brange, prcp_bmin;
	double srad_brange, dayl_brange, vp_brange, swe_brange;
	double srad_bmin, dayl_bmin, vp_bmin, swe_bmin;
	
	double tmax_scale_ratio, tmin_scale_ratio, tmax_bytrange, tmin_bytrange;
	double prcp_scale_ratio, prcp_bytrange;
	double srad_scale_ratio, dayl_scale_ratio, srad_bytrange, dayl_bytrange;
	double vp_scale_ratio, vp_bytrange;
	double swe_scale_ratio, swe_bytrange;
	double tmax_valmin, tmin_valmin, prcp_valmin;
	double srad_valmin, dayl_valmin, vp_valmin, swe_valmin;
	double sum_prcp, ann_prcp, effann_prcp;
	
	double bsg1, bsg2, bsg3;
	double pa;
	double sum_pet, ann_pet, *save_pet, *tday_ts;
	double tmink, pet, ratio, tdewk;
	double multiplier;
	long int good_mask;
	int dayl_out;
	int isloop;
	int snow_flag;

	/* assuming command-line parameter passing */
	/* copy MAP or LIST to iocase string */
	strcpy(iocase,argv[1]);
	if (!isalpha(ioc=toupper(iocase[0])))
	{
		printf("Error on command line...\n");
		printf("First parameter should be:\n");
		printf("MAP or LIST, to indicate the I/O mode ... exiting\n");
		exit(1);
	}
	if (ioc == 'M')
	{
		map = 1;
		add = 7;
		good_mask = atoi(argv[2]);
		strcpy(mask_nc_f.name,argv[3]);
		strcpy(dem_nc_f.name,argv[4]);
		strcpy(slp_nc_f.name,argv[5]);
		strcpy(asp_nc_f.name,argv[6]);
		strcpy(ehoriz_nc_f.name,argv[7]);
		strcpy(whoriz_nc_f.name,argv[8]);
		strcpy(avghoriz_nc_f.name,argv[9]);
	}
	else if (ioc == 'L')
	{
		map = 0;
		add = 0;
		/* scan next parameter as the list filename */
		strcpy(list_f.name,argv[2]);
		if (file_open(&list_f, 'r'))
		{
			printf("Error opening %s for binary read\n",list_f.name);
			exit(1);
		}
	}
	else
	{
		printf("Error on command line:\n");
		printf("First parameter should be\n");
		printf("MAP or LIST, to indicate the I/O mode ... exiting\n");
		exit(1);
	}
	/* the remaining parameters are the same for both MAP and LIST inputs,
	but the parameter numbers are offset for map input */
	strcpy(tmax_f.name,argv[3+add]);
	strcpy(tmin_f.name,argv[4+add]);
	strcpy(prcp_f.name,argv[5+add]);
	strcpy(outprefix,argv[6+add]);
	dt = atof(argv[7+add]);
	dayl_out = atoi(argv[8+add]);
	snow_flag = atoi(argv[9+add]);
	
/* -----------------------  netcdf calls for mask file --------------------------------------- */
/*                                                                                             */
        strcat(mask_nc_f.name,".nc");
        if( status = nc_open(mask_nc_f.name, 0, &ncid_mask) )
        {
                printf("Error opening %s for netcdf read, exiting\n", mask_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_mask, "ncols", &colsid_mask);
        status = nc_inq_dimlen(ncid_mask, colsid_mask, &lenp);
        mask_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_mask, "nrows", &rowsid_mask);
        status = nc_inq_dimlen(ncid_mask, rowsid_mask, &lenp);
        mask_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_mask, NC_GLOBAL, "tag", mask_nc_hdr.tag);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "ulx", &mask_nc_hdr.ulx);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "uly", &mask_nc_hdr.uly);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "cellsize", &mask_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_mask, NC_GLOBAL, "gctp_outsys", &mask_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_mask, NC_GLOBAL, "gctp_outzone", &mask_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_mask, NC_GLOBAL, "gctp_outdatum", &mask_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "gctp_outparm", mask_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_mask, NC_GLOBAL, "gctp_fn27", mask_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_mask, NC_GLOBAL, "gctp_fn83", mask_nc_hdr.gctp.fn83);
        if (strncmp(mask_nc_hdr.tag,IMGTAG,6))
        {
                printf("mask file not a registered image. Use register_image first.\n");
                printf("mask tag: %s\n",mask_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_mask, "image", &imageid_mask);
        if (!(mask_array = (long int*) malloc(mask_nc_hdr.nrows * mask_nc_hdr.ncols * sizeof(long int)))) exit(1);
        status = nc_get_var_long (ncid_mask, imageid_mask, mask_array);

        nc_close(ncid_mask);
/*                                                                                                 */
/* -----------------------  end netcdf calls for mask file --------------------------------------- */

/* -----------------------  netcdf calls for DEM file --------------------------------------- */
/*                                                                                             */
        strcat(dem_nc_f.name,".nc");
        if( status = nc_open(dem_nc_f.name, 0, &ncid_dem) )
        {
                printf("Error opening %s for netcdf read, exiting\n", dem_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_dem, "ncols", &colsid_dem);
        status = nc_inq_dimlen(ncid_dem, colsid_dem, &lenp);
        dem_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_dem, "nrows", &rowsid_dem);
        status = nc_inq_dimlen(ncid_dem, rowsid_dem, &lenp);
        dem_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_dem, NC_GLOBAL, "tag", dem_nc_hdr.tag);
        status = nc_get_att_double (ncid_dem, NC_GLOBAL, "ulx", &dem_nc_hdr.ulx);
        status = nc_get_att_double (ncid_dem, NC_GLOBAL, "uly", &dem_nc_hdr.uly);
        status = nc_get_att_double (ncid_dem, NC_GLOBAL, "cellsize", &dem_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_dem, NC_GLOBAL, "gctp_outsys", &dem_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_dem, NC_GLOBAL, "gctp_outzone", &dem_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_dem, NC_GLOBAL, "gctp_outdatum", &dem_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_dem, NC_GLOBAL, "gctp_outparm", dem_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_dem, NC_GLOBAL, "gctp_fn27", dem_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_dem, NC_GLOBAL, "gctp_fn83", dem_nc_hdr.gctp.fn83);
        if (strncmp(dem_nc_hdr.tag,IMGTAG,6))
        {
                printf("dem file not a registered image. Use register_image first.\n");
                printf("dem tag: %s\n",dem_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_dem, "image", &imageid_dem);
        /* if (!(dem_array = (int*) malloc(dem_nc_hdr.nrows * dem_nc_hdr.ncols * sizeof(int)))) exit(1); */
        /* status = nc_get_var_int (ncid_dem, imageid_dem, dem_array); */

        /* nc_close(ncid_dem); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for DEM file --------------------------------------- */

/* -----------------------  netcdf calls for slope file --------------------------------------- */
/*                                                                                             */
        strcat(slp_nc_f.name,".nc");
        if( status = nc_open(slp_nc_f.name, 0, &ncid_slp) )
        {
                printf("Error opening %s for netcdf read, exiting\n", slp_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_slp, "ncols", &colsid_slp);
        status = nc_inq_dimlen(ncid_slp, colsid_slp, &lenp);
        slp_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_slp, "nrows", &rowsid_slp);
        status = nc_inq_dimlen(ncid_slp, rowsid_slp, &lenp);
        slp_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_slp, NC_GLOBAL, "tag", slp_nc_hdr.tag);
        status = nc_get_att_double (ncid_slp, NC_GLOBAL, "ulx", &slp_nc_hdr.ulx);
        status = nc_get_att_double (ncid_slp, NC_GLOBAL, "uly", &slp_nc_hdr.uly);
        status = nc_get_att_double (ncid_slp, NC_GLOBAL, "cellsize", &slp_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_slp, NC_GLOBAL, "gctp_outsys", &slp_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_slp, NC_GLOBAL, "gctp_outzone", &slp_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_slp, NC_GLOBAL, "gctp_outdatum", &slp_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_slp, NC_GLOBAL, "gctp_outparm", slp_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_slp, NC_GLOBAL, "gctp_fn27", slp_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_slp, NC_GLOBAL, "gctp_fn83", slp_nc_hdr.gctp.fn83);
        if (strncmp(slp_nc_hdr.tag,IMGTAG,6))
        {
                printf("slp file not a registered image. Use register_image first.\n");
                printf("slp tag: %s\n",slp_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_slp, "image", &imageid_slp);
        /* if (!(slp_array = (float*) malloc(slp_nc_hdr.nrows * slp_nc_hdr.ncols * sizeof(float)))) exit(1); */
        /* status = nc_get_var_float (ncid_slp, imageid_slp, slp_array); */

        /* nc_close(ncid_slp); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for slope file --------------------------------------- */

/* -----------------------  netcdf calls for aspect file --------------------------------------- */
/*                                                                                             */
        strcat(asp_nc_f.name,".nc");
        if( status = nc_open(asp_nc_f.name, 0, &ncid_asp) )
        {
                printf("Error opening %s for netcdf read, exiting\n", asp_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_asp, "ncols", &colsid_asp);
        status = nc_inq_dimlen(ncid_asp, colsid_asp, &lenp);
        asp_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_asp, "nrows", &rowsid_asp);
        status = nc_inq_dimlen(ncid_asp, rowsid_asp, &lenp);
        asp_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_asp, NC_GLOBAL, "tag", asp_nc_hdr.tag);
        status = nc_get_att_double (ncid_asp, NC_GLOBAL, "ulx", &asp_nc_hdr.ulx);
        status = nc_get_att_double (ncid_asp, NC_GLOBAL, "uly", &asp_nc_hdr.uly);
        status = nc_get_att_double (ncid_asp, NC_GLOBAL, "cellsize", &asp_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_asp, NC_GLOBAL, "gctp_outsys", &asp_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_asp, NC_GLOBAL, "gctp_outzone", &asp_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_asp, NC_GLOBAL, "gctp_outdatum", &asp_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_asp, NC_GLOBAL, "gctp_outparm", asp_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_asp, NC_GLOBAL, "gctp_fn27", asp_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_asp, NC_GLOBAL, "gctp_fn83", asp_nc_hdr.gctp.fn83);
        if (strncmp(asp_nc_hdr.tag,IMGTAG,6))
        {
                printf("asp file not a registered image. Use register_image first.\n");
                printf("asp tag: %s\n",asp_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_asp, "image", &imageid_asp);
        /* if (!(asp_array = (float*) malloc(asp_nc_hdr.nrows * asp_nc_hdr.ncols * sizeof(float)))) exit(1); */
        /* status = nc_get_var_float (ncid_asp, imageid_asp, asp_array); */

        /* nc_close(ncid_asp); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for aspect file --------------------------------------- */

/* -----------------------  netcdf calls for ehoriz file --------------------------------------- */
/*                                                                                             */
        strcat(ehoriz_nc_f.name,".nc");
        if( status = nc_open(ehoriz_nc_f.name, 0, &ncid_ehoriz) )
        {
                printf("Error opening %s for netcdf read, exiting\n", ehoriz_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_ehoriz, "ncols", &colsid_ehoriz);
        status = nc_inq_dimlen(ncid_ehoriz, colsid_ehoriz, &lenp);
        ehoriz_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_ehoriz, "nrows", &rowsid_ehoriz);
        status = nc_inq_dimlen(ncid_ehoriz, rowsid_ehoriz, &lenp);
        ehoriz_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_ehoriz, NC_GLOBAL, "tag", ehoriz_nc_hdr.tag);
        status = nc_get_att_double (ncid_ehoriz, NC_GLOBAL, "ulx", &ehoriz_nc_hdr.ulx);
        status = nc_get_att_double (ncid_ehoriz, NC_GLOBAL, "uly", &ehoriz_nc_hdr.uly);
        status = nc_get_att_double (ncid_ehoriz, NC_GLOBAL, "cellsize", &ehoriz_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_ehoriz, NC_GLOBAL, "gctp_outsys", &ehoriz_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_ehoriz, NC_GLOBAL, "gctp_outzone", &ehoriz_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_ehoriz, NC_GLOBAL, "gctp_outdatum", &ehoriz_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_ehoriz, NC_GLOBAL, "gctp_outparm", ehoriz_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_ehoriz, NC_GLOBAL, "gctp_fn27", ehoriz_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_ehoriz, NC_GLOBAL, "gctp_fn83", ehoriz_nc_hdr.gctp.fn83);
        if (strncmp(ehoriz_nc_hdr.tag,IMGTAG,6))
        {
                printf("ehoriz file not a registered image. Use register_image first.\n");
                printf("ehoriz tag: %s\n",ehoriz_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_ehoriz, "image", &imageid_ehoriz);
        /* if (!(ehoriz_array = (float*) malloc(ehoriz_nc_hdr.nrows * ehoriz_nc_hdr.ncols * sizeof(float)))) exit(1); */
        /* status = nc_get_var_float (ncid_ehoriz, imageid_ehoriz, ehoriz_array); */

        /* nc_close(ncid_ehoriz); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for ehoriz file --------------------------------------- */

/* -----------------------  netcdf calls for whoriz file --------------------------------------- */
/*                                                                                             */
        strcat(whoriz_nc_f.name,".nc");
        if( status = nc_open(whoriz_nc_f.name, 0, &ncid_whoriz) )
        {
                printf("Error opening %s for netcdf read, exiting\n", whoriz_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_whoriz, "ncols", &colsid_whoriz);
        status = nc_inq_dimlen(ncid_whoriz, colsid_whoriz, &lenp);
        whoriz_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_whoriz, "nrows", &rowsid_whoriz);
        status = nc_inq_dimlen(ncid_whoriz, rowsid_whoriz, &lenp);
        whoriz_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_whoriz, NC_GLOBAL, "tag", whoriz_nc_hdr.tag);
        status = nc_get_att_double (ncid_whoriz, NC_GLOBAL, "ulx", &whoriz_nc_hdr.ulx);
        status = nc_get_att_double (ncid_whoriz, NC_GLOBAL, "uly", &whoriz_nc_hdr.uly);
        status = nc_get_att_double (ncid_whoriz, NC_GLOBAL, "cellsize", &whoriz_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_whoriz, NC_GLOBAL, "gctp_outsys", &whoriz_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_whoriz, NC_GLOBAL, "gctp_outzone", &whoriz_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_whoriz, NC_GLOBAL, "gctp_outdatum", &whoriz_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_whoriz, NC_GLOBAL, "gctp_outparm", whoriz_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_whoriz, NC_GLOBAL, "gctp_fn27", whoriz_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_whoriz, NC_GLOBAL, "gctp_fn83", whoriz_nc_hdr.gctp.fn83);
        if (strncmp(whoriz_nc_hdr.tag,IMGTAG,6))
        {
                printf("whoriz file not a registered image. Use register_image first.\n");
                printf("whoriz tag: %s\n",whoriz_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_whoriz, "image", &imageid_whoriz);
        /* if (!(whoriz_array = (float*) malloc(whoriz_nc_hdr.nrows * whoriz_nc_hdr.ncols * sizeof(float)))) exit(1); */
        /* status = nc_get_var_float (ncid_whoriz, imageid_whoriz, whoriz_array); */

        /* nc_close(ncid_whoriz); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for whoriz file --------------------------------------- */

/* -----------------------  netcdf calls for avghoriz file --------------------------------------- */
/*                                                                                             */
        strcat(avghoriz_nc_f.name,".nc");
        if( status = nc_open(avghoriz_nc_f.name, 0, &ncid_avghoriz) )
        {
                printf("Error opening %s for netcdf read, exiting\n", avghoriz_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_avghoriz, "ncols", &colsid_avghoriz);
        status = nc_inq_dimlen(ncid_avghoriz, colsid_avghoriz, &lenp);
        avghoriz_nc_hdr.ncols = lenp;
        status = nc_inq_dimid (ncid_avghoriz, "nrows", &rowsid_avghoriz);
        status = nc_inq_dimlen(ncid_avghoriz, rowsid_avghoriz, &lenp);
        avghoriz_nc_hdr.nrows = lenp;
        status = nc_get_att_text (ncid_avghoriz, NC_GLOBAL, "tag", avghoriz_nc_hdr.tag);
        status = nc_get_att_double (ncid_avghoriz, NC_GLOBAL, "ulx", &avghoriz_nc_hdr.ulx);
        status = nc_get_att_double (ncid_avghoriz, NC_GLOBAL, "uly", &avghoriz_nc_hdr.uly);
        status = nc_get_att_double (ncid_avghoriz, NC_GLOBAL, "cellsize", &avghoriz_nc_hdr.cellsize);
        status = nc_get_att_long (ncid_avghoriz, NC_GLOBAL, "gctp_outsys", &avghoriz_nc_hdr.gctp.outsys);
        status = nc_get_att_long (ncid_avghoriz, NC_GLOBAL, "gctp_outzone", &avghoriz_nc_hdr.gctp.outzone);
        status = nc_get_att_long (ncid_avghoriz, NC_GLOBAL, "gctp_outdatum", &avghoriz_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_avghoriz, NC_GLOBAL, "gctp_outparm", avghoriz_nc_hdr.gctp.outparm);
        status = nc_get_att_text (ncid_avghoriz, NC_GLOBAL, "gctp_fn27", avghoriz_nc_hdr.gctp.fn27);
        status = nc_get_att_text (ncid_avghoriz, NC_GLOBAL, "gctp_fn83", avghoriz_nc_hdr.gctp.fn83);
        if (strncmp(avghoriz_nc_hdr.tag,IMGTAG,6))
        {
                printf("avghoriz file not a registered image. Use register_image first.\n");
                printf("avghoriz tag: %s\n",avghoriz_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_avghoriz, "image", &imageid_avghoriz);
        /* if (!(avghoriz_array = (float*) malloc(avghoriz_nc_hdr.nrows * avghoriz_nc_hdr.ncols * sizeof(float)))) exit(1); */
        /* status = nc_get_var_float (ncid_avghoriz, imageid_avghoriz, avghoriz_array); */

        /* nc_close(ncid_avghoriz); */
/*                                                                                                 */
/* -----------------------  end netcdf calls for avghoriz file --------------------------------------- */

	/* open the input data files */
	printf("Tmax Binary: %s", tmax_f.name);
	if (file_open(&tmax_f, 'r'))
	{
		printf("Error opening %s for binary read\n",tmax_f.name);
		exit(1);
	}
	printf("Tmin Binary: %s", tmin_f.name);
	if (file_open(&tmin_f, 'r'))
	{
		printf("Error opening %s for binary read\n",tmin_f.name);
		exit(1);
	}
	printf("Prcp Binary: %s", prcp_f.name);
	if (file_open(&prcp_f, 'r'))
	{
		printf("Error opening %s for binary read\n",prcp_f.name);
		exit(1);
	}
		
	/* error checking on file headers */
	if (map) /* MAP input */
	{
		/* read mask file header and verify */
		if (strncmp(mask_nc_hdr.tag,IMGTAG,6))
		{
			printf("Mask not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* read dem file header and verify */
		if (strncmp(dem_nc_hdr.tag,IMGTAG,6))
		{
			printf("DEM not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* read slp file header and verify */
		if (strncmp(slp_nc_hdr.tag,IMGTAG,6))
		{
			printf("Slope not a registered image. Use register_image first.\n");
			exit(1);
		}

		
		/* read aspect file header and verify */
		if (strncmp(asp_nc_hdr.tag,IMGTAG,6))
		{
			printf("Aspect not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* read east horizon file header and verify */
		if (strncmp(ehoriz_nc_hdr.tag,IMGTAG,6))
		{
			printf("East-horizon not a registered image. Use register_image first.\n");
			exit(1);
		}

		
		/* read west horizon file header and verify */
		if (strncmp(whoriz_nc_hdr.tag,IMGTAG,6))
		{
			printf("west-horizon not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* read average horizon file header and verify */
		if (strncmp(avghoriz_nc_hdr.tag,IMGTAG,6))
		{
			printf("average-horizon not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* compare headers for DEM and mask */
		test = (mask_nc_hdr.ncols == dem_nc_hdr.ncols)*(mask_nc_hdr.nrows == dem_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == dem_nc_hdr.ulx)*(mask_nc_hdr.uly == dem_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == dem_nc_hdr.cellsize);
		if (!test) 
		{
			printf("Header information differs between mask and DEM...exiting\n");
			exit(1);
		}
		
		/* compare headers for slope and mask */
		test = (mask_nc_hdr.ncols == slp_nc_hdr.ncols)*(mask_nc_hdr.nrows == slp_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == slp_nc_hdr.ulx)*(mask_nc_hdr.uly == slp_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == slp_nc_hdr.cellsize);
		if (!test)
		{
			printf("Header information differs between mask and Slope...exiting\n");
			exit(1);
		}
		
		/* compare headers for aspect and mask */
		test = (mask_nc_hdr.ncols == asp_nc_hdr.ncols)*(mask_nc_hdr.nrows == asp_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == asp_nc_hdr.ulx)*(mask_nc_hdr.uly == asp_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == asp_nc_hdr.cellsize);
		if (!test)
		{
			printf("Header information differs between mask and Aspect...exiting\n");
			exit(1);
		}
	
		/* compare headers for east horizon and mask */
		test = (mask_nc_hdr.ncols == ehoriz_nc_hdr.ncols)*(mask_nc_hdr.nrows == ehoriz_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == ehoriz_nc_hdr.ulx)*(mask_nc_hdr.uly == ehoriz_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == ehoriz_nc_hdr.cellsize);
		if (!test)
		{
			printf("Header information differs between mask and ehoriz...exiting\n");
			exit(1);
		}
	
		/* compare headers for west horizon and mask */
		test = (mask_nc_hdr.ncols == whoriz_nc_hdr.ncols)*(mask_nc_hdr.nrows == whoriz_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == whoriz_nc_hdr.ulx)*(mask_nc_hdr.uly == whoriz_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == whoriz_nc_hdr.cellsize);
		if (!test)
		{
			printf("Header information differs between mask and whoriz...exiting\n");
			exit(1);
		}
	
		/* compare headers for average horizon and mask */
		test = (mask_nc_hdr.ncols == avghoriz_nc_hdr.ncols)*(mask_nc_hdr.nrows == avghoriz_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == avghoriz_nc_hdr.ulx)*(mask_nc_hdr.uly == avghoriz_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == avghoriz_nc_hdr.cellsize);
		if (!test)
		{
			printf("Header information differs between mask and avghoriz...exiting\n");
			exit(1);
		}
	
	}
	else /* LIST input */
	{
		/* read list file header and verify */
		fread(&listhdr, sizeof(lsthdr_struct), 1, list_f.ptr);
		if (strcmp(listhdr.tag,LSTTAG))
		{
			printf("List not registered. Use register_list first.\n");
			exit(1);
		}

		/* check that there are five elements per cell (x,y,z,slp,asp) */
		if (listhdr.nvals != 5)
		{
			printf("Expecting 5 values per cell in list file... exiting\n");
			exit(1);
		}
	}

	/* read headers from TMAX, TMIN, and PRCP files */
	if (fread(&tmaxhdr, sizeof(dayouthdr_struct), 1, tmax_f.ptr) != 1)
	{
		printf("Error reading header from Tmax daily output file\n");
		exit(1);
	}
	if (fread(&tminhdr, sizeof(dayouthdr_struct), 1, tmin_f.ptr) != 1)
	{
		printf("Error reading header from Tmin daily output file\n");
		exit(1);
	}
	if (fread(&prcphdr, sizeof(dayouthdr_struct), 1, prcp_f.ptr) != 1)
	{
		printf("Error reading header from Prcp daily output file\n");
		exit(1);
	}
	
	/* error checking on the tmax and tmin daily input files */
	if (strcmp(tmaxhdr.tag, TAIRTAG))
	{
		printf("%s is not a daily air temperature output file\n",tmax_f.name);
		exit(1);
	}
	if (strcmp(tminhdr.tag, TAIRTAG))
	{
		printf("%s is not a daily air temperature output file\n",tmin_f.name);
		exit(1);
	}
	if (strcmp(prcphdr.tag, PRCPTAG))
	{
		printf("%s is not a daily precipitation output file\n",prcp_f.name);
		exit(1);
	}

	test = ((tmaxhdr.start_yday == tminhdr.start_yday) && (tmaxhdr.start_yday == prcphdr.start_yday)) &&
		   ((tmaxhdr.start_year == tminhdr.start_year) && (tmaxhdr.start_year == prcphdr.start_year)) &&
		   ((tmaxhdr.ndays == tminhdr.ndays) && (tmaxhdr.ndays == prcphdr.ndays));
	if (!test)
	{
		printf("Error\n");
		printf("Tmax, Tmin, and/or Prcp daily input files have different start dates\n");
		printf("and/or different number of days\n");
		printf("Tmax start_yday = %d\n",tmaxhdr.start_yday);
		printf("Tmin start_yday = %d\n",tminhdr.start_yday);
		printf("Prcp start_yday = %d\n",prcphdr.start_yday);
		printf("Tmax start_year = %d\n",tmaxhdr.start_year);
		printf("Tmin start_year = %d\n",tminhdr.start_year);
		printf("Prcp start_year = %d\n",prcphdr.start_year);
		printf("Tmax ndays = %d\n",tmaxhdr.ndays);
		printf("Tmin ndays = %d\n",tminhdr.ndays);
		printf("Prcp ndays = %d\n",prcphdr.ndays);
		exit(1);
	}
        /* compare projection parameters between map (or list) and tmax/tmin/prcp */
        if (map)
        {
		/* set the local gctp structure */
		gctp = mask_nc_hdr.gctp;
	}
	else
	{
		/* list input */
		if (memcmp(&tmaxhdr.listhdr.gctp,&listhdr.gctp,sizeof(gctp_struct)))
		{
			printf("Different projection in list and tmax input.\n");
			exit(1);
		}
		if (memcmp(&tminhdr.listhdr.gctp,&listhdr.gctp,sizeof(gctp_struct)))
		{
			printf("Different projection in list and tmin input.\n");
			exit(1);
		}
		if (memcmp(&prcphdr.listhdr.gctp,&listhdr.gctp,sizeof(gctp_struct)))
		{
			printf("Different projection in list and prcp input.\n");
			exit(1);
		}
		
		/* set the local gctp structure */
		gctp = listhdr.gctp;
	}
	
	/* set internal variables */
	if (map)
	{
		ncols = mask_nc_hdr.ncols;
		nrows = mask_nc_hdr.nrows;
		ncells = ncols * nrows;
		map_ulx = mask_nc_hdr.ulx;
		map_uly = mask_nc_hdr.uly;
		cellsize = mask_nc_hdr.cellsize;
	}
	else
	{
		ncells = listhdr.ncells;
	}
	ndays = tmaxhdr.ndays;

	// * initialize byte scaling parameters for tmax, tmin, and prcp */
	// if (init_byte_scale(tmaxhdr.valmax,tmaxhdr.valmin,tmaxhdr.bytmax,
	// 	tmaxhdr.bytmin,&tmax_scale_ratio,&tmax_brange,&tmax_bmin))
	// {
	// 	printf("Error in init_byte_scale ... exiting\n");
	// 	exit(1);
	// }
	// tmax_bytrange = (double)tmax_brange;
	// tmax_valmin = tmaxhdr.valmin;
	
	// if (init_byte_scale(tminhdr.valmax,tminhdr.valmin,tminhdr.bytmax,
	// 	tminhdr.bytmin,&tmin_scale_ratio,&tmin_brange,&tmin_bmin))
	// {
	// 	printf("Error in init_byte_scale ... exiting\n");
	// 	exit(1);
	// }
	// tmin_bytrange = (double)tmin_brange;
	// tmin_valmin = tminhdr.valmin;
	
	// if (init_byte_scale(prcphdr.valmax,prcphdr.valmin,prcphdr.bytmax,
	// 	prcphdr.bytmin,&prcp_scale_ratio,&prcp_brange,&prcp_bmin))
	// {
	// 	printf("Error in init_byte_scale ... exiting\n");
	// 	exit(1);
	// }
	// prcp_bytrange = (double)prcp_brange;
	// prcp_valmin = prcphdr.valmin;
	
	// /* initialize byte scaling parameters for srad, dayl, vp, and swe */
	// if (init_byte_scale(SRAD_VALMAX, SRAD_VALMIN, SRAD_BYTMAX, SRAD_BYTMIN,
	// 	&srad_scale_ratio, &srad_brange, &srad_bmin))
	// {
	// 	printf("Error in init_byte_scale() for srad\n");
	// 	exit(1);
	// }
	// srad_bytrange = (double)srad_brange;
	// srad_valmin = SRAD_VALMIN;
	
	// if (init_byte_scale(DAYL_VALMAX, DAYL_VALMIN, DAYL_BYTMAX, DAYL_BYTMIN,
	// 	&dayl_scale_ratio, &dayl_brange, &dayl_bmin))
	// {
	// 	printf("Error in init_byte_scale() for dayl\n");
	// 	exit(1);
	// }
	// dayl_bytrange = (double)dayl_brange;
	// dayl_valmin = DAYL_VALMIN;
	
	// if (init_byte_scale(VP_VALMAX, VP_VALMIN, VP_BYTMAX, VP_BYTMIN,
	// 	&vp_scale_ratio, &vp_brange, &vp_bmin))
	// {
	// 	printf("Error in init_byte_scale() for vp\n");
	// 	exit(1);
	// }
	// vp_bytrange = (double)vp_brange;
	// vp_valmin = VP_VALMIN;
	
	if (snow_flag)
	{
		// if (init_byte_scale(SWE_VALMAX, SWE_VALMIN, SWE_BYTMAX, SWE_BYTMIN,
		// 	&swe_scale_ratio, &swe_brange, &swe_bmin))
		// {
		// 	printf("Error in init_byte_scale() for swe\n");
		// 	exit(1);
		// }
		// swe_bytrange = (double)swe_brange;
		// swe_valmin = SWE_VALMIN;
	}
	
	/* fill the daily output file header structures (srad, daylength, and vp) */
	sradhdr.start_yday = daylhdr.start_yday = vphdr.start_yday = swehdr.start_yday = tmaxhdr.start_yday;
	sradhdr.start_year = daylhdr.start_year = vphdr.start_year = swehdr.start_year = tmaxhdr.start_year;
	strcpy(sradhdr.tag, SRADTAG);
	strcpy(daylhdr.tag, DAYLTAG);
	strcpy(vphdr.tag, VPTAG);
	strcpy(swehdr.tag, SWETAG);
	sradhdr.valmax = SRAD_VALMAX;
	sradhdr.valmin = SRAD_VALMIN;
	sradhdr.bytmax = SRAD_BYTMAX;
	sradhdr.bytmin = SRAD_BYTMIN;
	daylhdr.valmax = DAYL_VALMAX;
	daylhdr.valmin = DAYL_VALMIN;
	daylhdr.bytmax = DAYL_BYTMAX;
	daylhdr.bytmin = DAYL_BYTMIN;
	vphdr.valmax = VP_VALMAX;
	vphdr.valmin = VP_VALMIN;
	vphdr.bytmax = VP_BYTMAX;
	vphdr.bytmin = VP_BYTMIN;
	swehdr.valmax = SWE_VALMAX;
	swehdr.valmin = SWE_VALMIN;
	swehdr.bytmax = SWE_BYTMAX;
	swehdr.bytmin = SWE_BYTMIN;
	sradhdr.ndays = daylhdr.ndays = vphdr.ndays = swehdr.ndays = ndays;
	sradhdr.ioflag = daylhdr.ioflag = vphdr.ioflag = swehdr.ioflag = map;
	if (map)
	{
		sradhdr.maskhdr = daylhdr.maskhdr = vphdr.maskhdr = swehdr.maskhdr = mask_nc_hdr;
	}
	else
	{
		sradhdr.listhdr = daylhdr.listhdr = vphdr.listhdr = swehdr.listhdr = listhdr;
	}

	/* create output filenames and open files for binary write */
	if (bin_out(&srad_out_f, outprefix, ".dmsrad_out")) exit(1);
	if (bin_out(&vp_out_f, outprefix, ".dmvp_out")) exit(1);
	if (bin_out(&srad_tavg_f, outprefix, ".dmsrad_tavg")) exit(1);
	if (bin_out(&vp_tavg_f, outprefix, ".dmvp_tavg")) exit(1);
	if (bin_out(&srad_savg_f, outprefix, ".dmsrad_savg")) exit(1);
	if (bin_out(&vp_savg_f, outprefix, ".dmvp_savg")) exit(1);
	if (snow_flag)
	{
		if (bin_out(&swe_out_f, outprefix, ".dmswe_out")) exit(1);
		if (bin_out(&swe_tavg_f, outprefix, ".dmswe_tavg")) exit(1);
		if (bin_out(&swe_savg_f, outprefix, ".dmswe_savg")) exit(1);
	}
	if (dayl_out)
	{
		if (bin_out(&dayl_out_f, outprefix, ".dmdayl_out")) exit(1);
		if (bin_out(&dayl_tavg_f, outprefix, ".dmdayl_tavg")) exit(1);
		if (bin_out(&dayl_savg_f, outprefix, ".dmdayl_savg")) exit(1);
	}
	
	/* write output headers to the daily output files */
	if (fwrite(&sradhdr, sizeof(dayouthdr_struct), 1, srad_out_f.ptr) != 1)
	{
		printf("Error writing header to srad outfile\n");
		exit(1);
	}
	if (fwrite(&vphdr, sizeof(dayouthdr_struct), 1, vp_out_f.ptr) != 1)
	{
		printf("Error writing header to vp outfile\n");
		exit(1);
	}
	if (snow_flag)
	{
		if (fwrite(&swehdr, sizeof(dayouthdr_struct), 1, swe_out_f.ptr) != 1)
		{
			printf("Error writing header to swe outfile\n");
			exit(1);
		}
	}
	if (dayl_out)
	{
		if (fwrite(&daylhdr, sizeof(dayouthdr_struct), 1, dayl_out_f.ptr) != 1)
		{
			printf("Error writing header to dayl outfile\n");
			exit(1);
		}
	}
	
	/* initialize some parameters */
	dh = dt / SECPERRAD; /* calculates hour-angle step */

	/* allocate memory for output */
	if (!(tmax_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(btmin_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(bprcp_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(yday_ts = (int*) malloc(ndays * sizeof(int)))) exit(1);
	if (!(tmin_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(tday_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(prcp_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(dtr_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(smooth_dtr_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(dayl_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(srad_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(vp_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(swe_ts = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(dayl_savg = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(srad_savg = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(vp_savg = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(swe_savg = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(parray = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(window = (double*) malloc((ndays+90)*sizeof(double)))) exit(1);
	if (!(t_fmax = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(tdew = (double*) malloc(ndays * sizeof(double)))) exit(1);
	if (!(save_pet = (double*) malloc(ndays * sizeof(double)))) exit(1);

	/* zero the time series output arrays */
	for (day=0 ; day<ndays ; day++)
	{
		dayl_ts[day] = srad_ts[day] = vp_ts[day] =swe_ts[day] = 0.0;
		dayl_savg[day] = srad_savg[day] = vp_savg[day] = swe_savg[day] = 0.0;
	}
	
	/* fill the yearday timeseries */
	start_yday = sradhdr.start_yday;
	start_year = sradhdr.start_year;
	for (day=0 ; day<ndays ; day++)
	{
		if (yearday(start_yday, start_year, day, &yday_ts[day], &junk))
		{
			printf("Error in yearday()\n");
			exit(1);
		}
		
		/* force Dec 31 = Dec 30 for a leap year */
		if (yday_ts[day] == 365)
		{
			yday_ts[day] = 364;
		}
	}

	/* fill the sin and cos arrays */
	/* generate declination arrays */
	for (day=0 ; day<365 ; day++)
	{
		angle = MINDECL * cos(((double)day + DAYSOFF) * RADPERDAY);
		sin_decl[day] = sin(angle);
		cos_decl[day] = cos(angle);
	}
	/* generate arrays for angles from 0 to 2pi, and -pi to pi */
	for (i=0 ; i<6284 ; i++)
	{
		angle = (double)i / 1000.0;
		sin_0to2pi[i] = sin(angle);
		cos_0to2pi[i] = cos(angle);
		angle = ((double)(i) - 3141.59) / 1000.0;
		sin_negpitopi[i] = sin(angle);
		cos_negpitopi[i] = cos(angle);
	}

	/* initialize inverse transformation function (x,y) --> (lon,lat) */
	inv_init(gctp.outsys, gctp.outzone, gctp.outparm, gctp.outdatum,
		gctp.fn27, gctp.fn83, &iflg, inv_trans);

	/* begin simulation, loop through ncells */
	n_in = 0;
	for (i=0 ; i<ncells ; i++)
	{
		srad_tavg = 0.0;
		dayl_tavg = 0.0;
		vp_tavg = 0.0;
		swe_tavg = 0.0;
		
		/* read mask value if MAP I/O */
		if (map)
		{
			mask     =          mask_array[i];

			if (mask != good_mask)
			{
				/* skip to next cell */
				continue;
			}
			
			/* this cell is in the simulation region */
			/* calculate column and row for map_x and map_y coords */
			n_in ++;
			
			index[0]= row = i/ncols;
			index[1]= col = i%ncols;
			map_x = ((double)col * cellsize) + map_ulx;
			map_y = map_uly - ((double)row * cellsize);
			/* read single values out of netcdf files for this point */
			status=nc_get_var1_int(ncid_dem, imageid_dem, index, &z);
			status=nc_get_var1_float(ncid_slp, imageid_slp, index, &slp);
			status=nc_get_var1_float(ncid_asp, imageid_asp, index, &asp);
			status=nc_get_var1_float(ncid_ehoriz, imageid_ehoriz, index, &ehoriz);
			status=nc_get_var1_float(ncid_whoriz, imageid_whoriz, index, &whoriz);
			status=nc_get_var1_float(ncid_avghoriz, imageid_avghoriz, index, &avghoriz);
			
			elev = (double)z;
			slope    = (double) slp;
			aspect   = (double) asp;
		}
		else
		{
			/* list input */
			if (fread(celldat, sizeof(double), 5, list_f.ptr) != 5)
			{
				printf("Error reading from list file... exiting\n");
				exit(1);
			}
			map_x = celldat[0];
			map_y = celldat[1];
			elev = celldat[2];
			slope = celldat[3];
			aspect = celldat[4];
		}
		
		/* convert map_x, map_y to lon,lat to get latitude */
		inv_trans[gctp.outsys](map_x, map_y, &longitude, &latitude);
		
		/* read ndays of data from tmax, tmin, and prcp input files */
		if (fread(tmax_ts, sizeof(double), ndays, tmax_f.ptr) != ndays)
		{
			printf("Error reading from tmax input file \n");
			exit(1);
		}
		if (fread(btmin_ts, sizeof(double), ndays, tmin_f.ptr) != ndays)
		{
			printf("Error reading from tmin input file\n");
			exit(1);
		}
		if (fread(bprcp_ts, sizeof(double), ndays, prcp_f.ptr) != ndays)
		{
			printf("Error reading from prcp input file\n");
			exit(1);
		}
		
		/* convert byte tmax, tmin to doubles and build DTR time series */
		avgdtr = 0.0;
		for (day=0 ; day<ndays ; day++)
		{
			// tmax = (((double)tmax_ts[day] - (double)tmax_bmin)/tmax_scale_ratio) + tmax_valmin;
			// tmin = (((double)btmin_ts[day] - (double)tmin_bmin)/tmin_scale_ratio) + tmin_valmin;
			tmax = (double)tmax_ts[day];
			tmin = (double)btmin_ts[day];
			tmean = (tmax + tmin)/2.0;
			dtr_ts[day] = tmax-tmin;
			if (dtr_ts[day] < 0.0) 
			{
				dtr_ts[day] = 0;
				tmin = tmax = tmean;
			}
			tmin_ts[day] = tmin;
			tday_ts[day] = ((tmax - tmean)*TDAYCOEF) + tmean;
			// prcp_ts[day] = (((double)bprcp_ts[day] - (double)prcp_bmin)/prcp_scale_ratio) + prcp_valmin;
			prcp_ts[day] = (double)bprcp_ts[day];
			avgdtr += dtr_ts[day];
		}
		avgdtr /= (double)ndays;
		
		/* smooth dtr array: After Bristow and Campbell, 1984 */
		if (ndays >= 30) /* use 30-day antecedent smoothing window */
		{
			if (pulled_boxcar(dtr_ts, smooth_dtr_ts, ndays, 30, 0))
			{
				printf("Error in boxcar smoothing, predict_daily_srad.c\n");
				exit(1);
			}
		}
		else /* smoothing window width = ndays */
		{
			if (pulled_boxcar(dtr_ts, smooth_dtr_ts, ndays, ndays, 0))
			{
				printf("Error in boxcar smoothing, predict_daily_srad.c\n");
				exit(1);
			}
		}
					
		/* check for flat gridcell, slope < 1.0 degrees */
		/* slope, aspect, and latitude are indexed by thousandths of
		radians	*/
		if (slope < RADPERDEG)
		{
			flat_flag = 1;
			slp_ind = 0;
			asp_ind = 0;
		}
		else
		{
			flat_flag = 0;
			slp_ind = (int)(slope*1000.0);  
			asp_ind = (int)(aspect*1000.0);
		}
		lat_ind = (int)(latitude*1000.0);
		
		/* calculate the annual total precip for decision between
		simple and arid-corrected humidity algorithm */
		sum_prcp = 0.0;
		for (day=0 ; day<ndays ; day++)
		{
			sum_prcp += prcp_ts[day];
		}
		ann_prcp = (sum_prcp/(double)ndays) * 365.25;
		if (ann_prcp == 0.0) ann_prcp = 1.0;

		/* Generate the effective annual precip, based on a 3-month
		moving-window. Requires some special case handling for the
		beginning of the record and for short records. */
		/* check if there are at least 90 days in this input file, if not,
		use a simple total scaled to effective annual precip */
		if (ndays < 90)
		{
			sum_prcp = 0.0;
			for (day=0 ; day<ndays ; day++)
			{
				sum_prcp += prcp_ts[day];
			}
			effann_prcp = (sum_prcp/(double)ndays) * 365.25;
			/* if the effective annual precip for this period
			is less than 8 cm, set the effective annual precip to 8 cm
			to reflect an arid condition, while avoiding possible
			division-by-zero errors and very large ratios (PET/Pann) */
			if (effann_prcp < 8.0) effann_prcp = 8.0;
			for (day=0 ; day<ndays ; day++)
			{
				parray[day] = effann_prcp;
			}
		}
		else
		{
			/* Check if the yeardays at beginning and the end of this input file
			match up. If so, use parts of the three months at the end
			of the input file to generate effective annual precip for
			the first 3-months. Otherwise, duplicate the first 90 days
			of the record. Assuming base-0 treatment for yeardays */
			start_yday = yday_ts[0];
			end_yday = yday_ts[ndays-1];
			if (start_yday != 0)
			{
				isloop = (end_yday == start_yday-1) ? 1 : 0;
			}
			else
			{
				isloop = (end_yday == 364 || end_yday == 365) ? 1 : 0;
			}

			/* fill the first 90 days of window */
			for (day=0 ; day<90 ; day++)
			{
				if (isloop) window[day] = prcp_ts[ndays-90+day];
				else window[day] = prcp_ts[day];
			}
			/* fill the rest of the window array */
			for (day=0 ; day<ndays ; day++)
			{
				window[day+90] = prcp_ts[day];
			}

			/* for each day, calculate the effective annual precip from 
			scaled 90-day total */
			for (day=0 ; day<ndays ; day++)
			{
				sum_prcp = 0.0;
				for (j=0 ; j<90 ; j++)
				{
					sum_prcp += window[day+j];
				}
				sum_prcp = (sum_prcp/90.0) * 365.25;
				/* if the effective annual precip for this 90-day period
				is less than 8 cm, set the effective annual precip to 8 cm
				to reflect an arid condition, while avoiding possible
				division-by-zero errors and very large ratios (PET/Pann) */
				parray[day] = (sum_prcp < 8.0) ? 8.0 : sum_prcp;
			}
		} /* end if ndays >= 90 */	
	
		if (snow_flag)
		{
			/* begin snowpack calculation for albedo effect */
			/* first pass to initialize SWE array */
			snowpack = 0.0;
			for (day=0 ; day<ndays ; day++)
			{
				newsnow = 0.0;
				snowmelt = 0.0;
				if (tmin_ts[day] <= SNOW_TCRIT) newsnow = prcp_ts[day];
				else snowmelt = SNOW_TRATE * (tmin_ts[day] - SNOW_TCRIT);
				snowpack += newsnow - snowmelt;
				if (snowpack < 0.0) snowpack = 0.0;
				swe_ts[day] = snowpack;
			}
			/* use the first pass to set the initial snowpack conditions for the
			first day of data */
			start_yday = yday_ts[0];
			if (start_yday == 0) prev_yday = 364;
			else prev_yday = start_yday-1;
			count = 0;
			sum = 0.0;
			for (day=1 ; day<ndays ; day++)
			{
				if (yday_ts[day] == start_yday || yday_ts[day] == prev_yday)
				{
					count ++;
					sum += swe_ts[day];
				}
			}
			/* Proceed with correction if there are valid days to reinitialize
			the snowpack estimates. Otherwise use the first-pass estimate. */
			if (count)
			{
				snowpack = sum/(double)count;
				for (day=0 ; day<ndays ; day++)
				{
					newsnow = 0.0;
					snowmelt = 0.0;
					if (tmin_ts[day] <= SNOW_TCRIT) newsnow = prcp_ts[day];
					else snowmelt = SNOW_TRATE * (tmin_ts[day] - SNOW_TCRIT);
					snowpack += newsnow - snowmelt;
					if (snowpack < 0.0) snowpack = 0.0;
					swe_ts[day] = snowpack;
				}
			}
		} /* end if snow_flag */

		/*****************************************
		 *                                       *
		 * start of the main radiation algorithm *
		 *                                       *
		 *****************************************/

		/* before starting the iterative algorithm between humidity and 
		radiation, calculate all the variables that don't depend on 
		humidity so they only get done once. */

		/* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
		t1 = 1.0 - (LR_STD * elev)/T_STD;
		t2 = G_STD / (LR_STD * (R/MA));
		pratio = pow(t1,t2);

		/* STEP (2) correct initial transmittance for elevation */ 
		trans1 = pow(TBASE,pratio);
		
		/* STEP (3) build 366-day array of ttmax0, potential rad, and daylength */
		/* initialize the parameters that don't change with yday */
		coslat = cos_negpitopi[lat_ind + 3142];
		sinlat = sin_negpitopi[lat_ind + 3142];
		cosslope = cos_0to2pi[slp_ind];
		sinslope = sin_0to2pi[slp_ind];
		coshalfslope = cos_0to2pi[slp_ind/2];
		cosasp = cos_0to2pi[asp_ind];
		sinasp = sin_0to2pi[asp_ind];
		coszeh_ind = (int)((HALFPI-ehoriz)*1000.0);
		coszwh_ind = (int)((HALFPI-whoriz)*1000.0);
		coszeh = cos_0to2pi[coszeh_ind];
		coszwh = cos_0to2pi[coszwh_ind];

		/* begin loop through yeardays */
		for (day=0 ; day<365 ; day++)
		{
			/* generate cosine and sine parameters from arrays */
			yday = yday_ts[day];
			cosdecl = cos_decl[yday];
			sindecl = sin_decl[yday];

			/* do some precalculations for beam-slope geometry (bsg) */
			bsg1 = -sinslope * sinasp * cosdecl;
			bsg2 = (-cosasp * sinslope * sinlat + cosslope * coslat) * cosdecl;
			bsg3 = (cosasp * sinslope * coslat + cosslope * sinlat) * sindecl;

			/* calculate daylength as a function of lat and decl */
			cosegeom = coslat * cosdecl;
			sinegeom = sinlat * sindecl;
			coshss = -(sinegeom) / cosegeom;
			if (coshss < -1.0) coshss = -1.0;  /* 24-hr daylight */
			if (coshss > 1.0) coshss = 1.0;    /* 0-hr daylight */
			hss = acos(coshss);                /* hour angle at sunset (radians) */
			/* daylength (seconds) */
			daylength[day] = 2.0 * hss * SECPERRAD;

			/* solar constant as a function of yearday (W m-2) */
			sc = 1368.0 + 45.5*sin((2.0*PI*(double)day/365.25) + 1.7);
			/* extraterrestrial radiation perpendicular to beam, total over
			the timestep (J m-2) */
			dir_beam_topa = sc * dt;

			sum_trans = 0.0;
			sum_flat_potrad = 0.0;
			sum_slope_potrad = 0.0;

			/* begin sub-daily hour-angle loop, from -hss to hss */
			for (h=-hss ; h<hss ; h+=dh)
			{
				/* precalculate cos and sin of hour angle */
				hh = (int)(h*1000.0);
				cosh = cos_negpitopi[hh + 3142];
				sinh = sin_negpitopi[hh + 3142];

				/* calculate cosine of solar zenith angle */
				cza = cosegeom * cosh + sinegeom;

				/* calculate cosine of beam-slope angle */
				cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;

				/* check if sun is above a flat horizon */
				if (cza > 0.0) 
				{
					/* when sun is above the ideal (flat) horizon, do all the
					flat-surface calculations to determine daily total
					transmittance, and save flat-surface potential radiation
					for later calculations of diffuse radiation */

					/* potential radiation for this time period, flat surface,
					top of atmosphere */
					dir_flat_topa = dir_beam_topa * cza;

					/* determine optical air mass */
					am = 1.0/(cza + 0.0000001);
					if (am > 2.9)
					{
						ami = (int)(acos(cza)/RADPERDEG) - 69;
						if (ami < 0) ami = 0;
						if (ami > 20) ami = 20;
						am = optam[ami];
					}

					/* correct instantaneous transmittance for this optical
					air mass */
					trans2 = pow(trans1,am);

					/* instantaneous transmittance is weighted by potential
					radiation for flat surface at top of atmosphere to get
					daily total transmittance */
					sum_trans += trans2 * dir_flat_topa;

					/* keep track of total potential radiation on a flat
					surface for ideal horizons */
					sum_flat_potrad += dir_flat_topa;

					/* keep track of whether this time step contributes to
					component 1 (direct on slope) */
					if ((h<0.0 && cza>coszeh && cbsa>0.0) ||
						(h>=0.0 && cza>coszwh && cbsa>0.0))
					{
						/* sun between east and west horizons, and direct on
						slope. this period contributes to component 1 */
						sum_slope_potrad += dir_beam_topa * cbsa;
					}

				} /* end if sun above ideal horizon */

			} /* end of sub-daily hour-angle loop */

			/* calculate maximum daily total transmittance and daylight average
			flux density for a flat surface and the slope */
			if (daylength[day])
			{
				ttmax0[day] = sum_trans / sum_flat_potrad;
				flat_potrad[day] = sum_flat_potrad / daylength[day];
				slope_potrad[day] = sum_slope_potrad / daylength[day];
			}
			else
			{
				ttmax0[day] = 0.0;
				flat_potrad[day] = 0.0;
				slope_potrad[day] = 0.0;
			}
		} /* end of 365 days loop */

		/* force yearday 366 = yearday 365 */
		ttmax0[365] = ttmax0[364];
		flat_potrad[365] = flat_potrad[364];
		slope_potrad[365] = slope_potrad[364];
		daylength[365] = daylength[364];

		/* STEP (4)  calculate the sky proportion for diffuse radiation */
		/* uses the product of spherical cap defined by average horizon angle
		and the great-circle truncation of a hemisphere. this factor does not
		vary by yearday. */
		horizon_scalar = 1.0 - sin(avghoriz);
		if (slope > avghoriz) slope_excess = slope - avghoriz;
		else slope_excess = 0.0;
		if (avghoriz > HALFPI) slope_scalar = 0.0;
		else
		{
			slope_scalar = 1.0 - (slope_excess/(PI - 2.0*avghoriz));
			if (slope_scalar < 0.0) slope_scalar = 0.0;
		}
		sky_prop = horizon_scalar * slope_scalar;

		//BWM debugging
	        //printf("sky_prop: %lf\n",sky_prop);		
		/* STEP (5)  some variables can be calculated outside the
		  iteration loop for srad and humidity ... */ 
		/* b parameter, and t_fmax not varying with Tdew, so these can be
		calculated once, outside the iteration between radiation and humidity
		estimates. Requires storing t_fmax in an array. */
		for (day=0 ; day<ndays ; day++)
		{	
			/* b parameter from 30-day average of DTR */
			b = B0 + B1 * exp(-B2 * smooth_dtr_ts[day]);

			/* proportion of daily maximum transmittance */
			t_fmax[day] = 1.0 - 0.9 * exp(-b * pow(dtr_ts[day],C));

			/* correct for precipitation if this is a rain day */
			if (prcp_ts[day] > 0.2) t_fmax[day] *= RAIN_SCALAR;
		}
		
		/* STEP (6) Begin iterative estimation of radiation and humidity */
		/* As a first approximation, calculate radiation assuming
		that Tdew = Tmin */
		for (day=0 ; day<ndays ; day++)
		{
			yday = yday_ts[day];
			tdew[day] = tmin_ts[day];
			pva = 610.7 * exp(17.38 * tdew[day] / (239.0 + tdew[day]));
			vp_ts[day] = pva;
			t_tmax = ttmax0[yday] + ABASE * pva;

			/* final daily total transmittance */
			t_final = t_tmax * t_fmax[day];

			if (t_final < 0) t_final = 0;
						
			/* estimate fraction of radiation that is diffuse, on an
			instantaneous basis, from relationship with daily total
			transmittance in Jones (Plants and Microclimate, 1992)
			Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
			Fig 6.14, p. 122. */
			pdif = -1.25*t_final + 1.25;
			if (pdif > 1.0) pdif = 1.0;
			if (pdif < 0.0) pdif = 0.0;

			/* estimate fraction of radiation that is direct, on an
			instantaneous basis */
			pdir = 1.0 - pdif;

			/* the daily total radiation is estimated as the sum of the
			following two components:
			1. The direct radiation arriving during the part of
			   the day when there is direct beam on the slope.
			2. The diffuse radiation arriving over the entire daylength
			   (when sun is above ideal horizon).
			*/

			/* component 1 */
			srad1 = slope_potrad[yday] * t_final * pdir;

			/* component 2 (diffuse) */
			/* includes the effect of surface albedo in raising the diffuse
			radiation for obstructed horizons */
			srad2 = flat_potrad[yday] * t_final * pdif * 
				(sky_prop + DIF_ALB*(1.0-sky_prop)); 

			/* snow pack influence on radiation */
			/* in the AgForMet paper, this correction is included as an
			absolute amount of radiation added when snow is present. This
			will not work well in high latitudes where daylength and
			potential radiation are low. I converted this to be a proportion
			of the daily radiation, so at low radiation, the correction is
			smaller. Based on the mean daily total radiation for the 24
			Austrian sites, and the original regression equation. The
			new relationship causes a correction from 11% to 37% of daily
			total radiation, for >0 to 30 cm SWE. */	
			if (swe_ts[day] > 0.0)
			{
				/* snow correction as a proportion */
				sc = 0.11 + 0.00867 * swe_ts[day];
				/* set a maximum correction of 37% */
				if (sc > 0.37) sc = 0.37;
			}
			else sc = 0.0;

			/* save daily radiation and daylength */
			srad_ts[day] = (srad1 + srad2) * (1.0+sc);

//BWM debugging
/*            if(180 == day)
            {
                printf("srad1: %lf 2: %lf ts: %lf\n", srad1, srad2, srad_ts[day]);
            }
*/
			dayl_ts[day] = daylength[yday];
			
		} /* end of first radiation iteration ndays loop */
	
		/* STEP (7)  Humidity estimation */
		/* the MTCLIM version of humidity correction causes a discontinuous
		estimation at the boundary between arid/non-arid. Changing slightly
		to fix this problem. */
		/* estimate air pressure at site */
		pa = pratio * P_STD;
		for (day=0 ; day<ndays ; day++)
		{
			pet = calc_pet(srad_ts[day],tday_ts[day],pa,dayl_ts[day]);
			/* if the effective annual PET, calculated assuming all
			days with this day's PET, is more than 2.5 times the
			effective annual precipitation (previously calculated as the
			30-day moving window of precipitation), then apply an aridity
			correction. The correction ramps up smoothly from a multiplier of
			1.0 at the ratio 2.5, to a multiplier 0.95 at the ratio 10.0. */
			ratio=(pet*365.25)/parray[day];
			if (ratio > 2.5)
			{
				tmink = tmin_ts[day] + 273.15;
				multiplier = 1.0166667 - 0.0066667*ratio;
				if (multiplier > 1.0) multiplier = 1.0;
				if (multiplier < 0.95) multiplier = 0.95;
				tdewk = tmink*multiplier;
				tdew[day] = tdewk - 273.15;
				
				/* Revise estimate of radiation using new Tdew, only for 
				days with the aridity correction */
				yday = yday_ts[day];
				pva = 610.7 * exp(17.38 * tdew[day] / (239.0 + tdew[day]));
				vp_ts[day] = pva;
				t_tmax = ttmax0[yday] + ABASE * pva;

				/* final daily total transmittance */
				t_final = t_tmax * t_fmax[day];

				if (t_final < 0) t_final = 0;

				/* estimate fraction of radiation that is diffuse, on an
				instantaneous basis, from relationship with daily total
				transmittance in Jones (Plants and Microclimate, 1992)
				Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
				Fig 6.14, p. 122. */
				pdif = -1.25*t_final + 1.25;
				if (pdif > 1.0) pdif = 1.0;
				if (pdif < 0.0) pdif = 0.0;

				/* estimate fraction of radiation that is direct, on an
				instantaneous basis */
				pdir = 1.0 - pdif;

				/* the daily total radiation is estimated as the sum of the
				following two components:
				1. The direct radiation arriving during the part of
				   the day when there is direct beam on the slope.
				2. The diffuse radiation arriving over the entire daylength
				   (when sun is above ideal horizon).
				*/

				/* component 1 */
				srad1 = slope_potrad[yday] * t_final * pdir;

				/* component 2 (diffuse) */
				/* includes the effect of surface albedo in raising the diffuse
				radiation for obstructed horizons */
				srad2 = flat_potrad[yday] * t_final * pdif * 
					(sky_prop + DIF_ALB*(1.0-sky_prop)); 

				/* snow pack influence on radiation */	
				if (swe_ts[day] > 0.0)
				{
					/* snow correction as a proportion */
					sc = 0.11 + 0.00867 * swe_ts[day];
					/* set a maximum correction of 37% */
					if (sc > 0.37) sc = 0.37;
				}
				else sc = 0.0;

				/* save updated daily radiation for this arid day */
				srad_ts[day] = (srad1 + srad2) * (1.0+sc);
				
				/* Final iteration step, use new radiaiton to do new 
				PET calculation, then correct tdew from tmin and calculate
				vapor pressure */
				pet = calc_pet(srad_ts[day],tday_ts[day],pa,dayl_ts[day]);
				ratio=(pet*365.25)/parray[day];
				multiplier = 1.0166667 - 0.0066667*ratio;
				if (multiplier > 1.0) multiplier = 1.0;
				if (multiplier < 0.95) multiplier = 0.95;
				tdewk = tmink*multiplier;
				tdew[day] = tdewk - 273.15;
				/* calculate vapor pressure from updated tdew */
				pva = 610.7 * exp(17.38 * tdew[day] / (239.0 + tdew[day]));
				vp_ts[day] = pva;
				
			} /* end arid correction */
			
		} /* end days loop testing for arid correction */

		/* byte scale and save in time series output arrays */
		for (day=0 ; day<ndays ; day++)
		{
			/* srad */
			srad_tavg += srad_ts[day];
			srad_savg[day] += srad_ts[day];

			/* vapor pressure */
			vp_tavg += vp_ts[day];
			vp_savg[day] += vp_ts[day];

			/* snow water equivalent */
			if (snow_flag)
			{
				swe_tavg += swe_ts[day];
				swe_savg[day] += swe_ts[day];
			}
			
			/* daylength output is optional */
			if (dayl_out)
			{
				dayl_tavg += dayl_ts[day];
				dayl_savg[day] += dayl_ts[day];
			}
		} /* end of day = ndays loop */
			
		/* generate average for image arrays and write to output arrays */
		srad_tavg /= (double)ndays;
		vp_tavg /= (double)ndays;
		fwrite(&srad_tavg, sizeof(double), 1, srad_tavg_f.ptr);
		fwrite(&vp_tavg, sizeof(double), 1, vp_tavg_f.ptr);
		if (snow_flag)
		{
			swe_tavg /= (double)ndays;
			fwrite(&swe_tavg, sizeof(double), 1, swe_tavg_f.ptr);
		}
		if (dayl_out)
		{
			dayl_tavg /= (double)ndays;
			fwrite(&dayl_tavg, sizeof(double), 1, dayl_tavg_f.ptr);
		}
		
		/* write daily output */
		fwrite(srad_ts, sizeof(double), ndays, srad_out_f.ptr); 
		fwrite(vp_ts, sizeof(double), ndays, vp_out_f.ptr); 
		if (snow_flag)
		{
			fwrite(swe_ts, sizeof(double), ndays, swe_out_f.ptr); 
		}
		if (dayl_out)
		{
			fwrite(dayl_ts, sizeof(double), ndays, dayl_out_f.ptr);
		}

	} /* end of i=ncells loop */
	
	/* close netcdf files */
	nc_close(ncid_dem);
	nc_close(ncid_slp);
	nc_close(ncid_asp);
	nc_close(ncid_ehoriz);
	nc_close(ncid_whoriz);
	nc_close(ncid_avghoriz);

	/* generate spatial averages for timeseries arrays */
	for (day=0 ; day<ndays ; day++)
	{
		srad_savg[day] /= (double) n_in;
		vp_savg[day] /= (double) n_in;
		if (snow_flag)
		{
			swe_savg[day] /= (double) n_in;
		}
		if (dayl_out)
		{
			dayl_savg[day] /= (double) n_in;
		}
	}
	
	/* write time series data to files */
/*    for(day = 0; i< ndays; i++)
    {
        printf("srad: %lf\n", srad_savg[i]);
    }
*/
	fwrite(srad_savg, sizeof(double), ndays, srad_savg_f.ptr);
	fwrite(vp_savg, sizeof(double), ndays, vp_savg_f.ptr);
	if (snow_flag)
	{
		fwrite(swe_savg, sizeof(double), ndays, swe_savg_f.ptr);
	}
	if (dayl_out)
	{
		fwrite(dayl_savg, sizeof(double), ndays, dayl_savg_f.ptr);
	}
	
	/* close files and exit */
	if (!map) fclose(list_f.ptr);
	fclose(tmax_f.ptr);
	fclose(tmin_f.ptr);
	fclose(srad_out_f.ptr);
	fclose(srad_tavg_f.ptr);
	fclose(srad_savg_f.ptr);
	fclose(vp_out_f.ptr);
	fclose(vp_tavg_f.ptr);
	fclose(vp_savg_f.ptr);
	if (snow_flag)
	{
		fclose(swe_out_f.ptr);
		fclose(swe_tavg_f.ptr);
		fclose(swe_savg_f.ptr);
	}
	if (dayl_out)
	{
		fclose(dayl_out_f.ptr);
		fclose(dayl_tavg_f.ptr);
		fclose(dayl_savg_f.ptr);
	}
	
        printf ("Execution of predict_daily_srad_vp_nc complete.\n");

	return(0);

} /* end of main */

	

