/*
predict_daily_tair_xval.c
Peter Thornton

3/31/00
Based on cmd_predict_daily_tair_xval_autopt.c, but instead of searching
for the optimal value for gsp, this routine just uses the input init_gsp
and does one xval calculation for that value. Output is the same format
as for autopt code.

An ammended version of xval1 that automates the optimal choice of gsp
given all the other parameters.

Generates the single-set cross-validation error statistics for air temperature
predictions using the interpolate and predict_daily_tair algorithms.

3/16/97
Changed the way that rad90 and lr stdev are being calculated, to remove the
temporal component from the variance calculation, as in the tastdev code.
12/9/09, PET: modified to use only the in_tile stations in calculation of mean
and stdev of tair, when inflag is set.
1/16/10, PET: converted output to netcdf, added daily obs and predictions
for each station in tile as ouput. Also added stnx,y,z to outfile.
9/20/2015, PET: Adding code to force the xval routine to use an exact number of stations
for each prediction point, in place of the iterative station density algorithm. 
Implemented using a cpp define EXACT_ANS. THis is one of the changes for Dayemt V3 processing.
*/

#include "metsrc2.h"
#include "netcdf.h"
#define EXACT_ANS

/* the following array contains the fractional filter radius needed to get
90% of the filter weight for filters having GSP in the set 0.5, 1.0, 1.5, ...,
25.0.  Used below to determine the mean radius of 90% influence for the
optimized gsp value */
double f90[50] = {
0.815295, 0.802175, 0.787455, 0.771103, 0.753171, 0.733808, 0.713272, 0.691912,
0.670136, 0.648369, 0.626999, 0.606345, 0.586638, 0.568018, 0.550549, 0.534235,
0.519040, 0.504903, 0.491747, 0.479494, 0.468063, 0.457380, 0.447377, 0.437989,
0.429159, 0.420838, 0.412978, 0.405541, 0.398491, 0.391795, 0.385425, 0.379356,
0.373564, 0.368030, 0.362734, 0.357661, 0.352794, 0.348122, 0.343630, 0.339307,
0.335144, 0.331130, 0.327257, 0.323517, 0.319902, 0.316405, 0.313021, 0.309744,
0.306567, 0.303486};

int main(int argc, char *argv[])
{
	/* variable declarations */
	interpolate_struct intstr;
	stnmetahdr_struct metahdr;
	stnmeta_struct metastr;
	daydathdr_struct datahdr;
	dat_struct datastr;
	daymet_tair_struct	tairstr;
	imghdr_struct maskhdr, demhdr;

	char *good;
	char *listgood;
	
	file ini;
	file stnmeta_f, stndata_f;
	file mask_f, dem_f;
	file out_f;

	char endkey[8];
	char round[32];
	char outprefix[80];
	int *dem_array;
	short count;

	int start_day, stop_day, nsimdays;
	int ndays,nstns,nregr;
	int i,j,k,l;
	int day;
	int x,y;
	int type;
	int nrows, ncols;
	int width, hw;
	int stn_filex, stn_filey;
	int test;
	int inflag, inval, *in, n_in;
	int smwidth;
	int ngoodpts;
	int drop,min1,done,first;
	int job;
	int nstns_nws, nstns_snot;
	int *snot;
	int ngoodlr,goodlr;
	
	long off;
	long int obs_off, obs_offset;
	long img_off, img_offset, stn_offset, koff, koffset;

	double pred, obs;
	double pred_avg, obs_avg;
	double rn;
	double dfirad;
	double ul_mapx, ul_mapy;
	double cellsize;
	double fzwidth;
	double *kernel;
	double stn_mapx, stn_mapy;
	double f_x, f_y;
	double value,sum_wt;
	double init_gsp;
	double err;
	double old_err, new_err, old_gsp, new_gsp,opt_err,opt_gsp;
	double delta;
	
	double b, ae, mae1, bias1;
	double mae2, bias2;
	double *ngooddays, *ngoodstns;
	double *xlrdaysum, *xlrdaysumsq;
	double *ylrdaysum, *ylrdaysumsq;
	double *zlrdaysum, *zlrdaysumsq;
	double xlrdaymean,ylrdaymean,zlrdaymean;
	double ngood;
	double nobs,n;
	double tamean, tastdv;
	double xlrsumsq, xlrmean, xlrstdv;
	double ylrsumsq, ylrmean,ylrstdv;
	double zlrsumsq, zlrmean,zlrstdv;
	double daymean, daystdv, daysumsq;
	double tasumsq;
	double rad90sumsq, rad90mean, rad90stdv;
	double frad,rad; 
	double nlr;
	double dz,sumdz; 
	double t1;

	size_t lenp;
        stnmetahdr_struct meta_nc_hdr;
        imghdr_struct mask_nc_hdr, dem_nc_hdr;
        long int *mask_array;
        file mask_nc_f, dem_nc_f, stnmeta_nc_f;
        file stndata_nc_f;
        int status, ncid;
        int nstnid, descid, yrdid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid, stnzid, goodid;
        int stnfxid, stnfyid, stnfzid;
        int stntypeid, stnuc1id, stnuc2id;
	int in_tileid;
        int tairid;
        int ncid_mask, colsid_mask, rowsid_mask, imageid_mask;
        int ncid_dem, colsid_dem, rowsid_dem, imageid_dem;
        int ncid_meta;
	int *mflag, *code;
	int mflagid, codeid;
        int year_day, descriptor_length;
	char *stnnames, *stnids;
	int *in_tile;
        double *elevs, *lats, *lons, *tair;
	double *stnx, *stny, *stnfx, *stnfy;
        int fill_flag = 1;
	/* new variables for xval netcdf output */
	int ncid_out;
	int outid_stnname, outid_stnid, outid_desc;
	int outid_tileid, outid_nstns, outid_nstns3x3, outid_nstndays, outid_rad90mean, outid_rad90stdv;
	int outid_daymae, outid_pormae, outid_bias;
	int outid_tamean;
	int outid_xlrmean, outid_xlrstdv, outid_ylrmean, outid_ylrstdv, outid_zlrmean, outid_zlrstdv;
	int nstndays;
	char att[128];
	int *out_good;
	double *out_obs, *out_pred;
	size_t len_stns, len_days;
	int outdid_stns, outdid_days;
	int outid_obs, outid_pred, outid_good;
	int out_did[2];
	int out_did1[1];
	int nstns_intile;
	int stndims[2];
	char *out_stnnames, *out_stnids;
	double *out_stnx, *out_stny, *out_stnz;
	int outid_stnx, outid_stny, outid_stnz;
	double miss = -9999.0;
	/* new variables for exact_ans bubble sort */
	short *bsi;
	short temp_short;
	int ans;
	double wtsum;
	
	
	/* check for appropriate number of command line arguments */
	if (argc != 14)
	{
		printf("Usage: \n");
		printf("<exec> <job#> <mask> <dem> <meta> <data> <zfilter> <dfirad> <init_gsp> <ans> <smooth> <inflag> <inval> <outpre>\n");
		exit(1);
	}

	/* scan command line arguments into parameters */
	job = atoi(argv[1]);
        strcpy(mask_nc_f.name, argv[2]);
	strcpy(dem_nc_f.name, argv[3]);
	strcpy(stnmeta_nc_f.name, argv[4]);
	strcpy(stndata_nc_f.name, argv[5]);
	fzwidth = atof(argv[6]);
	dfirad = atof(argv[7]);
	init_gsp = atof(argv[8]);
	intstr.ans = atof(argv[9]);
	smwidth = atoi(argv[10]);
	inflag = atoi(argv[11]);
	inval = atoi(argv[12]);
	strcpy(outprefix,argv[13]);
	ans=(int)intstr.ans;
	
	/* create netcdf xval output file */
	strcpy(out_f.name, outprefix);
	strcat(out_f.name,"_xval.nc");
	if (status = nc_create(out_f.name, 0, &ncid_out))
	{
		printf("Error creating %s as netcdf file, exiting\n", out_f.name);
		exit(1);
	}
	/* define scalar variables for xval output file */
	if (status = nc_def_var(ncid_out, "tileid", NC_INT, 0, 0, &outid_tileid))
	{
		printf("Error creating tileid variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"Daymet tile ID");
	nc_put_att_text(ncid_out, outid_tileid, "longname", strlen(att), att);
	strcpy(att,"index");
	nc_put_att_text(ncid_out, outid_tileid, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "nstns3x3", NC_INT, 0, 0, &outid_nstns3x3))
	{
		printf("Error creating nstns3x3 variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"number of stations in 3x3 tile region");
	nc_put_att_text(ncid_out, outid_nstns3x3, "longname", strlen(att), att);
	strcpy(att,"stations");
	nc_put_att_text(ncid_out, outid_nstns3x3, "units", strlen(att), att);
	
	if (status = nc_def_var(ncid_out, "nstns", NC_INT, 0, 0, &outid_nstns))
	{
		printf("Error creating nstns variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"number of stations evaluated (central tile)");
	nc_put_att_text(ncid_out, outid_nstns, "longname", strlen(att), att);
	strcpy(att,"stations");
	nc_put_att_text(ncid_out, outid_nstns, "units", strlen(att), att);
	
	if (status = nc_def_var(ncid_out, "nstndays", NC_INT, 0, 0, &outid_nstndays))
	{
		printf("Error creating nstndays variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"number of station-days evaluated");
	nc_put_att_text(ncid_out, outid_nstndays, "longname", strlen(att), att);
	strcpy(att,"days");
	nc_put_att_text(ncid_out, outid_nstndays, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "rad90mean", NC_DOUBLE, 0, 0, &outid_rad90mean))
	{
		printf("Error creating rad90mean variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"mean: radius capturing 90% of filter kernel weight");
	nc_put_att_text(ncid_out, outid_rad90mean, "longname", strlen(att), att);
	strcpy(att,"m");
	nc_put_att_text(ncid_out, outid_rad90mean, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "rad90stdv", NC_DOUBLE, 0, 0, &outid_rad90stdv))
	{
		printf("Error creating rad90stdv variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"std dev: radius capturing 90% of filter kernel weight");
	nc_put_att_text(ncid_out, outid_rad90stdv, "longname", strlen(att), att);
	strcpy(att,"m");
	nc_put_att_text(ncid_out, outid_rad90stdv, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "daymae", NC_DOUBLE, 0, 0, &outid_daymae))
	{
		printf("Error creating daymae variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"mean absolute error for single-day predictions");
	nc_put_att_text(ncid_out, outid_daymae, "longname", strlen(att), att);
	strcpy(att,"degrees C");
	nc_put_att_text(ncid_out, outid_daymae, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "pormae", NC_DOUBLE, 0, 0, &outid_pormae))
	{
		printf("Error creating annmae variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"mean absolute error for period-of-record predictions");
	nc_put_att_text(ncid_out, outid_pormae, "longname", strlen(att), att);
	strcpy(att,"degrees C");
	nc_put_att_text(ncid_out, outid_pormae, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "bias", NC_DOUBLE, 0, 0, &outid_bias))
	{
		printf("Error creating bias variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"mean prediction bias");
	nc_put_att_text(ncid_out, outid_bias, "longname", strlen(att), att);
	strcpy(att,"degrees C");
	nc_put_att_text(ncid_out, outid_bias, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "tamean", NC_DOUBLE, 0, 0, &outid_tamean))
	{
		printf("Error creating ppmean variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"mean observed temperature");
	nc_put_att_text(ncid_out, outid_tamean, "longname", strlen(att), att);
	strcpy(att,"degrees C");
	nc_put_att_text(ncid_out, outid_tamean, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "xlrmean", NC_DOUBLE, 0, 0, &outid_xlrmean))
	{
		printf("Error creating xlrmean variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: mean x-component");
	nc_put_att_text(ncid_out, outid_xlrmean, "longname", strlen(att), att);
	strcpy(att,"degrees C/m");
	nc_put_att_text(ncid_out, outid_xlrmean, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "xlrstdv", NC_DOUBLE, 0, 0, &outid_xlrstdv))
	{
		printf("Error creating xlrstdv variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: among-station std dev of x-component");
	nc_put_att_text(ncid_out, outid_xlrstdv, "longname", strlen(att), att);
	strcpy(att,"degrees C/m");
	nc_put_att_text(ncid_out, outid_xlrstdv, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "ylrmean", NC_DOUBLE, 0, 0, &outid_ylrmean))
	{
		printf("Error creating ylrmean variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: mean y-component");
	nc_put_att_text(ncid_out, outid_ylrmean, "longname", strlen(att), att);
	strcpy(att,"degrees C/m");
	nc_put_att_text(ncid_out, outid_ylrmean, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "ylrstdv", NC_DOUBLE, 0, 0, &outid_ylrstdv))
	{
		printf("Error creating ylrstdv variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: among-station std dev of y-component");
	nc_put_att_text(ncid_out, outid_ylrstdv, "longname", strlen(att), att);
	strcpy(att,"degrees C/m");
	nc_put_att_text(ncid_out, outid_ylrstdv, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "zlrmean", NC_DOUBLE, 0, 0, &outid_zlrmean))
	{
		printf("Error creating zlrmean variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: mean z-component");
	nc_put_att_text(ncid_out, outid_zlrmean, "longname", strlen(att), att);
	strcpy(att,"derees C/m");
	nc_put_att_text(ncid_out, outid_zlrmean, "units", strlen(att), att);

	if (status = nc_def_var(ncid_out, "zlrstdv", NC_DOUBLE, 0, 0, &outid_zlrstdv))
	{
		printf("Error creating zlrstdv variable for %s, exiting\n", out_f.name);
		exit(1);
	}
	strcpy(att,"3-d regression: among-station std dev of z-component");
	nc_put_att_text(ncid_out, outid_zlrstdv, "longname", strlen(att), att);
	strcpy(att,"degrees C/m");
	nc_put_att_text(ncid_out, outid_zlrstdv, "units", strlen(att), att);
	
	
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
        if (strncmp(mask_nc_hdr.tag,IMGTAG,7))
        {
                printf("Error: mask file not a registered image. Use register_image first.\n");
                exit(1);
        }

        status = nc_inq_varid (ncid_mask, "image", &imageid_mask);

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
        if (strncmp(dem_nc_hdr.tag,IMGTAG,7))
        {
                printf("dem file not a registered image. Use register_image first.\n");
                exit(1);
        }

        status = nc_inq_varid (ncid_dem, "image", &imageid_dem);

        nc_close(ncid_dem);
/*                                                                                                 */
/* -----------------------  end netcdf calls for DEM file --------------------------------------- */

/* -----------------------  netcdf calls for station meta file --------------------------------------- */
/*                                                                                                     */
        if( status = nc_open(stnmeta_nc_f.name, NC_WRITE, &ncid_meta) )
        {
                printf("Error opening %s for netcdf read, exiting\n", stnmeta_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_meta, "stations", &nstnid);
        status = nc_inq_dimlen(ncid_meta, nstnid, &lenp);
        meta_nc_hdr.nstns = lenp;
        status = nc_inq_dimid (ncid_meta, "year_day", &yrdid);
        status = nc_inq_dimlen(ncid_meta, yrdid, &lenp);
        year_day = lenp;
        status = nc_inq_dimid (ncid_meta, "descriptor_length", &descid);
        status = nc_inq_dimlen(ncid_meta, descid, &lenp);
        descriptor_length = lenp;

        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);

        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        if (!(stnnames = (char*) malloc(meta_nc_hdr.nstns * descriptor_length * sizeof(char)))) 
	{
		printf("Error allocating stnnames array.\n");
		exit(1);
	}
        status = nc_get_var_text (ncid_meta, stnnameid, stnnames);

        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);
        if (!(stnids = (char*) malloc(meta_nc_hdr.nstns * descriptor_length * sizeof(char))))
	{
		printf("Error allocating stnids array.\n");
		exit(1);
	}
        status = nc_get_var_text (ncid_meta, stnidid, stnids);

        status = nc_inq_varid (ncid_meta, "in_tile", &in_tileid);
        if (!(in_tile = (int*) malloc(meta_nc_hdr.nstns * sizeof(int))))
	{
                printf("Error allocating in_tile array.\n");
		exit(1);
	}
        status = nc_get_var_int (ncid_meta, in_tileid, in_tile);

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating elevs array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating lats array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating lons array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnlonid, lons);

        status = nc_inq_varid (ncid_meta, "stnx", &stnxid);
        if (!(stnx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating stnx array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnxid, stnx);

        status = nc_inq_varid (ncid_meta, "stny", &stnyid);
        if (!(stny = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating stny array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnyid, stny);

        status = nc_inq_varid (ncid_meta, "stnfx", &stnfxid);
        if (!(stnfx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating stnfx array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnfxid, stnfx);

        status = nc_inq_varid (ncid_meta, "stnfy", &stnfyid);
        if (!(stnfy = (double*) malloc(meta_nc_hdr.nstns * sizeof(double))))
	{
                printf("Error allocating stnfy array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, stnfyid, stnfy);

        status = nc_inq_varid (ncid_meta, "mflag", &mflagid);
        if (!(mflag = (int*) malloc(meta_nc_hdr.nstns * year_day * sizeof(int))))
	{
                printf("Error allocating mflag array.\n");
		exit(1);
	}
        status = nc_get_var_int (ncid_meta, mflagid, mflag);

        status = nc_inq_varid (ncid_meta, "code", &codeid);
        if (!(code = (int*) malloc(meta_nc_hdr.nstns * year_day * sizeof(int))))
	{
                printf("Error allocating code array.\n");
		exit(1);
	}
        status = nc_get_var_int (ncid_meta, codeid, code);

        status = nc_inq_varid (ncid_meta, "tmax", &tairid);
        if (status)
        {
                status = nc_inq_varid (ncid_meta, "tmin", &tairid);
        }
        if (status)
        {
                status = nc_inq_varid (ncid_meta, "tair", &tairid);
        }
        if (!(tair = (double*) malloc(meta_nc_hdr.nstns * year_day * sizeof(double))))
	{
                printf("Error allocating tair array.\n");
		exit(1);
	}
        status = nc_get_var_double (ncid_meta, tairid, tair);


        status = nc_get_att_int(ncid_meta, NC_GLOBAL, "xyloc_flag", &meta_nc_hdr.xyloc_flag);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */

	/* if z-filtering width other than 0.0, open the special ini file
	containing filter type and mask and DEM filenames */
	if (fzwidth)
	{
		strcpy(ini.name, "zfilter.ini");
		if (file_open(&ini, 'i'))
		{
			printf("Error opening %s as special zfilter ini file\n",ini.name);
			exit(1);
		}
		if (scan_value(ini, &type, 'i'))
		{
			printf("Error scanning type\n");
			exit(1);
		}
		if (scan_open(ini, &mask_f, 'r'))
		{
			printf("Error opening mask file\n");
			exit(1);
		}
		if (scan_open(ini, &dem_f, 'r'))
		{
			printf("Error opening dem file\n");
			exit(1);
		}
		fclose(ini.ptr);
	}
		
	/* perform error checking on input file headers */
	if (fzwidth)
	{
		/* read mask file header and verify */
		if (strncmp(mask_nc_hdr.tag,IMGTAG,6))
		{
			printf("Error: %s not a registered image. Use register_image first.\n",mask_f.name);
			exit(1);
		}
		
		/* read dem file header and verify */
		if (strncmp(dem_nc_hdr.tag,IMGTAG,6))
		{
			printf("Error: %s not a registered image. Use register_image first.\n",dem_f.name);
			exit(1);
		}
		
		/* compare headers for DEM and mask */
		test = (mask_nc_hdr.ncols     == dem_nc_hdr.ncols   ) * (mask_nc_hdr.nrows == dem_nc_hdr.nrows) *
			(mask_nc_hdr.ulx      == dem_nc_hdr.ulx     ) * (mask_nc_hdr.uly   == dem_nc_hdr.uly  ) *
			(mask_nc_hdr.cellsize == dem_nc_hdr.cellsize);
		if (!test)
		{
			printf("Error: Header information differs between mask and DEM...exiting\n");
			exit(1);
		}
	}
	
	/* set internal variables */
	if (fzwidth)
	{
		ncols = mask_nc_hdr.ncols;
		nrows = mask_nc_hdr.nrows;
		cellsize = dem_nc_hdr.cellsize;
		ul_mapx = dem_nc_hdr.ulx;
		ul_mapy = dem_nc_hdr.uly;
	}
	nstns = tairstr.nstns = intstr.nstns = meta_nc_hdr.nstns;
	ndays = tairstr.ndays = year_day;
	tairstr.noz = 0;
	/* the next two lines hardwire the xval period to the number of days
	in the data files */
	start_day = 0;
	nsimdays = ndays;

	if (fzwidth)
	{
		/* allocate zfilter memory */
		if (!(dem_array = (int*) malloc(nrows*ncols*sizeof(int))))
		{
                	printf("Error allocating dem array.\n");
			exit(1);
		}
		if (!(mask_array = (long int*) malloc(nrows*ncols*sizeof(long int))))
		{
                	printf("Error allocating mask array.\n");
			exit(1);
		}
	
		/* read and close dem and mask files */
		fread(dem_array,sizeof(int),nrows*ncols,dem_f.ptr);
		fread(mask_array,sizeof(long int),nrows*ncols,mask_f.ptr);
		fclose(dem_f.ptr);
		fclose(mask_f.ptr);
		
		/* set size of elevation filter window */
		width = (int) (fzwidth / cellsize);
		if (!(width % 2)) width++;
		hw = width/2;

		/* generate filter kernel */
		if (!(kernel = (double*) malloc(width * width * sizeof(double))))
		{
			printf("Error allocating for kernel... exiting\n");
			exit(1);
		}
		if (make_kernel(kernel, width, type))
		{
			printf("Error in call to make_kernel() ... exiting\n");
			exit(1);
		}
	}
	
	/* allocate memory */
	if (!(good = (char*) malloc(nstns * ndays * sizeof(char))))
	{
		printf("Error allocating good array... exiting\n");
		exit(1);
	}
	if (!(tairstr.stnobs = (double*) malloc(nstns * ndays * sizeof(double))))
	{
		printf("Error allocating stnsobs array... exiting\n");
		exit(1);
	}
	if (!(tairstr.stnsmobs = (double*) malloc(nstns * ndays * sizeof(double))))
	{
		printf("Error allocating stnsmobs array... exiting\n");
		exit(1);
	}
	if (!(tairstr.stnx = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating stnx array... exiting\n");
		exit(1);
	}
	if (!(tairstr.stny = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating stny array... exiting\n");
		exit(1);
	}
	if (!(tairstr.stnz = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating stnz array... exiting\n");
		exit(1);
	}
	if (!(intstr.stnx = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating stnx array... exiting\n");
		exit(1);
	}
	if (!(intstr.stny = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating stny array... exiting\n");
		exit(1);
	}
	if (!(intstr.sqdist = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating sqdist array... exiting\n");
		exit(1);
	}
	if (!(intstr.wt = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating wt array... exiting\n");
		exit(1);
	}
	if (!(intstr.id = (short*) malloc(nstns * sizeof(short))))
	{
		printf("Error allocating id array... exiting\n");
		exit(1);
	}
	/* index array needed for bubble sort algorithm */
	if (!(bsi = (short*) malloc(nstns * sizeof(short)))) {printf("Couild not malloc bsi\n"); exit(1);}
	/* istr and tstr share space for wt and id arrays */
	tairstr.wt = intstr.wt;
	tairstr.id = intstr.id;
	if (inflag)
	{
		if (!(in = (int*) malloc(nstns * sizeof(int))))
		{
			printf("Error allocating in array... exiting\n");
			exit(1);
		}
	}
	if (!(snot = (int*) malloc(nstns * sizeof(int))))
	{
		printf("Error allocating snot array... exiting\n");
		exit(1);
	}
	
	/* allocate memory for output arrays */
	if (!(ngooddays = (double*) malloc(nstns * sizeof(double))))
	{
		printf("Error allocating ngooddays array... exiting\n");
		exit(1);
	}
	if (!(ngoodstns = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating ngoodstns array... exiting\n");
		exit(1);
	}
	if (!(xlrdaysum = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating xlrdaysum array... exiting\n");
		exit(1);
	}
	if (!(xlrdaysumsq = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating xlrdaysumsq array... exiting\n");
		exit(1);
	}
	if (!(ylrdaysum = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating ylrdaysum array... exiting\n");
		exit(1);
	}
	if (!(ylrdaysumsq = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating ylrdaysumsq array... exiting\n");
		exit(1);
	}
	if (!(zlrdaysum = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating zlrdaysum array... exiting\n");
		exit(1);
	}
	if (!(zlrdaysumsq = (double*) malloc(ndays * sizeof(double))))
	{
		printf("Error allocating zlrdaysumsq array... exiting\n");
		exit(1);
	}
	
	/* read stnmeta data and fill appropriate arrays */
	nstns_nws = nstns_snot = 0;
	nstns_intile = 0;
	for (i=0 ; i<nstns ; i++)
	{
		/* station type: 0 = NWS  1 = SNOTEL */
		snot[i] = 0;
		if (snot[i])
		{
			nstns_snot++;
		}
		else
		{
			nstns_nws++;
		} 
		if (inflag)
		{
			in[i] = in_tile[i];
			if (in[i] == inval) nstns_intile++;
		}
		
		if (meta_nc_hdr.xyloc_flag)
		{
			stn_mapx = intstr.stnx[i] = metastr.fx = stnfx[i];
			stn_mapy = intstr.stny[i] = metastr.fy = stnfy[i];
			tairstr.stnx[i] = metastr.fx = stnfx[i];;
			tairstr.stny[i] = metastr.fy = stnfy[i];;
		}
		else
		{
			stn_mapx = intstr.stnx[i] = metastr.x = stnx[i];
			stn_mapy = intstr.stny[i] = metastr.y = stny[i];
			tairstr.stnx[i] = metastr.x = stnx[i];
			tairstr.stny[i] = metastr.y = stny[i];
		}

		if (fzwidth)
		{	
			/**********************/
			/* zfiltering routine */
			/**********************/
			
			/* calculate initial file (x,y) coordinates for station	*/
			f_x = (stn_mapx - ul_mapx)/cellsize;
			f_y = (ul_mapy - stn_mapy)/cellsize;
			sprintf(round,"%.0lf",f_x);
			stn_filex = atoi(round);
			sprintf(round,"%.0lf",f_y);
			stn_filey = atoi(round);
			stn_offset = stn_filey*ncols + stn_filex;

			/* test if station is within mask, and within DEM data region */
			/* if not, use original elevation , continue */
			if ((stn_filex < 0) || (stn_filex >= ncols) ||
				(stn_filey < 0) || (stn_filey >= nrows) ||
				(!mask_array[stn_offset]))
			{
				tairstr.stnz[i] = metastr.z = elevs[i];
				continue;
			}

			/* station is in mask, and inside DEM data region */
			/* begin loops through kernel to calculate filtered elevation */
			value = 0.0;
			sum_wt = 0.0;
			for (k=0 ; k<width ; k++)
			{
				/* check for bounds of image */
				y = stn_filey + k - hw;
				if ((y < 0) || (y >= nrows))
					continue;

				koff = k * width;
				img_off = y * ncols;

				for (l=0 ; l<width ; l++)
				{
					/* check for bounds of image */
					x = stn_filex + l - hw;
					if ((x < 0) || (x >= ncols))
						continue;

					koffset = koff + l;
					img_offset = img_off + x;

					/* finally, check mask for convolution overlap */
					if (!mask_array[img_offset])
						continue;

					/* accumulate convolved values */
					value += dem_array[img_offset] * kernel[koffset];
					sum_wt += kernel[koffset];
				}
			}

			/* check that at least some good values were found */
			if (!sum_wt)
			{
				tairstr.stnz[i] = metastr.z = elevs[i];
				continue;
			}
			tairstr.stnz[i] = value/sum_wt;
		}
		else /* fzwidth == 0 */
		{
			tairstr.stnz[i] = metastr.z = elevs[i];
		}
		
	} /* end of nstns loop */
		
	if (nstns_intile)
	{
		/* allocate space for xval ouput */
		if (!(out_stnnames = (char*)malloc(nstns_intile * descriptor_length * sizeof(char))))
		{
			printf("Error allocating for out_stnnames array.\n");
			exit(1);
		}
		if (!(out_stnids = (char*)malloc(nstns_intile * descriptor_length * sizeof(char))))
		{
			printf("Error allocating for out_stnids array.\n");
			exit(1);
		}
		if (!(out_stnx = (double*)malloc(nstns_intile * sizeof(double))))
		{
			printf("Error allocating for out_stnx array... exiting\n");
			exit(1);
		}
		if (!(out_stny = (double*)malloc(nstns_intile * sizeof(double))))
		{
			printf("Error allocating for out_stny array... exiting\n");
			exit(1);
		}
		if (!(out_stnz = (double*)malloc(nstns_intile * sizeof(double))))
		{
			printf("Error allocating for out_stnz array... exiting\n");
			exit(1);
		}
		if (!(out_obs = (double*)malloc(nstns_intile*ndays* sizeof(double))))
		{
			printf("Error allocating for out_obs array... exiting\n");
			exit(1);
		} 
		if (!(out_pred = (double*)malloc(nstns_intile*ndays* sizeof(double))))
		{
			printf("Error allocating for out_pred array... exiting\n");
			exit(1);
		} 
		if (!(out_good = (int*)malloc(nstns_intile*ndays* sizeof(int))))
		{
			printf("Error allocating for out_good array... exiting\n");
			exit(1);
		} 

		/* netcdf define dimensions and array variables for output arrays */
		len_stns = nstns_intile;
		len_days = ndays;
		nc_def_dim(ncid_out, "stns", len_stns, &outdid_stns);
        	nc_def_dim(ncid_out, "descriptor_length", descriptor_length, &outid_desc);
		nc_def_dim(ncid_out, "days", len_days, &outdid_days);

		stndims[0] = outdid_stns;
		stndims[1] = outid_desc;
		out_did[0] = outdid_stns;
		out_did[1] = outdid_days;
		out_did1[0] = outdid_stns;
		
        	nc_def_var (ncid_out, "station_name",      NC_CHAR,   2, stndims, &outid_stnname);
        	nc_def_var (ncid_out, "station_id",        NC_CHAR,   2, stndims, &outid_stnid);
		
		nc_def_var(ncid_out, "stnx", NC_DOUBLE, 1, out_did1, &outid_stnx);
		strcpy(att,"station projected x coordinate");
		nc_put_att_text(ncid_out, outid_stnx, "longname", strlen(att), att);
		strcpy(att,"m");
		nc_put_att_text(ncid_out, outid_stnx, "units", strlen(att), att);

		nc_def_var(ncid_out, "stny", NC_DOUBLE, 1, out_did1, &outid_stny);
		strcpy(att,"station projected y coordinate");
		nc_put_att_text(ncid_out, outid_stny, "longname", strlen(att), att);
		strcpy(att,"m");
		nc_put_att_text(ncid_out, outid_stny, "units", strlen(att), att);

		nc_def_var(ncid_out, "stnz", NC_DOUBLE, 1, out_did1, &outid_stnz);
		strcpy(att,"station elevation");
		nc_put_att_text(ncid_out, outid_stnz, "longname", strlen(att), att);
		strcpy(att,"m");
		nc_put_att_text(ncid_out, outid_stnz, "units", strlen(att), att);

		nc_def_var(ncid_out, "obs", NC_DOUBLE, 2, out_did, &outid_obs);
		strcpy(att,"observed temperature");
		nc_put_att_text(ncid_out, outid_obs, "longname", strlen(att), att);
		strcpy(att,"degrees C");
		nc_put_att_text(ncid_out, outid_obs, "units", strlen(att), att);

		nc_def_var(ncid_out, "pred", NC_DOUBLE, 2, out_did, &outid_pred);
		strcpy(att,"predicted temperature");
		nc_put_att_text(ncid_out, outid_pred, "longname", strlen(att), att);
		strcpy(att,"degrees C");
		nc_put_att_text(ncid_out, outid_pred, "units", strlen(att), att);

		nc_def_var(ncid_out, "good", NC_INT, 2, out_did, &outid_good);
		strcpy(att,"good value flag (0=missing obs, 1=good obs)");
		nc_put_att_text(ncid_out, outid_good, "longname", strlen(att), att);
		strcpy(att,"flag");
		nc_put_att_text(ncid_out, outid_good, "units", strlen(att), att);
	} /* end if nstns_intile */
	
	/* take output file out of define mode, ready for write */
	nc_enddef(ncid_out);
	
	/* read station daily data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			tairstr.stnobs[obs_offset] = tair[obs_offset];
			if (mflag[obs_offset] == 1) good[obs_offset] = 0;
                        else good[obs_offset] = 1;
		}
	}
	/* calculate mean and spatial standard deviation */
	tamean = tastdv = tasumsq = nobs = 0.0;
	for (j=0 ; j<ndays ; j++)
	{
		daymean = daystdv = daysumsq = n = 0.0;
		for (i=0 ; i<nstns ; i++)
		{
			/* test for inclusion/exclusion */
			if (inflag)
			{
				if (in[i] != inval) continue;
			}
			obs_offset = i*ndays + j;
			if (good[obs_offset])
			{
				n++;
				nobs++;
				value = tairstr.stnobs[obs_offset];
				tamean += value;
				daymean += value;
				daysumsq += value*value;
			}
		}
		/* calculate this day's standard deviation */
		daymean /= n;
		daystdv = sqrt((daysumsq/n)-(daymean*daymean));
		tastdv += daystdv;
	}
	tamean /= nobs;
	tastdv /= (double)ndays;

	/* smooth the station data arrays by station
	For cross-validation code, this smoothing algorithm also requires an
	array of characters indicating good and missing data (1=good, 0=miss) */
	for (i=0 ; i<nstns ; i++)
	{
		off = i * ndays;
		if (xv_boxcar_smooth(good+off,tairstr.stnobs+off,tairstr.stnsmobs+off,
			ndays,smwidth,1))
		{
			printf("Error in call to xv_boxcar_smooth()... exiting\n");
			exit(1);
		}
	}
	
	/* assign interpolation gsp */
	intstr.gsp = init_gsp;
	
	/* initialize the mean for rad90 and lr, find frad */
	rad90mean = rad90stdv = rad90sumsq = 0.0;
	for (day=0 ; day<nsimdays ; day++)
	{
		xlrdaysum[day] = 0.0;
		xlrdaysumsq[day] = 0.0;
		ylrdaysum[day] = 0.0;
		ylrdaysumsq[day] = 0.0;
		zlrdaysum[day] = 0.0;
		zlrdaysumsq[day] = 0.0;
		ngoodstns[day] = 0.0;
	}
	i = (int)((intstr.gsp * 2.0) - 1.0);
	/* force to the range 0:49 */
	if (i > 49) i = 49;
	frad = f90[i];

	/* calculate initial density filter parameters */
	intstr.dfsqrad = dfirad * dfirad;
	intstr.inv_dfsqrad = 1.0 / intstr.dfsqrad;
	intstr.dfarea = PI * intstr.dfsqrad;
	intstr.trunc = exp(-intstr.gsp);
	intstr.dftrunc = exp(-DFGSP);
	intstr.dfavgwt = ((1.0 - intstr.dftrunc)/DFGSP) - intstr.dftrunc;

	/* start loop through stations */
	mae1 = bias1 = 0.0;
	mae2 = bias2 = 0.0;
	tairstr.switch_init = 1;
	n_in = 0;
	ngood = 0.0;
	nlr = 0.0;
	ngoodlr = 0;
	for (i=0 ; i<nstns ; i++)
	{
		/* test for inclusion/exclusion */
		if (inflag)
		{
			if (in[i] != inval) continue;
		}

		/* copy station name and station ID from input arrays */
		for (j=0 ; j<descriptor_length ; j++)
		{
			out_stnnames[n_in*descriptor_length + j]=stnnames[i*descriptor_length + j];
			out_stnids[n_in*descriptor_length + j]=stnids[i*descriptor_length + j];
		}
		out_stnx[n_in]=tairstr.stnx[i];
		out_stny[n_in]=tairstr.stny[i];
		out_stnz[n_in]=tairstr.stnz[i];
		n_in ++;

		/* assign x,y,z for this station as the prediction point */
		intstr.x = intstr.stnx[i];
		intstr.y = intstr.stny[i];
		tairstr.mapx = tairstr.stnx[i];
		tairstr.mapy = tairstr.stny[i];
		tairstr.elev = tairstr.stnz[i];

		/* calculate the squared distance from this point to each station */
		if (calc_sqdist(&intstr))
		{
			printf("Error in call to calc_sqdist()... exiting\n");
			exit(1);
		}
#ifdef EXACT_ANS
		/* use a bubble-sort algorithm to get exactly ANS stations in the prediction
		list for each point. Since we don't want the prediction point included
		for xval, find nearest ANS+1, then reject the prediction point station
		at the end. */
		for (j=0 ; j<nstns ; j++) bsi[j]=j;
		for (j=0 ; j<ans+1 ; j++) {
		  for (k=nstns-1 ; k>j ; k--) {
		    if (intstr.sqdist[bsi[k]] < intstr.sqdist[bsi[k-1]]) {
		      temp_short = bsi[k];
		      bsi[k]=bsi[k-1];
		      bsi[k-1]=temp_short;
		    }
		  }
		}
		
		/* set the weighting kernel radius to 1m greater than the furthest included station */
		/* the extra 1.0m ensures that none of the stations in the list will have wt=0 */
		intstr.sqrad=intstr.sqdist[bsi[ans]]+1.0;
		
		/* fill the list of included stations, being sure to leave out the prediction point
		station */
		k=0;
		wtsum=0.0;
		for (j=0 ; j<ans+1 ; j++) {
		  if (bsi[j] == i) continue;
		  intstr.id[k]=bsi[j];
		  intstr.wt[k]=exp(-intstr.gsp * intstr.sqdist[bsi[j]] * (1.0/intstr.sqrad)) - intstr.trunc;
		  wtsum += intstr.wt[k];
		  k++;
		}
		/* at this point k should equal ANS - exit with error if not true */
		if (k != ans)
		{
			printf("Error in call to cmd_predict_daily_tair_xval3_nc, exiting\n");
			printf("Failed test that k == ANS after removing the prediction point.\n");
			printf("Part of EXACT_ANS code.\n");
			printf("Failed for station # %d\n",i);
			exit(1);
		}
		/* normalize weights */
		if (wtsum) {
		  for (j=0 ; j<ans ; j++) intstr.wt[j] /= wtsum;
		} else {
		  printf("Error in call to cmd_predict_daily_tair_xval3_nc, exiting.\n");
		  printf("wtsum = 0 after bubble sort.\n");
		  exit(1);
		}
		intstr.count = (short)ans;
#else

		/* iterative density algorithm to calculate squared search radius */
		if (calc_search_rad(&intstr))
		{
			printf("Error in call to calc_search_rad()... exiting\n");
			exit(1);
		}
		/* generate cross-validation weighted station list for this point */
		if (xv_weight_list(&intstr, i))
		{
			printf("Error in call to xv_weight_list()... exiting\n");
			exit(1);
		}
#endif

		/* based on search radius, calculate the f90 radius at this gsp */
		rad = sqrt(intstr.sqrad) * frad;
		rad90mean += rad;
		rad90sumsq += rad*rad;

		/* assign station count and n_regression_points */
		tairstr.count = count = intstr.count;
		tairstr.nregr = nregr = ((count * count) - count) / 2;

		/* allocate memory depending on count */
		if (!(listgood = (char*) malloc(count * ndays * sizeof(char))))
		{
			printf("Error allocating listgood array... exiting\n");
			exit(1);   
		}
		if (!(tairstr.listx = (double*) malloc(count * sizeof(double))))
		{
			printf("Error allocating listx array... exiting\n");
			exit(1);
		}
		if (!(tairstr.listy = (double*) malloc(count * sizeof(double))))
		{
			printf("Error allocating listy array... exiting\n");
			exit(1);
		}
		if (!(tairstr.listz = (double*) malloc(count * sizeof(double))))
		{
			printf("Error allocating listz array... exiting\n");
			exit(1);
		}
		if (!(tairstr.listobs = (double*) malloc(count*ndays*sizeof(double))))
		{
			printf("Error allocating listobs array... exiting\n");
			exit(1);
		}
		if (!(tairstr.listsmobs = (double*) malloc(count*ndays*sizeof(double))))
		{
			printf("Error allocating listmobs array... exiting\n");
			exit(1);
		}
		if (!(tairstr.regx = (double**) malloc(3 * sizeof(double*))))
		{
			printf("Error allocating regx array... exiting\n");
			exit(1);
		}
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.regx[j] = (double*)malloc(nregr*sizeof(double))))
			{
				printf("Error allocating regx array... exiting\n");
				exit(1);
			}
		}
		if (!(tairstr.regdx = (double**) malloc(3 * sizeof(double*))))
		{
			printf("Error allocating regdx array... exiting\n");
			exit(1);
		}
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.regdx[j] = (double*)malloc(nregr*sizeof(double))))
			{
				printf("Error allocating regdx array... exiting\n");
				exit(1);
			}
		}
		if (!(tairstr.regy = (double*) malloc(nregr * sizeof(double))))
		{
			printf("Error allocating regy array... exiting\n");
			exit(1);
		}
		if (!(tairstr.regwtall = (double*) malloc(nregr * sizeof(double))))
		{
			printf("Error allocating regwtall array... exiting\n");
			exit(1);
		}
		if (!(tairstr.regwt = (double*) malloc(nregr * sizeof(double))))
		{
			printf("Error allocating regwt array... exiting\n");
			exit(1);
		}
		if (!(tairstr.inv = (double**) malloc(3 * sizeof(double*))))
		{
			printf("Error allocating inv array... exiting\n");
			exit(1);
		}
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.inv[j] = (double*)malloc(3*sizeof(double))))
			{
				printf("Error allocating inv array... exiting\n");
				exit(1);
			}
		}

		/* fill the station-list arrays */
		if (xv_fill_list(good, listgood, &tairstr))
		{
			printf("Error filling station list arrays ... exiting\n");
			exit(1);
		}

		/* for now, using all lapse rates */
		goodlr = 1;

		/* generate regression arrays that don't vary between days */
		if (xv_tair_regr_xwt_2switch(&tairstr))
		{
			printf("Error in xv_tair_regr_xwt_2switch() ... exiting\n");
			exit(1);
		}

		obs_off = i * ndays + start_day;
		ngooddays[i] = 0.0;
		obs_avg = pred_avg = 0.0;

		/* begin loop through nsimdays */
		for (day=0 ; day<nsimdays ; day++)
		{
			/* first check if there is a good observation for the day */
			if (!good[obs_off+day])
			{
				out_obs[(n_in-1)*nsimdays+day]=-9999.0;
				out_pred[(n_in-1)*nsimdays+day]=-9999.0;
				out_good[(n_in-1)*nsimdays+day]=0;
				continue;
			}

			/* initialize daily prediction variables */
			tairstr.ta = 0.0;

			/* generate offset into the station list observations */
			tairstr.dam_off = (day + start_day) * count;

			/* test for enough good data to generate regression */
			/* remove this test, since this comes after call to fill missing values 
			ngoodpts = 0;
			for (j=0 ; j<count ; j++)
			{
				if (listgood[tairstr.dam_off + j])
				{
					ngoodpts ++;
				}
			}

			if (ngoodpts < 3)
			{
				continue;
			}
			*/
			ngoodpts = count;

			/* generate regression arrays that change by day */
			if (xv_tair_regr_y_2switch(listgood, &tairstr))
			{
				printf("Error in xv_tair_regr_y_switch() ... exiting\n");
				exit(1);
			}

			/* weighted least squares slope */
			if (xv_tair_mlr(&tairstr))
			{
				printf("Error in wt_regr() ... exiting\n");
				exit(1);
			}

			/* update mean and stdv for slope (lapse rate) */
			if (goodlr)
			{
				ngoodstns[day]++;
				xlrdaysum[day] += tairstr.coef[0];
				xlrdaysumsq[day] += tairstr.coef[0] * tairstr.coef[0];
				ylrdaysum[day] += tairstr.coef[1];
				ylrdaysumsq[day] += tairstr.coef[1] * tairstr.coef[1];
				zlrdaysum[day] += tairstr.coef[2];
				zlrdaysumsq[day] += tairstr.coef[2] * tairstr.coef[2];
			}

			/* predict tair for the day */
			if (xv_predict_tair_noint(listgood, &tairstr))
			{
				printf("Error in xv_predict_tair_noint() ... exiting\n");
				exit(1);
			}

			/* calculate the error and bias for this station-day */
			pred = tairstr.ta;
			obs = tairstr.stnobs[obs_off+day];

			out_obs[(n_in-1)*nsimdays+day]=obs;
			out_pred[(n_in-1)*nsimdays+day]=pred;
			out_good[(n_in-1)*nsimdays+day]=1;
			
			b = pred - obs;
			ae = fabs(b);

			/* this day has good observation and good predictor data */
			ngooddays[i]++;
			ngood++;

			bias1 += b;
			mae1 += ae;

			pred_avg += pred;
			obs_avg += obs;

		} /* end of nsimdays loop */

		/* free memory depending on count */
		free(listgood);
		free(tairstr.listobs);
		free(tairstr.listsmobs);
		free(tairstr.listx);
		free(tairstr.listy);
		free(tairstr.listz);
		for (j=0 ; j<3 ; j++)
		{
			free(tairstr.regx[j]);
		}
		free(tairstr.regx);
		for (j=0 ; j<3 ; j++)
		{
			free(tairstr.regdx[j]);
		}
		free(tairstr.regdx);
		free(tairstr.regy);
		free(tairstr.regwtall);
		free(tairstr.regwt);
		for (j=0 ; j<3 ; j++)
		{
			free(tairstr.inv[j]);
		}
		free(tairstr.inv);

		/* calculate period-of-record averages */
		if (ngooddays[i])
		{
			rn = ngooddays[i];
			pred_avg /= rn;
			obs_avg /= rn;
			b = pred_avg - obs_avg;
			ae = fabs(b);
			mae2 += rn * ae;
			bias2 += rn * b;
		}

		/* switch initial sign for next station */
		tairstr.switch_init = -tairstr.switch_init;

	} /* end of nstations loop */

	if (nstns_intile)
	{
		/* write netcdf output arrays */
		nc_put_var_text(ncid_out, outid_stnname, out_stnnames);
		nc_put_var_text(ncid_out, outid_stnid, out_stnids);
		nc_put_var_double(ncid_out, outid_stnx, out_stnx);
		nc_put_var_double(ncid_out, outid_stny, out_stny);
		nc_put_var_double(ncid_out, outid_stnz, out_stnz);
		nc_put_var_int(ncid_out, outid_good, out_good);
		nc_put_var_double(ncid_out, outid_obs, out_obs);
		nc_put_var_double(ncid_out, outid_pred, out_pred);
	}

	/* calculate final xval stats */
	xlrstdv = xlrmean = 0.0;
	ylrstdv = ylrmean = 0.0;
	zlrstdv = zlrmean = 0.0;
	if (ngood)
	{
		rn = ngood;
		mae1 /= rn;
		bias1 /= rn;
		mae2 /= rn;
		bias2 /= rn;

		/* calculate means and stdevs for lr */
		for (day=0 ; day<nsimdays ; day++)
		{
			xlrdaymean = xlrdaysum[day]/ngoodstns[day];
			xlrmean += xlrdaymean;
			xlrstdv += sqrt((xlrdaysumsq[day]/ngoodstns[day]) - (xlrdaymean*xlrdaymean));
			ylrdaymean = ylrdaysum[day]/ngoodstns[day];
			ylrmean += ylrdaymean;
			ylrstdv += sqrt((ylrdaysumsq[day]/ngoodstns[day]) - (ylrdaymean*ylrdaymean));
			zlrdaymean = zlrdaysum[day]/ngoodstns[day];
			zlrmean += zlrdaymean;
			t1=(zlrdaysumsq[day]/ngoodstns[day]) - (zlrdaymean*zlrdaymean);
			if (t1 < 0.0) t1=0.0;
			zlrstdv += sqrt(t1);
		}
		xlrmean /= (double)nsimdays;
		xlrstdv /= (double)nsimdays;
		ylrmean /= (double)nsimdays;
		ylrstdv /= (double)nsimdays;
		zlrmean /= (double)nsimdays;
		zlrstdv /= (double)nsimdays;

		/* calculate the final means and standard deviations for rad90 */
		rad90mean /= (double)n_in;
		rad90stdv = sqrt((rad90sumsq/(double)n_in) - (rad90mean*rad90mean));

		/* the optimized error is defined as the average of the
		daily mean absolute error and the period-of-record
		mean absolute error */
		opt_err = (mae1 + mae2) / 2.0;

	}
	else 
	{
		printf("No good station days in tile, can't calculate xval stats for this tile.\n");
		rad90mean = miss;
		rad90stdv = miss;
		mae1 = miss;
		mae2 = miss;
		bias1 = miss;
		tamean = miss;
		xlrmean = miss;
		xlrstdv = miss;
		ylrmean = miss;
		ylrstdv = miss;
		zlrmean = miss;
		zlrstdv = miss;
	}

	/* netcdf write */
	nstndays = (int)ngood;
	nc_put_var1_int(ncid_out, outid_tileid,		0, &job);
	nc_put_var1_int(ncid_out, outid_nstns3x3, 	0, &nstns);
	nc_put_var1_int(ncid_out, outid_nstns, 		0, &nstns_intile);
	nc_put_var1_int(ncid_out, outid_nstndays, 	0, &nstndays);
	nc_put_var1_double(ncid_out, outid_rad90mean,	0, &rad90mean);
	nc_put_var1_double(ncid_out, outid_rad90stdv,	0, &rad90stdv);
	nc_put_var1_double(ncid_out, outid_daymae,	0, &mae1);
	nc_put_var1_double(ncid_out, outid_pormae,	0, &mae2);
	nc_put_var1_double(ncid_out, outid_bias,	0, &bias1);
	nc_put_var1_double(ncid_out, outid_tamean,	0, &tamean);
	nc_put_var1_double(ncid_out, outid_xlrmean,	0, &xlrmean);
	nc_put_var1_double(ncid_out, outid_xlrstdv,	0, &xlrstdv);
	nc_put_var1_double(ncid_out, outid_ylrmean,	0, &ylrmean);
	nc_put_var1_double(ncid_out, outid_ylrstdv,	0, &ylrstdv);
	nc_put_var1_double(ncid_out, outid_zlrmean,	0, &zlrmean);
	nc_put_var1_double(ncid_out, outid_zlrstdv,	0, &zlrstdv);

	/* close netcdf output file */
	nc_close(ncid_out);
		
	/* free memory not dependent on count */
	if (fzwidth)
	{
		free(dem_array);
		free(mask_array);
		free(kernel);
	}
	if (inflag) free(in);
	free(good);
	free(tairstr.stnobs);
	free(tairstr.stnsmobs);
	free(tairstr.stnz);
	free(intstr.stnx);
	free(intstr.stny);
	free(intstr.sqdist);
	free(intstr.wt);
	free(intstr.id);
	free(snot);
	free(ngooddays);
	free(ngoodstns);
	free(xlrdaysum);
	free(xlrdaysumsq);	
	free(ylrdaysum);
	free(ylrdaysumsq);	
	free(zlrdaysum);
	free(zlrdaysumsq);	

	return(0);
	
} /* end of main */
