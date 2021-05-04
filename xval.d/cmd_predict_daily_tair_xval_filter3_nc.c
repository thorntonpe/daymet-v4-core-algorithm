/*
cmd_predict_daily_tair_xval_filter3.c
Peter Thornton
2/13/97

A revision of original xval_filter, but uses a commend line instead of
a *.ini file

Generates the single-set cross-validation error statistics for air temperature
predictions using the interpolate and predict_daily_tair algorithms.

Based on predict_daily_tair_xval1, but now includes a second loop through
all stations in which days from the previous loop with absolute errors 
greater than 10 C are now set to missing values.  Output is appended
to the end of the original output files.

Updated:
9/6/00, PET: from cmd_predict_daily_tair_xval_filter.c, added support for
  3-d regressions
  
Sept 05 CSH: added netcdf code
*/

#include "metsrc2.h"
#include "netcdf.h"

int main(int argc, char *argv[])
{
	/* variable declarations */
	interpolate_struct intstr;
	stnmetahdr_struct metahdr, metahdrjunk;
	stnmeta_struct metastr;
	daydathdr_struct datahdr, datahdrjunk;
	dat_struct datastr;
	daymet_tair_struct	tairstr;

	char *good, *good2;
	char *listgood;
	
	file ini;
	file stnmeta_f, stndata_f;
	file filtmeta_f, filtdata_f;

	char endkey[8];
	char round[32];
	char outprefix[80];

	short count;

	int ndays,nstns,nregr;
	int i,j,k,l;
	int day;
	int x,y;
	int type;
	int nrows, ncols;
	int width, hw;
	int stn_filex, stn_filey;
	int test;
	int smwidth;
	int ngoodpts;
	int *keep;
	int nkeep,nfiltgood;
	int cmd;
	
	long off;
	long int obs_off, obs_offset;
	long stn_offset;

	double pred, obs;
	double rn;
	double dfirad;
	double ul_mapx, ul_mapy;
	double cellsize;
	double stn_mapx, stn_mapy;
	double f_x, f_y;
	double value,sum_wt;
	
	double b, ae;
	double critp;

	size_t lenp;
        stnmetahdr_struct meta_nc_hdr;
        file stnmeta_nc_f;
        file stndata_nc_f;
        int status, ncid;
        int nstnid, descid, yrdid;
	int no_dataid;
	int in_tileid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid, stnzid, goodid;
        int stnfxid, stnfyid, stnfzid;
        int stntypeid, stnuc1id, stnuc2id;
        int tairid;
	int tair_name;
        int ncid_meta;
	char ayear[4];
	int year;
        int year_day, descriptor_length;
	char *stnnames, *stnids;
        double *elevs, *lats, *lons, *tair;
	double *stnx, *stnfx, *stny, *stnfy;
	double *no_data;
	int *in_tile;
	char *stnnames_filt, *stnids_filt;
        double *elevs_filt, *lats_filt, *lons_filt, *tair_filt;
	double *stnx_filt, *stnfx_filt, *stny_filt, *stnfy_filt;
	double *no_data_filt;
	int *in_tile_filt;
	int *mflag, *mflag_filt;
	int mflagid, mflagdims[2];
	int *code, *code_filt;
	int codeid, codedims[2];
	int ii, tairdims[2], stndims[2], gooddims[2];
	int error_check;

	error_check=1;

	/* check for command line arguments */
	if (argc != 9)
	{
		printf("Usage: <exec> <meta_f> <data_f> <cript> <dfirad> <gsp> <ans> <smwidth> <outpre>\n");
		exit(1);
	}
	/* get initialization data from command line */
	strcpy(stnmeta_nc_f.name,argv[1]);
	strcpy(stndata_nc_f.name,argv[2]);
	critp = atof(argv[3]);
	dfirad = atof(argv[4]);
	intstr.gsp = atof(argv[5]);
	intstr.ans = atof(argv[6]);
	smwidth = atoi(argv[7]);
	strcpy(outprefix,argv[8]);
	

/* -----------------------  netcdf calls for station meta file --------------------------------------- */
/*                                                                                                     */
        if( status = nc_open(stnmeta_nc_f.name, NC_SHARE, &ncid_meta) )
        {
                printf("Error opening %s for netcdf read, exiting\n", stnmeta_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_meta, "stations", &nstnid);
        status = nc_inq_dimlen(ncid_meta, nstnid, &lenp);
        meta_nc_hdr.nstns = lenp;
	nstns = meta_nc_hdr.nstns;
        status = nc_inq_dimid (ncid_meta, "year_day", &yrdid);
        status = nc_inq_dimlen(ncid_meta, yrdid, &lenp);
        year_day = lenp;
        status = nc_inq_dimid (ncid_meta, "descriptor_length", &descid);
        status = nc_inq_dimlen(ncid_meta, descid, &lenp);
        descriptor_length = lenp;

        status = nc_get_att_text (ncid_meta, NC_GLOBAL, "year", ayear);
        year = atoi(ayear);
        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "days_per_year", &year_day);

        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        if (!(stnnames = (char*) malloc(nstns * descriptor_length * sizeof(char)))) exit(1);
        status = nc_get_var_text (ncid_meta, stnnameid, stnnames);

        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);
        if (!(stnids = (char*) malloc(nstns * descriptor_length * sizeof(char)))) exit(1);
        status = nc_get_var_text (ncid_meta, stnidid, stnids);

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlonid, lons);

        status = nc_inq_varid (ncid_meta, "no_data", &no_dataid);
        if (!(no_data = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, no_dataid, no_data);

        status = nc_inq_varid (ncid_meta, "in_tile", &in_tileid);
        if (!(in_tile = (int*) malloc(nstns * sizeof(int)))) exit(1);
        status = nc_get_var_int (ncid_meta, in_tileid, in_tile);

        status = nc_inq_varid (ncid_meta, "stnx", &stnxid);
        if (!(stnx = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnxid, stnx);

        status = nc_inq_varid (ncid_meta, "stny", &stnyid);
        if (!(stny = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnyid, stny);

        status = nc_inq_varid (ncid_meta, "stnfx", &stnfxid);
        if (!(stnfx = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnfxid, stnfx);

        status = nc_inq_varid (ncid_meta, "stnfy", &stnfyid);
        if (!(stnfy = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnfyid, stnfy);

	tair_name=0;
        status = nc_inq_varid (ncid_meta, "tmax", &tairid);
        if (status)
        {
		tair_name=1;
                status = nc_inq_varid (ncid_meta, "tmin", &tairid);
        }
        if (!(tair = (double*) malloc(nstns * year_day * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, tairid, tair);

        status = nc_get_att_int(ncid_meta, NC_GLOBAL, "xyloc_flag", &meta_nc_hdr.xyloc_flag);
	
	printf("xyloc_flag = %d\n",meta_nc_hdr.xyloc_flag);
	
	nc_close(ncid_meta);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */

	/* set internal variables */
	nstns = tairstr.nstns = intstr.nstns = meta_nc_hdr.nstns;
	ndays = tairstr.ndays = year_day;
	tairstr.noz = 0;

    //BWM
    //printf("nstns: %d ndays: %d\n", nstns, ndays);
	/* allocate memory */
	if (!(good  = (char*) malloc(nstns * ndays * sizeof(char)))) exit(1);
	if (!(good2 = (char*) malloc(nstns * ndays * sizeof(char)))) exit(1);
	if (!(mflag =  (int*) malloc(nstns * ndays * sizeof(int)))) exit(1);
	if (!(code  =  (int*) malloc(nstns * ndays * sizeof(int)))) exit(1);
	if (!(keep  =  (int*) malloc(nstns * sizeof(int)))) exit(1);
	if (!(tairstr.stnobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(tairstr.stnsmobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(tairstr.stnx = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(tairstr.stny = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(tairstr.stnz = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(intstr.stnx = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(intstr.stny = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(intstr.sqdist = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(intstr.wt = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(intstr.id = (short*) malloc(nstns * sizeof(short)))) exit(1);
	/* istr and tstr share space for wt and id arrays */
	tairstr.wt = intstr.wt;
	tairstr.id = intstr.id;
	
	/* read stnmeta data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		if (meta_nc_hdr.xyloc_flag)
		{
			stn_mapx = intstr.stnx[i] = stnfx[i];
			stn_mapy = intstr.stny[i] = stnfy[i];
			tairstr.stnx[i] = stnfx[i];
			tairstr.stny[i] = stnfy[i];
		}
		else
		{
			stn_mapx = intstr.stnx[i] = stnx[i];
			stn_mapy = intstr.stny[i] = stny[i];
			tairstr.stnx[i] = stnx[i];
			tairstr.stny[i] = stny[i];
		}

		tairstr.stnz[i] = elevs[i];
	} /* end of nstns loop */

	/* read station daily data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			tairstr.stnobs[obs_offset] = tair[obs_offset];
			if (tair[obs_offset] == no_data[i]) {
				good[obs_offset] = 0;
				mflag[obs_offset] = 1;
			}
			else {
				good[obs_offset] = 1;
				mflag[obs_offset] = 0;
			}
			code[obs_offset] = 0;
		}
	}
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

	/* calculate initial density filter parameters */
	intstr.dfsqrad = dfirad * dfirad;
	intstr.inv_dfsqrad = 1.0 / intstr.dfsqrad;
	intstr.dfarea = PI * intstr.dfsqrad;
	intstr.trunc = exp(-intstr.gsp);
	intstr.dftrunc = exp(-DFGSP);
	intstr.dfavgwt = ((1.0 - intstr.dftrunc)/DFGSP) - intstr.dftrunc;
	
	/* start loop through stations */
	tairstr.switch_init = 1;
	for (i=0 ; i<nstns ; i++)
	{
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
		
		/* assign station count to pstr (id and wt already pointed) */
		tairstr.count = count = intstr.count;
		tairstr.nregr = nregr = ((count * count) - count) / 2;

		//BWM printf("station #%d count=%d ndays: %d\n",i,count, ndays);
	
		/* allocate memory depending on count */
		if (!(listgood = (char*) malloc(count * ndays * sizeof(char)))) exit(1);   
		if (!(tairstr.listx = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(tairstr.listy = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(tairstr.listz = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(tairstr.listobs = (double*) malloc(count*ndays*sizeof(double)))) exit(1);
		if (!(tairstr.listsmobs = (double*) malloc(count*ndays*sizeof(double)))) exit(1);
		if (!(tairstr.regx = (double**) malloc(3 * sizeof(double*)))) exit(1);
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.regx[j] = (double*)malloc(nregr*sizeof(double)))) exit(1);
		}
		if (!(tairstr.regy = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(tairstr.regwtall = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(tairstr.regwt = (double*) malloc(nregr * sizeof(double)))) exit(1);

		/* fill the station-list arrays */
		if (xv_fill_list(good, listgood, &tairstr))
		{
			printf("Error filling station list arrays ... exiting\n");
			exit(1);
		}

		/* generate regression arrays that don't vary between days */
		if (xv_tair_regr_xwt_2switch(&tairstr))
		{
			printf("Error in xv_tair_regr_xwt_2switch() ... exiting\n");
			exit(1);
		}

		obs_off = i * ndays;

		/* begin loop through ndays */
		for (day=0 ; day<ndays ; day++)
		{
			/* first check if there is a good observation for the day */
			good2[obs_off + day] = 1;
			if (!good[obs_off+day])
			{
				good2[obs_off + day] = 0;
				continue;
			}

			/* initialize daily prediction variables */
			tairstr.ta = 0.0;
			
			/* generate offset into the station list observations */
			tairstr.dam_off = day * count;
			
			/* test for enough good data to generate regression */
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
				printf("less than 3 good stations for the day %d of %d\n", day, ndays);
				continue;
			}
			
			/* generate regression arrays that change by day */
			if (xv_tair_regr_y_2switch(listgood, &tairstr))
			{
				printf("Error in xv_tair_regr_y_switch() ... exiting\n");
				exit(1);
			}

			/* weighted least squares slope and intercept */
			if (xv_tair_mlr(&tairstr))
			{
				printf("Error in wt_regr() ... exiting\n");
				exit(1);
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
			
			b = pred - obs;
			ae = fabs(b);
			
			/* set good2 array for filtering */
			if (ae >= 10.0) good2[obs_off + day] = 0;

		} /* end of ndays loop */
			
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
		free(tairstr.regy);
		free(tairstr.regwtall);
		free(tairstr.regwt);
		
		/* switch initial sign for next station */
		tairstr.switch_init = -tairstr.switch_init;
		
	} /* end of nstations loop */

	/* loop through stations and cast out all those which exceed the 
	missing days criteria after the xval filtering */
	nkeep = 0;
	for (i=0 ; i<nstns ; i++)
	{
		nfiltgood = 0;
		keep[i] = 0;
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			if (good2[obs_offset]) nfiltgood++;
		}
		if ((double)nfiltgood/(double)ndays >= critp)
		{
			/* this station passes the filtering criteria, and will be
			retained in the filtered stnmeta and stndata files */
			keep[i] = 1;
			nkeep++;
		}
	}
	printf("nkeep = %d\n",nkeep);
	error_check=2;
		
	/* allocate output memory */
        if (!(stnnames_filt =   (char*) malloc(nkeep * descriptor_length * sizeof(char)))) exit(1);
        if (!(stnids_filt   =   (char*) malloc(nkeep * descriptor_length * sizeof(char)))) exit(1);
        if (!(elevs_filt    = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(lats_filt     = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(lons_filt     = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(no_data_filt  = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(in_tile_filt  =    (int*) malloc(nkeep * sizeof(int)))) exit(1);
        if (!(stnx_filt     = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(stny_filt     = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(stnfx_filt    = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(stnfy_filt    = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        if (!(tair_filt     = (double*) malloc(nkeep * year_day * sizeof(double)))) exit(1);
        if (!(mflag_filt    =    (int*) malloc(nkeep * year_day * sizeof(int)))) exit(1);
        if (!(code_filt     =    (int*) malloc(nkeep * year_day * sizeof(int)))) exit(1);

	ii = 0;

	/* read through the meta and daily data, and write out the station records that 
	are being kept to the new meta and data output files */
	for (i=0 ; i<nstns ; i++)
	{
		if (keep[i])
		{
			for (j=0 ; j<descriptor_length ; j++)
			{
				stnnames_filt[ii*descriptor_length+j] = stnnames[i*descriptor_length+j];
				  stnids_filt[ii*descriptor_length+j] =   stnids[i*descriptor_length+j];
			}
			elevs_filt[ii]   = elevs[i];
			lats_filt[ii]    = lats[i];
			lons_filt[ii]    = lons[i];
			no_data_filt[ii] = no_data[i];
			in_tile_filt[ii] = in_tile[i];
			stnx_filt[ii]    = stnx[i];
			stny_filt[ii]    = stny[i];
			stnfx_filt[ii]   = stnfx[i];
			stnfy_filt[ii]   = stnfy[i];
			obs_off = i*ndays;
			for (j=0 ; j<ndays ; j++)
			{
				obs_offset = obs_off + j;
				tair_filt[ii*ndays+j] = tair[obs_offset];
				mflag_filt[ii*ndays+j] = mflag[obs_offset];
				code_filt[ii*ndays+j] = 0;  
				if (!mflag[obs_offset] && !good2[obs_offset])
				{
					mflag_filt[ii*ndays+j] = 1;
					code_filt[ii*ndays+j] = 9;  /* filtered value */
				}
			}
			ii++;
		}
	}
	
/* -------------------  netcdf calls to write new station meta file ---------------------------------- */
/*                                                                                                     */
        if( status = nc_create (stnmeta_nc_f.name, NC_SHARE, &ncid_meta) )
        {
                printf("Error opening temp file for netcdf output, exiting\n");
                exit(1);
        }
        status = nc_def_dim (ncid_meta, "stations",          nkeep,             &nstnid);
        status = nc_def_dim (ncid_meta, "descriptor_length", descriptor_length, &descid);
        status = nc_def_dim (ncid_meta, "year_day",          year_day,           &yrdid);

	stndims[0] = nstnid;
	stndims[1] = descid;
        status = nc_def_var (ncid_meta, "station_name",      NC_CHAR,   2, stndims, &stnnameid);
        status = nc_def_var (ncid_meta, "station_id",        NC_CHAR,   2, stndims, &stnidid);
        status = nc_def_var (ncid_meta, "station_elevation", NC_DOUBLE, 1, &nstnid, &stnelevid);
        status = nc_def_var (ncid_meta, "station_latitude",  NC_DOUBLE, 1, &nstnid, &stnlatid);
        status = nc_def_var (ncid_meta, "station_longitude", NC_DOUBLE, 1, &nstnid, &stnlonid);
        status = nc_def_var (ncid_meta, "no_data",           NC_DOUBLE, 1, &nstnid, &no_dataid);
        status = nc_def_var (ncid_meta, "in_tile",           NC_INT,    1, &nstnid, &in_tileid);
        status = nc_def_var (ncid_meta, "stnx",              NC_DOUBLE, 1, &nstnid, &stnxid);
        status = nc_def_var (ncid_meta, "stny",              NC_DOUBLE, 1, &nstnid, &stnyid);
        status = nc_def_var (ncid_meta, "stnfx",             NC_DOUBLE, 1, &nstnid, &stnfxid);
        status = nc_def_var (ncid_meta, "stnfy",             NC_DOUBLE, 1, &nstnid, &stnfyid);
	codedims[0] = nstnid;
	codedims[1] = yrdid;
        status = nc_def_var (ncid_meta, "mflag",             NC_INT,    2, codedims, &mflagid);
        status = nc_def_var (ncid_meta, "code",              NC_INT,    2, codedims, &codeid);

	tairdims[0] = nstnid;
	tairdims[1] = yrdid;
	if(tair_name==0) status = nc_def_var (ncid_meta, "tmax", NC_DOUBLE, 2, tairdims, &tairid);
	if(tair_name==1) status = nc_def_var (ncid_meta, "tmin", NC_DOUBLE, 2, tairdims, &tairid);

        status = nc_put_att_int(ncid_meta,  NC_GLOBAL, "xyloc_flag", NC_INT, 1, &meta_nc_hdr.xyloc_flag);
        status = nc_put_att_text(ncid_meta, NC_GLOBAL, "year", 4, ayear);
        status = nc_put_att_int(ncid_meta,  NC_GLOBAL, "days_per_year", NC_INT, 1, &year_day);

	status = nc_enddef (ncid_meta);

        status = nc_put_var_text   (ncid_meta, stnnameid, stnnames_filt);
        status = nc_put_var_text   (ncid_meta, stnidid, stnids_filt);
        status = nc_put_var_double (ncid_meta, stnelevid, elevs_filt);
        status = nc_put_var_double (ncid_meta, stnlatid, lats_filt);
        status = nc_put_var_double (ncid_meta, stnlonid, lons_filt);
        status = nc_put_var_double (ncid_meta, no_dataid, no_data_filt);
        status = nc_put_var_int    (ncid_meta, in_tileid, in_tile_filt);
        status = nc_put_var_double (ncid_meta, stnxid, stnx_filt);
        status = nc_put_var_double (ncid_meta, stnyid, stny_filt);
        status = nc_put_var_double (ncid_meta, stnfxid, stnfx_filt);
        status = nc_put_var_double (ncid_meta, stnfyid, stnfy_filt);
        status = nc_put_var_int    (ncid_meta, mflagid, mflag_filt);
        status = nc_put_var_int    (ncid_meta, codeid, code_filt);
        status = nc_put_var_double (ncid_meta, tairid, tair_filt);
	nc_sync(ncid_meta);
	nc_close(ncid_meta);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */	

	/* free memory not depending on count */
	free(stnnames);
	free(stnids);
	free(elevs);
	free(lats);
	free(lons);
	free(no_data);
	free(in_tile);
	free(stnx);
	free(stny);
	free(stnfx);
	free(stnfy);
	free(tair);
	free(good);
	free(good2);
	free(mflag);
	free(code);
	free(keep);
	free(tairstr.stnobs);
	free(tairstr.stnsmobs);
	free(tairstr.stnx);
	free(tairstr.stny);
	free(tairstr.stnz);
	free(intstr.stnx);
	free(intstr.stny);
	free(intstr.sqdist);
	free(intstr.wt);
	free(intstr.id);
	
	return(0);
	
} /* end of main */
