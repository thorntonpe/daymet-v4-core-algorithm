/*
 cmd_prcp_ndays_filter_nc.c
 Peter Thornton
 10/12/09
 
 Since station-years are not currently being filtered as they are
 extracted from the portal for backend processing, a new routine is
 introduced here that checks the number of missing days for each station-year
 and throws out the stations that don't meet the user-specfied criteria.
 
 */

#include "metsrc2.h"
#include "netcdf.h"

int main(int argc, char *argv[])
{
	/* variable declarations */
 
        file prcp_nc_f;
 
	char *good, *good2;
	char *listgood;

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
	int xyloc_flag;
        int status, ncid;
        int nstnid, descid, yrdid;
	int no_dataid;
	int in_tileid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid, stnzid, goodid;
        int stnfxid, stnfyid, stnfzid;
        int stntypeid, stnuc1id, stnuc2id;
        int prcpid;
	int prcp_name;
        int ncid_meta;
	int year;
        int year_day, descriptor_length;
	char *stnnames, *stnids;
        double *elevs, *lats, *lons, *prcp;
	double *stnx, *stnfx, *stny, *stnfy;
	double *no_data;
	int *in_tile;
	char *stnnames_filt, *stnids_filt;
        double *elevs_filt, *lats_filt, *lons_filt, *prcp_filt;
	double *stnx_filt, *stnfx_filt, *stny_filt, *stnfy_filt;
	double *no_data_filt;
	int *in_tile_filt;
	int *mflag, *mflag_filt;
	int mflagid, mflagdims[2];
	int ii, prcpdims[2], stndims[2], gooddims[2];
	int codeid, codedims[2];

	/* check for command line arguments */
	if (argc != 3)
	{
		printf("Usage: <exec> <meta_f> <cript>\n");
		printf("argc = %d\n",argc);
		exit(1);
	}
	/* get initialization data from command line */
	strcpy(prcp_nc_f.name,argv[1]);
	critp = atof(argv[2]);
/* -----------------------  netcdf calls for station meta file --------------------------------------- */
/*                                                                                                     */
        if( status = nc_open(prcp_nc_f.name, NC_WRITE, &ncid_meta) )
        {
                printf("Error opening %s for netcdf read, exiting\n", prcp_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_meta, "stations", &nstnid);
        status = nc_inq_dimlen(ncid_meta, nstnid, &lenp);
	nstns = lenp;
        status = nc_inq_dimid (ncid_meta, "year_day", &yrdid);
        status = nc_inq_dimlen(ncid_meta, yrdid, &lenp);
        year_day = lenp;
        status = nc_inq_dimid (ncid_meta, "descriptor_length", &descid);
        status = nc_inq_dimlen(ncid_meta, descid, &lenp);
        descriptor_length = lenp;
        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "year", &year);
        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "days_per_year", &year_day);
	ndays = year_day;

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

        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "xyloc_flag", &xyloc_flag);
        if (xyloc_flag == 1)
        {
                status = nc_inq_varid (ncid_meta, "stnfx", &stnfxid);
                if (!(stnfx = (double*) malloc(nstns * sizeof(double)))) exit(1);
                status = nc_get_var_double (ncid_meta, stnfxid, stnfx);

                status = nc_inq_varid (ncid_meta, "stnfy", &stnfyid);
                if (!(stnfy = (double*) malloc(nstns * sizeof(double)))) exit(1);
                status = nc_get_var_double (ncid_meta, stnfyid, stnfy);
        }

        status = nc_inq_varid (ncid_meta, "precip", &prcpid);
        if (!(prcp = (double*) malloc(nstns * year_day * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, prcpid, prcp);

	nc_close(ncid_meta);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */
	/* allocate memory */
	if (!(good  = (char*) malloc(nstns * ndays * sizeof(char)))) exit(1);
	if (!(mflag =  (int*) malloc(nstns * ndays * sizeof(int)))) exit(1);
	if (!(keep  =  (int*) malloc(nstns * sizeof(int)))) exit(1);
	/* read station daily data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			if (prcp[obs_offset] == no_data[i]) {
				good[obs_offset] = 0;
				mflag[obs_offset] = 1;
			}
			else
			{
				good[obs_offset] = 1;
				mflag[obs_offset] = 0;
			}
		}
	}
	/* loop through stations and cast out all those which exceed the 
	missing days criteria  */
	nkeep = 0;
	for (i=0 ; i<nstns ; i++)
	{
		nfiltgood = 0;
		keep[i] = 0;
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			if (good[obs_offset]) nfiltgood++;
		}
		if ((double)nfiltgood/(double)ndays >= critp)
		{
			/* this station passes the filtering criteria, and will be
			retained in the filtered stnmeta and stndata files */
			keep[i] = 1;
			nkeep++;
		}
	}
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
        if (xyloc_flag)
	{
		if (!(stnfx_filt    = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        	if (!(stnfy_filt    = (double*) malloc(nkeep * sizeof(double)))) exit(1);
        }
	if (!(prcp_filt     = (double*) malloc(nkeep * year_day * sizeof(double)))) exit(1);
        if (!(mflag_filt    =    (int*) malloc(nkeep * year_day * sizeof(int)))) exit(1);

	/* read through the meta and daily data, and write out the station records that 
	are being kept to the new meta and data output files */
	ii = 0;
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
			if (xyloc_flag)
			{
				stnfx_filt[ii]   = stnfx[i];
				stnfy_filt[ii]   = stnfy[i];
			}
			obs_off = i*ndays;
			for (j=0 ; j<ndays ; j++)
			{
				obs_offset = obs_off + j;
				prcp_filt[ii*ndays+j] = prcp[obs_offset];
				mflag_filt[ii*ndays+j] = mflag[obs_offset];
			}
			ii++;
		}
	}

/* -------------------  netcdf calls to write new station meta file ---------------------------------- */
/*                                                                                                     */
        if( status = nc_create (prcp_nc_f.name, 0, &ncid_meta) )
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
        if (xyloc_flag)
	{
		status = nc_def_var (ncid_meta, "stnfx",             NC_DOUBLE, 1, &nstnid, &stnfxid);
        	status = nc_def_var (ncid_meta, "stnfy",             NC_DOUBLE, 1, &nstnid, &stnfyid);
	}
	codedims[0] = nstnid;
	codedims[1] = yrdid;
        status = nc_def_var (ncid_meta, "mflag",             NC_INT,    2, codedims, &mflagid);

	prcpdims[0] = nstnid;
	prcpdims[1] = yrdid;
	status = nc_def_var (ncid_meta, "precip", NC_DOUBLE, 2, prcpdims, &prcpid);

        status = nc_put_att_int(ncid_meta,  NC_GLOBAL, "xyloc_flag", NC_INT, 1, &xyloc_flag);
        status = nc_put_att_int(ncid_meta,  NC_GLOBAL, "year",       NC_INT, 1, &year);
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
	if (xyloc_flag)
	{
        	status = nc_put_var_double (ncid_meta, stnfxid, stnfx_filt);
        	status = nc_put_var_double (ncid_meta, stnfyid, stnfy_filt);
	}
        status = nc_put_var_int    (ncid_meta, mflagid, mflag_filt);
        status = nc_put_var_double (ncid_meta, prcpid, prcp_filt);

	nc_close(ncid_meta);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */	

	/* free memory */
	free(stnnames);
	free(stnids);
	free(elevs);
	free(lats);
	free(lons);
	free(no_data);
	free(in_tile);
	free(stnx);
	free(stny);
	if (xyloc_flag == 1)
	{
		free(stnfx);
		free(stnfy);
	}
	free(prcp);
	free(good);
	free(mflag);
	free(keep);
		
	free(stnnames_filt);
	free(stnids_filt);
	free(elevs_filt);
	free(lats_filt);
	free(lons_filt);
	free(no_data_filt);
	free(in_tile_filt);
	free(stnx_filt);
	free(stny_filt);
	if (xyloc_flag == 1)
	{
		free(stnfx_filt);
		free(stnfy_filt);
	}
	free(prcp_filt);
	free(mflag_filt);
	
	return(0);
} /* end of main */
	
