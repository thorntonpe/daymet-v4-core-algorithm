/* 
predict_daily_tair3_nc.c
Peter Thronton

Revision history:
July 2005:  netcdf version - CSH
6/27/95:  Revision of mtclim3d_tair.c.
7/9/96 :  Revision of daymet_tair.c (/scratch5/daymet.d/orig.d)
7/10/96:  Revision of daymet_tair_test.c added general data structure for
function calls
7/12/96:  Revision of daymet_tair_test2.c clean-up for final version
11/11/96: Revision of daymet_tair.c, to include new metsrc structures.
1/28/97:  Changed error checking routine for list input to look for a file with
		  5 elements (x,y,z,slp,asp) instead of the older version which 
		  had 3 (x,y,z)
4/2/97: Added an option for passing all required parameters on the 
		command-line
4/2/97: altered treatment of start_day and stop_day: if both are set to 0,
 then all days are used
8/15/00: Adapted from predict_daily_tair.c. Changed tavg output so that only
 good_value mask region gets a value, which means that these images need to be
 filtered through a mask for display. Also added the good_value parameter for
 US grid. Removed LIST option.
9/1/00: Adapted from predict_daily_tair2.c, adding three-variable multiple
 linear regression model.
*/

#include "metsrc2.h"
#include "netcdf.h"

int main(int argc, char* argv[])
{
	/* variable declarations */
	stnmetahdr_struct metahdr;
	stnmeta_struct metastr;
	daymet_tair_struct tairstr;
	imghdr_struct maskhdr, demhdr;
	lsthdr_struct listhdr;
	dayouthdr_struct outhdr;

	file ini;
	file mask_f;
	file dem_f;
	file list_f;
	file ct_f, id_f, wt_f;
	file t_out_f;
	file t_tavg_f;
	file t_savg_f;
	file tnoz_out_f;
	file tnoz_tavg_f;
	file tnoz_savg_f;
	file mx_tavg_f;
	file mx_savg_f;
	file my_tavg_f;
	file my_savg_f;
	file mz_tavg_f;
	file mz_savg_f;

	//char endkey[8];
	//char iocase[16];
	char round[16];
	char ioc,map;
	char outprefix[80];
	double *t_ts;
	double *tnoz_ts;
	unsigned char brange,bmin;
	long int mask;

	short count;
	int z;

	int start_day, stop_day, nsimdays;
	int ncols, nrows, ncells;
	int ndays,nstns,nregr;
	int i,j;
	int day;
	int n_in;
	int bytmax, bytmin;
	int noz;
	int smwidth;
	int test;
	long int good_value;

	long int off;
	long int obs_off, obs_offset;

	int imapx,imapy;
	double ulx,uly,cellsize;
	double t, tnoz;
	double t_tavg, tnoz_tavg;
	double mx_tavg;
	double my_tavg;
	double mz_tavg;
	double *t_savg_ts, *tnoz_savg_ts;
	double *mx_savg_ts;
	double *my_savg_ts;
	double *mz_savg_ts;
	double valmax,valmin,scale_ratio;
	double bytrange;
	double rn;
	//double celldat[5];

	//char userin[80];

	size_t lenp;
	int xyloc_flag;
	imghdr_struct mask_nc_hdr, dem_nc_hdr;
        stnmetahdr_struct meta_nc_hdr;
	dat_struct data_nc_str;
	file mask_nc_f, dem_nc_f, stnmeta_nc_f, stndata_nc_f;
	long int *mask_array;
	int *dem_array;
	int ncid_mask, colsid_mask, rowsid_mask, imageid_mask;
	int ncid_dem, colsid_dem, rowsid_dem, imageid_dem;
      	int nstnid, descid, yrdid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid, stnzid;
        int stnfxid, stnfyid, stnfzid;
        int stntypeid, stnuc1id, stnuc2id;
        int tairid;
	int ncid_meta;
	char ayear[5];
        int year, year_days, descriptor_length;
        double *elevs, *lats, *lons, *stnx, *stny, *tair;
	double *stnfx, *stnfy;
	int status;
	size_t index[2];
	
	/* assuming command-line parameter passing */
	map = 1;
	strcpy(mask_nc_f.name,argv[1]);
	good_value = atoi(argv[2]);
	strcpy(dem_nc_f.name,argv[3]);
	strcpy(stnmeta_nc_f.name,argv[4]);
	strcpy(stndata_nc_f.name,argv[5]);
	strcpy(ct_f.name,argv[6]);
	strcpy(id_f.name,argv[7]);
	strcpy(wt_f.name,argv[8]);
	smwidth = atoi(argv[9]);
	noz = atoi(argv[10]);
	start_day = atoi(argv[11]);
	stop_day = atoi(argv[12]);
	valmax = atof(argv[13]);
	valmin = atof(argv[14]);
	bytmax = atoi(argv[15]);
	bytmin = atoi(argv[16]);
	strcpy(outprefix,argv[17]);

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
        if (strncmp(dem_nc_hdr.tag,IMGTAG,7))
        {
                printf("dem file not a registered image. Use register_image first.\n");
                printf("dem tag: %s\n",dem_nc_hdr.tag);
                exit(1);
        }

        status = nc_inq_varid (ncid_dem, "image", &imageid_dem);
/*                                                                                                 */
/* -----------------------  end netcdf calls for DEM file --------------------------------------- */

/* -----------------------  netcdf calls for station meta file --------------------------------------- */
/*                                                                                                     */
        if( status = nc_open(stnmeta_nc_f.name, 0, &ncid_meta) )
        {
                printf("Error opening %s for netcdf read, exiting\n", stnmeta_nc_f.name);
                exit(1);
        }
        status = nc_inq_dimid (ncid_meta, "stations", &nstnid);
        status = nc_inq_dimlen(ncid_meta, nstnid, &lenp);
        meta_nc_hdr.nstns = lenp;
        status = nc_inq_dimid (ncid_meta, "descriptor_length", &descid);
        status = nc_inq_dimlen(ncid_meta, descid, &lenp);
        descriptor_length = lenp;
        status = nc_inq_dimid (ncid_meta, "year_day", &yrdid);
        status = nc_inq_dimlen(ncid_meta, yrdid, &lenp);
        year_days = lenp;

        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "year", &year);
        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlonid, lons);

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "stnx", &stnxid);
        if (!(stnx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnxid, stnx);

        status = nc_inq_varid (ncid_meta, "stny", &stnyid);
        if (!(stny = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnyid, stny);

        status = nc_get_att_int (ncid_meta, NC_GLOBAL, "xyloc_flag", &xyloc_flag);
	if (xyloc_flag == 1)
        {       
                status = nc_inq_varid (ncid_meta, "stnfx", &stnfxid);
                if (!(stnfx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
                status = nc_get_var_double (ncid_meta, stnfxid, stnfx);

                status = nc_inq_varid (ncid_meta, "stnfy", &stnfyid);
                if (!(stnfy = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
                status = nc_get_var_double (ncid_meta, stnfyid, stnfy);
        }

        status = nc_inq_varid (ncid_meta, "tmax", &tairid);
        if (status)
        {
                status = nc_inq_varid (ncid_meta, "tmin", &tairid);
        }
        if (!(tair = (double*) malloc(meta_nc_hdr.nstns * year_days * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, tairid, tair);

        nc_close(ncid_meta);

	meta_nc_hdr.zfilt_flag = 0;
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */

	/* open the input data files */
	if (file_open(&ct_f,'r'))
	{
		printf("Error opening %s for binary read\n",ct_f.name);
		exit(1);
	}
	if (file_open(&id_f,'r'))
	{
		printf("Error opening %s for binary read\n",id_f.name);
		exit(1);
	}
	if (file_open(&wt_f,'r'))
	{
		printf("Error opening %s for binary read\n",wt_f.name);
		exit(1);
	}
		
	/* error checking on file headers */
	if (map)
	{
		/* read mask file header and verify */
		if (strncmp(mask_nc_hdr.tag,IMGTAG,7))
		{
			printf("Mask not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* read dem file header and verify */
		if (strncmp(dem_nc_hdr.tag,IMGTAG,7))
		{
			printf("DEM not a registered image. Use register_image first.\n");
			exit(1);
		}

		/* compare headers for DEM and mask */
		test = (mask_nc_hdr.ncols == dem_nc_hdr.ncols)*(mask_nc_hdr.nrows == dem_nc_hdr.nrows)*
			(mask_nc_hdr.ulx == dem_nc_hdr.ulx)*(mask_nc_hdr.uly == dem_nc_hdr.uly)*
			(mask_nc_hdr.cellsize == dem_nc_hdr.cellsize);
	}
	else
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

	/* set internal variables */
	if (map)
	{
		ncols = mask_nc_hdr.ncols;
		nrows = mask_nc_hdr.nrows;
		ulx = mask_nc_hdr.ulx;
		uly = mask_nc_hdr.uly;
		cellsize = mask_nc_hdr.cellsize;
		ncells = ncols * nrows;
	}
	else
	{
		ncells = listhdr.ncells;
	}

	nstns = tairstr.nstns = meta_nc_hdr.nstns;
	ndays = tairstr.ndays = year_days;

	tairstr.noz = (char)noz;
	
	/* some other basic error checking */
	if (start_day > stop_day)
	{
		printf("start_day > stop_day ... exiting\n");
		exit(1);
	}
	if (stop_day >= ndays)
	{
		printf("stop_day >= ndays ... exiting\n");
		printf("remember to reference days from base zero (0,1,2...)\n");
		exit(1);
	}
	/* test for special start and stop day combo:
	when both are 0, use all ndays */
	if (start_day == 0 && stop_day == 0)
	{
		nsimdays = ndays;
	}
	else
	{
		nsimdays = (stop_day - start_day) + 1;
	}
	rn = (double)nsimdays;
	
	/* initialize byte scaling parameters */
	// if (init_byte_scale(valmax,valmin,bytmax,bytmin,&scale_ratio,&brange,&bmin))
	// {
	// 	printf("Error in init_byte_scale ... exiting\n");
	// 	exit(1);
	// }
	// bytrange = (double)brange;

	printf("outhdr.start_year: %i year: %i\n", outhdr.start_year, year);	
	/* fill the daily output file header structure */
	/* determine the yearday and year for the start of output */
	if (yearday(0, year, start_day,
		&outhdr.start_yday, &outhdr.start_year)) 
	{
		printf("Error in call to yearday()...exiting\n");
		exit(1);
	}
	strcpy(outhdr.tag, TAIRTAG);
	outhdr.valmax = valmax;
	outhdr.valmin = valmin;
	outhdr.bytmax = bytmax;
	outhdr.bytmin = bytmin;
	outhdr.ndays = ndays;
	outhdr.ioflag = map;
	if (map)
	{
		outhdr.maskhdr = mask_nc_hdr;
	}
	else
	{
		outhdr.listhdr = listhdr;
	}

	/* create output filenames and open files for binary write */
	if (bin_out(&t_out_f, outprefix, ".dmtair_out")) exit(1);
	if (bin_out(&t_tavg_f, outprefix, ".dmtair_tavg")) exit(1);
	if (bin_out(&t_savg_f, outprefix, ".dmtair_savg")) exit(1);
	if (bin_out(&mx_tavg_f, outprefix, ".dmtair_mxtavg")) exit(1);
	if (bin_out(&mx_savg_f, outprefix, ".dmtair_mxsavg")) exit(1);
	if (bin_out(&my_tavg_f, outprefix, ".dmtair_mytavg")) exit(1);
	if (bin_out(&my_savg_f, outprefix, ".dmtair_mysavg")) exit(1);
	if (bin_out(&mz_tavg_f, outprefix, ".dmtair_mztavg")) exit(1);
	if (bin_out(&mz_savg_f, outprefix, ".dmtair_mzsavg")) exit(1);
	if (noz)
	{
		if (bin_out(&tnoz_out_f, outprefix, ".dmtair_nozout")) exit(1);
		if (bin_out(&tnoz_tavg_f, outprefix, ".dmtair_noztavg")) exit(1);
		if (bin_out(&tnoz_savg_f, outprefix, ".dmtair_nozsavg")) exit(1);
	}


	printf("predict.d/predict_daily_tair3_nc.c outhdr.start_year: %s %i\n", outprefix, outhdr.start_year);	
	/* write headers to daily data output files */
	outhdr.noz = 0;
	if (fwrite(&outhdr, sizeof(dayouthdr_struct), 1, t_out_f.ptr) != 1)
	{
		printf("Error writing header to outfile... exiting\n");
		exit(1);
	}
	if (noz)
	{
		outhdr.noz = 1;
		if (fwrite(&outhdr, sizeof(dayouthdr_struct), 1, tnoz_out_f.ptr) != 1)
		{
			printf("Error writing header to noz outfile... exiting\n");
			exit(1);
		}
	}
		
	/* allocate array memory */
	if (!(tairstr.stnobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(tairstr.stnsmobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(tairstr.stnx = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(tairstr.stny = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(tairstr.stnz = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(t_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(t_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(mx_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(my_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(mz_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (noz)
	{
		if (!(tnoz_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
		if (!(tnoz_savg_ts = (double*) malloc(nsimdays * sizeof(double))))	exit(1);
	}
	
	/* initialize time series arrays */
	for (i=0 ; i<nsimdays ; i++)
	{
		t_ts[i] =  0;
		t_savg_ts[i] = mx_savg_ts[i] = my_savg_ts[i] = mz_savg_ts[i] = 0.0;
		if (noz)
		{
			tnoz_ts[i] = 0;
			tnoz_savg_ts[i] = 0.0;
		}
	}
	
	/* read stnmeta data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
/*
		if (meta_nc_hdr.zfilt_flag) tairstr.stnz[i] = metastr.fz;
		else tairstr.stnz[i] = metastr.z;
*/

		if (xyloc_flag) 
		{
			tairstr.stnx[i] = metastr.fx = stnfx[i];
			tairstr.stny[i] = metastr.fy = stnfy[i];
		}
		else
		{
			tairstr.stnx[i] = stnx[i];
			tairstr.stny[i] = stny[i];
		}
		tairstr.stnz[i] = elevs[i];
	}

	/* read station daily data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			tairstr.stnobs[obs_offset] = tair[obs_offset];
		}
	}
	
	/* smooth the station data arrays by station (event smoothing algorithm) */
	for (i=0 ; i<nstns ; i++)
	{
		off = i * ndays;
		if (boxcar_smooth(tairstr.stnobs+off,tairstr.stnsmobs+off,
			ndays,smwidth,1))
		{
			printf("Error in call to boxcar_smooth()... exiting\n");
			exit(1);
		}
	}
	
	/* begin simulation, loop through ncells */
	n_in = 0;
	tairstr.switch_init = 1;
	
	for (i=0 ; i<ncells ; i++)
	{
		/* read mask value if MAP I/O */
		mask = mask_array[i];
		if (mask != good_value)
		{
			/* skip to next cell */
			continue;
		}

		/* this cell is in the simulation region */
		n_in ++;
		imapy = i/ncols;
		imapx = i - (imapy*ncols);
		index[0]=imapy;
		index[1]=imapx;
		status=nc_get_var1_int(ncid_dem, imageid_dem, index, &z);
		tairstr.mapy = uly - (double)imapy*cellsize;
		tairstr.mapx = ulx + (double)imapx*cellsize;
		tairstr.elev = (double)z;
		
		/* initialize temporal variables for this cell */
		t_tavg = 0.0;
		mx_tavg = 0.0;
		my_tavg = 0.0;
		mz_tavg = 0.0;
		
		if (noz)
		{
			tnoz_tavg = 0.0;
		}
		
		/* read station count for this cell */
		fread(&count, sizeof(short), 1, ct_f.ptr);
		tairstr.count = count;

		/* allocate memory depending on count */
		tairstr.nregr = nregr = ((count * count) - count) / 2;
		if (!(tairstr.id = (short*) malloc(count * sizeof(short)))) exit(1);
		if (!(tairstr.wt = (double*) malloc(count * sizeof(double)))) exit(1);
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
		if (!(tairstr.regdx = (double**) malloc(3 * sizeof(double*)))) exit(1);
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.regdx[j] = (double*)malloc(nregr*sizeof(double)))) exit(1);
		}
		if (!(tairstr.regy = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(tairstr.regwt = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(tairstr.inv = (double**) malloc(3 * sizeof(double*)))) exit(1);
		for (j=0 ; j<3 ; j++)
		{
			if (!(tairstr.inv[j] = (double*)malloc(3*sizeof(double)))) exit(1);
		}

		/* read station id and wt arrays */
		fread(tairstr.id, sizeof(short), count, id_f.ptr);
		fread(tairstr.wt, sizeof(double), count, wt_f.ptr);
		
		/* fill the station-list arrays */
		if (fill_list(&tairstr))
		{
			printf("Error filling station list arrays ... exiting\n");
			exit(1);
		}
		
		/* generate regression arrays that don't vary between days */
		if (tair_regr_xwt_2switch(&tairstr))
		{
			printf("Error in tair_regr_xwt_switch() ... exiting\n");
			exit(1);
		}

		/* begin loop through nsimdays */
		for (day=0 ; day<nsimdays ; day++)
		{
			/* generate offset into the station list observations */
			tairstr.dam_off = (day + start_day) * count;
			
			/* generate regression arrays that change by day */
			if (tair_regr_y_2switch(&tairstr))
			{
				printf("Error in tair_regr_y_2switch() ... exiting\n");
				exit(1);
			}

			/* weighted multiple least squares regression to get slopes and intercept */
			if (mlr(&tairstr))
			{
				printf("Error in mlr() ... exiting\n");
				exit(1);
			}
			
			/* predict tair for the day */
			if (predict_tair_noint(&tairstr))
			{
				printf("Error in predict_tair_noint() ... exiting\n");
				exit(1);
			}
			
			/* update spatial and temporal averages */
			t = tairstr.ta;
			tnoz = tairstr.tanoz;
			t_tavg += t;
			t_savg_ts[day] += t;
			mx_tavg += tairstr.coef[0];
			mx_savg_ts[day] += tairstr.coef[0];
			my_tavg += tairstr.coef[1];
			my_savg_ts[day] += tairstr.coef[1];
			mz_tavg += tairstr.coef[2];
			mz_savg_ts[day] += tairstr.coef[2];
			if (noz)
			{
				tnoz_tavg += tnoz;
				tnoz_savg_ts[day] += tnoz;
			}

			/* byte scale daily predictions */
			// if ((t = (t - valmin) * scale_ratio) < 0.0) t = 0.0;
			// if (t > bytrange) t = bytrange;
			// sprintf(round, "%.0lf", t);
			// t_ts[day] = (char)atoi(round) + bmin;
			t_ts[day] = t;
			
			if (noz)
			{
				// if ((tnoz = (tnoz - valmin) * scale_ratio) < 0.0) tnoz = 0.0;
				// if (tnoz > bytrange) tnoz = bytrange;
				// sprintf(round, "%.0lf", tnoz);
				// tnoz_ts[day] = (char)atoi(round) + bmin;
				tnoz_ts[day] = tnoz;
			}
			
		} /* end of nsimdays loop */
		
		/* free memory depending on count */
		free(tairstr.id);
		free(tairstr.wt);
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
		free(tairstr.regwt);
		for (j=0 ; j<3 ; j++)
		{
			free(tairstr.inv[j]);
		}
		free(tairstr.inv);
		
		/* write daily data to outfiles */
		fwrite(t_ts, sizeof(double), nsimdays, t_out_f.ptr);
		if (noz)
		{
			fwrite(tnoz_ts, sizeof(double), nsimdays, tnoz_out_f.ptr);
		}
		
		/* calculate temporal average */
		t_tavg /= rn;
		mx_tavg /= rn;
		my_tavg /= rn;
		mz_tavg /= rn;
		if (noz)
		{
			tnoz_tavg /= rn;
		}
		
		/* write temporal averages to outfiles */
		fwrite(&t_tavg, sizeof(double), 1, t_tavg_f.ptr);
		fwrite(&mx_tavg, sizeof(double), 1, mx_tavg_f.ptr);
		fwrite(&my_tavg, sizeof(double), 1, my_tavg_f.ptr);
		fwrite(&mz_tavg, sizeof(double), 1, mz_tavg_f.ptr);
		if (noz)
		{
			fwrite(&tnoz_tavg, sizeof(double), 1, tnoz_tavg_f.ptr);
		}
		
		/* switch initial sign for next cell */
		tairstr.switch_init = -tairstr.switch_init;
		
	} /* end of ncells loop */
	
	nc_close(ncid_dem);
			
	/* calculate spatial averages for time series arrays */
	if (map) rn = (double)n_in;
	if (!map) rn = (double)ncells;
	
	for (i=0 ; i<nsimdays ; i++)
	{
		t_savg_ts[i] /= rn;
		mx_savg_ts[i] /= rn;
		my_savg_ts[i] /= rn;
		mz_savg_ts[i] /= rn;
		if (noz)
		{
			tnoz_savg_ts[i] /= rn;
		}
	}

	/* write time series arrays to outfiles */
	fwrite(t_savg_ts, sizeof(double), nsimdays, t_savg_f.ptr);
	fwrite(mx_savg_ts, sizeof(double), nsimdays, mx_savg_f.ptr);
	fwrite(my_savg_ts, sizeof(double), nsimdays, my_savg_f.ptr);
	fwrite(mz_savg_ts, sizeof(double), nsimdays, mz_savg_f.ptr);
	if (noz)
	{
		fwrite(tnoz_savg_ts, sizeof(double), nsimdays, tnoz_savg_f.ptr);
	}

	/* close output files */
	fclose(t_out_f.ptr);
	fclose(t_tavg_f.ptr);
	fclose(t_savg_f.ptr);
	fclose(mx_tavg_f.ptr);
	fclose(mx_savg_f.ptr);
	fclose(my_tavg_f.ptr);
	fclose(my_savg_f.ptr);
	fclose(mz_tavg_f.ptr);
	fclose(mz_savg_f.ptr);
	if (noz)
	{
		fclose(tnoz_out_f.ptr);
		fclose(tnoz_tavg_f.ptr);
		fclose(tnoz_savg_f.ptr);
	}
	fclose(ct_f.ptr);
	fclose(id_f.ptr);
	fclose(wt_f.ptr);
	
	/* free remaining memory */
	free(tairstr.stnobs);
	free(tairstr.stnsmobs);
	free(tairstr.stnz);
	free(t_ts);
	free(t_savg_ts);
	free(mx_savg_ts);
	free(my_savg_ts);
	free(mz_savg_ts);
	if (noz)
	{
		free(tnoz_ts);
		free(tnoz_savg_ts);
	}
	
	return(0);
	
} /* end of main */
