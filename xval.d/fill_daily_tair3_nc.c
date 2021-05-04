/*
fill_daily_tair3_nc.c
Peter Thornton
11/11/96

Fills in the missing data values in a station data file using a specified 
set of interpolation and prediction parameters.
The missing data flags remain in the data file, and so these missing
values can be filled in again (overwritten) using other interpolation and/or
prediction parameters at a later time if desired.

Updated:
9/5/00, PET: from fill_daily_tair.c, added 3-variable regression
July 2005 : netcdf version - CSH
4/12/2020, PET: reduce number of stations when not enough good data (short station lists like PR)
*/

#include "metsrc2.h"
#include "dailywx.h"
#include "cproj.h"
#include "proj.h"
#include "netcdf.h"

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
        gctp_struct gctp;
        long iflg;

	char *good;
	int gooddims[2];
	char *listgood;
	
	file ini;
	file stnmeta_f, stndata_f;
	file mask_f, dem_f;
	file tmp_f;

	char endkey[8];
	char round[32];
	char outprefix[80];
	char jobstr[80];

	short *dem_array;
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
	int sstn;
	int n_fill;
	int job;
	
	long off;
	long int obs_off, obs_offset;
	long img_off, img_offset, stn_offset, koff, koffset;

	double dfirad;
	double ul_mapx, ul_mapy;
	double cellsize;
	double fzwidth;
	double *kernel;
	double stn_mapx, stn_mapy;
	double f_x, f_y;
	double value,sum_wt;
	double t1;
	
	size_t lenp;
	int xyloc_flag;
	stnmetahdr_struct meta_nc_hdr;
	imghdr_struct mask_nc_hdr, dem_nc_hdr;
        long int *mask_array;
        file mask_nc_f, stnmeta_nc_f;
        file stndata_nc_f;
        int status, ncid;
        int nstnid, descid, yrdid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid, stnzid, goodid;
        int stnfxid, stnfyid, stnfzid;
        int stntypeid, stnuc1id, stnuc2id;
        int tairid;
        int ncid_mask, colsid_mask, rowsid_mask, imageid_mask;
        int ncid_meta;
	int *mflag, *code;
	int mflagid, codeid;
        int year_day, descriptor_length;
        double *elevs, *lats, *lons, *tair;
	double *stnx, *stny, *stnfx, *stnfy;
	int fill_flag = 1;
	
	short *bsi;
	short temp_short;
	double wtsum;
    int ngood_today;

	/* check for appropriate number of command line arguments */
	if (argc != 9)
	{
		printf("Usage: \n");
		printf("<exec> <mask> <meta> <data> <zfilter> <dfirad> <init_gsp> <ans> <smwidth>\n");
		exit(1);
	}
	/* assuming command-line parameter passing */
	/* scan command line arguments into parameters */
       /* scan next parameter as the mask filename */
        strcpy(mask_nc_f.name,argv[1]);
	strcpy(stnmeta_nc_f.name, argv[2]);
	strcpy(stndata_nc_f.name, argv[3]);
	fzwidth = atof(argv[4]);
	dfirad = atof(argv[5]);
	intstr.gsp = atof(argv[6]);
	intstr.ans = atof(argv[7]);
	smwidth = atoi(argv[8]);
	
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
        status = nc_get_att_text   (ncid_mask, NC_GLOBAL, "tag",            mask_nc_hdr.tag);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "ulx",           &mask_nc_hdr.ulx);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "uly",           &mask_nc_hdr.uly);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "cellsize",      &mask_nc_hdr.cellsize);
        status = nc_get_att_long   (ncid_mask, NC_GLOBAL, "gctp_outsys",   &mask_nc_hdr.gctp.outsys);
        status = nc_get_att_long   (ncid_mask, NC_GLOBAL, "gctp_outzone",  &mask_nc_hdr.gctp.outzone);
        status = nc_get_att_long   (ncid_mask, NC_GLOBAL, "gctp_outdatum", &mask_nc_hdr.gctp.outdatum);
        status = nc_get_att_double (ncid_mask, NC_GLOBAL, "gctp_outparm",   mask_nc_hdr.gctp.outparm);
        status = nc_get_att_text   (ncid_mask, NC_GLOBAL, "gctp_fn27",      mask_nc_hdr.gctp.fn27);
        status = nc_get_att_text   (ncid_mask, NC_GLOBAL, "gctp_fn83",      mask_nc_hdr.gctp.fn83);
        if (strncmp(mask_nc_hdr.tag,IMGTAG,7))
        {
                printf("mask file not a registered image. Use register_image first.\n");
                exit(1);
        }

        status = nc_inq_varid (ncid_mask, "image", &imageid_mask);
        if (!(mask_array = (long int*) malloc(mask_nc_hdr.nrows * mask_nc_hdr.ncols * sizeof(long int)))) exit(1);
        status = nc_get_var_long (ncid_mask, imageid_mask, mask_array);

        nc_close(ncid_mask);
/*                                                                                                 */
/* -----------------------  end netcdf calls for mask file --------------------------------------- */

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

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlonid, lons);

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

        status = nc_inq_varid (ncid_meta, "mflag", &mflagid);
        if (!(mflag = (int*) malloc(meta_nc_hdr.nstns * year_day * sizeof(int))))
	{
                printf("Error allocating mflag array.\n");
		exit(1);
	}
        status = nc_get_var_int (ncid_meta, mflagid, mflag);

        status = nc_inq_varid (ncid_meta, "tmax", &tairid);
	if (status)
	{
		status = nc_inq_varid (ncid_meta, "tmin", &tairid);
	}
        if (!(tair = (double*) malloc(meta_nc_hdr.nstns * year_day * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, tairid, tair);

        meta_nc_hdr.zfilt_flag = 0;
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
		if (scan_value(ini, &type, 'i')) exit(1);
		if (scan_open(ini, &mask_f, 'r')) exit(1);
		if (scan_open(ini, &dem_f, 'r')) exit(1);
		fclose(ini.ptr);
	}

	/* perform error checking on input file headers */
	if (fzwidth)
	{
		/* read mask file header and verify */
		fread(&maskhdr,sizeof(imghdr_struct), 1, mask_f.ptr);
		if (strcmp(maskhdr.tag,IMGTAG))
		{
			printf("%s not a registered image. Use register_image first.\n",mask_f.name);
			exit(1);
		}
		
		/* read dem file header and verify */
		fread(&demhdr,sizeof(imghdr_struct), 1, dem_f.ptr);
		if (strcmp(demhdr.tag,IMGTAG))
		{
			printf("%s not a registered image. Use register_image first.\n",dem_f.name);
			exit(1);
		}
		/* check that dem is sizeof(short) */
		
		/* compare headers for DEM and mask */
		test = (maskhdr.ncols == demhdr.ncols)*(maskhdr.nrows == demhdr.nrows)*
			(maskhdr.ulx == demhdr.ulx)*(maskhdr.uly == demhdr.uly)*
			(maskhdr.cellsize == demhdr.cellsize);
		if ((!test) || memcmp(&maskhdr.gctp,&demhdr.gctp,sizeof(gctp_struct)))
		{
			printf("Header information differs between mask and DEM...exiting\n");
			exit(1);
		}
	}
	
	if (fzwidth)
	{
		/* compare gctp structures between images and stnmeta */
		if (memcmp(&maskhdr.gctp,&metahdr.gctp,sizeof(gctp_struct)))
		{
			printf("Different projection parameters in images and stnmeta.\n");
			exit(1);
		}
	}

	/* set internal variables */
	if (fzwidth)
	{
		ncols = maskhdr.ncols;
		nrows = maskhdr.nrows;
		cellsize = demhdr.cellsize;
		ul_mapx = demhdr.ulx;
		ul_mapy = demhdr.uly;
	}
	nstns = tairstr.nstns = intstr.nstns = meta_nc_hdr.nstns;
	printf("Number stations: %d\n", nstns); 
	ndays = tairstr.ndays = year_day;
	tairstr.noz = 0;

	if (fzwidth)
	{
		/* allocate zfilter memory */
		if (!(dem_array = (short*) malloc(nrows*ncols*sizeof(short)))) exit(1);
		if (!(mask_array = (long int*) malloc(nrows*ncols*sizeof(long int)))) exit(1);
	
		/* read and close dem and mask files */
		fread(dem_array,sizeof(short),nrows*ncols,dem_f.ptr);
		fread(mask_array,sizeof(char),nrows*ncols,mask_f.ptr);
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
	if (!(good = (char*) malloc(nstns * ndays * sizeof(char)))) exit(1);
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

	if (!(bsi = (short*) malloc(nstns * sizeof(short)))) exit(1);
	
	/* istr and tstr share space for wt and id arrays */
	tairstr.wt = intstr.wt;
	tairstr.id = intstr.id;
	
	/* read stnmeta data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		if (xyloc_flag)
		{
			tairstr.stnx[i] = stn_mapx = intstr.stnx[i] = stnfx[i];
			tairstr.stny[i] = stn_mapy = intstr.stny[i] = stnfy[i];
		}
		else
		{
			tairstr.stnx[i] = stn_mapx = intstr.stnx[i] = stnx[i];
			tairstr.stny[i] = stn_mapy = intstr.stny[i] = stny[i];
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
				tairstr.stnz[i] = elevs[i];
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
				tairstr.stnz[i] = elevs[i];
				continue;
			}
			tairstr.stnz[i] = value/sum_wt;
		}
		else /* fzwidth == 0 */
		{
			tairstr.stnz[i] = elevs[i];
		}
		
	} /* end of nstns loop */

	/* read station daily data and fill appropriate arrays */
    for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			tairstr.stnobs[obs_offset] = tair[obs_offset];
			if (mflag[obs_offset] == 1) good[obs_offset] = 0;
			else {
                		good[obs_offset] = 1;
			}
		}
	}
	
	/* smooth the station data arrays by station
	For missing-value-fill code, this smoothing algorithm also requires an
	array of characters indicating good and missing data (1=good, 0=miss) */
	for (i=0 ; i<nstns ; i++)
	{
		off = i * ndays;
		if (fill_boxcar_smooth(good+off,tairstr.stnobs+off,tairstr.stnsmobs+off,
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
	n_fill = 0;
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
		
		obs_off = i * ndays;

		/* begin loop through nsimdays */
		for (day=0 ; day<ndays ; day++)
		{
            /* only proceed with data filling if this day has a missing value */
			/* otherwise, copy the good value to temporary file and continue */
			if (good[obs_off+day])
			{
				continue;
			}
			
			n_fill++;
			
			/* do the bubble-sort to find the ANS nearest stations with good data */
			for (j=0 ; j<nstns ; j++) bsi[j]=j;
			for (j=0 ; j<intstr.ans ; j++)
			{
				for (k=nstns-1 ; k>j ; k--)
				{
					if (good[bsi[k]*ndays+day] && ((intstr.sqdist[bsi[k]] < intstr.sqdist[bsi[k-1]]) || !good[bsi[k-1]*ndays+day])) 
					{
						/* switch these two stations in bsi */
						temp_short = bsi[k];
						bsi[k]=bsi[k-1];
						bsi[k-1]=temp_short;
					}
				}
			}
            
            /* PET (4/12/2020): Count the number of good stations for the day */
            /* up to a max of ANS */
            ngood_today = 0;
            for (j=0 ; j<intstr.ans ; j++)
            {
                if (!good[bsi[j]*ndays+day])
                {
                    break;
                }
                ngood_today++;
            }
            /* If no good stations for the day, exit with error */
            if (ngood_today < 2)
            {
                printf("Error: no good stations for station %d, day %d. Exiting.\n",i,day);
                exit(1);
            }
            
            /* PET (4/12/2020): Moved the allocation for count and nregr inside the days loop */
            /* This prevents using stations that don't have good data, when the total number */
            /* of stations in the station list is low (like for the Puerto Rico runs) */
            intstr.count = ngood_today; 
            tairstr.count = count = intstr.count;
            tairstr.nregr = nregr = ((count * count) - count) / 2;
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
            
			/* set intstr.sqrad as the distance to the ngood_today station in bsi list */
			/* this forces the truncated filter weight to 0 just past the ANS station in list */
			intstr.sqrad = intstr.sqdist[bsi[intstr.count-1]] + 1.0;
			
			/* generate the weight and id lists */
			wtsum = 0.0;
			for (j=0 ; j<count ; j++)
			{
				intstr.id[j]=bsi[j];
				intstr.wt[j]=exp(-intstr.gsp*intstr.sqdist[bsi[j]]*(1.0/intstr.sqrad))-intstr.trunc;
				wtsum += intstr.wt[j];
			}
			for (j=0 ; j<count ; j++)
			{
				intstr.wt[j] /= wtsum;
			}
			
			/* fill the station-list arrays */
			/* at the moment this is filling the arrays for all days, which is unnecessary, but won't hurt anything */
			if (xv_fill_list(good, listgood, &tairstr))
			{
				printf("Error filling station list arrays ... exiting\n");
				exit(1);
			}

			/* generate regression arrays */
			if (xv_tair_regr_xwt_2switch(&tairstr))
			{
				printf("Error in xv_tair_regr_xwt_2switch() ... exiting\n");
				exit(1);
			}

			/* initialize daily prediction */
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
			
			/* A filled data value is required, by hook or by crook, so 
			most desirable option is DayMet routines, but if there isn't
			enough good local data for the station-day in question, then
			resort to the nearest neighbor with good data and use a default
			lapse rate. */
			if (ngoodpts >= 3)
			{
				// /* generate regression arrays that change by day */
				// if (xv_tair_regr_y_2switch(listgood, &tairstr))
				// {
				// 	printf("Error in xv_tair_regr_y_switch() ... exiting\n");
				// 	exit(1);
				// }

				// /* weighted least squares slope and intercept */
				// if (xv_tair_mlr(&tairstr))
				// {
				// 	printf("Error in wt_regr() ... exiting\n");
				// 	exit(1);
				// }

				tairstr.coef[0] = tairstr.coef[1] = tairstr.coef[3] = 0.0;
                tairstr.coef[2] = DEFAULT_TVSZ;
				
				/* predict tair for the day */
				if (xv_predict_tair_noint(listgood, &tairstr))
				{
					printf("Error in xv_predict_tair_noint() ... exiting\n");
					exit(1);
				}
			} /* end if ngoodpts >= 3 */
			else if (ngoodpts) /* 1 or 2 good stations */
			{
				/* find the closest of these stations */
				t1=0.0;
				for (j=0 ; j<count ; j++)
				{
					if (listgood[tairstr.dam_off + j] && tairstr.wt[j] > t1)
					{
						sstn = j;
						t1 = tairstr.wt[j];
					}
				}
				/* if a good station is found (should always be true) */
				if (t1)
				{
					tairstr.ta = tairstr.listobs[tairstr.dam_off + sstn] + 
					(DEFAULT_TVSZ * (tairstr.elev - tairstr.listz[sstn]));
				}
				else
				{
					printf("Error: ngoodpts > 0, but no good stations found\n");
					printf("This should never happen!\n");
					printf("exiting...\n");
					exit(1);
				}
			}
			else /* no good satations */
			{
				printf("Error: No good stations found\n");
				printf("Suggest using a larger value for average number of stations\n");
				printf("exiting...\n");
                		for (j=0 ; j<count ; j++)
                		{
                    			printf("stnin %d\t%d\t%0.4lf\t%0.4lf\t%d\n",j,tairstr.id[j],tairstr.wt[j],tairstr.listobs[tairstr.dam_off + j],listgood[tairstr.dam_off +j]);
                		}
				exit(1);
			}
			
			/* by now there is either a filled data value or an exit on error,
			so update the data structure and write to temporary file */
			tair[obs_off+day] = tairstr.ta;
            
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
					
		} /* end of nsimdays loop */
		
		/* switch initial sign for next station */
		tairstr.switch_init = -tairstr.switch_init;
		
	} /* end of nstations loop */
	
/* -----------------------  netcdf calls for station meta file update  --------------------------------------- */
/*                                                                                                             */
                status = nc_redef(ncid_meta);

		gooddims[0] = nstnid;
		gooddims[1] = yrdid;
                status = nc_def_var(ncid_meta, "good", NC_BYTE, 2, gooddims, &goodid);

                status = nc_put_att_int(ncid_meta, NC_GLOBAL, "fill_flag", NC_INT, 1, &fill_flag);

                status = nc_enddef(ncid_meta);

                status = nc_put_var_double(ncid_meta, tairid, tair);
                status = nc_put_var_uchar(ncid_meta, goodid, good);
/*                                                                                                             */
/* -----------------------  end netcdf calls for station meta file update  ----------------------------------- */
	
	/* close files and exit */
        nc_close(ncid_meta);
	
	/* free remaining memory */
	if (fzwidth)
	{
		free(dem_array);
		free(mask_array);
		free(kernel);
	}
	free(good);
	free(tairstr.stnobs);
	free(tairstr.stnsmobs);
	free(tairstr.stnz);
	free(intstr.stnx);
	free(intstr.stny);
	free(intstr.sqdist);
	free(intstr.wt);
	free(intstr.id);

	exit(0);
	
} /* end of main */
