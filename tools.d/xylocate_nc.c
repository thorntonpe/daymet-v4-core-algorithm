/* 
xylocate.c
Peter Thornton

8/29/96

Given a station meta data file (*.stnmeta), this program searches a DEM
within a user-specified radius around the official location for each
station, and generates a corrected (x,y) location for the station by
minimizing the error between the official station elevation and the
DEM elevations in the specified neighborhood.
The station meta file is updated.

If using an initialization file, enter that filename on the command line.
Otherwise, the user is prompted for all required inputs.

5/27/97 Added a command-line option, detected by the number of command-line
parameters. Format:
<exec> <mask> <dem> <stnmeta> <search_width> <max_error>

12/9/09, PET: switched to using the DEM no_data value to determine which 
stations should not be corrected. For now the no_data value is hardwired to -9999.
This was introduced during the ORCA testing.

3/24/10, PET: multiple clean-up changes to get rid of unused code left in since
netcdf conversion. Also removed use of mask file - was opened and read (large memory
allocation) but never used.

*/

#include "metsrc.h"
#include "dailywx.h"
#include "cproj.h"
#include "proj.h"
#include "netcdf.h"

int main(int argc, char *argv[])
{
	double search_width;
	double ul_mapx, ul_mapy;
	double cellsize;
	double max_error, min_dz;
	double stn_mapx, stn_mapy, stn_realz;
	double f_x, f_y;
	double dz, demz;
	double d,nd;
	double mean_z;
	
	int n_rows, n_cols, nstns;
	int width, hw;
	int stn_filex, stn_filey, sx, sy;
	int nx, ny;
	int i,j,k,l;
	int in;
	int pout;
	
	size_t lenp;
	size_t index[2];
	long dem_cell, dem_stn;

	char outprefix[80];
	char round[32];

	long iflg;
	long (*for_trans[MAXPROJ + 1])();

	gctp_struct gctp;
        imghdr_struct dem_nc_hdr;
        file dem_nc_f, stnmeta_nc_f;
        int status, ncid;
        int nstnid, descid, yrdid;
        int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
        int stnxid, stnyid;
        int stnfxid, stnfyid;
        int ncid_dem, colsid_dem, rowsid_dem, imageid_dem;
        int ncid_meta;
        int year, year_day, descriptor_length;
        double *elevs, *lats, *lons;
	double *stnx, *stny;
	double *stnfx, *stnfy;
	int xyloc_flag = 0;
	/* hardwire the dem no_data value */ 
	int no_data = -9999;
	
	/* check for number of command line arguments */
	if (argc == 6)
	{
		strcpy(dem_nc_f.name, argv[2]);
		strcpy(stnmeta_nc_f.name, argv[3]);
		search_width = atof(argv[4]);
		max_error = atof(argv[5]);
		pout = 0;
	}
	else
	{
		printf("usage: exe <mask.nc> <dem.nc> <stnmeta.nc> <search radius (m)> <max_elev error (m)>\nexiting\n");
		exit(1);
	}


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
/*                                                                                             */
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
        nstns = lenp;
        status = nc_inq_dimid (ncid_meta, "year_day", &yrdid);
        status = nc_inq_dimlen(ncid_meta, yrdid, &lenp);
        year_day = lenp;
        status = nc_inq_dimid (ncid_meta, "descriptor_length", &descid);
        status = nc_inq_dimlen(ncid_meta, descid, &lenp);
        descriptor_length = lenp;

        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(nstns * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, stnlonid, lons);
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */

        /* using gctp structure just defined, initialize the transformation
        function for (lon,lat) --> (x,y) */
        gctp = dem_nc_hdr.gctp;
        for_init(gctp.outsys,gctp.outzone,gctp.outparm,gctp.outdatum,gctp.fn27,
                gctp.fn83,&iflg,for_trans);
	
	/* set internal variables */
	n_cols = dem_nc_hdr.ncols;
	n_rows = dem_nc_hdr.nrows;
	cellsize = dem_nc_hdr.cellsize;
	ul_mapx = dem_nc_hdr.ulx;
	ul_mapy = dem_nc_hdr.uly;

	/* allocate memory */
	if (!(stnx = (double*) malloc(nstns*sizeof(double)))) exit(1);
	if (!(stny = (double*) malloc(nstns*sizeof(double)))) exit(1);
	if (!(stnfx = (double*) malloc(nstns*sizeof(double)))) exit(1);
	if (!(stnfy = (double*) malloc(nstns*sizeof(double)))) exit(1);
	
	/* set size of location search window */
	width = (int) (search_width / cellsize);
	if (!(width % 2)) width++;
	hw = width/2;

	/* read station meta data one station at a time */
	for (i=0 ; i<nstns ; i++)
	{
                for_trans[gctp.outsys](lons[i]*D2R,lats[i]*D2R,&stn_mapx,&stn_mapy);
		stnx[i] = stn_mapx;
		stny[i] = stn_mapy;
		stn_realz = elevs[i];

		/* calculate initial file (x,y) coordinates for station	*/
		f_x = (stn_mapx - ul_mapx)/cellsize;
		f_y = (ul_mapy - stn_mapy)/cellsize;
		sprintf(round,"%.0lf",f_x);
		stn_filex = atoi(round);
		sprintf(round,"%.0lf",f_y);
		stn_filey = atoi(round);
		
		/* set netcdf index for single-value extraction */
		index[0]=stn_filey;
		index[1]=stn_filex;
		status=nc_get_var1_long(ncid_dem, imageid_dem, index, &dem_stn);

		/* test if station is within mask, and within DEM data region */
		/* if not, use original location , continue */
		/* if so, scan DEM for best estimate of station location */
		if ((stn_filex < 0) || (stn_filex >= n_cols) ||
			(stn_filey < 0) || (stn_filey >= n_rows) ||
			(dem_stn == no_data))
		{
			stnfx[i] = stnx[i];
			stnfy[i] = stny[i];
			
			/* send message to terminal that this station is outside mask */
			if (pout) printf("Station #%5d: No location correction\n",i); 
			continue;
		}
		
		// Use elevation from the DEM if station elevation from metadata < -999 m.
		// That catches the missing values in station data (-999.9), and also the
		// possiblity of a bad (low) value.
		// Also use DEM value if metadata elevation is > 9000 m.
		// If bad station elevation and missing data in DEM, assume coast.
		if ((stn_realz < -999.0) || (stn_realz > 9000.0))
		{
			stnfx[i] = stnx[i];
			stnfy[i] = stny[i];
			if (dem_stn == no_data)
			{
			    elevs[i] = 0.0;
			}
			else
			{
			    elevs[i] = dem_stn;
			}
			continue;
		}
		
		/* station is in mask, and inside DEM data region */
		/* search for DEM cell that minimizes elevation error */
		min_dz = max_error;
		nd = 1e6;
		nx = stn_filex;
		ny = stn_filey;
		for (j=-hw ; j<=hw ; j++)
		{
			sy = stn_filey + j;
			/* test if row is within mask */
			if ((sy < 0) || (sy >= n_rows))
				continue;

			for (k=-hw ; k<=hw ; k++)
			{
				sx = stn_filex + k;
				/* test if column is within mask */
				if ((sx < 0) || (sx >= n_cols))
					continue;

				/* set netcdf index for single-value extraction */
				index[0]=sy;
				index[1]=sx;
				status=nc_get_var1_long(ncid_dem, imageid_dem, index, &dem_cell);

				/* test if scan-cell is within DEM data region */
				if (dem_cell == no_data)
					continue;

				demz = (double) dem_cell;
				dz = fabs(demz - stn_realz);
				if (dz < min_dz)
				{
					min_dz = dz;
					nx = sx;
					ny = sy;
					nd = sqrt(j*j + k*k);
				}
				
				/* in the case of a tie, the cell nearer the original center
				wins */
				if (dz == min_dz)
				{
					d = sqrt(j*j + k*k);
					if (d < nd)
					{
						nx = sx;
						ny = sy;
						nd = d;
					}
				}
			}
		}
		
		/* print if no gridcell error less than initial min_dz is found */
		if (min_dz == max_error)
		{
			if (pout) printf("Station #%5d: No location correction\n",i);
		}
			
		/* convert file coord back to map coord for output */
		stn_filex = nx;
		stn_filey = ny;
		stnfx[i] = ul_mapx + (double)stn_filex * cellsize;
		stnfy[i] = ul_mapy - (double)stn_filey * cellsize;
				
	} /* end of nstns loop */
	
	nc_close(ncid_dem);

	/* testing - calculate mean stnz */
	mean_z = 0.0;
	for (i=0 ; i<nstns ; i++)
	{
		mean_z += elevs[i];
	}
	mean_z = mean_z/(double)nstns;
	printf("mean_z = %lf\t%s\n",mean_z,stnmeta_nc_f.name);
	

/* -----------------------  netcdf calls for station meta file update  --------------------------------------- */
/*                                                                                                             */
        status = nc_redef(ncid_meta);
	if (status != NC_NOERR)
	{
		printf("Error entering redef mode for stnmeta file\n");
		exit(1);
	}

        status = nc_def_var(ncid_meta, "stnx",  NC_DOUBLE, 1, &nstnid, &stnxid);
	if (status != NC_NOERR)
	{
		printf("Error defvar stnx stnmeta file\n");
		exit(1);
	}
        status = nc_def_var(ncid_meta, "stny",  NC_DOUBLE, 1, &nstnid, &stnyid);
	if (status != NC_NOERR)
	{
		printf("Error defvar stny stnmeta file\n");
		exit(1);
	}
        status = nc_def_var(ncid_meta, "stnfx", NC_DOUBLE, 1, &nstnid, &stnfxid);
	if (status != NC_NOERR)
	{
		printf("Error defvar stnfx stnmeta file\n");
		exit(1);
	}
        status = nc_def_var(ncid_meta, "stnfy", NC_DOUBLE, 1, &nstnid, &stnfyid);
	if (status != NC_NOERR)
	{
		printf("Error defvar stnfy stnmeta file\n");
		exit(1);
	}

	xyloc_flag = 1;
        status = nc_put_att_int(ncid_meta, NC_GLOBAL, "xyloc_flag", NC_INT, 1, &xyloc_flag);
	if (status != NC_NOERR)
	{
		printf("Error att_put xyloc_flag stnmeta file\n");
		exit(1);
	}

        status = nc_enddef(ncid_meta);
	if (status != NC_NOERR)
	{
		printf("Error enddef stnmeta file\n");
		exit(1);
	}

        status = nc_put_var_double(ncid_meta, stnxid,  stnx);
	if (status != NC_NOERR)
	{
		printf("Error put stnx stnmeta file\n");
		exit(1);
	}
        status = nc_put_var_double(ncid_meta, stnyid,  stny);
	if (status != NC_NOERR)
	{
		printf("Error put stny stnmeta file\n");
		exit(1);
	}
        status = nc_put_var_double(ncid_meta, stnfxid, stnfx);
	if (status != NC_NOERR)
	{
		printf("Error put stnfx stnmeta file\n");
		exit(1);
	}
        status = nc_put_var_double(ncid_meta, stnfyid, stnfy);
	if (status != NC_NOERR)
	{
		printf("Error put stnfy stnmeta file\n");
		exit(1);
	}
        status = nc_put_var_double(ncid_meta, stnelevid, elevs);
	if (status != NC_NOERR)
	{
		printf("Error put elevs stnmeta file\n");
		exit(1);
	}

	status = nc_sync(ncid_meta);
	if (status != NC_NOERR)
	{
		printf("Error sync stnmeta file\n");
		exit(1);
	}
	status = nc_close(ncid_meta);
	if (status != NC_NOERR)
	{
		printf("Error close stnmeta file\n");
		exit(1);
	}
/*                                                                                                             */
/* -----------------------  end netcdf calls for station meta file update  ----------------------------------- */
	
	free(elevs);
	free(lats);
	free(lons);
	free(stnx);
	free(stny);
	free(stnfx);
	free(stnfy);
	
	return(0);
}

	
