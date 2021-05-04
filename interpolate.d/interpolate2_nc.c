/* 
interpolate2_nc.c
Peter Thornton

Revision history:
July 2005  changed to netcdf I/O - CSH
8/6/96   Revision of daymet_make_cart.c, accept either MAP or LIST I/O
9/6/96   adapted from daymet_interpolate.c, now general to metsrc
4/2/97   added command-line parameter passing option
8/14/00  adapted from interpolate.c, added a new parameter that specifies 
		the mask value for inclusion (for US daymet runs with multiple 
		2x2 degree regions). Removed LIST support.
11/14/2014 Revised to use a bubble-sort algorithm to find exactly ANS stations for the prediction at each point. Also doing some cleanup on debug code. Added a preproc flag EXACT_ANS.
The addition of bubble sort was tested for tile 12293, year 1998, and found to remove all signs of linear tracks or striping in the output files. It was also demopnstrated that every gridcell get exactly ANS input stations, and that all stations receive non-zero weight.
*/

#include "metsrc.h"
#include "dailywx.h"
#include "cproj.h"
#include "proj.h"
#include "netcdf.h"
#define EXACT_ANS


int main(int argc, char *argv[])
{
	/* variable declarations */
	interpolate_struct str;
	lsthdr_struct list_hdr;

	file ini;
	file list_f;
	file ct_f, id_f, wt_f;
	
	double ulmapx, ulmapy, cellsize;
	double dfirad;
	double cellx, celly;
	double celldat[5];
	
	char endkey[8];
	char iocase[16];
	char c,ioc,map;
	char outprefix[80];
	long int mask;
	char userin[80];

	int ncols, nrows, ncells, nstns;
	int i, j, k, row, col;
	int err;
	long int good_value;

	size_t lenp;
	int xyloc_flag;
	imghdr_struct mask_nc_hdr;
	stnmetahdr_struct meta_nc_hdr;
	long int *mask_array;
	file mask_nc_f, stnmeta_nc_f;
	int status, ncid;
	int nstnid, descid, yrdid;
	int stnnameid, stnidid, stnelevid, stnlatid, stnlonid;
	int stnxid, stnyid, stnzid;
	int stnfxid, stnfyid, stnfzid;
	int stntypeid, stnuc1id, stnuc2id;
	int tmaxid;
	int ncid_mask, colsid_mask, rowsid_mask, imageid_mask;
	int ncid_meta;
	int ans;
	//int iyear = 2002;
	//int iyrdays = 366;
	int year_days, descriptor_length;
	double *elevs, *lats, *lons, *stnx, *stny;
	double *stnfx, *stnfy;
	short *bsi;
	short temp_short;
	double wtsum;

	/* check for number of command line arguments */
	/* assume that parameters are being passed on the command line */
	map = 1;
	/* scan next parameter as the mask filename */
	strcpy(mask_nc_f.name,argv[1]);
	/* scan the next parameter as the good mask value */
	good_value = atoi(argv[2]);
	/* read command-line parameters that don't depend on I/O type */
	strcpy(stnmeta_nc_f.name,argv[3]);
	dfirad = atof(argv[4]);
	str.gsp = atof(argv[5]);
	str.ans = atof(argv[6]);
	strcpy(outprefix,argv[7]);
	ans=(int)str.ans;
	
/* -----------------------  netcdf calls for mask file --------------------------------------- */
/*                                                                                             */
	//printf("Read in NetCDF\n");
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
                exit(1);
        }

        status = nc_inq_varid (ncid_mask, "image", &imageid_mask);
        if (!(mask_array = (long int*) malloc(mask_nc_hdr.nrows * mask_nc_hdr.ncols * sizeof(long int)))) { printf("Could no malloc mask_array"); exit(1);}
        status = nc_get_var_long (ncid_mask, imageid_mask, mask_array);

	nc_close(ncid_mask);
/*                                                                                                 */
/* -----------------------  end netcdf calls for mask file --------------------------------------- */

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
  
        status = nc_inq_varid (ncid_meta, "station_name", &stnnameid);
        status = nc_inq_varid (ncid_meta, "station_id", &stnidid);

        status = nc_inq_varid (ncid_meta, "station_elevation", &stnelevid);
        if (!(elevs = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc elevs\n"); exit(1); }
        status = nc_get_var_double (ncid_meta, stnelevid, elevs);

        status = nc_inq_varid (ncid_meta, "station_latitude", &stnlatid);
        if (!(lats = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc lats\n"); exit(1); }
        status = nc_get_var_double (ncid_meta, stnlatid, lats);

        status = nc_inq_varid (ncid_meta, "station_longitude", &stnlonid);
        if (!(lons = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc lons\n"); exit(1); }
        status = nc_get_var_double (ncid_meta, stnlonid, lons);

        status = nc_inq_varid (ncid_meta, "stnx", &stnxid);
        if (!(stnx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc stnx\n"); exit(1); }
        status = nc_get_var_double (ncid_meta, stnxid, stnx);

        status = nc_inq_varid (ncid_meta, "stny", &stnyid);
        if (!(stny = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc stny\n"); exit(1); }
        status = nc_get_var_double (ncid_meta, stnyid, stny);

        status = nc_get_att_int (ncid_mask, NC_GLOBAL, "xyloc_flag", &xyloc_flag);
        if (xyloc_flag == 1)
        {
                status = nc_inq_varid (ncid_meta, "stnfx", &stnfxid);
                if (!(stnfx = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc stnfx\n"); exit(1); }
                status = nc_get_var_double (ncid_meta, stnfxid, stnfx);

                status = nc_inq_varid (ncid_meta, "stnfy", &stnfyid);
                if (!(stnfy = (double*) malloc(meta_nc_hdr.nstns * sizeof(double)))) { printf("Could no malloc stnfy\n"); exit(1); }
                status = nc_get_var_double (ncid_meta, stnfyid, stnfy);
        }
/*                                                                                                     */
/* -------------------  end netcdf calls for station meta file --------------------------------------- */

	/* error checking on file headers */
	if (strncmp(mask_nc_hdr.tag,IMGTAG,7))
	{
		printf("Mask not a registered image. Use register_image first.\n");
		exit(1);
	}

	/* create output filenames and open for binary write */
	if (bin_out(&ct_f, outprefix, ".interp_ct")) { printf("bin_out error ct\n"); exit(1); }
	if (bin_out(&id_f, outprefix, ".interp_id")) { printf("bin_out error id\n"); exit(1); }
	if (bin_out(&wt_f, outprefix, ".interp_wt")) { printf("bin_out error wt\n"); exit(1); }

	/* set internal variables */
	ncols = mask_nc_hdr.ncols;
	nrows = mask_nc_hdr.nrows;
	cellsize = mask_nc_hdr.cellsize;
	ulmapx = mask_nc_hdr.ulx;
	ulmapy = mask_nc_hdr.uly;
	ncells = ncols * nrows;
	nstns = str.nstns = meta_nc_hdr.nstns;
	/* allocate memory */
	if (!(str.stnx = (double*) malloc(nstns * sizeof(double)))) { printf("Could no malloc str.stnx\n"); exit(1); }
	if (!(str.stny = (double*) malloc(nstns * sizeof(double)))) { printf("Could no malloc str.stny\n"); exit(1); }
	if (!(str.sqdist = (double*) malloc(nstns * sizeof(double)))) { printf("Could no malloc str.sqdist\n"); exit(1); }
	if (!(str.wt = (double*) malloc(nstns * sizeof(double)))) { printf("Could no malloc str.wt\n"); exit(1); }
	if (!(str.id = (short*) malloc(nstns * sizeof(short)))) { printf("Could no malloc str.id\n"); exit(1); }
	/* index array needed for bubble sort algorithm */
	if (!(bsi = (short*) malloc(nstns * sizeof(short)))) {printf("Couild not malloc bsi\n"); exit(1);}
	/* transfer station (x,y) data from netcdf arrays */
	for (i=0 ; i<nstns ; i++)
    {
            if (xyloc_flag)
            {
                    str.stnx[i] = stnfx[i];
                    str.stny[i] = stnfy[i];
            }
            else
            {
                    str.stnx[i] = stnx[i];
                    str.stny[i] = stny[i];
            }
    }
	
	/* calculate initial density filter parameters */
	str.dfsqrad = dfirad * dfirad;
	str.inv_dfsqrad = 1.0 / str.dfsqrad;
	str.dfarea = PI * str.dfsqrad;
	
	/* calculate parameters for the Gaussian filter */
	str.trunc = exp(-str.gsp);
	str.dftrunc = exp(-DFGSP);
	str.dfavgwt = ((1.0 - str.dftrunc)/DFGSP) - str.dftrunc;

	/* begin main loop through ncells */
	for (i=0 ; i<ncells ; i++)
	{
		/* MAP I/O */
		mask = mask_array[i];
		if (mask != good_value)
		{
			/* cell not in interpolation region, continue */
			continue;
		}

		/* cell is in interpolation region */
		/* calculate cell (x,y) coordinates */
		row = i/ncols;
		col = i%ncols;
		str.x = ulmapx + ((double)col * cellsize);
		str.y = ulmapy - ((double)row * cellsize);
	        
		/* calculate the squared distance from this point to each station */
		if (calc_sqdist(&str))
		{
			printf("Error in call to calc_sqdist()... exiting\n");
			exit(1);
		}

#ifdef EXACT_ANS
		/* use a bubble-sort algorithm to get exactly ANS stations in the prediction list for each point */
		for (j=0 ; j<nstns ; j++) bsi[j]=j;
		for (j=0 ; j<ans ; j++) {
		  for (k=nstns-1 ; k>j ; k--) {
		    if (str.sqdist[bsi[k]] < str.sqdist[bsi[k-1]]) {
		      temp_short = bsi[k];
		      bsi[k]=bsi[k-1];
		      bsi[k-1]=temp_short;
		    }
		  }
		}
		/* set the weighting kernel radius to 1m greater than the furthest included station */
		/* the extra 1.0m ensures that none of the stations in the list will have wt=0 */
		str.sqrad=str.sqdist[bsi[ans-1]]+1.0;
		str.count = (short)ans;		
		wtsum = 0.0;
		for (j=0 ; j<ans ; j++) {
		  str.id[j]=bsi[j];
		  str.wt[j]=exp(-str.gsp * str.sqdist[bsi[j]] * (1.0/str.sqrad)) - str.trunc;
		  wtsum += str.wt[j];
		}
		/* normalize weights */
		if (wtsum) {
		  for (j=0 ; j<ans ; j++) str.wt[j] /= wtsum;
		} else {
		  printf("wtsum = 0... exiting\n");
		  exit(1);
		}
#else

		if (initial_radius(&str))
                {
                    printf("Error in call to initial_radius()... exiting\n");
                    exit(1);
                }

		/* iterative density algorithm to calculate squared search radius */
		if (calc_search_rad(&str))
		{
			printf("Error in call to calc_search_rad()... exiting\n");
			exit(1);
		}
	
		/* generate weighted station list for this point */
		if (weight_list(&str))
		{
			printf("Error in call to weight_list()... exiting\n");
			exit(1);
		}
#endif
		
		/* write count, id, and wt to output files */
        size_t temp;
		temp = fwrite(&str.count, sizeof(short), 1, ct_f.ptr);
        if(temp != 1) { printf("Error writing str.count to ct_f.ptr\n"); exit(1); } 
 
		temp = fwrite(str.id, sizeof(short), str.count, id_f.ptr);
		if(temp != str.count) {printf("Error writing str.id to id_f.ptr\n"); exit(1);}

        temp = fwrite(str.wt, sizeof(double), str.count, wt_f.ptr);
		if(temp != str.count)  {printf("Error writing str.wt to wt_f.prt\n"); exit(1); }
	} /* end of ncells loop */

	/* free memory */
	free(str.stnx);
	free(str.stny);
	free(str.sqdist);
	free(str.id);
	free(str.wt);
	
	/* close files */
	fclose(ct_f.ptr);
	fclose(id_f.ptr);
	fclose(wt_f.ptr);

	nc_close(ncid_meta);
    
    printf("Successful exit\n");
	exit(0);

} /* end of main */
		

