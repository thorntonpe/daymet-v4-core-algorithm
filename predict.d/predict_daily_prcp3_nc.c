/*
predict_daily_prcp3_nc.c
Peter Thornton

Revision history:
July 2005:   netcdf version - CSH
6/29/95:  Revision of daymet_prcp.c.
7/26/95:  Revision of daymet_prcp.c, added NDPI
8/1/96:  Revision of crb_daymet_prcp.c, general data structure for f()'s
9/20/96:  Revision to include the new metsrc.h structures
4/2/97: Added an option for passing all required parameters on the 
		command-line
4/2/97: altered treatment of start_day and stop_day: if both are set to 0,
 then all days are used
8/15/00: Adapted from predict_daily_prcp.c. Changed tavg output so that only
 good_value mask region gets a value, which means that these images need to be
 filtered through a mask for display. Also added the good_value parameter for
 US grid. Removed LIST option.
9/2/00: Adapted from predict_daily_prcp2.c, adding three-variable multiple
 linear regression model.
9/16/09: PET - fixed bug in section that updates stats. Added test for noz.
1/11/10: PET - fixed bug in read of DEM, introduced in netcdf code.



Output file notes:
*ttot_f  =  total precip over nsimdays 
*tavg_f  =  average daily amount for wet days
*tfrq_f  =  temporal frequency, fraction of nsimdays with precip

*stot_f  =  average daily precip over ncells
*savg_f  =  average daily precip over only wet cells
*sfrq_f  =  spatial frequency, fraction of ncells with precip

slope (m) files based on # days/cells with prcp 

*/

#include "metsrc2.h"
#include "netcdf.h"

int main(int argc, char* argv[])
{
	/* variable declarations */
	stnmetahdr_struct metahdr;
	stnmeta_struct metastr;
	daymet_prcp_struct prcpstr;
	imghdr_struct maskhdr, demhdr;
	lsthdr_struct listhdr;
	dayouthdr_struct outhdr;

	file ini;
	file mask_f;
	file dem_f;
	file list_f;
	file ct_f, id_f, wt_f;
	file p_out_f;
	file p_ttot_f;
	file p_tavg_f;
	file p_tfrq_f;
	file p_stot_f;
	file p_savg_f;
	file p_sfrq_f;
	file pnoz_out_f;
	file pnoz_ttot_f;
	file pnoz_tavg_f;
	file pnoz_tfrq_f;
	file pnoz_stot_f;
	file pnoz_savg_f;
	file pnoz_sfrq_f;
	file mx_tavg_f;
	file mx_savg_f;
	file my_tavg_f;
	file my_savg_f;
	file mz_tavg_f;
	file mz_savg_f;
	
	char endkey[8];
	char iocase[16];
	char round[16];
	char ioc,map;
	char outprefix[80];
	double *p_ts;
	double *pnoz_ts;
	// unsigned char brange,bmin;
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
	int nwetdays, nnozwetdays;
	int nwetstns, wet;
	int test;
	long int good_value;

	long int off;
	long int obs_off, obs_offset;

	int imapx,imapy;
	double ulx,uly,cellsize;
	double p, pnoz;
	double p_ttot, p_tavg, p_tfrq;
	double pnoz_ttot, pnoz_tavg, pnoz_tfrq;
	double mx_tavg;
	double my_tavg;
	double mz_tavg;
	double *p_stot_ts, *p_savg_ts, *p_sfrq_ts;
	double *pnoz_stot_ts, *pnoz_savg_ts, *pnoz_sfrq_ts;
	double *mx_savg_ts;
	double *my_savg_ts;
	double *mz_savg_ts;
	double valmax,valmin,scale_ratio;
	// double bytrange;
	double rn, rn2;
	double pop, pop_crit;
	double *pop_ptr;
	double celldat[5];
	
	daydathdr_struct datahdr;
	dat_struct datastr;
	file stnmeta_f, stndata_f;

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
        int prcpid;
        int ncid_meta;
        int year, year_days, descriptor_length;
        double *elevs, *lats, *lons, *stnx, *stny, *prcp;
	double *stnfx, *stnfy;
        int status;
	size_t index[2];

	/* assuming command-line parameter passing */
	map = 1;
	/* scan next two parameters as the mask and dem filenames */
	strcpy(mask_nc_f.name,argv[1]);
	good_value = atoi(argv[2]);
	strcpy(dem_nc_f.name,argv[3]);
	strcpy(stnmeta_nc_f.name,argv[4]);
	strcpy(stndata_nc_f.name,argv[5]);
	strcpy(ct_f.name,argv[6]);
	strcpy(id_f.name,argv[7]);
	strcpy(wt_f.name,argv[8]);
	smwidth = atoi(argv[9]);
	pop_crit = atof(argv[10]);
	prcpstr.f_max = atof(argv[11]);
	noz = atoi(argv[12]);
	start_day = atoi(argv[13]);
	stop_day = atoi(argv[14]);
	valmax = atof(argv[15]);
	valmin = atof(argv[16]);
	bytmax = atoi(argv[17]);
	bytmin = atoi(argv[18]);
	strcpy(outprefix,argv[19]);

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

        status = nc_inq_varid (ncid_meta, "precip", &prcpid);
        if (!(prcp = (double*) malloc(meta_nc_hdr.nstns * year_days * sizeof(double)))) exit(1);
        status = nc_get_var_double (ncid_meta, prcpid, prcp);

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
               	if (!test) 
                {
                        printf("Header information differs between mask and DEM...exiting\n");
                        exit(1);
                }
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
	nstns = prcpstr.nstns = meta_nc_hdr.nstns;
	ndays = prcpstr.ndays = year_days;

	prcpstr.noz = (char)noz;

	/* fill the daily output file header structrue */
	/* determine the yearday and year for the start of output */
	if (yearday(0, year, start_day,
		&outhdr.start_yday, &outhdr.start_year)) 
	{
		printf("Error in call to yearday()...exiting\n");
		exit(1);
	}
	strcpy(outhdr.tag, PRCPTAG);
	outhdr.valmax = valmax;
	outhdr.valmin = valmin;
	outhdr.bytmax = bytmax;
	outhdr.bytmin = bytmin;
	outhdr.ndays = ndays;
	outhdr.ioflag = map;
	if (map)
	{
		outhdr.maskhdr = mask_nc_hdr;
		//populate listhdr with consistant information
		strcpy(listhdr.tag,"0123456\0");
		listhdr.ncells=0;
		listhdr.nvals=0;
		listhdr.typesize=0;
		gctp_struct gctp;
		gctp.outsys=0.0;
		gctp.outzone=0.0;
		gctp.outdatum=0.0;
		double array[15]={0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,\
			0.0,0.0,0.0,0.0,0.0,0.0};
		memccpy(gctp.outparm, array, 15, sizeof(double));
		strcpy(gctp.fn27,"false\0");
		strcpy(gctp.fn83,"false\0");
		listhdr.gctp = gctp;
		outhdr.listhdr = listhdr;
	}
	else
	{
		outhdr.listhdr = listhdr;
		//populate maskhdr with consistant information
		istrcpy(maskhdr.tag,"0123456\0");
		maskhdr.ncols=0;
		maskhdr.nrows=0;
		maskhdr.typesize=0.0;
		maskhdr.ulx=0.0;
		maskhdr.uly=0.0;
		maskhdr.cellsize=0.0;
		gctp_struct gctp;
                gctp.outsys=0.0;
                gctp.outzone=0.0;
                gctp.outdatum=0.0;
		double array[15]={0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,\
			0.0,0.0,0.0,0.0,0.0,0.0};
                memccpy(gctp.outparm,array, 15, sizeof(double));
                strcpy(gctp.fn27,"false\0");
                strcpy(gctp.fn83,"false\0");
		maskhdr.gctp = gctp;
		outhdr.maskhdr = maskhdr;
	}
	
	
	/* some other basic error checking */
	if (start_day > stop_day)
	{
		printf("start_day > stop_day ... exiting\n");
		exit(1);
	}
	if (stop_day >= ndays)
	{
		printf("stop_day >= ndays ... exiting\n");
		printf("remember to reference yeardays from base zero (0,1,2...)\n");
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
	
	/* create output filenames and open files for binary write */
	printf("bin_out %s prefix %s suffix %s\n", &p_out_f, outprefix, &p_out_f.ptr);
	if (bin_out(&p_out_f, outprefix, ".dmprcp_out")) exit(1);
	if (bin_out(&p_ttot_f, outprefix, ".dmprcp_ttot")) exit(1);
	if (bin_out(&p_tavg_f, outprefix, ".dmprcp_tavg")) exit(1);
	if (bin_out(&p_tfrq_f, outprefix, ".dmprcp_tfrq")) exit(1);
	if (bin_out(&p_stot_f, outprefix, ".dmprcp_stot")) exit(1);
	if (bin_out(&p_savg_f, outprefix, ".dmprcp_savg")) exit(1);
	if (bin_out(&p_sfrq_f, outprefix, ".dmprcp_sfrq")) exit(1);
	if (bin_out(&mx_tavg_f, outprefix, ".dmprcp_mxtavg")) exit(1);
	if (bin_out(&mx_savg_f, outprefix, ".dmprcp_mxsavg")) exit(1);
	if (bin_out(&my_tavg_f, outprefix, ".dmprcp_mytavg")) exit(1);
	if (bin_out(&my_savg_f, outprefix, ".dmprcp_mysavg")) exit(1);
	if (bin_out(&mz_tavg_f, outprefix, ".dmprcp_mztavg")) exit(1);
	if (bin_out(&mz_savg_f, outprefix, ".dmprcp_mzsavg")) exit(1);
	if (noz)
	{
		if (bin_out(&pnoz_out_f, outprefix, ".dmprcp_nozout")) exit(1);
		if (bin_out(&pnoz_ttot_f, outprefix, ".dmprcp_nozttot")) exit(1);
		if (bin_out(&pnoz_tavg_f, outprefix, ".dmprcp_noztavg")) exit(1);
		if (bin_out(&pnoz_tfrq_f, outprefix, ".dmprcp_noztfrq")) exit(1);
		if (bin_out(&pnoz_stot_f, outprefix, ".dmprcp_nozstot")) exit(1);
		if (bin_out(&pnoz_savg_f, outprefix, ".dmprcp_nozsavg")) exit(1);
		if (bin_out(&pnoz_sfrq_f, outprefix, ".dmprcp_nozsfrq")) exit(1);
	}
		
	/* write headers to daily data output files */
	outhdr.noz = 0;
	if (fwrite(&outhdr, sizeof(dayouthdr_struct), 1, p_out_f.ptr) != 1)
	{
		printf("Error writing header to outfile... exiting\n");
		exit(1);
	}
	if (noz)
	{
		outhdr.noz = 1;
		if (fwrite(&outhdr, sizeof(dayouthdr_struct), 1, pnoz_out_f.ptr) != 1)
		{
			printf("Error writing header to noz outfile... exiting\n");
			exit(1);
		}
	}

	/* allocate array memory */
	if (!(prcpstr.stnobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(prcpstr.stnsmobs = (double*) malloc(nstns * ndays * sizeof(double)))) exit(1);
	if (!(prcpstr.stnx = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(prcpstr.stny = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(prcpstr.stnz = (double*) malloc(nstns * sizeof(double)))) exit(1);
	if (!(p_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(p_stot_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(p_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(p_sfrq_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(mx_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(my_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (!(mz_savg_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	if (noz)
	{
		if (!(pnoz_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
		if (!(pnoz_stot_ts = (double*) malloc(nsimdays * sizeof(double))))	exit(1);
		if (!(pnoz_savg_ts = (double*) malloc(nsimdays * sizeof(double))))	exit(1);
		if (!(pnoz_sfrq_ts = (double*) malloc(nsimdays * sizeof(double)))) exit(1);
	}
	
	/* initialize time series arrays */
	for (i=0 ; i<nsimdays ; i++)
	{
		p_ts[i] =  0;
		p_stot_ts[i] = p_sfrq_ts[i] = 0.0;
		p_savg_ts[i] = mx_savg_ts[i] = my_savg_ts[i] = mz_savg_ts[i] = 0.0;
		if (noz)
		{
			pnoz_ts[i] = 0;
			pnoz_stot_ts[i] = pnoz_savg_ts[i] = pnoz_sfrq_ts[i] = 0.0;
		}
	}
	
	/* read stnmeta data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		if (xyloc_flag) 
		{
			prcpstr.stnx[i] = metastr.fx = stnfx[i];
			prcpstr.stny[i] = metastr.fy = stnfy[i];
		}
		else
		{
			prcpstr.stnx[i] = stnx[i];
			prcpstr.stny[i] = stny[i];
		}
		prcpstr.stnz[i] = elevs[i];
	}

	/* read station daily data and fill appropriate arrays */
	for (i=0 ; i<nstns ; i++)
	{
		obs_off = i*ndays;
		for (j=0 ; j<ndays ; j++)
		{
			obs_offset = obs_off + j;
			prcpstr.stnobs[obs_offset] = prcp[obs_offset];
		}
	}
	
	/* smooth the station data arrays by station (event smoothing algorithm) */
	for (i=0 ; i<nstns ; i++)
	{
		off = i * ndays;
		if (prcp_boxcar_smooth(prcpstr.stnobs+off,prcpstr.stnsmobs+off,
			ndays,smwidth,1))
		{
			printf("Error in call to prcp_boxcar_smooth()... exiting\n");
			exit(1);
		}
	}
	
	/* begin simulation, loop through ncells */
	n_in = 0;
	prcpstr.switch_init = 1;
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
		prcpstr.mapy = uly - (double)imapy*cellsize;
		prcpstr.mapx = ulx + (double)imapx*cellsize;
		prcpstr.elev = (double)z;
		
		/* initialize temporal variables for this cell */
		nwetdays = 0;
		p_ttot = p_tavg = p_tfrq = 0.0;
		mx_tavg = 0.0;
		my_tavg = 0.0;
		mz_tavg = 0.0;
		
		if (noz)
		{
			nnozwetdays = 0;
			pnoz_ttot = pnoz_tavg = pnoz_tfrq = 0.0;
		}
		
		/* read station count for this cell */
		fread(&count, sizeof(short), 1, ct_f.ptr);
		prcpstr.count = count;

		/* allocate memory depending on count */
		prcpstr.nregr = nregr = ((count * count) - count) / 2;
		if (!(prcpstr.id = (short*) malloc(count * sizeof(short)))) exit(1);
		if (!(prcpstr.wt = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(prcpstr.listx = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(prcpstr.listy = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(prcpstr.listz = (double*) malloc(count * sizeof(double)))) exit(1);
		if (!(prcpstr.listobs = (double*) malloc(count*ndays*sizeof(double)))) exit(1);
		if (!(prcpstr.listsmobs = (double*) malloc(count*ndays*sizeof(double)))) exit(1);
		if (!(prcpstr.regx = (double**) malloc(3 * sizeof(double*)))) exit(1);
		for (j=0 ; j<3 ; j++)
		{
			if (!(prcpstr.regx[j] = (double*)malloc(nregr*sizeof(double)))) exit(1);
		}
		if (!(prcpstr.regy = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(prcpstr.regwtall = (double*) malloc(nregr * sizeof(double)))) exit(1);
		if (!(prcpstr.regwt = (double*) malloc(nregr * sizeof(double)))) exit(1);

		/* read station id and wt arrays */
		fread(prcpstr.id, sizeof(short), count, id_f.ptr);
		fread(prcpstr.wt, sizeof(double), count, wt_f.ptr);

		/* fill the station-list arrays */
		if (fill_list(&prcpstr))
		{
			printf("Error filling station list arrays ... exiting\n");
			exit(1);
		}

		/* generate regression arrays that don't vary between days */
		if (prcp_regr_xwt_2switch(&prcpstr))
		{
			printf("Error in prcp_regr_xwt_2switch() ... exiting\n");
			exit(1);
		}

		/* begin loop through nsimdays */
		for (day=0 ; day<nsimdays ; day++)
		{

			/* initialize daily prediction variables */
			p = pnoz = 0.0;
			
			/* generate offset into the station list observations */
			prcpstr.dam_off = (day + start_day) * count;
			
			/* generate precipitation occurrence prediction */
			pop_ptr = prcpstr.listobs+prcpstr.dam_off;
			pop = 0.0;
			nwetstns = 0;
			for (j=0 ; j<count ; j++)
			{
				pop += prcpstr.wt[j] * (wet = (pop_ptr[j] > 0.0));
				nwetstns += wet;
			}
			
			/* if precipitation occurrence threshold exceeded, predict 
			precipitation amount */
			if (pop > pop_crit)
			{
				/* only generate regression stats when nwetstns > 3 */
				if (nwetstns > 3)
				{
					/* generate regression arrays that change by day */
					if (prcp_regr_y_2switch(&prcpstr))
					{
						printf("Error in prcp_regr_y_2switch() ... exiting\n");
						exit(1);
					}
					
					/* weighted multiple least squares fit */
					if (prcp_mlr(&prcpstr))
					{
						printf("Error in wt_regr() ... exiting\n");
						exit(1);
					}
					
				} /* end if nwet > 3 */
				
				else
				{
					prcpstr.coef[0] = prcpstr.coef[1] = prcpstr.coef[2] = 
						prcpstr.coef[3] = 0.0;
				}
			
				/* predict prcp for the day */
				if (predict_prcp_noint(&prcpstr))
				{
					printf("Error in predict_prcp() ... exiting\n");
					exit(1);
				}
				
				if (p = prcpstr.pp)
				{
					nwetdays++;
					p_sfrq_ts[day]++;
					p_stot_ts[day] += p;
					p_ttot += p;
					mx_tavg += prcpstr.coef[0];
					mx_savg_ts[day] += prcpstr.coef[0];
					my_tavg += prcpstr.coef[1];
					my_savg_ts[day] += prcpstr.coef[1];
					mz_tavg += prcpstr.coef[2];
					mz_savg_ts[day] += prcpstr.coef[2];
				}
					
				/* bug fixed 9/16/09, PET: added test for noz in following if block */
				if (noz && (pnoz = prcpstr.ppnoz))
				{
					nnozwetdays++;
					pnoz_sfrq_ts[day]++;
					pnoz_stot_ts[day] += pnoz;
					pnoz_ttot += pnoz;
				}
				
			} /* end if pop > pop_crit */

			/* byte scale daily predictions */
			// if ((p < 0.0) p = 0.0;
			// if (p > bytrange) p = bytrange;
			// sprintf(round, "%.0lf", p);
			// sprintf(round, "%.2f", p);
			// p_ts[day] = (char)atoi(round) + bmin;
			p_ts[day] = p;

			if (noz)
			{
				// if ((pnoz < 0.0) pnoz = 0.0;
				// if (pnoz > bytrange) pnoz = bytrange;
				// sprintf(round, "%.0lf", pnoz);
				// sprintf(round, "%.2f", pnoz);
				// pnoz_ts[day] = (char)atoi(round) + bmin;
				// pnoz_ts[day] = (char)atoi(round);
				pnoz_ts[day] = pnoz;
			}
				
		} /* end of nsimdays loop */
						
		/* free memory depending on count */
		free(prcpstr.id);
		free(prcpstr.wt);
		free(prcpstr.listobs);
		free(prcpstr.listsmobs);
		free(prcpstr.listx);
		free(prcpstr.listy);
		free(prcpstr.listz);
		for (j=0 ; j<3 ; j++)
		{
			free(prcpstr.regx[j]);
		}
		free(prcpstr.regx);
		free(prcpstr.regy);
		free(prcpstr.regwtall);
		free(prcpstr.regwt);
		
		/* write daily data to outfiles */
		fwrite(p_ts, sizeof(double), nsimdays, p_out_f.ptr);
		if (noz)
		{
			fwrite(pnoz_ts, sizeof(double), nsimdays, pnoz_out_f.ptr);
		}
		
		/* finish temporal variable calculations */
		rn2 = (double)nwetdays;
		p_tfrq = rn2 / rn;
		if (rn2)
		{
			p_tavg = p_ttot / rn2;
			mx_tavg /= rn2;
			my_tavg /= rn2;
			mz_tavg /= rn2;
		}
		else
		{
			p_tavg = mx_tavg = my_tavg = mz_tavg = 0.0;
		}
		
		if (noz)
		{
			rn2 = (double)nnozwetdays;
			pnoz_tavg = pnoz_ttot / rn2;
			pnoz_tfrq = rn2 / rn;
		}
		
		/* write temporal averages to outfiles */
		fwrite(&p_ttot, sizeof(double), 1, p_ttot_f.ptr);
		fwrite(&p_tavg, sizeof(double), 1, p_tavg_f.ptr);
		fwrite(&p_tfrq, sizeof(double), 1, p_tfrq_f.ptr);
		fwrite(&mx_tavg, sizeof(double), 1, mx_tavg_f.ptr);
		fwrite(&my_tavg, sizeof(double), 1, my_tavg_f.ptr);
		fwrite(&mz_tavg, sizeof(double), 1, mz_tavg_f.ptr);
		if (noz)
		{
			fwrite(&pnoz_ttot, sizeof(double), 1, pnoz_ttot_f.ptr);
			fwrite(&pnoz_tavg, sizeof(double), 1, pnoz_tavg_f.ptr);
			fwrite(&pnoz_tfrq, sizeof(double), 1, pnoz_tfrq_f.ptr);
		}
		
		/* switch initial sign for next cell */
		prcpstr.switch_init = -prcpstr.switch_init;
		
	} /* end of ncells loop */

	/* calculate spatial averages for time series arrays */
	if (map) rn = (double)n_in;
	if (!map) rn = (double)ncells;
	
	for (i=0 ; i<nsimdays ; i++)
	{
		rn2 = p_sfrq_ts[i];
		if (rn2)
		{
			p_savg_ts[i] = p_stot_ts[i] / rn2;
			mx_savg_ts[i] /= rn2;
			my_savg_ts[i] /= rn2;
			mz_savg_ts[i] /= rn2;
		}
		else
		{
			p_savg_ts[i] = mx_savg_ts[i] = my_savg_ts[i] = mz_savg_ts[i] = 0.0;
		}
		p_sfrq_ts[i] /= rn;
		p_stot_ts[i] /= rn;
		
		if (noz)
		{
			rn2 = pnoz_sfrq_ts[i];
			if (rn2)
			{
				pnoz_savg_ts[i] = pnoz_stot_ts[i] / rn2;
			}
			else
			{
				pnoz_savg_ts[i] = 0.0;
			}
			pnoz_stot_ts[i] /= rn;
			pnoz_sfrq_ts[i] /= rn;
		}
	}
		
	/* write time series arrays to outfiles */
	fwrite(p_stot_ts, sizeof(double), nsimdays, p_stot_f.ptr);
	fwrite(p_savg_ts, sizeof(double), nsimdays, p_savg_f.ptr);
	fwrite(p_sfrq_ts, sizeof(double), nsimdays, p_sfrq_f.ptr);
	fwrite(mx_savg_ts, sizeof(double), nsimdays, mx_savg_f.ptr);
	fwrite(my_savg_ts, sizeof(double), nsimdays, my_savg_f.ptr);
	fwrite(mz_savg_ts, sizeof(double), nsimdays, mz_savg_f.ptr);
	if (noz)
	{
		fwrite(pnoz_stot_ts, sizeof(double), nsimdays, pnoz_stot_f.ptr);
		fwrite(pnoz_savg_ts, sizeof(double), nsimdays, pnoz_savg_f.ptr);
		fwrite(pnoz_sfrq_ts, sizeof(double), nsimdays, pnoz_sfrq_f.ptr);
	}
	
	/* close output files */
	fclose(p_out_f.ptr);
	fclose(p_ttot_f.ptr);
	fclose(p_tavg_f.ptr);
	fclose(p_tfrq_f.ptr);
	fclose(p_stot_f.ptr);
	fclose(p_savg_f.ptr);
	fclose(p_sfrq_f.ptr);
	fclose(mx_tavg_f.ptr);
	fclose(mx_savg_f.ptr);
	fclose(my_tavg_f.ptr);
	fclose(my_savg_f.ptr);
	fclose(mz_tavg_f.ptr);
	fclose(mz_savg_f.ptr);
	if (noz)
	{
		fclose(pnoz_out_f.ptr);
		fclose(pnoz_ttot_f.ptr);
		fclose(pnoz_tavg_f.ptr);
		fclose(pnoz_tfrq_f.ptr);
		fclose(pnoz_stot_f.ptr);
		fclose(pnoz_savg_f.ptr);
		fclose(pnoz_sfrq_f.ptr);
	}
	fclose(ct_f.ptr);
	fclose(id_f.ptr);
	fclose(wt_f.ptr);
	
	/* free remaining memory */
	free(prcpstr.stnobs);
	free(prcpstr.stnsmobs);
	free(prcpstr.stnz);
	free(p_ts);
	free(p_stot_ts);
	free(p_savg_ts);
	free(p_sfrq_ts);
	free(mx_savg_ts);
	free(my_savg_ts);
	free(mz_savg_ts);
	if (noz)
	{
		free(pnoz_ts);
		free(pnoz_stot_ts);
		free(pnoz_savg_ts);
		free(pnoz_sfrq_ts);
	}

	return(0);
	
} /* end of main */
