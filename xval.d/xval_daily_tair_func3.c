/*
xval_daily_tair_func3.c
Peter Thornton
11/11/96
*/

#ifdef _AIX
#include <fp.h>
#endif

#include "metsrc2.h"

/* xv_boxcar_smooth() performs a windowed smoothing on the input array, returns
result in output array. Both arrays must be doubles. n=array length,
w = windowing width, w_flag (0=flat boxcar, 1=ramped boxcar, e.g. [1 2 3 2 1])
For x-validation code, this function includes a check for missing data, and
requires an additional array as input, consisting of character data, 
1 for good values, 0 for missing values, corresponding to the actual input
data value array.

Updated:
9/2/00, PET: added support for three-dimensional regression.
*/
int xv_boxcar_smooth(char* good, double* input, double* output, int n, int w, int w_flag);
int xv_boxcar_smooth(char* good, double* input, double* output, int n, int w, int w_flag)
{
	int ok=1;
    int tail,i,j;
    int* wt;
    double total,sum;
    int always_good = 1;

    if (ok && (w > n/2))
    {
        printf("Boxcar window longer than 1/2 array length...\n");
        printf("Resize window and try again\n");
        ok=0;
    }

    /* establish the lengths of the boxcar tails */
    if (ok)
    {
	    if (!(w % 2))
	        w += 1;
	    tail = w/2;
	}
	
    if (ok && (!(wt = (int*) malloc(w * sizeof(int)))))
    {
    	printf("Allocation error in tair_boxcar... Exiting\n");
    	ok=0;
    }

    /* when w_flag != 0, use linear ramp to weight tails, 
    otherwise use constant weight */
	if (ok)
	{
	    if (w_flag)
	    {
	        for (i=0 ; i<tail ; i++)
	            wt[i] = i+1;
	        for (i=0 ; i<= tail ; i++)
	            wt[i+tail] = tail + 1 - i;
	    }
	    else
	        for (i=0 ; i<w ; i++)
	            wt[i] = 1;

	   
	    for (i=0 ; i<n ; i++)
	    {
	        total = 0.0;
	        sum = 0.0;
	        if (i < tail)
	       	{
	            for (j=tail - i ; j<w ; j++)
	            {
	            	if (always_good)
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if ((i >= tail) && (i < n-tail))
	       	{
	            for (j=0 ; j<w ; j++)
	            {
	            	if (always_good)
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if (i >= n-tail)
	        {
	            for (j=0 ; j<tail+n-i ; j++)
	            {
	            	if (always_good)
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
		        }
	        }
	        output[i] = total/sum;
	        
	    } /* end for i=nelements */
	    
	    free(wt);
	    
	} /* end if ok */
	
    return (!ok);
}   

int fill_boxcar_smooth(char* good, double* input, double* output, int n, int w, int w_flag);
int fill_boxcar_smooth(char* good, double* input, double* output, int n, int w, int w_flag)
{
	int ok=1;
    int tail,i,j;
    int* wt;
    double total,sum;

    if (ok && (w > n/2))
    {
        printf("Boxcar window longer than 1/2 array length...\n");
        printf("Resize window and try again\n");
        ok=0;
    }

    /* establish the lengths of the boxcar tails */
    if (ok)
    {
	    if (!(w % 2))
	        w += 1;
	    tail = w/2;
	}
	
    if (ok && (!(wt = (int*) malloc(w * sizeof(int)))))
    {
    	printf("Allocation error in tair_boxcar... Exiting\n");
    	ok=0;
    }

    /* when w_flag != 0, use linear ramp to weight tails, 
    otherwise use constant weight */
	if (ok)
	{
	    if (w_flag)
	    {
	        for (i=0 ; i<tail ; i++)
	            wt[i] = i+1;
	        for (i=0 ; i<= tail ; i++)
	            wt[i+tail] = tail + 1 - i;
	    }
	    else
	        for (i=0 ; i<w ; i++)
	            wt[i] = 1;

	   
	    for (i=0 ; i<n ; i++)
	    {
	        total = 0.0;
	        sum = 0.0;
	        if (i < tail)
	       	{
	            for (j=tail - i ; j<w ; j++)
	            {
	            	if (good[i+j-tail])
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if ((i >= tail) && (i < n-tail))
	       	{
	            for (j=0 ; j<w ; j++)
	            {
	            	if (good[i+j-tail])
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if (i >= n-tail)
	        {
	            for (j=0 ; j<tail+n-i ; j++)
	            {
	            	if (good[i+j-tail])
	            	{
		                total += input[i+j-tail] * wt[j];
		                sum += (double) wt[j];
		            }
		        }
	        }
	        output[i] = total/sum;
	        
	    } /* end for i=nelements */
	    
	    free(wt);
	    
	} /* end if ok */
	
    return (!ok);
}   

/* generates weighted station list for cross-validation routines */
int xv_weight_list(interpolate_struct* str, int curstn);
int xv_weight_list(interpolate_struct* str, int curstn)
{
	int ok=1;
	short i, nstns, count;
	double sqd, sqrad, inv_sqrad, trunc, gsp;
	double wtsum;	
	double *sqdist;
	
	nstns = (short)str->nstns;
	sqdist = str->sqdist;
	sqrad = str->sqrad;
	gsp = str->gsp;
	trunc = str->trunc;
	inv_sqrad = 1.0/sqrad;	

	count = 0;
	wtsum = 0.0;
	
	for (i=0 ; i<nstns ; i++)
	{
		if ((i != curstn) && (sqd = sqdist[i]) < sqrad)
		{
			str->id[count] = i;
			str->wt[count] = exp(-gsp * sqd * inv_sqrad) - trunc;
			wtsum += str->wt[count];
			count++;
		}
	}
	
	/* normalize weights */
	if (count > 1)
	{
		for (i=0 ; i<count ; i++)
		{
			str->wt[i] /= wtsum;
		}
		str->count = count;
	}
	else
	{
		printf("less than 2 stations in list\n");
		ok = 0;
	}
	
	return (!ok);
}

int xv_fill_list(char* stmgood, char* listgood, daymet_tair_struct* str);
int xv_fill_list(char* stmgood, char* listgood, daymet_tair_struct* str)
/* fills the station-major list arrays (station-major -> day-major) 
includes an array of good value flags for the cross-validation routine */
{
	int ok=1;
	int i,j;
	int ndays,nstns;
	long int stm_off;
	
	ndays = str->ndays;
	nstns = str->count;
	
	for (i=0 ; i<nstns ; i++)
	{
		stm_off = str->id[i] * ndays;
		for (j=0 ; j<ndays ; j++)
		{
			str->listobs[j*nstns + i] = str->stnobs[stm_off + j];
			str->listsmobs[j*nstns + i] = str->stnsmobs[stm_off + j];
			listgood[j*nstns + i] = stmgood[stm_off + j];
		}
		str->listx[i] = str->stnx[str->id[i]];
		str->listy[i] = str->stny[str->id[i]];
		str->listz[i] = str->stnz[str->id[i]];
	}
	
	return (!ok);
}

int xv_tair_regr_xwt_2switch(daymet_tair_struct* str);
int xv_tair_regr_xwt_2switch(daymet_tair_struct* str)
/* make time invariant regression arrays */
/* generates the x and wt variables once that are constant between days */
/* this version employs a switching mechanism to keep differencing unbiased */
/* this version uses a double switching mechanism, alternates initial sign */
{
	int ok=1;
	int nstns;
	int sign_init, sign;
	int stni, stnj, n=0;
	double xi,yi,elevi,wti;
	
	nstns = str->count;
	sign_init = str->switch_init;
	
	for (stni=0 ; stni<nstns ; stni++)
	{
		sign = sign_init;
		xi = str->listx[stni];
		yi = str->listy[stni];
		elevi = str->listz[stni];
		wti = str->wt[stni];
		for (stnj=stni+1 ; stnj<nstns ; stnj++)
		{
			str->regx[0][n] = sign * (xi - str->listx[stnj]);
			str->regx[1][n] = sign * (yi - str->listy[stnj]);
			str->regx[2][n] = sign * (elevi - str->listz[stnj]);
			str->regwtall[n] = wti * str->wt[stnj];
			n++;
			sign = -sign;
		}
	}
	
	return (!ok);
}

int xv_tair_regr_y_2switch(char* listgood, daymet_tair_struct* str);
int xv_tair_regr_y_2switch(char* listgood, daymet_tair_struct* str)
/* make time variant regression array */
/* generates the y variable for daily regressions */
/* this version employs a switching mechanism to keep differencing unbiased */
/* this version uses a double switching mechanism, alternates initial sign */
/* 9/20/2015, PET: replacing the check of good with always good, since this is called
   after the fill_missing_values call.
*/
{
	int ok=1;
	int nstns;
	int sign_init, sign;
	int stni, stnj, n=0;
	double* tair;
	double tair1, tair2;
	char* good;
	int always_good = 1;
	
	nstns = str->count;
	sign_init = str->switch_init;
	tair = str->listsmobs+str->dam_off;
	good = listgood+str->dam_off;
	for (stni=0 ; stni<nstns-1 ; stni++)
	{
		sign = sign_init;
		tair1 = tair[stni];
		if (always_good)
		{
			for (stnj=stni+1 ; stnj<nstns ; stnj++)
			{
				tair2 = tair[stnj];
				if (always_good)
				{
					str->regy[n] = sign * (tair1-tair2);
					str->regwt[n] = str->regwtall[n];
				}
				else /* missing data */
				{
					str->regy[n] = str->regwt[n] = 0.0;
				}
				n++;
				sign = -sign;
				
			} /* end for stnj */
			
		} /* end if (tair1) */
		
		else /* missing data */
		{
			for (stnj=stni+1 ; stnj<nstns ; stnj++)
			{
				str->regy[n] = str->regwt[n] = 0.0;
				n++;
			}
		}
		
	} /* end for stni */
	
	return (!ok);
}

int xv_tair_mlr(daymet_tair_struct* str);
int xv_tair_mlr(daymet_tair_struct* str)
{
	int ok=1;
	
	double **x;
	double *y;
	double *w;
	int npts;
	double *coef;
	double x1m,x2m,x3m,ym,wm;
	double sx1,sx2,sx3,ssx1,ssx2,ssx3,ssy;
	double sx1x2,sx1x3,sx2x3;
	double sx1y,sx2y,sx3y;
	
	double dx1,dx2,dx3,dy;
	double n,wt;
		
	double r[3];
	double a,b,c,aa,bb,cc,det;
	double inv[3][3];
	double sigx[3], sigy;
	double t1;
	
	int i,j;
	
	npts = str->nregr;
	x=str->regx;
	y=str->regy;
	w=str->regwt;
	coef=str->coef;
	
	n = (double)npts;
	x1m=x2m=x3m=ym=wm=0.0;
	sx1=sx2=sx3=ssx1=ssx2=ssx3=ssy=0.0;
	sx1x2=sx1x3=sx2x3=0.0;
	sx1y=sx2y=sx3y=0.0;
	
	/* step 1 : calculate means for x1,x2,x3,y,w */
	for (i=0 ; i<npts ; i++)
	{
		x1m += w[i]*x[0][i];
		x2m += w[i]*x[1][i];
		x3m += w[i]*x[2][i];
		ym += w[i]*y[i];
		wm += w[i];
	}
	x1m /= wm;
	x2m /= wm;
	x3m /= wm;
	ym /= wm;
	wm /= n;
/*
	printf("tair_mlr: %f %f %f %f\n", x1m, x2m, x3m, ym);
*/
	
	/* step 2 : calculate sums */
	for (i=0 ; i<npts ; i++)
	{
		wt=w[i]/wm;
		dx1=x[0][i]-x1m;
		dx2=x[1][i]-x2m;
		dx3=x[2][i]-x3m;
		dy=y[i]-ym;
		
		sx1 += wt*dx1;
		ssx1 += wt*dx1*dx1;
		sx2 += wt*dx2;
		ssx2 += wt*dx2*dx2;
		sx3 += wt*dx3;
		ssx3 += wt*dx3*dx3;
		ssy += wt*dy*dy;
		sx1x2 += wt*dx1*dx2;
		sx1x3 += wt*dx1*dx3;
		sx2x3 += wt*dx2*dx3;
		sx1y += wt*dx1*dy;
		sx2y += wt*dx2*dy;
		sx3y += wt*dx3*dy;
	}	
	
	/* step 3 : calculate R */
	r[0]=sx1y/sqrt(ssx1*ssy);
	r[1]=sx2y/sqrt(ssx2*ssy);
	r[2]=sx3y/sqrt(ssx3*ssy);
	
	/* step 4 : calculate a b and c in diagonal ARRAY */
	a=sx1x2/sqrt(ssx2*ssx1);
	b=sx1x3/sqrt(ssx1*ssx3);
	c=sx2x3/sqrt(ssx2*ssx3);
	aa = a*a;
	bb = b*b;
	cc = c*c;
	det=aa+bb+cc-2.0*a*b*c-1.0;
	
	/* step 5 : calculate inverse of ARRAY */
	/* first the diagonal elements */
	inv[0][0]=(cc-1.0)/det;
	inv[1][1]=(bb-1.0)/det;
	inv[2][2]=(aa-1.0)/det;
	/* next the off-diagonals */
	inv[0][1] = inv[1][0] = (a-c*b)/det;
	inv[0][2] = inv[2][0] = (b-a*c)/det;
	inv[1][2] = inv[2][1] = (c-a*b)/det;
	
	/* step 6 : calculate SIGMAY */
	sigy = sqrt(ssy/(wm*(n-1)));
	
	/* step 7 : calculate SIGMAX */
	sigx[0]=sqrt(ssx1/(wm*(n-1)));
	sigx[1]=sqrt(ssx2/(wm*(n-1)));
	sigx[2]=sqrt(ssx3/(wm*(n-1)));
	
	/* step 8 : calculate coef = R # (inverse of ARRAY) * (SIGMAY/SIGMAX) */
	for (i=0 ; i<3 ; i++)
	{
		coef[i]=0.0;
		for (j=0 ; j<3 ; j++)
		{
			coef[i] += r[j]*inv[i][j];
		}
		coef[i] *= sigy/sigx[i];
	}
	
	/* step 9 : calculate the constant term (returned in coef[3]) */
	coef[3] = ym - (coef[0]*x1m + coef[1]*x2m + coef[2]*x3m);

	/* test for NANQ in coefficients */
	for (i=0 ; i<3 ; i++)
	{
#ifdef _AIX
		if (IS_NAN(coef[i])) 
#else
		if (isnan(coef[i])) 
#endif
		{
			for (j=0 ; j<3 ; j++)
			{
				coef[j] = 0.0;
			}
		}
	}
	
	/* set some constraints on z coefficient (elevation) */
	if (coef[2] > 0.001) coef[2] = 0.001;
	if (coef[2] < -0.012) coef[2] = -0.012;
	
	return(!ok);
}


int xv_predict_tair_noint(char* listgood, daymet_tair_struct* str);
int xv_predict_tair_noint(char* listgood, daymet_tair_struct* str)
/* generates predicted tair, with and without elevation correction 
no regression intercept (assumed = 0.0) 
Uses a list of good data indicators for cross-validation treatment of
missing data */
/* 9/20/2015, PET: Modified to consider all station data to be good, since it has been filled. */
{
	int ok=1;
	int nstns;
	int stni;
	double tc, t, sumwt;
	double* tair;
	double mx,my,mz,x,y,z;
	char *good;
	int always_good = 1;
	double max_stn_day_tair;
	double tair_factor = 10.0;
	
	mx = str->coef[0];
	my = str->coef[1];
	mz = str->coef[2];
	x = str->mapx;
	y = str->mapy;
	z = str->elev;
	nstns = str->count;
	tair = str->listobs+str->dam_off;
	good = listgood+str->dam_off;
	tc = t = sumwt = 0.0;
	
	max_stn_day_tair = 0.0;
	for (stni=0 ; stni<nstns ; stni++)
	{
		if (always_good)
		{
			tc += str->wt[stni]*(tair[stni] + mx*(x-str->listx[stni]) + 
				my*(y-str->listy[stni]) + mz*(z-str->listz[stni])); 
			sumwt += str->wt[stni];

			if (tair[stni] > max_stn_day_tair) max_stn_day_tair = tair[stni];
		}
	}
	
	tc /= sumwt;
	if (tc > max_stn_day_tair + tair_factor) tc = max_stn_day_tair + tair_factor;

	str->ta = tc;
	if (str->noz)
	{
		t /= sumwt;
		if (t > max_stn_day_tair + tair_factor) t = max_stn_day_tair + tair_factor;
		str->tanoz = t;
	}
	
	return (!ok);
}
