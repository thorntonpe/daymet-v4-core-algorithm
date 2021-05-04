/*
xval_daily_prcp_func3.c
Peter Thornton
11/4/96

Holds functions called from prcp cross-validation code.  Typically these
are functions that have to include checking for the exclusion of the
cross-validation station, or for exclusion of missing data values.
Updated:
9/2/00, PET: added support for three-dimensional regression.
*/

#include "metsrc2.h"

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
		printf("fewer than 2 stations in list\n");
		ok = 0;
	}
	
	return (!ok);
}

int xv_fill_list(char* stmgood, char* listgood, daymet_prcp_struct* str);
int xv_fill_list(char* stmgood, char* listgood, daymet_prcp_struct* str)
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


int xv_prcp_regr_y_2switch(char* listgood, daymet_prcp_struct* str);
int xv_prcp_regr_y_2switch(char* listgood, daymet_prcp_struct* str)
/* make time variant regression array */
/* generates the y variable for daily regressions */
/* this version employs a switching mechanism to keep differencing unbiased */
/* this version uses a double switching mechanism, alternates initial sign */
/* PET (10/18/2015: adding always_good option since this is now being called with only good data */
{
	int ok=1;
	int nstns;
	int sign_init, sign;
	int stni, stnj, n=0;
	double* prcp;
	double prcp1, prcp2;
	char* good;
	int always_good = 1;
	
	nstns = str->count;
	sign_init = str->switch_init;
	prcp = str->listsmobs+str->dam_off;
	good = listgood+str->dam_off;
	
	for (stni=0 ; stni<nstns-1 ; stni++)
	{
		sign = sign_init;
		prcp1 = prcp[stni];
		if (prcp1 && always_good)
		{
			for (stnj=stni+1 ; stnj<nstns ; stnj++)
			{
				prcp2 = prcp[stnj];
				if (prcp2 && always_good)
				{
					/* str->regy[n] = sign * (prcp1-prcp2)/(prcp1+prcp2); */
					str->regy[n] = sign * (prcp1-prcp2);
					str->regwt[n] = str->regwtall[n];
				}
				else /* !prcp2 or missing data */
				{
					str->regy[n] = str->regwt[n] = 0.0;
				}
				n++;
				sign = -sign;
				
			} /* end for stnj */
			
		} /* end if (prcp1) */
		
		else /* !prcp1 or missing data */
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


int xv_predict_prcp_noint(char* listgood, daymet_prcp_struct* str)
/* generates predicted prcp, with and without elevation correction 
assumes the intercept from regression is 0.0 
Uses a list of good data indicators for cross-validation treatment of
missing data */
/* PET (10/18/2015): added always_good option since this now is called with good data always */
{
	int ok=1;
	int nstns;
	int stni;
	double p1, pc, p, sumwt;
	double* prcp;
	double mx,my,mz,x,y,z;
	double f, f_max;
	double max_stn_day_prcp;
	double max_factor = 2.0;
	char *good;
	int always_good = 1;
	
	mx = str->coef[0];
	my = str->coef[1];
	mz = str->coef[2];
	x = str->mapx;
	y = str->mapy;
	z = str->elev;
	f_max = str->f_max;
	nstns = str->count;
	prcp = str->listobs+str->dam_off;
	good = listgood+str->dam_off;
	pc = p = sumwt = 0.0;
	
	max_stn_day_prcp = 0.0;
	for (stni=0 ; stni<nstns ; stni++)
	{
		if ((p1 = prcp[stni]) && always_good)
		{
			pc += str->wt[stni]*(p1 + mx*(x-str->listx[stni]) + 
				my*(y-str->listy[stni]) + mz*(z-str->listz[stni])); 
			sumwt += str->wt[stni];
			/*
			f = mx*(x-str->listx[stni]) + my*(y-str->listy[stni]) + mz*(z - str->listz[stni]);
			if (f > f_max) f = f_max;
			if (f < -f_max) f = -f_max;
			pc += str->wt[stni] * p1 * (1.0 + f) / (1.0 - f);
			sumwt += str->wt[stni];
			*/
			
			if (prcp[stni] > max_stn_day_prcp) max_stn_day_prcp = prcp[stni];
		}
	}
	
	if (sumwt)
	{
		pc /= sumwt;
	}
	else
	{
		printf("predict_prcp() sumwt = 0.0 ... exiting\n");
		ok=0;
	}
	
	/* set some constraints on daily precip */
	/* PET (10/18/2015): reset the upper limit from 20 cm/day to 40 cm/day */
	if (pc < 0.0) pc = 0.0;
// 	if (pc > 40.0) pc = 40.0;
    
    // constrain to range of obs values, times a fixed factor
    if (pc > max_stn_day_prcp * max_factor) pc = max_stn_day_prcp * max_factor;

	str->pp = pc;
	
	return (!ok);
}


/* xv_prcp_boxcar_smooth()  generates smoothed values
from non-zero observations in the smoothing window, and returns a value which
reflects the average value per non-zero element.
For x-validation code, this function includes a check for missing data, and
requires an additional array as input, consisting of character data, 
1 for good values, 0 for missing values, corresponding to the actual input
data value array.
10/25/2015, PET: Switching to always_good treatment, since for V3 code the 
arrays are all filled with good values before xval is called.
*/
   
int xv_prcp_boxcar_smooth(char* good, double* input,double* output,int n,int w,int w_flag);
int xv_prcp_boxcar_smooth(char* good, double* input,double* output,int n,int w,int w_flag)
{
	int ok=1;
    int tail,i,j;
    int *wt;
    double p,total,sum;
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
    	printf("Allocation error in prcp_boxcar... Exiting\n");
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
	            	if ((p = input[i+j-tail]) && always_good)
	            	{
		                total += p * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if ((i >= tail) && (i < n-tail))
	        {
	            for (j=0 ; j<w ; j++)
	            {
	            	if ((p = input[i+j-tail]) && always_good)
	            	{
	            		total += p * wt[j];
	                	sum += (double) wt[j];
	                }
	            }
	        }
	        else if (i >= n-tail)
	       	{
	            for (j=0 ; j<tail+n-i ; j++)
	            {
	            	if ((p = input[i+j-tail]) && always_good)
	            	{
	                	total += p * wt[j];
	                	sum += (double) wt[j];
	                }
	            }
	        }
	        if (sum)
	        	output[i] = total/sum;
	        else
	        	output[i] = 0.0;
	        	
	    } /* end for i = n elements */
	
		free(wt);
		    
	} /* end if (ok) */
	
	
    return (!ok);
} 

int fill_prcp_boxcar_smooth(char* good, double* input,double* output,int n,int w,int w_flag);
int fill_prcp_boxcar_smooth(char* good, double* input,double* output,int n,int w,int w_flag)
{
	int ok=1;
    int tail,i,j;
    int *wt;
    double p,total,sum;
    
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
    	printf("Allocation error in prcp_boxcar... Exiting\n");
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
	            	if ((p = input[i+j-tail]) && good[i+j-tail])
	            	{
		                total += p * wt[j];
		                sum += (double) wt[j];
		            }
	            }
	        }
	        else if ((i >= tail) && (i < n-tail))
	        {
	            for (j=0 ; j<w ; j++)
	            {
	            	if ((p = input[i+j-tail]) && good[i+j-tail])
	            	{
	            		total += p * wt[j];
	                	sum += (double) wt[j];
	                }
	            }
	        }
	        else if (i >= n-tail)
	       	{
	            for (j=0 ; j<tail+n-i ; j++)
	            {
	            	if ((p = input[i+j-tail]) && good[i+j-tail])
	            	{
	                	total += p * wt[j];
	                	sum += (double) wt[j];
	                }
	            }
	        }
	        if (sum)
	        	output[i] = total/sum;
	        else
	        	output[i] = 0.0;
	        	
	    } /* end for i = n elements */
	
		free(wt);
		    
	} /* end if (ok) */
	
	
    return (!ok);
} 
