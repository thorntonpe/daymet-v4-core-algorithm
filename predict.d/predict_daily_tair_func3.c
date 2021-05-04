/* 
predict_daily_tair_func.c
Peter Thornton

Part of the metsrc code family. This file contains functions called from
predict_daily_tair.c

Revision history:
11/11/96: revised from daymet_tair_func.c, to include the metsrc structures
9/1/00: revised from predict_daiy_tair_func.c, added support for 3-variable
 multiple linear regression model.
*/
/* $$$NCW Removed for non-AIX compatitibility  */
#ifdef _AIX
#include <fp.h>
#endif
/*#include "../library.d/fileio_tools.h"  */
#include "metsrc2.h"
/*#include "predict_daily_tair_func3.h"  */

// int init_byte_scale(double valmax, double valmin, int bytmax, int bytmin, 
// 	double* scale_ratio, unsigned char* brange, unsigned char* bmin);
// int init_byte_scale(double valmax, double valmin, int bytmax, int bytmin, 
// 	double* scale_ratio, unsigned char* brange, unsigned char* bmin)
// /* initializes byte scaling parameters and checks for range errors */
// {
// 	int ok=1;
	
// 	if (ok && (valmax <= valmin))
// 	{
// 		printf("Error: max value <= min value ... exiting\n");
// 		ok=0;
// 	}
	
// 	if (ok && (bytmax <= bytmin))
// 	{
// 		printf("Error: max byte value <= min byte value ... exiting\n");
// 		ok=0;
// 	}
	
// 	if (ok && (bytmin < 0))
// 	{
// 		printf("Error: min byte value < 0 ... exiting\n");
// 		ok=0;
// 	}
	
// 	if (ok && (bytmax > 255))
// 	{
// 		printf("Error: max byte value > 255 ... exiting\n");
// 		ok=0;
// 	}
	
// 	*brange = (unsigned char)(bytmax - bytmin);
// 	*scale_ratio = (double)(*brange)/(valmax-valmin);
// 	*bmin = (unsigned char)bytmin;
	
// 	return (!ok);
// }

/* boxcar_smooth() performs a windowed smoothing on the input array, returns
result in output array. Both arrays must be doubles. n=array length,
w = windowing width, w_flag (0=flat boxcar, 1=ramped boxcar, e.g. [1 2 3 2 1])
*/

// int boxcar_smooth(double* input, double* output, int n, int w, int w_flag);
int boxcar_smooth(double* input, double* output, int n, int w, int w_flag)
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
    	printf("Allocation error in boxcar_smooth... Exiting\n");
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
	                total += input[i+j-tail] * wt[j];
	                sum += (double) wt[j];
	            }
	        }
	        else if ((i >= tail) && (i < n-tail))
	       	{
	            for (j=0 ; j<w ; j++)
	            {
	                total += input[i+j-tail] * wt[j];
	                sum += (double) wt[j];
	            }
	        }
	        else if (i >= n-tail)
	        {
	            for (j=0 ; j<tail+n-i ; j++)
	            {
	                total += input[i+j-tail] * wt[j];
	                sum += (double) wt[j];
	            }
	        }
	        output[i] = total/sum;
	        
	    } /* end for i=nelements */
	    
		free(wt);
		
	} /* end if ok */
	
	return (!ok);
}   

int fill_list(daymet_tair_struct* str);
int fill_list(daymet_tair_struct* str)
/* fills the station-major list arrays (station-major -> day-major) */
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
		}
		str->listx[i] = str->stnx[str->id[i]];
		str->listy[i] = str->stny[str->id[i]];
		str->listz[i] = str->stnz[str->id[i]];
	}
	
	return (!ok);
}

int tair_regr_xwt_2switch(daymet_tair_struct* str);
int tair_regr_xwt_2switch(daymet_tair_struct* str)
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
	
	double **x;
	double *w;
	double x1m,x2m,x3m,ym,wm;
	double sx1,sx2,sx3,ssx1,ssx2,ssx3;
	double sx1x2,sx1x3,sx2x3;
	
	double dx1,dx2,dx3;
	double dn,wt;
		
	double r[3];
	double a,b,c,aa,bb,cc,det;
	double inv[3][3];
	double sigx[3];
	double t1;
	
	int i,j,npts;
	
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
			str->regwt[n] = wti * str->wt[stnj];
			n++;
			sign = -sign;
		}
	}
	
	/* many of the regression steps can be done in the time-invariant step */
	x=str->regx;
	w=str->regwt;
	dn = (double)str->nregr;
	npts = str->nregr;
	x1m=x2m=x3m=wm=0.0;
	sx1=sx2=sx3=ssx1=ssx2=ssx3=0.0;
	sx1x2=sx1x3=sx2x3=0.0;
	
	/* step 1 : calculate means for x1,x2,x3,w */
	for (i=0 ; i<npts ; i++)
	{
		x1m += w[i]*x[0][i];
		x2m += w[i]*x[1][i];
		x3m += w[i]*x[2][i];
		wm += w[i];
	}
	x1m /= wm;
	x2m /= wm;
	x3m /= wm;
	str->x1m=x1m;
	str->x2m=x2m;
	str->x3m=x3m;
	str->wm1=wm;
	wm /= dn;
	str->wm=wm;
	
	/* step 2 : calculate sums */
	for (i=0 ; i<npts ; i++)
	{
		wt=w[i]/wm;
		str->regdx[0][i] = dx1=x[0][i]-x1m;
		str->regdx[1][i] = dx2=x[1][i]-x2m;
		str->regdx[2][i] = dx3=x[2][i]-x3m;
		
		sx1 += wt*dx1;
		ssx1 += wt*dx1*dx1;
		sx2 += wt*dx2;
		ssx2 += wt*dx2*dx2;
		sx3 += wt*dx3;
		ssx3 += wt*dx3*dx3;
		sx1x2 += wt*dx1*dx2;
		sx1x3 += wt*dx1*dx3;
		sx2x3 += wt*dx2*dx3;
	}	
	str->ssx1=ssx1;
	str->ssx2=ssx2;
	str->ssx3=ssx3;
	
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
	str->inv[0][0]=(cc-1.0)/det;
	str->inv[1][1]=(bb-1.0)/det;
	str->inv[2][2]=(aa-1.0)/det;
	/* next the off-diagonals */
	str->inv[0][1] = str->inv[1][0] = (a-c*b)/det;
	str->inv[0][2] = str->inv[2][0] = (b-a*c)/det;
	str->inv[1][2] = str->inv[2][1] = (c-a*b)/det;
	
	/* step 7 : calculate SIGMAX */
	str->sigx[0]=sqrt(ssx1/(wm*(n-1)));
	str->sigx[1]=sqrt(ssx2/(wm*(n-1)));
	str->sigx[2]=sqrt(ssx3/(wm*(n-1)));
	
	return (!ok);
}

int tair_regr_y_2switch(daymet_tair_struct* str);
int tair_regr_y_2switch(daymet_tair_struct* str)
/* make time variant regression array */
/* generates the y variable for daily regressions */
/* this version employs a switching mechanism to keep differencing unbiased */
/* this version uses a double switching mechanism, alternates initial sign */
{
	int ok=1;
	int nstns;
	int sign_init, sign;
	int stni, stnj, n=0;
	double* tair;
	double tairi;
	
	nstns = str->count;
	sign_init = str->switch_init;
	tair = str->listsmobs+str->dam_off;
	
	for (stni=0 ; stni<nstns ; stni++)
	{
		sign = sign_init;
		tairi = tair[stni];
		for (stnj=stni+1 ; stnj<nstns ; stnj++)
		{
			str->regy[n] = sign * (tairi - tair[stnj]);
			n++;
			sign = -sign;
		}
	}
	
	return (!ok);
}

int mlr(daymet_tair_struct* str);
int mlr(daymet_tair_struct* str)
{
	int ok=1;
	double **x;
	double *y;
	double *w;
	double ym,wm;
	double ssy;
	double sx1y,sx2y,sx3y;
	
	double dx1,dx2,dx3,dy;
	double n,wt;
		
	double r[3];
	double coef[4];
	double **inv;
	double *sigx, sigy;
	double t1;
	
	int i,j,npts;
	
	x=str->regx;
	y=str->regy;
	w=str->regwt;
	inv=str->inv;
	sigx=str->sigx;
	n = (double)str->nregr;
	npts = str->nregr;
	ym=wm=0.0;
	ssy=0.0;
	sx1y=sx2y=sx3y=0.0;
	
	/* step 1 : calculate mean for y */
	for (i=0 ; i<npts ; i++)
	{
		ym += w[i]*y[i];
	}
	ym /= str->wm1;
	wm = str->wm;
	
	/* step 2 : calculate sums */
	for (i=0 ; i<npts ; i++)
	{
		wt=w[i]/wm;
		dx1=str->regdx[0][i];
		dx2=str->regdx[1][i];
		dx3=str->regdx[2][i];
		dy=y[i]-ym;
		
		ssy += wt*dy*dy;
		sx1y += wt*dx1*dy;
		sx2y += wt*dx2*dy;
		sx3y += wt*dx3*dy;
	}	
	
	/* step 3 : calculate R */
	r[0]=sx1y/sqrt(str->ssx1*ssy);
	r[1]=sx2y/sqrt(str->ssx2*ssy);
	r[2]=sx3y/sqrt(str->ssx3*ssy);
	
	/* step 6 : calculate SIGMAY */
	sigy = sqrt(ssy/(wm*(n-1)));
	
	/* step 8 : calculate coef = R # (inverse of ARRAY) * (SIGMAY/SIGMAX) */
	for (i=0 ; i<3 ; i++)
	{
		coef[i]=0.0;
		for (j=0 ; j<3 ; j++)
		{
			coef[i] += r[j]*inv[i][j];
		}
		coef[i] *= sigy/sigx[i];
		str->coef[i] = coef[i];
	}
	
	/* step 9 : calculate the constant term (returned in coef[3]) */
	str->coef[3] = ym - (coef[0]*str->x1m + coef[1]*str->x2m + coef[2]*str->x3m);

	/* test for NANQ in coefficients */
	for (i=0 ; i<3 ; i++)
	{
#ifdef _AIX
        	if (IS_NAN(str->coef[i]))
#else
        // $$$NCW removed IS_NAN for non-AIX compatibility
                if (isnan(str->coef[i]))
#endif
		{
			for (j=0 ; j<3 ; j++)
			{
				str->coef[j] = 0.0;
			}
		}
	}

	/* set some constraints on z coefficient (elevation) */
	if (str->coef[2] > 0.001) str->coef[2] = 0.001;
	if (str->coef[2] < -0.012) str->coef[2] = -0.012;
	
	
	return (!ok);
}

int predict_tair_noint(daymet_tair_struct* str);
int predict_tair_noint(daymet_tair_struct* str)
/* generates predicted tair, with and without regression corrections 
   no regression intercept (assumed = 0.0) */
{
	int ok=1;
	int nstns;
	int stni;
	double tc, t, sumwt;
	double* tair;
	double mx,my,mz,x,y,z;
	double max_stn_day_tair;
	double tair_factor = 10.0;
	
        //BWM
	mx = str->coef[0];
	my = str->coef[1];
	mz = str->coef[2];

	x = str->mapx;
	y = str->mapy;
	z = str->elev;
	nstns = str->count;
	tair = str->listobs+str->dam_off;
	tc = t = sumwt = 0.0;
	
	max_stn_day_tair = 0.0;
	for (stni=0 ; stni<nstns ; stni++)
	{
		tc += str->wt[stni]*(tair[stni] + mx*(x-str->listx[stni]) + 
			my*(y-str->listy[stni]) + mz*(z-str->listz[stni])); 
		if (str->noz)
		{
			t += str->wt[stni]*tair[stni];
		}
		sumwt += str->wt[stni];

		if (tair[stni] > max_stn_day_tair) max_stn_day_tair = tair[stni];
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
