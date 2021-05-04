/*
predict_daily_prcp_func3.c
Peter Thornton
8/1/96

Stand-alone functions called from daymet_prcp.c
9/2/00: revised from predict_daiy_prcp_func.c, added support for 3-variable
 multiple linear regression model.
9/21/00: corrected discontinuity in lr check in prcp_mlr.
*/

/* $$$NCW - Removed for compatibility with non-AIX systems */
/*  Review if this should be conditionaly compiled */
#ifdef _AIX
#include <fp.h>
#endif

/* #include "predict_daily_prcp_func3.h" */
#include "metsrc2.h"



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


/* prcp_boxcar_smooth()  generates smoothed values
from non-zero observations in the smoothing window, and returns a value which
reflects the average value per non-zero element.
*/


int prcp_boxcar_smooth(double* input,double* output,int n,int w,int w_flag)
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
	            	if (p = input[i+j-tail])
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
	            	if (p = input[i+j-tail])
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
	            	if (p = input[i+j-tail])
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


int fill_list(daymet_prcp_struct* str)
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


int prcp_regr_xwt_2switch(daymet_prcp_struct* str)
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

	for (stni=0 ; stni<nstns-1 ; stni++)
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


int prcp_regr_y_2switch(daymet_prcp_struct* str)
/* make time variant regression array */
/* generates the y variable for daily regressions */
/* this version employs a switching mechanism to keep differencing unbiased */
/* this version uses a double switching mechanism, alternates initial sign */
{
	int ok=1;
	int nstns;
	int sign_init, sign;
	int stni, stnj, n=0;
	double* prcp;
	double prcp1, prcp2;

	nstns = str->count;
	sign_init = str->switch_init;
	prcp = str->listsmobs+str->dam_off;

	for (stni=0 ; stni<nstns-1 ; stni++)
	{
		sign = sign_init;
		prcp1 = prcp[stni];
		if (prcp1)
		{
			for (stnj=stni+1 ; stnj<nstns ; stnj++)
			{
				prcp2 = prcp[stnj];
				if (prcp2)
				{
					/* str->regy[n] = sign * (prcp1-prcp2)/(prcp1+prcp2); */
					str->regy[n] = sign * (prcp1-prcp2);
					str->regwt[n] = str->regwtall[n];
				}
				else /* !prcp2 */
				{
					str->regy[n] = str->regwt[n] = 0.0;
				}
				n++;
				sign = -sign;

			} /* end for stnj */

		} /* end if (prcp1) */

		else /* !prcp1 */
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


int prcp_mlr(daymet_prcp_struct* str)
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
	/* return with zeros for regression slopes in the case of a singular array */
	/* $$$NCW removed IS_NAN() for non-AIX compatibility... */
#ifdef _AIX
    if (IS_NAN(det))
#else
    if (isnan(det))
#endif
	{
		coef[0] = coef[1] = coef[2] = coef[3] = 0.0;
		return (!ok);
	}

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
                /* $$$NCW removed IS_NAN() for non-AIX compatibility... */
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
	if (coef[2] < 0.0) coef[2] = 0.0;
	if (coef[2] > 0.002) coef[2] = 0.002;

	/* set some constraints on x,y coefficients */
	if (coef[0] < -0.0001) coef[0] = -0.0001;
	if (coef[0] > 0.0001) coef[0] = 0.0001;
	if (coef[1] < -0.0001) coef[1] = -0.0001;
	if (coef[1] > 0.0001) coef[1] = 0.0001;

	return(!ok);
}


int predict_prcp_noint(daymet_prcp_struct* str)
/* generates predicted prcp, with and without elevation correction
assumes the intercept from regression is 0.0 */
{
	int ok=1;
	int nstns;
	int stni;
	char noz;
	double p1, pc, p, sumwt;
	double* prcp;
	double mx,my,mz,x,y,z;
	double f, f_max;
	double max_stn_day_prcp;
	double max_factor = 2.0;

	noz = str->noz;
	mx = str->coef[0];
	my = str->coef[1];
	mz = str->coef[2];
	x = str->mapx;
	y = str->mapy;
	z = str->elev;
	f_max = str->f_max;
	nstns = str->count;
	prcp = str->listobs+str->dam_off;
	pc = p = sumwt = 0.0;

	max_stn_day_prcp = 0.0;
	for (stni=0 ; stni<nstns ; stni++)
	{
		if (p1 = prcp[stni])
		{
			if (noz)
			{
				p += str->wt[stni] * p1;
			}

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
		if (noz)
		{
			p /= sumwt;
		}
	}
	else
	{
		printf("predict_prcp() sumwt = 0.0 ... exiting\n");
		ok=0;
	}

	/* set some constraints on daily precip */
	if (pc < 0.0) pc = 0.0;
	// if (pc > 20.0) pc = 20.0;

	// constrain to range of obs values, times a fixed factor
    if (pc > max_stn_day_prcp * max_factor) pc = max_stn_day_prcp * max_factor;
	
	str->pp = pc;

	if (noz)
	{
		if (p < 0.0) p = 0.0;
		// if (p > 20.0) p = 20.0;

		// constrain to range of obs values, times a fixed factor
    	if (pc > max_stn_day_prcp * max_factor) pc = max_stn_day_prcp * max_factor;
		str->ppnoz = p;
	}

	return (!ok);
}
