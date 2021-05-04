/* 
predict_daily_srad_vp_func.c
Peter Thornton
2/13/01

functions called from predict_daily_srad_vp.c
*/

#include "metsrc2.h"
#include "cproj.h"
#include "proj.h"

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define R        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
#define EPS      0.62196351      /* (MW/MA) unitless ratio of molec weights */

/* calc_pet() calculates the potential evapotranspiration for aridity 
corrections in calc_vpd(), according to Kimball et al., 1997 */
double calc_pet(double rad, double ta, double pa, double dayl)
{
	/* input parameters and units :
	double rad      (W/m2)  daylight average incident shortwave radiation
	double ta       (deg C) daylight average air temperature
	double pa       (Pa)    air pressure
	double dayl     (s)     daylength 
	*/
	
	double rnet;       /* (W m-2) absorbed shortwave radiation avail. for ET */
	double lhvap;      /* (J kg-1) latent heat of vaporization of water */ 
	double gamma;      /* (Pa K-1) psychrometer parameter */
	double dt = 0.2;   /* offset for saturation vapor pressure calculation */
	double t1, t2;     /* (deg C) air temperatures */
	double pvs1, pvs2; /* (Pa)   saturated vapor pressures */
	double pet;        /* (kg m-2 day-1) potential evapotranspiration */
	double s;          /* (Pa K-1) slope of saturated vapor pressure curve */

	/* calculate absorbed radiation, assuming albedo = 0.2  and ground
	heat flux = 10% of absorbed radiation during daylight */
	rnet = rad * 0.72;
		
    /* calculate latent heat of vaporization as a function of ta */
    lhvap = 2.5023e6 - 2430.54 * ta;
    
    /* calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
    where:
    cp       (J/kg K)   specific heat of air
    epsilon  (unitless) ratio of molecular weights of water and air
    */
    gamma = CP * pa / (lhvap * EPS);
    
    /* estimate the slope of the saturation vapor pressure curve at ta */
    /* temperature offsets for slope estimate */
    t1 = ta+dt;
    t2 = ta-dt;
    
    /* calculate saturation vapor pressures at t1 and t2, using formula from 
	Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
	Meteorol. Mag., 114:49-56.
	*/
    pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1));
    pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2));

    /* calculate slope of pvs vs. T curve near ta */
    s = (pvs1-pvs2) / (t1-t2);
    
    /* calculate PET using Priestly-Taylor approximation, with coefficient
    set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day */
	pet = (1.26 * (s/(s+gamma)) * rnet * dayl)/lhvap;
	
	/* return a value in centimeters/day, because this value is used in a ratio
	to annual total precip, and precip units are centimeters */
	return (pet/10.0);
}    
		
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

/* pulled_boxcar() calculates a moving average of antecedent values in an
array, using either a ramped (w_flag=1) or a flat (w_flag=0) weighting */	
int pulled_boxcar(double *input,double *output,int n,int w,int w_flag)
{
	int ok=1;
    int i,j;
    double *wt;
    double total,sum_wt;

    if (w > n) {
        printf("Boxcar longer than array...\n");
        printf("Resize boxcar and try again\n");
        ok=0;
    }
    
    if (ok && !(wt = (double*) malloc(w * sizeof(double))))
    {
    	printf("Allocation error in boxcar()\n");
    	ok=0;
    }
    
    if (ok)
    {
	    /* when w_flag != 0, use linear ramp to weight tails,
	    otherwise use constant weight */
	    sum_wt = 0.0;
	    if (w_flag)
	    {
	        for (i=0 ; i<w ; i++)
	       	{
	            wt[i] = (double)(i+1);
	            sum_wt += wt[i];
	        }
	    }
	    else
	    {
	        for (i=0 ; i<w ; i++)
	        { 	
	            wt[i] = 1.0;
	            sum_wt += wt[i];
	        }
	    }
	    
	    /* fill the output array, starting with the point where a full
	    boxcar can be calculated */
	    for (i=w-1 ; i<n ; i++)
	    {
	        total = 0.0;
	        for (j=0 ; j<w ; j++)
	        {
	            total += input[i-w+j+1] * wt[j];
	        }
	        output[i] = total/sum_wt;
	    }
	    
	    /* fill the first w elements of the output array with the value from
	    the first full boxcar */
	    for (i=0 ; i<w-1 ; i++)
	    {
	    	output[i] = output[w-1];
	    }
	    
	    free(wt);
	    
	} /* end if ok */
	
    return (!ok);
}
/* end of pulled_boxcar() */  
