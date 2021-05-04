/*
interpolate_func.c
Peter Thornton

Revision history:
8/6/96   written to support daymet_interpolate.c
9/6/96   adapted from daymet_interpolate_func.c, now general to metsrc
*/

#include "metsrc.h"
#include <float.h>

/* calculates the squared distance from point (x,y) to all stations */
int calc_sqdist(interpolate_struct* str)
{
	int ok=1;
	int i, nstns;
	double dx,dy,x,y;
	double *stnx, *stny;
	
	nstns = str->nstns;
	x = str->x;
	y = str->y;
	stnx = str->stnx;
	stny = str->stny;
	
	for (i=0 ; i<nstns ; i++)
	{
		dx = x - stnx[i];
		dy = y - stny[i];
		str->sqdist[i] = (dx * dx) + (dy * dy);
	}
	
	return (!ok);
}

/* Find initial search radius in attempt to stabilize calc_search_rad() */
int initial_radius(interpolate_struct* str)
{
    int ans = str->ans;
    int nstns = str->nstns;
    int index;
    double distance;
    double max=DBL_MAX;

    int i,j;

    index = -1;
    distance = DBL_MAX;

    //Find the smallest station
    for(i=0; i<nstns; i++)
    {
        if(str->sqdist[i] < distance)
        {
            //printf("Initial distance: %lf\n", str->sqdist[i]);
            distance = str->sqdist[i];
        }
    }

    printf("ans: %d\n", ans);
    //find the ans next smallest stations
    for(j=0; j<ans-1; j++)
    {
        double current = DBL_MAX;
        for(i=0; i<nstns; i++)
        {
            double s = str->sqdist[i];
            //Find next smallest station
            if( (s < current) && (s > distance) && (s > 0))
            {
                //printf("New distance: %lf\n", s);
                current = s;
                index = i;
            }
        }
        //printf("Current: %lf\n", current);
        distance = current;
    }
    //printf("inn before str->dfsqrad %lf\n", str->dfsqrad);
    //plus 1?
    str->dfsqrad = distance+1.0;
    str->inv_dfsqrad = 1.0 / str->dfsqrad;
    str->dfarea = PI * str->dfsqrad;

    printf("new initial str->dfsqrad %lf\n", str->dfsqrad);
    return (0);
}

/* iterative density algorithm to calculate squared search radius */
int calc_search_rad(interpolate_struct* str)
{
	int ok=1;
	int i,j,nstns;
	double dfsqrad, inv_dfsqrad, dfarea;
	double density;
	double w,ans,trunc,avgwt;
	int stncnt;
	double sqd, *sqdist;
	
	dfsqrad = str->dfsqrad;
        inv_dfsqrad = str->inv_dfsqrad;
	dfarea = str->dfarea;
	trunc = str->dftrunc;
	avgwt = str->dfavgwt;
	nstns = str->nstns;
	ans = str->ans;
	sqdist = str->sqdist;

	/* begin density iteration loop */
	for (i=0 ; (i<DFNIT && ok) ; i++)
	{
		w = 0.0;
		stncnt = 0;
		for (j=0 ; j<nstns ; j++)
		{
			/* test if station is within filter */
			if ((sqd = sqdist[j]) < dfsqrad)
			{
                                //printf("DFGSP: %lf\t sqd: %lf\t inv_dfsqrad: %lf\t trunc: %lf\t avgwt: %lf\n",DFGSP, sqd, inv_dfsqrad, trunc, avgwt);
				stncnt ++;
                                //printf("w: %lf\t expression: %lf\n", w, (exp(-DFGSP * sqd * inv_dfsqrad) - trunc) / avgwt);
				w += 1.0;
                                //w += (exp(-DFGSP * sqd * inv_dfsqrad) - trunc) / avgwt;
			}
		}
	        //printf("Stncnt: %d dfsqrad: %f density: %f\n", stncnt, dfsqrad, density);

		/* test for empty density filter */
		if (!stncnt)
		{
			printf("CALC_SEARCH_RAD: no stations found. Resetting radius size. \n");
			printf("inside zero station: i = %d\t dfsqrad = %lf\n",i,dfsqrad);

			dfsqrad *= 1.1;
			inv_dfsqrad = 1.0 / dfsqrad;
			dfarea = PI * dfsqrad;
			i--;
			/* test code */
			//printf("inside zero station: i = %d\t dfsqrad = %lf\n",i,dfsqrad);
			if(isinf(dfsqrad))
                	{
                        	printf("dfsqrad became infinite\n");
                        	exit(1);
                	}
		}
		else
		{
			density = w / dfarea;
                        //printf("density: %lf\t w: %lf\t dfarea: %lf\ti: %d\n", density, w, dfarea, i);

			/* calculate new density parameters for next iteration */
			if (i < DFNIT-1)
			{
				dfsqrad = DFNSM * ans / (density * PI);
				inv_dfsqrad = 1.0 / dfsqrad;
				dfarea = PI * dfsqrad;
                                //printf("dfsqrad: %lf\t DFNSM: %lf\t ans: %lf\t density: %lf\n", dfsqrad, DFNSM, ans, density);
                                
			}
		}
	    //printf("Stncnt: %d dfsqrad: %f density: %f\n", stncnt, dfsqrad, density);
	} /* end of density iteration loop */
	
	/* calculate squared radius for search filter */
	str->sqrad = ans / (density * PI);
	
	/* guarantee that there are at least 3 stations per gridcell */ 
	do 
	{
		stncnt = 0;
		for (j=0 ; j<nstns ; j++)
		{
			/* test if station is within filter */
			if (sqdist[j] < str->sqrad)
			{
				stncnt ++;
			}
		}
		if (stncnt < 3)
		{
			printf("stncnt = %d\tsqrad = %lf\n",stncnt,str->sqrad);
			str->sqrad *= 1.01;
		}
		if(isinf(str->sqrad))
		{
			printf("sqrad became infinite\n");
			exit(1);
		}
	} while (stncnt < 3);
	
	
	return (!ok);
}
	
/* generates weighted station list given interpolation parameters */
int weight_list(interpolate_struct* str)
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
		if ((sqd = sqdist[i]) < sqrad)
		{
			str->id[count] = i;
			str->wt[count] = exp(-gsp * sqd * inv_sqrad) - trunc;
			wtsum += str->wt[count];
			count++;
		}
	}
        printf("Final count: %d\t sqrad %lf\n", count, sqrad);
	/* normalize weights */
	if (wtsum)
	{
		for (i=0 ; i<count ; i++)
		{
			str->wt[i] /= wtsum;
		}
		str->count = count;
	}
	else
	{
		printf("no stations in list\n");
		ok = 0;
	}
	
	return (!ok);
}
