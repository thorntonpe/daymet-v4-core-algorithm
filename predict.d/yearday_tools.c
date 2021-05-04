/* yearday_tools.c */

int yearday(int first_yday, int first_year, int dayoff, int* ret_day, int* ret_year);
int yearday(int first_yday, int first_year, int dayoff, int* ret_day, int* ret_year)
{
	/* this function contains an option to always ignore leapyears.
	To do so, the user must set the input first_year == 0, in which
	case the daysinyear variable is always 365.
	For all other values of first_year, the program determines if it is a 
	leap year and assigns the number of days in the year accordingly */
	
	int ok=1;
	int nday;
	int year;
	int daysinyear;
	
	year = first_year;
	nday = first_yday + dayoff;
	if (!year || year%4)
	{
		/* not leap year */
		daysinyear = 365;
	}
	else
	{
		/* leap year */
		daysinyear = 366;
	}

	while (nday >= daysinyear)
	{
		/* offset places day in another year */
		nday -= daysinyear;
		if (year) year ++;

		/* is this a leap year? */
		if (!year || year%4)
		{
			/* not leap year */
			daysinyear = 365;
		}
		else
		{
			/* leap year */
			daysinyear = 366;
		}
	}
	
	*ret_day = nday;
	*ret_year = year;
	
	return (!ok);
}

int mday_to_yday(int year, int month, int day, int* ret_yday);
int mday_to_yday(int year, int month, int day, int* ret_yday)
{
	int ok=1;
	int mondays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	int leap,i;
	int yday = 0;
	
	leap = year%4 ? 0:1;
	
	for (i=0 ; i<month-1 ; i++)
	{
		yday += mondays[i];
		if (leap && (i == 1)) yday--;
	}
	
	yday += day-1;
	*ret_yday = yday;

	return(!ok);
}
