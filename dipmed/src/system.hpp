#ifndef INC_system_hpp
#define INC_system_hpp

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

/******************************************************************************/
inline int allocate_1D_array_of_doubles(double **parray, const int size, 
                                         const char *tag)
{
    *parray = (double *)calloc(size,sizeof(double)); 
    assert(*parray);
    if(*parray==NULL) printf("Can't allocate %s\n",tag), abort();
    else return size*sizeof(double);
}

/******************************************************************************/
inline void allocate_2D(double ***m, int d1, int d2) 
{
     *m = new double* [d1];
     for (int i=0; i<d1; i++)
     {
          (*m)[i] = new double [d2];
          for (int j=0; j<d2; j++) (*m)[i][j] = 0.0;
     }
}

/******************************************************************************/
inline void delete_2D(double **m, int d1) 
{
     for (int i=0; i<d1; i++) delete [] m[i];
     delete [] m;
}

/******************************************************************************/
inline void flush_output_buffers(void)
{
    fflush(stdout);
    fflush(stderr);
}

/******************************************************************************/
inline void format_seconds(char *timestr, const int strsize, const int totsecs)
{
     const int minchar = 35;
     const char *tmax = "999 days, 23:59:59 (86399999 s)";
     if(strsize < minchar) 
     {
          printf("\n");
          printf("format_seconds() error : timestr too small!\n");
          printf("The provided string should be at least %d characters wide.\n",minchar);
          printf("[max time = %s]\n",tmax);
          printf("\n");
          abort();
     }
          
     int secs = 0;
     int mins = 0;
     int hours = 0;
     int days = 0;
     
     const int sec_per_min = 60;
     const int min_per_hour = 60;
     const int hour_per_day = 24;
     
     secs  = totsecs;
     mins  = secs/sec_per_min;
     secs -= mins*sec_per_min;
     hours = mins/min_per_hour;
     mins -= hours*min_per_hour;
     days  = hours/hour_per_day;
     hours-= days*hour_per_day;
     
     sprintf(timestr,"%d days, %.02d:%.02d:%.02d (%d s)\n",days,hours,mins,secs,totsecs);
}

/******************************************************************************/
inline int free_array_of_doubles(double *parray, const int size)
{
    if(parray!=NULL)
    {
        free(parray);
        return size*sizeof(double);
    }
    else return 0;
}

/******************************************************************************/
inline void please_sleep(const int ms)
{
    double goal;
    //1 Second = 1000 milliseconds
    goal = (double)ms + (double)clock()/(CLOCKS_PER_SEC/1000.0);
    //Do nothing until goal is greater than the
    //current elapsed time
    while( goal >= clock()/(CLOCKS_PER_SEC/1000)){};
}

/******************************************************************************/
inline void print_character(const char character, const int n)
// Prints n times the provided character, plus \n.
{
    for(int i = n; i--;) printf("%c",character);
    printf("\n");
}

/******************************************************************************/
inline void print_progress(const char character, const int n, const int max)
{
     const int period = (int)max/100;
     if(n != 0)
     {
          if(n%period==0) printf("%c",character);
          if(n%(50*period)==0) printf(" %d%c\n",n/period,'%');
     }
     if(n == max) printf("\n");
     fflush(stdout);
}

/******************************************************************************/
inline void ShowRunTime(const time_t start, const time_t end)
{
     char timestr[35];
     format_seconds(timestr,sizeof(timestr),(int)difftime(end,start));
     puts("");
     print_character('*',60);
     puts("Run time");
     puts(timestr);
}

#endif // #ifndef INC_system_hpp
