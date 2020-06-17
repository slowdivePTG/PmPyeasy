#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "errmess.h"

#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
        int     i, imax=0, j, k;
        double  big, dum, sum, temp2, *vv;

  if (!(vv=(double *)malloc((n+1)*sizeof(double)))) errmess("malloc(vv)");

  *d=1.0;
  for (i=1; i<=n; i++)
  {
    big=0.0;
    for (j=1; j<=n; j++)
      if ((temp2=fabs(a[i][j])) > big) big=temp2;

    if (big == 0.0)
    {
      printf("ERROR! Singular matrix in routine LUDCMP (aga)\n");
      exit(-2);
    }
    vv[i]=1.0/big;
  }

  for (j=1; j<=n; j++)
  {
    for (i=1; i<j; i++)
    {
      sum=a[i][j];
      for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];

      a[i][j]=sum;
    }

    big=0.0;
    for (i=j; i<=n; i++)
    {
      sum=a[i][j];
      for (k=1; k<j; k++) sum -= a[i][k]*a[k][j];

      a[i][j]=sum;

      if ((dum=vv[i]*fabs(sum)) >= big)
      {
        big=dum;
        imax=i;
      }
    }

    if (j != imax)
    {
      for (k=1; k<=n; k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }

    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n)
    {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }

  free(vv);

  return;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 1!57. */
