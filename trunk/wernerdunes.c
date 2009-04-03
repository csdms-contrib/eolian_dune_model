#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

int **height,*iup,*jup,*idown,*jdown,lattice_size_x,lattice_size_y;
float **mask,thresh;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
       m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
       m[i] -= ncl;
       }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
     m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
     m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


void setupgridneighborsperiodic()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=lattice_size_x;
     iup[lattice_size_x]=1;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=lattice_size_y;
     jup[lattice_size_y]=1;
}

void avalanchedown(i,j)
int i,j;
{    float hshadow;
     int i2,violate,min,mini,minj;

     violate=0;min=height[i][j];mini=i;minj=j;
     if (height[i][j]-height[iup[i]][j]>thresh)
      {violate=1;
       if (height[iup[i]][j]<min)
        {min=height[iup[i]][j];mini=iup[i];minj=j;}}
     if (height[i][j]-height[idown[i]][j]>thresh)
      {violate=1;
       if (height[idown[i]][j]<min)
        {min=height[idown[i]][j];mini=idown[i];minj=j;}}
     if (height[i][j]-height[i][jup[j]]>thresh)
      {violate=1;
       if (height[i][jup[j]]<min)
        {min=height[i][jup[j]];mini=i;minj=jup[j];}}
     if (height[i][j]-height[i][jdown[j]]>thresh)
      {violate=1;
       if (height[i][jdown[j]]<min)
        {min=height[i][jdown[j]];mini=i;minj=jdown[j];}}
     if (violate==1)
      {height[i][j]--;
       height[mini][minj]++;
       avalanchedown(mini,minj);}
     else
      {hshadow=height[i][j]-0.5;
       i2=iup[i];
       while (height[i2][j]<hshadow)
        {if (mask[i2][j]<hshadow) mask[i2][j]=hshadow;
         i2=iup[i2];hshadow-=0.5;}}
}

void avalancheup(i,j)
int i,j;
{      float hshadow;
       int i2,violate,min,mini,minj;

       violate=0;min=height[i][j];mini=i;minj=j;
       if (height[iup[i]][j]-height[i][j]>thresh)
        {violate=1;
         if (height[iup[i]][j]>min)
          {min=height[iup[i]][j];mini=iup[i];minj=j;}}
       if (height[idown[i]][j]-height[i][j]>thresh)
        {violate=1;
         if (height[idown[i]][j]>min)
          {min=height[idown[i]][j];mini=idown[i];minj=j;}}
       if (height[i][jup[j]]-height[i][j]>thresh)
        {violate=1;
         if (height[i][jup[j]]>min)
          {min=height[i][jup[j]];mini=i;minj=jup[j];}}
       if (height[i][jdown[j]]-height[i][j]>thresh)
        {violate=1;
         if (height[i][jdown[j]]>min)
          {min=height[i][jdown[j]];mini=i;minj=jdown[j];}}
       if (violate==1)
        {height[i][j]++;
         height[mini][minj]--;
         avalancheup(mini,minj);}
       else
        {hshadow=height[i][j]-0.5;
         i2=iup[i];
         while ((hshadow<mask[i2][j])&&(hshadow>0))
          {if (height[i2][j]>=hshadow) mask[i2][j]=0; else mask[i2][j]=hshadow;
           i2=iup[i2];hshadow-=0.5;}}
}

main()
{    FILE *fp1,*fp2;
     int idum,t,l,ijump,jjump,i,j,duration;
     float psand,pbed,p;

     fp1=fopen("wernermodeltopo","w");
     fp2=fopen("wernermodelshadowmask","w");
     idum=-56;
     psand=0.6;
     pbed=0.4;
     thresh=2;
     l=5;
     lattice_size_x=300;
     lattice_size_y=300;
     duration=100;
     setupgridneighborsperiodic();
     height=imatrix(1,lattice_size_x,1,lattice_size_y);
     mask=matrix(1,lattice_size_x,1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      for (j=1;j<=lattice_size_y;j++)
       {height[i][j]=3;
        mask[i][j]=0;}
     for (t=1;t<=duration*lattice_size_x*lattice_size_y;t++)
      {if (t%(lattice_size_x*lattice_size_y)==0) printf("%d\n",t);
       ijump=(int)(ran3(&idum)*lattice_size_x)+1;
       while (ijump>lattice_size_x) ijump=(int)(ran3(&idum)*lattice_size_x)+1;
       jjump=(int)(ran3(&idum)*lattice_size_y)+1;
       while (jjump>lattice_size_y) jjump=(int)(ran3(&idum)*lattice_size_y)+1;
       while ((height[ijump][jjump]==0)||(mask[ijump][jjump]>0.1))
        {ijump=(int)(ran3(&idum)*lattice_size_x)+1;
         while (ijump>lattice_size_x) ijump=(int)(ran3(&idum)*lattice_size_x)+1;
         jjump=(int)(ran3(&idum)*lattice_size_y)+1;
         while (jjump>lattice_size_y)
          jjump=(int)(ran3(&idum)*lattice_size_y)+1;}
       height[ijump][jjump]--;
       avalancheup(ijump,jjump);
       ijump=(ijump+l)%lattice_size_x+1;
       if (mask[ijump][jjump]>0.1) p=1;
         else if (height[ijump][jjump]>0) p=psand; else p=pbed;
       while (ran3(&idum)>p)
        {ijump=(ijump+l)%lattice_size_x+1;
         if (height[ijump][jjump]>0) p=psand; else p=pbed;}
       height[ijump][jjump]++;
       avalanchedown(ijump,jjump);}
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fprintf(fp1,"%d\n",height[i][j]);
        fprintf(fp2,"%f\n",mask[i][j]);}
     fclose(fp1);
     fclose(fp2);
}
