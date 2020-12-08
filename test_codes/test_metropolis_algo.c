#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define N 1000000
#define NN 20000

double ran1(long *idum);
double gasdev(long *idum);
double q(double x, double y);
int metropolis(int numIter, double *pointsx, double *pointsy, double *pointsq);
double frequencychecker(int numiter, double inc, double *ptsfx, double *ptsfy, double *pfx, double *pfy, double *pfz);

double ptsx[N] = { 0 };
double ptsy[N] = { 0 };
double ptsq[N] = { 0 };
double px[NN] = { 0 };
double py[NN] = { 0 };
double pz[NN] = { 0 };

int main(){
	int j1=0;
    int j2=metropolis(N, ptsx, ptsy, ptsq);
    FILE *fp1;
    fp1=fopen("test_metropolis_algo_metropolis2d.txt", "w");
    for (j1=0; j1<j2; j1++){
        fprintf(fp1, "%f	%f\n", ptsx[j1], ptsy[j1]);
    }
    fclose(fp1);
    double inc1=0.1;
    int maxiy;
    maxiy=frequencychecker(j2, inc1, ptsx, ptsy, px, py, pz);
    printf("maxiy = %d",  maxiy);
    printf("\n\nDONE!");
    return EXIT_SUCCESS;
}

double ran1(long *idum){
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0 || !iy){
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--){
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

double gasdev(long *idum){
	double ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	if (*idum < 0) iset=0;
	if (iset == 0){
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else{
		iset=0;
		return gset;
	}
}

double q(double x, double y){
	return (double) ((exp(-(pow(x-1,2.0))/(2.0)))*(exp(-(pow(y-2,2.0))/(4.0))));
}

int metropolis(int numIter, double *pointsx, double *pointsy, double *pointsq){
	double xg = 0.0;
	double yg = 0.0;
	double s = 0.2;
	double fxyg = q(xg, yg);
	double xdash, ydash, fxygn, alpha, u, accratio;
	long seed1=-10, seed2=-20, i=0, k=0;
    while(i < numIter){
    	xdash = xg + s*gasdev(&seed2);
    	ydash = yg + s*gasdev(&seed2);
        fxygn = q(xdash, ydash);
        alpha = fxygn/fxyg;
        if(alpha>=1){
        	xg=xdash;
        	yg=ydash;
        	fxyg = fxygn;
        	pointsx[k] = xg;
        	pointsy[k] = yg;
        	pointsq[k] = fxygn;
        	k=k+1;
		}
		else{
			u = ran1(&seed1);
        	if (u<alpha){
        		xg=xdash;
        		yg=ydash;
        		fxyg = fxygn;
        		pointsx[k] = xg;
        		pointsy[k] = yg;
        		pointsq[k] = fxygn;
        		k=k+1;
			}
		}
        i++;
    }
    accratio=(double)k/N;
    printf("Number of accepted values = %d\nAcceptance Ratio = %f\n\n", k, accratio);
    return k;
}

double frequencychecker(int numiter, double inc, double *ptsfx, double *ptsfy, double *pfx, double *pfy, double *pfz){
	FILE *fp4;
	double xvalue=-4;
	double yvalue=-6;
	int ka=0, ia=0;
	int iy=0;
	fp4=fopen("test_metropolis_algo_metropolis2dfrequency.txt", "w");
	while (xvalue<=7.0){
		while(yvalue<=9.0){
			for(ia=0; ia<numiter; ia++){
				if(xvalue<=ptsfx[ia] && ptsfx[ia]<(xvalue+inc) && yvalue<=ptsfy[ia] && ptsfy[ia]<(yvalue+inc)){
					ka=ka+1;
				}
			}
			pfx[iy]=xvalue+(inc/2);
			pfy[iy]=yvalue+(inc/2);
			pfz[iy]=ka;
			fprintf(fp4, "%f	%f	%f\n", pfx[iy], pfy[iy], pfz[iy]);
			yvalue=yvalue+inc;
			iy=iy+1;
			ka=0;
		}
		fprintf(fp4, "\n");
		yvalue=-4;
		xvalue=xvalue+inc;
	}
	fclose(fp4);
	printf("Frequencies for 2D data created!\n\niy = %d\n\n", iy);
	return iy;
}
