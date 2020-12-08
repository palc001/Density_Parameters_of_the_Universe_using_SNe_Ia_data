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

#define N 100000
#define NN 30000

float ran1(long *idum);
float gasdev(long *idum);
double q(double x, double y);
int metropolis(int numIter, double *pointsx, double *pointsy, double *pointsq);
double frequencychecker(int numiter, double inc, double *ptsfx, double *ptsfy, double *pfx, double *pfy, double *pfz);
double chisqcompute(double mex, double dey);
double mzval(double z, double me, double de);

double ptsx[N] = { 0 };
double ptsy[N] = { 0 };
double ptsq[N] = { 0 };
double px[NN] = { 0 };
double py[NN] = { 0 };
double pz[NN] = { 0 };

int main(){
	int j1=0;
	int j2=metropolis(N, ptsx, ptsy, ptsq);
    FILE *fp1, *fp2, *fp3;
    fp2=fopen("test_MCMC_MH/mcmcomegade.txt", "w");
    fp1=fopen("test_MCMC_MH/mcmcomegam.txt", "w");
    fp3=fopen("test_MCMC_MH/mcmc.txt", "w");
    for (j1=0; j1<j2; j1++){
        fprintf(fp3, "%f	%f	%.22f\n", ptsx[j1], ptsy[j1], ptsq[j1]);
        fprintf(fp1, "%f\n", ptsx[j1]);
        fprintf(fp2, "%f\n", ptsy[j1]);
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    double inc1=0.005;
    int maxiy;
    maxiy=frequencychecker(j2, inc1, ptsx, ptsy, px, py, pz);
    printf("maxiy = %d",  maxiy);
    printf("\n\nDONE!");
    return EXIT_SUCCESS;
}

float ran1(long *idum){
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

float gasdev(long *idum){
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	if (*idum < 0) iset=0;
	if (iset == 0){
		do{
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		}while (rsq >= 1.0 || rsq == 0.0);
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

double chisqcompute(double mex, double dey){
	int n=0, i=0, j=0;
	double me1=mex, de1=dey;
	FILE *fdata;
	char ch;
	fdata=fopen("data/data1.txt", "r");
	while(1){
		ch=fgetc(fdata);
		if (ch==EOF) break;
		if (ch=='\n') n++;
	}
	fclose(fdata);
	int nn=n*3;
	double lscanner[nn];
	FILE *fdata1 = fopen("data/data1.txt", "r");
	for (i=0; i<nn; i++){
		fscanf(fdata1, "%lf", &lscanner[i]);
	}
	fclose(fdata1);
	double sigma[n][n];
	i=2;
	while(i<nn){
		sigma[j][j]=lscanner[i];
		i=i+3;
		j=j+1;
	}
	i=j=0;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(i!=j){
				sigma[i][j]=0;
			}
		}
	}
	double sigmainv[n][n];
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(i!=j){
				sigmainv[i][j]=0;
			}
			else{
				sigmainv[i][j]=(double)1/(sigma[i][j]);
			}
		}
	}
	double mobs[n][1];
	double mthe[n][1];
	double M[n][1];
	i=1;
	j=0;
	while(i<nn){
		mobs[j][0]=lscanner[i];
		i=i+3;
		j=j+1;
	}
	double zval[n];
	i=0;
	j=0;
	while(i<nn){
		zval[j]=lscanner[i];
		i=i+3;
		j=j+1;
	}
	i=j=0;
	for(i=0; i<n; i++){
		mthe[i][0]=mzval(zval[i], me1, de1);
	}
	i=j=0;
	for(i=0; i<n; i++){
		M[i][0]=mobs[i][0]-mthe[i][0];
	}
	i=0;
	double MT[1][n];
	for(i=0; i<n; i++){
		MT[0][i]=M[i][0];
	}
	double matA[1][n];
	for(i=0; i<n; i++){
		matA[0][i]=MT[0][i]*sigmainv[i][i];
	}
	i=0;
	double product, chisquare=0;
	for(i=0; i<n; i++){
		product=matA[0][i]*M[i][0];
		chisquare=chisquare+product;
	}
	return chisquare;
}

double mzval(double z, double me, double de){
	int i, j, n=1000;
	double s, integral, DL, logDL, mz;
	double cl=299792458, MM=-18.46;
	double a[n+1];
	double h=z/(double)n;
	a[0]=1;
	for (i=1; i<=n; i++){
		a[i]=i*h;
	}
	double c=0, d=0, z1, z2, f1, f2, f3;
	for (j=1; j<=n/2-1; j++){
		z1=a[2*j];
		f1=1/sqrt((1+z1)*(1+z1)*(1+me*z1)-z1*(2+z1)*de);
		c=c+f1;
	}
	for (j=1; j<=n/2; j++){
		z2=a[2*j-1];
		f2=1/sqrt((1+z2)*(1+z2)*(1+me*z2)-z2*(2+z2)*de);
		d=d+f2;
	}
	f3=1/sqrt((1+z)*(1+z)*(1+me*z)-z*(2+z)*de);
	s=a[0]+2*c+4*d+f3;
	integral=h*s/3;
	DL=cl*(1+z)*integral;
	logDL=log(DL)/log(10.0);
	mz=5*logDL+MM;
	return mz;
}

double q(double mx, double dy){
	double chisqval=chisqcompute(mx, dy);
	double qval=exp(-chisqval);
	return qval;
}

int metropolis(int numIter, double *pointsx, double *pointsy, double *pointsq){
	double xg = 1.0;
	double yg = 1.0-xg;
	double s=0.2;
	double fxyg = q(xg, yg);
	double xdash, ydash, fxygn, alpha, accratio;
	float u;
	long seed1=-10, seed2=-20, i=0, k=0;
    while(i < numIter){
    	xdash = xg + s*gasdev(&seed2);
    	ydash = 1.0-xdash;
        fxygn = q(xdash, ydash);
        alpha = fxygn/fxyg;
        if(alpha>=1){
        	xg=xdash;
        	yg=ydash;
        	fxyg=fxygn;
        	pointsx[k] = xg;
        	pointsy[k] = yg;
        	pointsq[k] = fxyg;
        	k=k+1;
		}
		else{
			u = ran1(&seed1);
        	if (u<alpha){
        		xg=xdash;
        		yg=ydash;
        		fxyg=fxygn;
        		pointsx[k] = xg;
        		pointsy[k] = yg;
        		pointsq[k] = fxyg;
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
	double xvalue=-0.1;
	double yvalue=0.4;
	int ka=0, ia=0;
	int iy=0;
	fp4=fopen("test_MCMC_MH/mcmcfrequency.txt", "w");
	while (xvalue<=0.6){
		while(yvalue<=1.1){
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
		yvalue=0.4;
		xvalue=xvalue+inc;
	}
	fclose(fp4);
	printf("Frequencies for 2D data created!\n\niy = %d\n\n", iy);
	return iy;
}
