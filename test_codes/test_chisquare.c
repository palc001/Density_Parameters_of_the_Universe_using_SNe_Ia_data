#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

double mzval(double z);

int main(){
	int n=0, i=0, j=0;
	FILE *fdata;
	char ch;
	fdata=fopen("E:/C Programming/data/data1.txt", "r");
	while(1){
		ch=fgetc(fdata);
		if (ch==EOF) break;
		if (ch=='\n') n++;
	}
	fclose(fdata);
	int nn=n*3;
	double lscanner[nn];
	FILE *fdata1 = fopen("E:/C Programming/data1.txt", "r");
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
		mthe[i][0]=mzval(zval[i]);
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
	printf("%f\n", chisquare);
	return 0;
}

double mzval(double z){
	int i, j, n=1000;
	double s, integral, DL, logDL, mz;
	double cl=299792458, M=-18.46;
	double m=0.198;
	double de=1.0-m;
	double a[n+1];
	double h=z/(double)n;
	a[0]=1;
	for (i=1; i<=n; i++){
		a[i]=i*h;
	}
	double c=0, d=0, z1, z2, f1, f2, f3;
	for (j=1; j<=n/2-1; j++){
		z1=a[2*j];
		f1=1/sqrt((1+z1)*(1+z1)*(1+m*z1)-z1*(2+z1)*de);
		c=c+f1;
	}
	for (j=1; j<=n/2; j++){
		z2=a[2*j-1];
		f2=1/sqrt((1+z2)*(1+z2)*(1+m*z2)-z2*(2+z2)*de);
		d=d+f2;
	}
	f3=1/sqrt((1+z)*(1+z)*(1+m*z)-z*(2+z)*de);
	s=a[0]+2*c+4*d+f3;
	integral=h*s/3;
	DL=cl*(1+z)*integral;
	logDL=log(DL)/log(10.0);
	mz=5*logDL+M;
	return mz;
}
