#include <stdio.h>
#include <math.h>
#include <conio.h>

int main(){
	int i, j, n;
	double s, integral, error;
	float tz, tm, tde;
	printf("Estimation of Integral using Simpson's One-third Rule\n\nEnter Values of OmegaM and OmegaDE respectively\n");
	scanf("%f %f", &tm, &tde);
	printf("\nEnter the value of z for integration\nz = ");
	scanf("%f", &tz);
	n=1000;
	double z=(double)tz;
	double m=(double)tm;
	double de=(double)tde;
	double a[n+1], h=z/n;
	a[0]=1;
	double c=0, d=0, z1, z2, f1, f2, f3;
	for (i=1; i<=n; i++){
		a[i]=i*h;
	}
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
	printf("\n\nThe value of the integral is = %.15f +/- %.30f\n", integral, error);
	return 0;
}
