#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <stdlib.h>

int main(){
	int i, j, n=1000;
	double s, integral, DL, logDL, mz;
	float tm, tde;
	printf("Generating plot\n\nEnter Values of OmegaM and OmegaDE respectively.\n");
	scanf("%f %f", &tm, &tde);
	double z=0, dz=0.01, cl=299792458, M=-18.46;
	double m=(double)tm;
	double de=(double)tde;
	printf("%f %f\n\n", m, de);
	FILE *f = fopen("test_luminosity_distance_error_generated_data.txt", "w");
	if (f == NULL){
       perror("test_luminosity_distance_error_generated_data.txt");
       exit(1);
	}
	while (z<=1){
		double a[n+1], h=z/n;
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
		
		double y=0, ez1=0, ez, error1, error2, error3, error;
		while (y<=z){
			ez=(3*(35*(y+1)*(y+1)*(y+1)*(y+1)*(-2*de+3*m*y+m+2)*(-2*de+3*m*y+m+2)*(-2*de+3*m*y+m+2)*(-2*de+3*m*y+m+2)-60*(y+1)*(y+1)*(-2*de+m*(6*y+4)+2)*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*(-2*de+3*m*y+m+2)*(-2*de+3*m*y+m+2)+96*m*(y+1)*(de*y*(y+2)-(y+1)*(y+1)*(m*y+1))*(de*y*(y+2)-(y+1)*(y+1)*(m*y+1))*(-2*de+3*m*y+m+2)+12*(-2*de+m*(6*y+4)+2)*(-2*de+m*(6*y+4)+2)*(de*y*(y+2)-(y+1)*(y+1)*(m*y+1)*(m*y+1))*(de*y*(y+2)-(y+1)*(y+1)*(m*y+1)*(m*y+1))))/(16*sqrt(((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))*((y+1)*(y+1)*(m*y+1)-de*y*(y+2))));
			if (ez1<ez){
				ez1=ez;
			}
			y=y+0.0001;
		}
		error1=abs(h*h*h*h*z*ez1/180);
		error2=DL*sqrt((0.01/z)*(0.01/z)+(error1/integral)*(error1/integral));
		error3=5*0.434*(error2/DL);
		error=mz*sqrt((0.03/3.17)*(0.03/3.17)+(error3/(5*logDL))*(error3/(5*logDL)));
		
		fprintf(f, "%f   %f   %le\n", z, mz, error1);
		z=z+dz;
	}
	fclose(f);
	return 0;
}
