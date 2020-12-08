#include <stdio.h>
#include <math.h>
#include <conio.h>

int main(){
	int i, j, n;
	double b=2*M_PI, s, integral;
	double error;
	printf("Estimation of Definite Integral of 'sin(x/50) from 0 to 2pi' using Simpson's One-third Rule\n\n");
	printf("Enter the number of subintervals you want. (More sub-intervals = More accuracy but slower processing speed)\n");
	printf("n = ");
	scanf("%d", &n);
	printf("\nThe number of sub-intervals you chose = ");
	printf("%d", n);
	
	if (n % 2 == 0){
		double a[n+1];
		double h=b/n;
		printf("\nValue of pi that I use is %.11f\n", M_PI);
		a[0]=0;
		double c, d, f1, f2;
		c=0;
		d=0;
		for (i=1; i<=n; i++){
			a[i]=i*h;
		}
		for (j=1; j<=n/2-1; j++){
			f1=sin(a[2*j]/50);
			c=c+f1;
		}
		for (j=1; j<=n/2; j++){
			f2=sin(a[2*j-1]/50);
			d=d+f2;
		}
		
		s=a[0]+2*c+4*d+sin(b/50);
		integral=h*s/3;
		error=b*b*b*b*b*sin(b/50)/(1.8*6.25*n*n*n*n*100000000);
		
		printf("\n\nThe value of the integral is = %.15f +/- %.30f\n", integral, error);	
	}
    
	else{
		printf("\nn has to be even for this to work. Please restart program and input an even n.");
	}
	
	return 0;
}
