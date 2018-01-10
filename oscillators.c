#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#define N 2
#define TMAX 100
#define dt 0.01
#define PI 3.14159
#define w0 2.5
#define w1 5

double InternalCoupling(double epsilon,double phi1,double phi2){
	double x;
	x=epsilon*sin(phi1-phi2);
	return x;
}

double ExternalCoupling(double epsilon,double phi1,double phi2){
	double x;
	x=epsilon*sin(phi1-phi2);
	return x;
}

double* Coupling(double* v,double* w,double** matrix,double epsilon1,double epsilon2,double phi){

	double* x;
    	x = (double*)malloc(N*sizeof(double));
	int i,j;

	for(i=0;i<N;i++){
		x[i]=w[i];
		for(j=0;j<N;j++){
			x[i]+=matrix[i][j]*InternalCoupling(epsilon1,v[j],v[i]);		
		}
		x[i]+=ExternalCoupling(epsilon2,phi,v[i]);  
	}	
	
	return x;
}

double** allocMatrix(int n,int m){

	double** matrix=(double**)malloc(n * sizeof(double*));
	int i;
	for(i=0;i<m;i++){
	    matrix[i] = (double*)malloc(m * sizeof(double));		
	}
	return matrix;
} 
void freeMatrix(double** matrix,int m){
	int i;
	for(i=0;i<m;i++){
	    free(matrix[i]);		
	}
	free(matrix);
} 

void RK2(double *v,double *w,double **matrix,double epsilon1,double epsilon2,double phi){

    double *vaux;
    double *K1,*K2;
    int i;

    vaux = (double*)malloc(N*sizeof(double));
    K1 = (double*)malloc(N*sizeof(double));
    K2 = (double*)malloc(N*sizeof(double));

    for(i=0;i<N;i++)
       	vaux[i]=v[i];
    K1=Coupling(vaux,w,matrix,epsilon1,epsilon2,phi);
   
    for(i=0;i<N;i++)
       	vaux[i]=v[i]+dt*K1[i];
    K2=Coupling(vaux,w,matrix,epsilon1,epsilon2,phi);
    
    for(i=0;i<N;i++){
       	v[i]=v[i]+0.5*(dt*K1[i]+dt*K2[i]);
	if (v[i]>2*PI){
		v[i]=v[i]-2*PI;
		}
	}

    free(vaux);
    free(K1);
    free(K2);
}

void Connectivity_1(double** matrix){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			matrix[i][j]=0;
			if(j==i-1 || j==i+1){
				matrix[i][j]=1;			
			} 
		}
	}
}

void ConnectivityAlltoAll(double** matrix){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			matrix[i][j]=1;			
		}
	}
}

void OscillatorsFreqProfile(double* w){

	int i;
	double Aw=1.0*(w1-w0)/N;
	for(i=0;i<N;i++){
		w[i]=w0+Aw*i;	
	}
}

void InitialConditions(double *v){

	int i;
	for(i=0;i<N;i++){
		v[i]=2*PI*(1.0*rand()/RAND_MAX);
	}
}

void PrintMatrix(double** matrix,int n,int m){

	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<m;j++)
			printf("%lf\t",matrix[i][j]);
		printf("\n");
	}
}

int main(){

	FILE* file;
	if((file = fopen("two_oscillators.dat", "w")) == NULL){
		printf("No puedo abrir el archivo.\n");
		exit(1);
	}

	srand(time(NULL));
	double* neuron=(double*)malloc(N*sizeof(double));
	double** conn=allocMatrix(N,N);
	double* freq=(double*)malloc(N*sizeof(double));
	double eps1=10,eps2=0;
	double phase=0;
	double t;
	int i;

	Connectivity_1(conn);
	OscillatorsFreqProfile(freq);
	InitialConditions(neuron);

	t=0;	
	while(t<TMAX){
		RK2(neuron,freq,conn,eps1,eps2,phase);
		fprintf(file, "%lf\t",t);
		for(i=0;i<N;i++){
			fprintf(file, "%lf\t",neuron[i]);
		}
		fprintf(file, "\n");
		t+=dt;
	}
	fclose(file);

	return;
}
