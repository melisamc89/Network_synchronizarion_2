#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#define N 100
#define TMAX 1000
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
			x[i]+=matrix[i][j]*InternalCoupling(epsilon1,v[i],v[j]);		
		}
		x[i]+=ExternalCoupling(epsilon2,v[i],phi);  
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
	if((file = fopen("oscillator_10.dat", "r")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}
	FILE* file2;
	if((file2 = fopen("frecuency_10.dat", "w")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}

	double* max=(double*)malloc(N*sizeof(double));
	int* counter=(int*)malloc(N*sizeof(int));
	double* freq_final=(double*)malloc(N*sizeof(double));

	int lines=TMAX/dt;
	double* data=(double*)malloc(N*sizeof(double));
	int i,j;
	double t;
	
	for(i=0;i<N;i++){
		max[i]=0;
		counter[i]=0;
	}


	for(j=0;j<lines;j++){
		fscanf(file, "%lf\t",&t);
		for(i=0;i<N;i++){
			fscanf(file,"%lf\t",&data[i]);
			if(data[i]>max[i]){
				max[i]=data[i];
			}
			if(data[i]<max[i]){
				counter[i]+=1;
				max[i]=data[i];				
				} 
		}
		fscanf(file,"\n");
		
	}

	for(i=0;i<N;i++){
		freq_final[i]=2.0*PI*counter[i]/TMAX;			
		fprintf(file2, "%lf\n",freq_final[i]);
		}


	free(max);
	free(counter);
	fclose(file);
	fclose(file2);

	return;
}
