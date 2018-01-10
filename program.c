#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#define N 100
#define TMAX 700
#define dt 0.01
#define PI 3.14159

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
			if(i==j)
				matrix[i][j]=0;		
		}
	}
}

void OscillatorsFreqProfile(double* w,double w1,double w2){

	int i;
	double Aw=1.0*(w2-w1)/N;
	for(i=0;i<N;i++){
		w[i]=w1+Aw*i;	
	}
}

double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}


void OscillatorsGaussianProfile(double* w){

	int i;
	double Aw=1.0*(w1-w0)/N;
	for(i=0;i<N;i++){
		w[i]=randn(w0,N*Aw/10);	
	}
}


void InitialConditions(double *v){

	int i;
	for(i=0;i<N;i++){
		//v[i]=2*PI*(1.0*rand()/RAND_MAX);
		v[i]=(0.1*rand()/RAND_MAX);
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
	if((file = fopen("oscillators_unif20.dat", "w")) == NULL){
		printf("No puedo abrir el archivo.\n");
		exit(1);
	}

	FILE* file2;
	if((file2 = fopen("sinchronization_unif20.dat", "w")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}

	FILE* file3;
	if((file3 = fopen("frecuency_unif20.dat", "w")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}

	srand(time(NULL));
	double* neuron=(double*)malloc(N*sizeof(double));
	double** conn=allocMatrix(N,N);
	double* freq=(double*)malloc(N*sizeof(double));
	double* max=(double*)malloc(N*sizeof(double));
	int* counter=(int*)malloc(N*sizeof(int));
	double* freq_final=(double*)malloc(N*sizeof(double));
	double eps1=20,eps2=0;
	double phase=0;
	double t;
	int i,j;

	for(i=0;i<N;i++){
		max[i]=0;
		counter[i]=0;
	}

	Connectivity_1(conn);
	//ConnectivityAlltoAll(conn);data6
	//PrintMatrix(conn,N,N);
	//OscillatorsGaussianProfile(freq);
	OscillatorsFreqProfile(freq);
	InitialConditions(neuron);

	t=0;	
	while(t<TMAX){
		RK2(neuron,freq,conn,eps1,eps2,phase);
		fprintf(file, "%lf\t",t);
		double x=0,y=0;
		double sinchr=0;
		for(i=0;i<N;i++){
			fprintf(file, "%lf\t",neuron[i]);
			x+=cos(neuron[i]);
			y+=sin(neuron[i]);
			if(neuron[i]>max[i]){
				max[i]=neuron[i];
			}
			if(neuron[i]<max[i]){
				counter[i]+=1;
				max[i]=neuron[i];				
				} 
		}
		fprintf(file, "\n");
		sinchr=1.0*sqrt(x*x+y*y)/N;
		fprintf(file2,"%lf\t%lf\n",t,sinchr);
		t+=dt;
	}

	for(i=0;i<N;i++){
		freq_final[i]=2.0*PI*counter[i]/TMAX;			
		fprintf(file3, "%lf\n",freq_final[i]);
		}


	free(max);
	free(counter);
	free(neuron);
	freeMatrix(conn,N);
	free(freq);
	fclose(file);
	fclose(file2);
	fclose(file3);

	return;
}
