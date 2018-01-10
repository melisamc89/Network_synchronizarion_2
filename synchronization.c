#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>

#define N 100
#define TMAX 100
#define dt 0.01
#define PI 3.14159

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

void PrintMatrix(double** matrix,int n,int m){

	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<m;j++)
			printf("%lf\t",matrix[i][j]);
		printf("\n");
	}
}

double* OscillatorsFrequency(double** PhaseEvol,int n,int rows){

	double* max=(double*)malloc(n*sizeof(double));
	int* counter=(int*)malloc(n*sizeof(int));
	double* freq=(double*)malloc(n*sizeof(double));
	int i,j;

	for(i=0;i<N;i++){
		max[i]=0;
		counter[i]=0;
	}

	for(j=0;j<rows;j++){
		for(i=0;i<n;i++){
			if(PhaseEvol[j][i]>max[i]){
				max[i]=PhaseEvol[j][i];
				}
			if(PhaseEvol[j][i]<max[i]){
				counter[i]+=1;
				max[i]=PhaseEvol[j][i];				
				} 
			}			
		}

	for(i=0;i<N;i++){
		freq[i]=1.0*counter[i]/TMAX;
		}

	free(max);
	free(counter);
	return freq;
	}

int main(){

	FILE* file;
	if((file = fopen("oscillator_10.dat", "r")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}
	FILE* file2;
	if((file2 = fopen("sinchronization_10.dat", "w")) == NULL){
		printf("No puedo abrir el archivo de lectura.\n");
		exit(1);
	}

	int lines=TMAX/dt;
	double* data=(double*)malloc(N*sizeof(double));
	printf("%i\t%i\n",lines,N);
	double* freq;
	int i,j;
	double t;


	for(j=0;j<lines;j++){
		fscanf(file, "%lf\t",&t);
		double x=0,y=0;
		double sinchr=0;
		for(i=0;i<N;i++){
			fscanf(file,"%lf\t",&data[i]);
			x+=cos(data[i]);
			y+=sin(data[i]);
		}
		fscanf(file,"\n");
		sinchr=sqrt(x*x+y*y)/N;
		fprintf(file2,"%lf\t%lf\n",t,sinchr);
		
	}
	fclose(file);
	fclose(file2);
	return;
}
