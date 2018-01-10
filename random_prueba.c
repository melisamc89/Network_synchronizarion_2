#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <time.h>


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

int main(){

	FILE *file1;
	if((file1 = fopen("random.dat", "w")) == NULL){
		printf("No puedo abrir el archivo.\n");
		exit(1);
	}

	srand(time(NULL));
	int i;
	double num;

	for(i=0;i<10000;i++){
		num=randn(5,2);
		fprintf(file1,"%lf\n",num);
	}

	return 0;
}
