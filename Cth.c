//#include "Cth.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

double g(double x, int k, double m){
  return cos(k*x)*sqrt(m*m+2*(1-cos(x)))/(4*M_PI);
}

double f(double x, int k, double m){
  return cos(k*x)/(4*M_PI*sqrt(m*m+2*(1-cos(x))));
}

double integrar_f(double *x, int len_x, int k, double m){
  double I = 0;
  double dx = x[1]-x[0];
  for(int i = 0; i < len_x-1; i++){
    I += (dx/6)*(f(x[i], k, m) + 4*f(0.5*(x[i]+x[i+1]), k, m) + f(x[i+1], k, m));
  }
  return I;
}

double integrar_g(double *x, int len_x, int k, double m){
  double I = 0;
  double dx = x[1]-x[0];
  for(int i = 0; i < len_x-1; i++){
    I += (dx/6)*(g(x[i], k, m) + 4*g(0.5*(x[i]+x[i+1]), k, m) + g(x[i+1], k, m));
  }
  return I;
}

int main(int argc, char* argv[]){
  double m = 0.0001;
  int R = 300;

  int t0 = time(NULL);

  if (argc > 1){
    sscanf(argv[1], "%lf", &m);
    if (argc > 2) sscanf(argv[2], "%d", &R);
  }

  double *res1 = (double *) malloc(R*sizeof(double));
  double *res2 = (double *) malloc(R*sizeof(double));

  int N_inter = 10000;
  double *x = (double *) malloc(N_inter*sizeof(double));
  x[0] = -M_PI;
  for(int i = 1; i < N_inter; i++){
    x[i] = x[i-1] + 2*M_PI/N_inter;
  }
  for(int k = 0; k < R; k++){
    res1[k] = integrar_f(x, N_inter, k, m);
    res2[k] = integrar_g(x, N_inter, k, m);

  }

  free(x);
  double *X = (double *) malloc(R*R*sizeof(double));
  double *P = (double *) malloc(R*R*sizeof(double));

  for(int i = 0; i < R; i++){
    for(int j = i; j < R; j++){
      X[i*R+j] = res1[j-i];
      X[j*R+i] = X[i*R+j];
      P[i*R+j] = res2[j-i];
      P[j*R+i] = P[i*R+j];
    }
  }

  free(res1);
  free(res2);

  char filename[255];
  sprintf(filename, "X_m=%1.4f.txt", m);
  FILE* fp = fopen(filename, "w");
  for(int i = 0; i < R; i++){
    for(int j = 0; j < R; j++){
      fprintf(fp, "%f ", X[i*R+j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  sprintf(filename, "P_m=%1.4f.txt", m);
  fp = fopen(filename, "w");
  for(int i = 0; i < R; i++){
    for(int j = 0; j < R; j++){
      fprintf(fp, "%f ", P[i*R+j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  free(X);
  free(P);
  printf("Tardo %d segundos\n", (int) time(NULL)-t0);
  return 0;
}
