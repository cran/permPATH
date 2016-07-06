// Copyright (C) 2016 Ivo D. Shterev

#include <R.h>

#include "ranker.h"

using namespace std;

extern "C"{
void YpermI(int *Y, const int N)
{
  for (int n = 1; n < N; n++){
    int idx = n*unif_rand();   // uniform over [0,n]
    int val = Y[n];
    Y[n] = Y[idx];
    Y[idx] = val;
  }
}

void YpermD(double *Y, const int N)
{
  for (int n = 1; n < N; n++){
    int idx = n*unif_rand();   // uniform over [0,n]
    double val = Y[n];
    Y[n] = Y[idx];
    Y[idx] = val;
  }
}

void Index0(const int *Y, int *N0, const int N)
{
  int count = 0;
  for (int n = 0; n < N; n++){
    if (Y[n] == 0){
      N0[count] = n;
      count++;
    }
  }
}

void SumSumSq(const double *X, double *Sum, double *SumSq, const int N, const int K){
  for (int k = 0; k < K; k++){
    // compute start of X segment
    int iter = k*N;
    for (int n = 0; n < N; n++){
      Sum[k] += X[iter+n];
      SumSq[k] += X[iter+n] * X[iter+n];
    }
  }
}

void Tstat(const double *X, double *T, const double *Sum, const double *SumSq, const int *N0, const int N, const int K, const int n0, const int n1)
{
  double s0;
  double s1;

  double x0;
  double x1;

  float xj;
  
  for(int k = 0; k < K; k++){
    // compute start of X segment
    int iter = k*N;

    s0 = 0.0;
    x0 = 0.0;
    for (int n = 0; n < n0; n++){
      xj = X[iter + N0[n]];
      x0 += xj;
      s0 += xj * xj;
    }
    x0 /= n0;
    x1  = (Sum[k] - n0*x0) / n1;
     
    s1  = (SumSq[k]-s0-n1*x1*x1) / (n1-1);
    s0 = (s0-n0*x0*x0) / (n0-1);

    T[k] = (x0-x1) / sqrt(s0/n0+s1/n1);
  }
}

void PearsonStat(const double *Y, const double *X, double *T, const double *meanX, const double *sumquadX, const int N, const int K, const double meanY, const double sumquadY)
{
  double stat;
  for (int k = 0; k < K; k++){
    // compute start of X segment
    int iter = k*N;

    stat = 0.0;

    // compute test statistic
    for (int n = 0; n < N; n++)   
      stat += X[iter+n] * Y[n];

    stat = (stat-N*meanY*meanX[k])/sqrt(sumquadX[k]*sumquadY);
    T[k] = stat*sqrt(N-2.0)/sqrt(1.0-stat*stat);
  }
}

void SpearmanStat(const double *Y, const double *X, double *T, const int N, const int K)
{
  double stat;
  for (int k = 0; k < K; k++){
    // compute start of X segment
    int iter = k*N;

    stat = 0.0;

    // compute test statistic
    for (int n = 0; n < N; n++)   
      stat += pow(X[iter+n]-Y[n], 2.0);

    T[k] = 1.0 - 6.0*stat/(N*N*N-N);
  }
}

void SpearmanStatNew(const double *Y, const double *X, double *T, const int N, const int K)
{
  for (int k = 0; k < K; k++){
    // compute start of X segment
    int iter = k*N;

    double meanx = 0.0;
    double meany = 0.0;
    for (int n = 0; n < N; n++){
      meanx += X[iter+n];
      meany += Y[n];
    }
    meanx = meanx / N;
    meany = meany / N;

    // compute test statistic
    double nom = 0.0;
    double xterm = 0.0;
    double yterm = 0.0;
    for (int n = 0; n < N; n++){   
      nom += (X[iter+n]-meanx)*(Y[n]-meany);
      xterm += (X[iter+n]-meanx)*(X[iter+n]-meanx);
      yterm += (Y[n]-meany)*(Y[n]-meany);
    }

    T[k] = nom/sqrt(xterm*yterm);
  }
}

void JT_test_norm(const double *X, const int *Y, const int C, const int K, const int N, double *T, double *mean, double *var)
{
  int *count;
  count = (int*)malloc(C*sizeof(int));

  int *ind;
  ind = (int*)malloc(N*sizeof(int));
  
  int cit = 0;
  for (int c = 0; c < C; c++){
    count[c] = 0;
    for (int n = 0; n < N; n++){
      if (Y[n] == c){
        ind[cit+count[c]] = n;
        count[c]++;
      }
    }
    cit += count[c];
  }

  int flag = 0;
  for (int c1 = 0; c1 < C-1; c1++){
    for (int c2 = c1+1; c2 < C; c2++){
      if (count[c1]==0 && count[c2]==0){
        flag = 1;
        break;
      }
    }
  }

  if (flag == 1){
    *mean = NA_REAL;
    *var = NA_REAL;
  }
  else{
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    for (int c = 0; c < C; c++){
      double val = count[c];
      term1 += val;
      term2 += pow(val, 2.0);
      term3 += pow(val,2.0)*(2.0*val+3.0);
    }
    *mean = (term1*term1 - term2) / 4.0;
    *var = (term1*term1*(2.0*term1+3.0)-term3) / 72.0;
  }
      
  for (int k = 0; k < K; k++){
    int xit = k*N;

    double valin = 0.0;
    double valeq = 0.0;   
    int cit = 0;

    // compare pairs
    for (int c1 = 0; c1 < C-1; c1++){
      for (int t1 = 0; t1 < count[c1]; t1++){
        double x1 = X[xit+ind[cit+t1]];
        for (int n=cit+count[c1]; n < N; n++){
          double x2 = X[xit+ind[n]];
          if (x1 < x2)
            valin++;
          if (x1 == x2)
            valeq++;
        }
      }
      cit += count[c1];
    }

    double stat = (valin + 0.5*valeq - *mean) / sqrt(*var);
    if(isnan(stat))
      stat = NA_REAL;
    if (stat == 0.0){
      if (flag == 1)
	stat = NA_REAL;
    }
    T[k] = stat;
  }

  free(count);
  free(ind);
}

void permpairw(const double *dif, const double *rdif, const int *n, const int *k, const int *b, double *T) // Wilcoxon test
{ 
  int N = *n;
  float Nf = N;
  int K = *k;
  int B = *b;

  // compute test statistics
  double *pr;
  pr = (double *)malloc(K*N*sizeof(double));   
  for (int k = 0; k < K; k++){
    for (int n = 0; n < N; n++){
      if(dif[k*N+n] > 0.0)
        pr[k*N+n] = rdif[k*N+n];
      else
        pr[k*N+n] = 0.0;
    }
  }

  // normalize test statistics
  double m = N*(N+1)/4;
  double s = sqrt(Nf*(Nf+1.0)*(2.0*Nf+1.0)/24.0);
 
  for (int k = 0; k < K; k++){
    double t = 0.0;
    for (int n = 0; n < N; n++){
      t += pr[k*N+n];
    }
    T[k] = (t-m) / s;
  }
    
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  int *yi;
  yi = (int *)malloc(N*sizeof(int));
  for (int b = 1; b <= B; b++){
    // update yi
    for (int n = 0; n < N; n++){
      if(unif_rand() < 0.5)
        yi[n] = -1;
      else
        yi[n] = 1;
    }

    for (int k = 0; k < K; k++){
      for (int n = 0; n < N; n++){
        if(dif[k*N+n]*yi[n] > 0.0)
          pr[k*N+n] = rdif[k*N+n];
        else
          pr[k*N+n] = 0.0;
      }
    }

    for (int k = 0; k < K; k++){
      double t = 0.0;
      for (int n = 0; n < N; n++){
        t += pr[k*N+n];
      }
      T[b*K+k] = (t-m) / s;
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(pr);
  free(yi);
}

void permpairt(const double *dif, const int *n, const int *k, const int *b, double *T)
{ 
  int N = *n;
  int K = *k;
  int B = *b;
  float Nf = N;

  // compute test statistics
  double *mean;
  mean = (double *)malloc(K*sizeof(double)); 
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += dif[k*N+n];
    mean[k] = val / N;
  }

  double *sd;
  sd = (double *)malloc(K*sizeof(double));   
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += pow(dif[k*N+n]-mean[k], 2.0);
    sd[k] = sqrt(val / (N-1.0));
  }

  for (int k = 0; k < K; k++)
    T[k] = sqrt(Nf)*mean[k]/sd[k];
 
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  int *yi;
  yi = (int *)malloc(N*sizeof(int));
  for (int b = 1; b <= B; b++){
    // update yi
    for (int n = 0; n < N; n++){
      if(unif_rand() < 0.5)
        yi[n] = -1;
      else
        yi[n] = 1;
    }
 
    for (int k = 0; k < K; k++){
      double val = 0.0;
      for (int n = 0; n < N; n++)
        val += yi[n]*dif[k*N+n];
      mean[k] = val / N;
    }
  
    for (int k = 0; k < K; k++){
      double val = 0.0;
      for (int n = 0; n < N; n++)
        val += pow(yi[n]*dif[k*N+n]-mean[k], 2.0);
      sd[k] = sqrt(val / (N-1.0));
    }

    for (int k = 0; k < K; k++)
      T[b*K+k] = sqrt(Nf)*mean[k]/sd[k];
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(mean);
  free(sd);
  free(yi);
}

void permpairw_path(const double *dif, const double *rdif, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, double *score_T) // Wilcoxon test
{ 
  int N = *n;
  float Nf = N;
  int K = *k;
  int B = *b;
  int paths = *PATHS;

  // compute test statistics
  double *pr;
  pr = (double *)malloc(K*N*sizeof(double));   
  for (int k = 0; k < K; k++){
    for (int n = 0; n < N; n++){
      if(dif[k*N+n] > 0.0)
        pr[k*N+n] = rdif[k*N+n];
      else
        pr[k*N+n] = 0.0;
    }
  }

  // normalize test statistics
  double m = N*(N+1)/4;
  double s = sqrt(Nf*(Nf+1.0)*(2*Nf+1.0)/24.0);
 
  for (int k = 0; k < K; k++){
    double t = 0.0;
    for (int n = 0; n < N; n++){
      t += pr[k*N+n];
    }
    score_T[(B+1)*paths+k] = (t-m) / s;
  }

  // compute scores
  string method = "average";  // Can also be "min" or "max" or "default"
  vector<double> a(K);
  for (int k = 0; k < K; k++)
    a[k] = fabs(score_T[(B+1)*paths+k]);
    
  vector<double> ranks;
  rankk(a, ranks, method);

  int c = 0;
  for (int p = 0; p < paths; p++){
    double val = 0.0;
    for (int n = 0; n < num_path[p]; n++){
      val += ranks[probeindex[c]-1];
      c++;
    }
    score_T[p] = val;
  }
    
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  int *yi;
  yi = (int *)malloc(N*sizeof(int));
  for (int b = 1; b <= B; b++){
    // update yi
    for (int n = 0; n < N; n++){
      if(unif_rand() < 0.5)
        yi[n] = -1;
      else
        yi[n] = 1;
    }

    for (int k = 0; k < K; k++){
      for (int n = 0; n < N; n++){
        if(dif[k*N+n]*yi[n] > 0.0)
          pr[k*N+n] = rdif[k*N+n];
        else
          pr[k*N+n] = 0.0;
      }
    }

    for (int k = 0; k < K; k++){
      double t = 0.0;
      for (int n = 0; n < N; n++){
        t += pr[k*N+n];
      }
      a[k] = fabs((t-m)/s);
    }

    vector<double> ranks;
    rankk(a, ranks, method);

    int c = 0;
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += ranks[probeindex[c]-1];
        c++;
      }
      score_T[b*paths+p] = val;
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(pr);
  free(yi);
}

void perm_path_w(const double *rdat, int *y, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, const char *gt[], double *score_T)
{ 
  int N = *n;
  int K = *k;
  int B = *b;
  int paths = *PATHS;
  string GT = *gt;

  // compute test statistics
  double *pr;
  pr = (double *)malloc(K*N*sizeof(double));   
  for (int k = 0; k < K; k++){
    for (int n = 0; n < N; n++){
      if(y[n] == 1)
        pr[k*N+n] = rdat[k*N+n];
      else
        pr[k*N+n] = 0.0;
    }
  }

  double n0 = 0.0;
  for (int n = 0; n < N; n++){
    if (y[n] == 0.0)
      n0++;
  }
  double n1 = N - n0;

  // normalize test statistics
  double m = n1 * (n0+n1+1.0) / 2.0; // summing over n1
  double s = sqrt(n0*n1*(n0+n1+1.0)/12.0);
 
  for (int k = 0; k < K; k++){
    double t = 0.0;
    for (int n = 0; n < N; n++){
      t += pr[k*N+n];
    }
    score_T[(B+1)*paths+k] = (t-m) / s;
  }

  // compute scores, mean of test statistics
  vector<double> a(K);
  for (int k = 0; k < K; k++)
    a[k] = score_T[(B+1)*paths+k];

  int c = 0;
  if (GT == "mean"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += a[probeindex[c]-1];
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (GT == "meanabs"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += fabs(a[probeindex[c]-1]);
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (GT == "maxmean"){
    int c = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = a[probeindex[c]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        c++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[p] = posval;
      else
        score_T[p] = negval;
    }
  }
   
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  for (int b = 1; b <= B; b++){
    YpermI(y, N);
    for (int k = 0; k < K; k++){
      for (int n = 0; n < N; n++){
        if(y[n] == 1)
          pr[k*N+n] = rdat[k*N+n];
        else
          pr[k*N+n] = 0.0;
      }
    }

    for (int k = 0; k < K; k++){
      double t = 0.0;
      for (int n = 0; n < N; n++){
        t += pr[k*N+n];
      }
      a[k] = (t-m)/s;
    }

    int c = 0;
    if (GT == "mean"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += a[probeindex[c]-1];
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (GT == "meanabs"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += fabs(a[probeindex[c]-1]);
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (GT == "maxmean"){
      int c = 0;
      for (int p = 0; p < paths; p++){
        double negval = 0.0;
        double posval = 0.0;
        double val = 0.0;
        int neg = 0;
        int pos = 0;

        for (int n = 0; n < num_path[p]; n++){
          val = a[probeindex[c]-1];
          if (val < 0.0){
            negval += val;
            neg++;
          }
          else{
            posval += val;
            pos++;
          }
          c++;
        }

        if (pos != 0)
          posval /= pos;
        if (neg != 0)
          negval /= neg;

        if (posval > fabs(negval))
          score_T[b*paths+p] = posval;
        else
          score_T[b*paths+p] = negval;
      }
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(pr);
}

void perm_path_t(const double *X, int *y, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, const char *a[], double *score_T)
{ 
  int N = *n;
  int K = *k;
  int B = *b;
  int paths = *PATHS;
  string A = *a;

  double *T;
  T = (double *)malloc(K*sizeof(double));
  double *Sum;
  Sum = (double *)malloc(K*sizeof(double));
  double *SumSq;
  SumSq = (double*)malloc(K*sizeof(double));

  for (int k = 0; k < K; k++){
    T[k] = 0.0;
    Sum[k] = 0.0;
    SumSq[k] = 0.0;
  }

  // compute n0 and n1
  int n0 = 0;
  for (int n = 0; n < N; n++){
    if (y[n] == 0)
      n0++;
  }
  int n1 = N - n0;

  // compute N0
  int *N0;
  N0 = (int *)malloc(n0*sizeof(int));
  Index0(y, N0, N);

  SumSumSq(X, Sum, SumSq, N, K);

  // compute T
  Tstat(X, T, Sum, SumSq, N0, N, K, n0, n1);
  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = T[k];

  // compute scores
  int c = 0;
  if (A == "mean"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += T[probeindex[c]-1];
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "meanabs"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += fabs(T[probeindex[c]-1]);
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "maxmean"){
    int c = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = T[probeindex[c]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        c++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[p] = posval;
      else
        score_T[p] = negval;
    }
  }
  
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  for (int b = 1; b <= B; b++){
    YpermI(y, N);
    Index0(y, N0, N);
    Tstat(X, T, Sum, SumSq, N0, N, K, n0, n1);

    int c = 0;
    if (A == "mean"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += T[probeindex[c]-1];
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "meanabs"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += fabs(T[probeindex[c]-1]);
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "maxmean"){
      int c = 0;
      for (int p = 0; p < paths; p++){
        double negval = 0.0;
        double posval = 0.0;
        double val = 0.0;
        int neg = 0;
        int pos = 0;

        for (int n = 0; n < num_path[p]; n++){
          val = T[probeindex[c]-1];
          if (val < 0.0){
            negval += val;
            neg++;
          }
          else{
            posval += val;
            pos++;
          }
          c++;
        }

        if (pos != 0)
          posval /= pos;
        if (neg != 0)
          negval /= neg;

        if (posval > fabs(negval))
          score_T[b*paths+p] = posval;
        else
          score_T[b*paths+p] = negval;
      }
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(T);
  free(Sum);
  free(SumSq);
  free(N0);
}

void perm_path_pearson(const double *X, double *y, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, const char *a[], double *score_T)
{ 
  int N = *n;
  int K = *k;
  int B = *b;
  int paths = *PATHS;
  string A = *a;

  double *T;
  T = (double *)malloc(K*sizeof(double));
  double *meanX;
  meanX = (double *)malloc(K*sizeof(double));
  double *sumquadX;
  sumquadX = (double *)malloc(K*sizeof(double));

  // compute mean of Y
  double meanY = 0.0;
  for (int n = 0; n < N; n++)
    meanY += y[n];
  meanY /= N;

  // compute meanX
  for (int k = 0; k < K; k++){
    meanX[k] = 0.0;
    int iter = k*N;
    for (int n = 0; n < N; n++)
      meanX[k] += X[iter+n];
    meanX[k] /= N;
  }

  // compute sumquadX
  for (int k = 0; k < K; k++){
    sumquadX[k] = 0.0;
    int iter = k*N;
    for (int n = 0; n < N; n++)
      sumquadX[k] += (X[iter+n]-meanX[k]) * (X[iter+n]-meanX[k]);
  }

  // compute sumquadY
  double sumquadY = 0.0;
  for (int n = 0; n < N; n++)
    sumquadY += (y[n]-meanY) * (y[n]-meanY);                                               

  // compute T
  PearsonStat(y, X, T, meanX, sumquadX, N, K, meanY, sumquadY);
  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = T[k];

  // compute scores
  int c = 0;
  if (A == "mean"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += T[probeindex[c]-1];
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "meanabs"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += fabs(T[probeindex[c]-1]);
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "maxmean"){
    int c = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = T[probeindex[c]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        c++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[p] = posval;
      else
        score_T[p] = negval;
    }
  }
  
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  for (int b = 1; b <= B; b++){
    YpermD(y, N);
    PearsonStat(y, X, T, meanX, sumquadX, N, K, meanY, sumquadY);

    int c = 0;
    if (A == "mean"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += T[probeindex[c]-1];
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "meanabs"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += fabs(T[probeindex[c]-1]);
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "maxmean"){
      int c = 0;
      for (int p = 0; p < paths; p++){
        double negval = 0.0;
        double posval = 0.0;
        double val = 0.0;
        int neg = 0;
        int pos = 0;

        for (int n = 0; n < num_path[p]; n++){
          val = T[probeindex[c]-1];
          if (val < 0.0){
            negval += val;
            neg++;
          }
          else{
            posval += val;
            pos++;
          }
          c++;
        }

        if (pos != 0)
          posval /= pos;
        if (neg != 0)
          negval /= neg;

        if (posval > fabs(negval))
          score_T[b*paths+p] = posval;
        else
          score_T[b*paths+p] = negval;
      }
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(T);
  free(meanX);
  free(sumquadX);
}

void perm_path_spearman(const double *X, double *y, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, const char *a[], double *score_T)
{ 
  int N = *n;
  int K = *k;
  int B = *b;
  int paths = *PATHS;
  string A = *a;

  double *T;
  T = (double *)malloc(K*sizeof(double));                                            

  // compute T
  SpearmanStatNew(y, X, T, N, K);
  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = T[k];

  // compute scores
  int c = 0;
  if (A == "mean"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += T[probeindex[c]-1];
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "meanabs"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += fabs(T[probeindex[c]-1]);
        c++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "maxmean"){
    int c = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = T[probeindex[c]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        c++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[p] = posval;
      else
        score_T[p] = negval;
    }
  }
  
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  for (int b = 1; b <= B; b++){
    YpermD(y, N);
    SpearmanStatNew(y, X, T, N, K);

    int c = 0;
    if (A == "mean"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += T[probeindex[c]-1];
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "meanabs"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += fabs(T[probeindex[c]-1]);
          c++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "maxmean"){
      int c = 0;
      for (int p = 0; p < paths; p++){
        double negval = 0.0;
        double posval = 0.0;
        double val = 0.0;
        int neg = 0;
        int pos = 0;

        for (int n = 0; n < num_path[p]; n++){
          val = T[probeindex[c]-1];
          if (val < 0.0){
            negval += val;
            neg++;
          }
          else{
            posval += val;
            pos++;
          }
          c++;
        }

        if (pos != 0)
          posval /= pos;
        if (neg != 0)
          negval /= neg;

        if (posval > fabs(negval))
          score_T[b*paths+p] = posval;
        else
          score_T[b*paths+p] = negval;
      }
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(T);
}

void perm_path_jt(const double *X, int *y, const int *c, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, double *mean, double *var, const int *b, const char *a[], double *score_T)
{ 
  const int N = *n;
  const int K = *k;
  const int B = *b;
  const int C = *c;
  const int paths = *PATHS;
  const string A = *a;

  double *T;
  T = (double *)malloc(K*sizeof(double));                                            

  // compute T
  JT_test_norm(X, y, C, K, N, T, mean, var);
  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = T[k];

  // compute scores
  int ci = 0;
  if (A == "mean"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += T[probeindex[ci]-1];
        ci++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "meanabs"){
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += fabs(T[probeindex[ci]-1]);
        ci++;
      }
      score_T[p] = val / num_path[p]; 
    }
  }
  if (A == "maxmean"){
    int ci = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = T[probeindex[ci]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        ci++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[p] = posval;
      else
        score_T[p] = negval;
    }
  }
  
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  for (int b = 1; b <= B; b++){
    YpermI(y, N);
    JT_test_norm(X, y, C, K, N, T, mean, var);

    int ci = 0;
    if (A == "mean"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += T[probeindex[ci]-1];
          ci++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "meanabs"){
      for (int p = 0; p < paths; p++){
        double val = 0.0;
        for (int n = 0; n < num_path[p]; n++){
          val += fabs(T[probeindex[ci]-1]);
          ci++;
        }
        score_T[b*paths+p] = val / num_path[p];
      }
    }
    if (A == "maxmean"){
      int ci = 0;
      for (int p = 0; p < paths; p++){
        double negval = 0.0;
        double posval = 0.0;
        double val = 0.0;
        int neg = 0;
        int pos = 0;

        for (int n = 0; n < num_path[p]; n++){
          val = T[probeindex[ci]-1];
          if (val < 0.0){
            negval += val;
            neg++;
          }
          else{
            posval += val;
            pos++;
          }
          ci++;
        }

        if (pos != 0)
          posval /= pos;
        if (neg != 0)
          negval /= neg;

        if (posval > fabs(negval))
          score_T[b*paths+p] = posval;
        else
          score_T[b*paths+p] = negval;
      }
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(T);
}

void permpairt_path(const double *dif, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, double *score_T) // Wilcoxon test
{ 
  int N = *n;
  float Nf = N;
  int K = *k;
  int B = *b;
  int paths = *PATHS;

 // compute test statistics
  double *mean;
  mean = (double *)malloc(K*sizeof(double)); 
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += dif[k*N+n];
    mean[k] = val / N;
  }

  double *sd;
  sd = (double *)malloc(K*sizeof(double));   
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += pow(dif[k*N+n]-mean[k], 2.0);
    sd[k] = sqrt(val / (N-1.0));
  }

  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = sqrt(Nf)*mean[k]/sd[k];

  // compute wilcoxon test statistic of probes in each set
  string method = "average";  // Can also be "min" or "max" or "default"
  vector<double> a(K);
  for (int k = 0; k < K; k++)
    a[k] = fabs(score_T[(B+1)*paths+k]);
    
  vector<double> ranks;
  rankk(a, ranks, method);

  int c = 0;
  for (int p = 0; p < paths; p++){
    double val = 0.0;
    for (int n = 0; n < num_path[p]; n++){
      val += ranks[probeindex[c]-1];
      c++;
    }
    score_T[p] = val;
  }
    
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  int *yi;
  yi = (int *)malloc(N*sizeof(int));
  for (int b = 1; b <= B; b++){
    // update yi
    for (int n = 0; n < N; n++){
      if(unif_rand() < 0.5)
        yi[n] = -1;
      else
        yi[n] = 1;
    }

    for (int k = 0; k < K; k++){
     double val = 0.0;
     for (int n = 0; n < N; n++)
       val += yi[n]*dif[k*N+n];
     mean[k] = val / N;
    }
  
    for (int k = 0; k < K; k++){
      double val = 0.0;
      for (int n = 0; n < N; n++)
        val += pow(yi[n]*dif[k*N+n]-mean[k], 2.0);
      sd[k] = sqrt(val / (N-1.0));
    }

    for (int k = 0; k < K; k++)
      a[k] = fabs(sqrt(Nf)*mean[k]/sd[k]);

    vector<double> ranks;
    rankk(a, ranks, method);

    int c = 0;
    for (int p = 0; p < paths; p++){
      double val = 0.0;
      for (int n = 0; n < num_path[p]; n++){
        val += ranks[probeindex[c]-1];
        c++;
      }
      score_T[b*paths+p] = val;
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(mean);
  free(sd);
  free(yi);
}

void permpair_maxmean_path(const double *dif, const int *probeindex, const int *num_path, const int *PATHS, const int *n, const int *k, const int *b, double *score_T) // t-test
{ 
  int N = *n;
  float Nf = N;
  int K = *k;
  int B = *b;
  int paths = *PATHS;

 // compute test statistics
  double *mean;
  mean = (double *)malloc(K*sizeof(double)); 
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += dif[k*N+n];
    mean[k] = val / N;
  }

  double *sd;
  sd = (double *)malloc(K*sizeof(double));   
  for (int k = 0; k < K; k++){
    double val = 0.0;
    for (int n = 0; n < N; n++)
      val += pow(dif[k*N+n]-mean[k], 2.0);
    sd[k] = sqrt(val / (N-1.0));
  }

  for (int k = 0; k < K; k++)
    score_T[(B+1)*paths+k] = sqrt(Nf)*mean[k]/sd[k];

  // compute maxmean score for each set
  vector<double> a(K);
  for (int k = 0; k < K; k++)
    a[k] = score_T[(B+1)*paths+k];

  int c = 0;
  for (int p = 0; p < paths; p++){
    double negval = 0.0;
    double posval = 0.0;
    double val = 0.0;
    int neg = 0;
    int pos = 0;

    for (int n = 0; n < num_path[p]; n++){
      val = a[probeindex[c]-1];
      if (val < 0.0){
        negval += val;
        neg++;
      }
      else{
        posval += val;
        pos++;
      }
      c++;
    }

    if (pos != 0)
      posval /= pos;
    if (neg != 0)
      negval /= neg;

    if (posval > fabs(negval))
      score_T[p] = posval;
    else
      score_T[p] = negval;
  }
    
  // compute permuted replicates
  // random number seed
  GetRNGstate();

  int *yi;
  yi = (int *)malloc(N*sizeof(int));
  for (int b = 1; b <= B; b++){
    // update yi
    for (int n = 0; n < N; n++){
      if(unif_rand() < 0.5)
        yi[n] = -1;
      else
        yi[n] = 1;
    }

    for (int k = 0; k < K; k++){
     double val = 0.0;
     for (int n = 0; n < N; n++)
       val += yi[n]*dif[k*N+n];
     mean[k] = val / N;
    }
 
    for (int k = 0; k < K; k++){
      double val = 0.0;
      for (int n = 0; n < N; n++)
        val += pow(yi[n]*dif[k*N+n]-mean[k], 2.0);
      sd[k] = sqrt(val / (N-1.0));
    }

    for (int k = 0; k < K; k++)
      a[k] = sqrt(Nf)*mean[k]/sd[k];

    int c = 0;
    for (int p = 0; p < paths; p++){
      double negval = 0.0;
      double posval = 0.0;
      double val = 0.0;
      int neg = 0;
      int pos = 0;

      for (int n = 0; n < num_path[p]; n++){
        val = a[probeindex[c]-1];
        if (val < 0.0){
          negval += val;
          neg++;
        }
        else{
          posval += val;
          pos++;
        }
        c++;
      }

      if (pos != 0)
        posval /= pos;
      if (neg != 0)
        negval /= neg;

      if (posval > fabs(negval))
        score_T[b*paths+p] = posval;
      else
        score_T[b*paths+p] = negval;
    }
  }

  // random seed
  PutRNGstate();

  // cleanup
  free(mean);
  free(sd);
  free(yi);
}
}

