/* WARM alpha=1 on the EVEN cycle Z_L: checks the martingale M = sum_e (-1)^e H(N_e - 1).
   Runs R independent replicas to time T; for each prints M, its quadratic variation QV,
   and the final alternating amplitude u (from L*artanh(u) = sum (-1)^e log X_e).
   usage: warm2 L T R seed */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint64_t s[2];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nx(void){uint64_t s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
static inline double u01(void){return (nx()>>11)*0x1.0p-53;}
int main(int c,char**v){
  long L=atol(v[1]); double T=atof(v[2]); long R=atol(v[3]); uint64_t seed=strtoull(v[4],0,10);
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)nx();
  double *N=malloc(sizeof(double)*L);
  printf("# L=%ld T=%g R=%ld seed=%llu\n# rep M QV sumlogratio u_from_logs u_from_altmean\n",L,T,R,(unsigned long long)seed);
  for(long r=0;r<R;r++){
    for(long e=0;e<L;e++)N[e]=1.0;
    double t=0,M=0,QV=0;
    while(t<T){
      t += -log(1.0-u01())/(double)L;
      long vt=(long)(u01()*L); long le=(vt-1+L)%L, re=vt; double A=N[le],B=N[re],S=A+B;
      QV += (1.0/S)*(1.0/A+1.0/B);
      if(u01()*S < A){ M += ((le&1)?-1.0:1.0)/A; N[le]+=1; }
      else           { M += ((re&1)?-1.0:1.0)/B; N[re]+=1; }
    }
    double sl=0,am=0; for(long e=0;e<L;e++){ sl += ((e&1)?-1.0:1.0)*log(N[e]); am += ((e&1)?-1.0:1.0)*N[e]/t; }
    printf("%ld %.6f %.6f %.6f %.6f %.6f\n",r,M,QV,sl,tanh(sl/L),am/L); fflush(stdout);
  }
  return 0;
}
