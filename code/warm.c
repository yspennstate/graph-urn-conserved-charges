/* WARM, alpha=1, on the cycle Z_L (edge e joins vertices e and e+1 mod L).
   Every vertex fires at rate 1; vertex v reinforces edge v-1 (left) or v (right)
   with probability proportional to the current counts. Initial counts 1.
   Prints, at log-spaced times t: mean, sd, min, max of N_e/t, and the
   alternating-mode amplitude A = (1/L) sum (-1)^e N_e / t (even L only),
   plus the lag-1 correlation of N_e/t - mean.
   usage: warm L Tmax seed */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint64_t s[2];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nx(void){uint64_t s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
static inline double u01(void){return (nx()>>11)*0x1.0p-53;}
int main(int c,char**v){
  long L=atol(v[1]); double T=atof(v[2]); uint64_t seed=strtoull(v[3],0,10);
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)nx();
  double *N=malloc(sizeof(double)*L); for(long e=0;e<L;e++)N[e]=1.0;
  double t=0, next=1.0; long long steps=0;
  printf("# L=%ld T=%g seed=%llu\n# t mean sd min max altamp lag1corr\n",L,T,(unsigned long long)seed);
  while(t<T){
    t += -log(1.0-u01())/(double)L;  /* superposition of L rate-1 clocks */
    long vtx = (long)(u01()*L); long le=(vtx-1+L)%L, re=vtx;
    if(u01()*(N[le]+N[re]) < N[le]) N[le]+=1; else N[re]+=1;
    steps++;
    if(t>=next){
      double m=0,q=0,mn=1e300,mx=-1e300,alt=0,lag=0;
      for(long e=0;e<L;e++){double x=N[e]/t;m+=x;q+=x*x;if(x<mn)mn=x;if(x>mx)mx=x;alt+=(e&1?-x:x);}
      m/=L; double var=q/L-m*m;
      for(long e=0;e<L;e++)lag+=(N[e]/t-m)*(N[(e+1)%L]/t-m);
      lag/=L;
      printf("%.4g %.6f %.6f %.6f %.6f %.6f %.4f\n",t,m,sqrt(var>0?var:0),mn,mx,alt/L,var>0?lag/var:0);
      fflush(stdout); next*=1.25;
    }
  }
  return 0;
}
