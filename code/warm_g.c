/* WARM alpha=1 on a small graph given as an edge list; R replicas to time T, rates 1, initial counts 1.
   Prints, per replica, X_e = N_e / t at the end and the max |vertex sum - 2|.
   usage: warm_g T R seed nV u0 w0 u1 w1 ...  */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint64_t s[2];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nx(void){uint64_t s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
static inline double u01(void){return (nx()>>11)*0x1.0p-53;}
int main(int c,char**v){
  double T=atof(v[1]); long R=atol(v[2]); uint64_t seed=strtoull(v[3],0,10); int nV=atoi(v[4]);
  int nE=(c-5)/2; int eu[64],ew[64]; for(int i=0;i<nE;i++){eu[i]=atoi(v[5+2*i]);ew[i]=atoi(v[6+2*i]);}
  int inc[16][16],deg[16]={0}; for(int i=0;i<nE;i++){inc[eu[i]][deg[eu[i]]++]=i; inc[ew[i]][deg[ew[i]]++]=i;}
  s[0]=seed*0x9E3779B97F4A7C15ULL+1; s[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)nx();
  double N[64];
  for(long r=0;r<R;r++){
    for(int i=0;i<nE;i++)N[i]=1; double t=0;
    while(t<T){
      t += -log(1.0-u01())/nV; int vt=(int)(u01()*nV);
      double S=0; for(int k=0;k<deg[vt];k++)S+=N[inc[vt][k]];
      double x=u01()*S; int k=0; while(k<deg[vt]-1 && x>=N[inc[vt][k]]){x-=N[inc[vt][k]];k++;}
      N[inc[vt][k]]+=1;
    }
    double dev=0; for(int a=0;a<nV;a++){double S=0; for(int k=0;k<deg[a];k++)S+=N[inc[a][k]]/t; if(fabs(S-2)>dev)dev=fabs(S-2);}
    printf("%ld",r); for(int i=0;i<nE;i++)printf(" %.4f",N[i]/t); printf(" dev=%.4f\n",dev);
  }
  return 0;
}
