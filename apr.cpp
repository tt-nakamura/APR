// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZX.h>
#include "Cycl.h"
using namespace NTL;

long IsPrime(long);
void CRT(long& x, long& m, long a, long p);

int comp(const void *a, const void *b) {
    return *(int*)a - *(int*)b;
}

long apr(const ZZ& n)
// primality proving by Adleman-Pomerance-Rumely method
// return 1 or 0 if n is prime or not, respectively
// assume n>=2
// reference: R. Crandall and C. Pomerance
//   "Prime Numbers: A Computational Perspective"
//     2nd edition, Algorithm 4.4.5
{
    long i,j,k,h,m,s,r,M,K,L;
    ZZ t,d,N;
    SqrRoot(t,n);
    K = long(log(log(n))/log(2) - 1);
    L = NumBits(t);
    if(K<1) K=1;
    long p[K<<1], q[L];// small primes
    PrimeSeq ps;
    for(i=0; i<K; i++) p[i] = ps.next();
    for(i=m=0, k=(1<<(K-1));; k<<=1, K++) {
        for(; i<k; i++) {
            for(h=i, j=1, r=2; h; h>>=1, j++)
                if(h&1) r *= p[j];
            r++;
            if(!(r>>31) && IsPrime(r)) q[m++] = r;
        }
        qsort(q, m, sizeof(long), comp);
        for(L=0, N=1; L<m;)
            if((N *= q[L++]) > t) goto a;
        p[K] = ps.next();
    }
a:  ;
    for(i=0; i<K; i++) if(n==p[i]) return 1;
    for(i=0; i<L; i++) if(n==q[i]) return 1;
    for(i=0; i<K; i++) if(divide(n,p[i])) return 0;
    GCD(d,n,N); if(d>1) return 0;
    bool c[K], pq[K][L];// table of divisibility
    for(i=0; i<K; i++) {
        c[i] = false;
        for(j=0; j<L; j++) {
            if((q[j]-1)%p[i]) pq[i][j] = false;
            else c[i] = pq[i][j] = true;
        }
    }
    Vec<long> ind[L];
    long g[L];
    for(j=0; j<L; j++) {
        h = q[j]-1;
        for(g[j]=2;; g[j]++) {// generator mod q
            for(i=0; i<K; i++)
                if(pq[i][j] && 
                   PowerMod(g[j], h/p[i], q[j]) == 1) break;
            if(i==K) break;
        }
        ind[j].SetLength(q[j]);
        for(k=0, r=1; k<h; k++) {// index mod q
            ind[j][r] = k;
            r *= g[j];
            r %= q[j];
        }
    }
    ZZ_p::init(n);
    Cycl G,H[L];
    ZZX f;
    long w[L], l[K][L];
    for(i=0; i<K; i++) {
        if(!c[i]) continue;
        Cycl::init(p[i]);
        h = p[i]-1;
        power(t,n,h); t--;
        for(s=0; divide(t,t,p[i]); s++);
        for(j=m=0; j<L; j++) {
            if(!pq[i][j]) continue;
            gsum(G, q[j], ind[j].data());// Gauss sum
            power(G,G,t);
            for(k=1; k<s; k++) {
                if((l[i][j] = G.IsRoot1()) >= 0) break;
                H[j] = G;
                power(G,G,p[i]);
            }
            if(k==s && (l[i][j] = G.IsRoot1()) < 0) return 0;
            if((w[j] = k) > m) m = k;
        }
        if(m==1) continue;
        for(j=0, k=L; j<L; j++) {
            if(!pq[i][j]) continue;
            if(w[j] < m) l[i][j] = 0;
            else if(l[i][j]) k=-1;
            else if(k==L) k=j;
        }
        if(k<0) continue;
        conv(f, rep(H[k]));
        if((m=deg(f)+1) < h) m--;
        for(j=0; j<=m; j++) {
            if(j) f[j-1]++;
            if(j<h) f[j]--;
            else for(k=0; k<h; k++) f[k]++;
            content(t,f);
            GCD(d,t,n);
            if(d>1) return 0;
        }
    }
    for(j=0, t=0, N=1; j<L; j++) {
        for(i=s=0, m=1; i<K; i++)
            if(pq[i][j]) CRT(s, m, l[i][j], p[i]);
        m = PowerMod(g[j], s, q[j]);
        CRT(t,N,m,q[j]);
    }
    if(sign(t) < 0) t+=N;
    for(i=1, M=2; i<K; i++) if(c[i]) M *= p[i];
    for(r=1, d=1; r<M; r++) {
        MulMod(d,d,t,N);
        if(divide(n,d) && d>1 && d<n) return 0;
    }
    return 1;
}