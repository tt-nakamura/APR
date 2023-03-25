// uses NTL
//   http://www.shoup.net/ntl

#include "Cycl.h"
using namespace NTL;

long Cycl::p;

void Cycl::set(const long *c)
// c = coefficients of polynomial representation
{
    a.SetLength(p-1);
    for(long i=p-2; i>=0; i--) a[i] = c[i] - c[p-1];
    a.normalize();
}

void Cycl::operator*=(const Cycl& x)
// multiplication modulo 1+X+...+X^{p-1}
{
    long i,d;
    if(&x == this) sqr(a,a);
    else a *= x.a;
    d = deg(a);
    if(d<p-1) return;
    for(; d>=p; d--) a[d-p] += a[d];
    for(i=0; i<d; i++) a[i] -= a[d];
    a.SetLength(d);
    a.normalize();
}

long Cycl::IsRoot1()
// return k if *this == e^{2pi ik/p} for some k (0<=k<p)
// return -1 otherwise
{
    long s(0),t(0),j,k;
    for(j=p-2; j>=0; j--) {
        if(IsOne(a[j])) { s++; k=j; }
        else if(a[j] == -1) t++;
        else if(!IsZero(a[j])) return -1;
    }
    if(t==p-1) return t;
    else if(s==1 && t==0) return k;
    else return -1;
}

void power(Cycl& x, const Cycl& a, const ZZ& n)
// x = a^n (mod 1+X+...+X^{p-1})
{
    long i;
    if(IsZero(n)) x=1;
    else if(&x==&a) {
        Cycl b(a);
        for(i=NumBits(n)-2; i>=0; i--) {
            x*=x;
            if(bit(n,i)) x*=b;
        }
    }
    else {
        x=a;
        for(i=NumBits(n)-2; i>=0; i--) {
            x*=x;
            if(bit(n,i)) x*=a;
        }
    }
}

void power(Cycl& x, const Cycl& a, long n)
// x = a^n (mod 1+X+...+X^{p-1})
{
    if(n==0) { x=1; return; }
    long m(1<<(NumBits(n)-1));
    if(&x==&a) {
        Cycl b(a);
        for(m>>=1; m; m>>=1) {
            x*=x;
            if(n&m) x*=b;
        }
    }
    else {
        x=a;
        for(m>>=1; m; m>>=1) {
            x*=x;
            if(n&m) x*=a;
        }
    }
}

void jsum1(Cycl& J, long q, const long *ind, long a)
// Jacobi sum
// J = sum_{j=2}^{q-1} chi_j^a chi_{1-j}
//   where chi_j = e^{2pi i ind[j]/p}
// q = prime number such that p divides q-1
// ind[j] = index of j (mod q) i.e.,
//   g^ind[j] = j (mod q) where g is generator
{
    long i,j(q-1),u,c[Cycl::p];
    for(i=0; i<Cycl::p; i++) c[i]=0;
    for(i=2; i<q; i++, j--) {
        u =ind[i]; u*=a;
        u+=ind[j]; u%=Cycl::p;
        c[u]++;
    }
    J.set(c);
}

void jsum2(Cycl& J, long q, const long *ind, long a)
// Jacobi sum
// J = sum_{j=2}^{q-1} chi_j^a chi_{1-j}^a
{
    long i,j(q-1),k((q+1)>>1),u,c[Cycl::p];
    for(i=0; i<Cycl::p; i++) c[i]=0;
    for(i=2; i<k; i++, j--) {
        u=ind[i]+ind[j]; u*=a; u%=Cycl::p;
        c[u]++;
    }
    for(i=0; i<Cycl::p; i++) c[i] <<= 1;
    u=ind[k]; u*=a; u<<=1; u%=Cycl::p;
    c[u]++;
    J.set(c);
}

void gsum(Cycl& G, long q, const long *ind, long k)
// Gauss sum
// if k==p, set G = G_1^p
// if 1<=k<p, set G = G_1^k/G_k
//   where G_k = sum_{j=1}^{q-1} chi_j^k e^{2pi ij/q}
{
    if(k<=1) G=1;
    else if(Cycl::p==2) G=-q;
    else {
        Cycl t;
        long j(k-1);
        while((j&1)==0) j>>=1;
        gsum(G,q,ind,j);
        while(j+1 < k) {
            jsum2(t,q,ind,j);
            G*=G;
            G*=t;
            j<<=1;
        }
        if(k==Cycl::p) G*=q;
        else {
            jsum1(t,q,ind,j);
            G*=t;
        }
    }
}