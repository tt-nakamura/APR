// uses NTL
//   http://www.shoup.net/ntl

#include<NTL/ZZ.h>
using namespace NTL;

long IsPrime(long n)
// primality proving by trial division
// return 1 or 0 if n is prime or not, respectively
// assume n>=2
{
    long p,m;
    if(n<2) return 0;
    m = SqrRoot(n);
    PrimeSeq ps;
    while((p = ps.next()) && p<=m)
        if(n%p == 0) return 0;
    if(p) return 1;
    else Error("too large n");
}

void CRT(long& x, long& m, long a, long p)
// incremental chinese remaindering
// compute y such that y=x (mod m), y=a (mod p), 0<=y<mp
//   and set x:=y, m:=mp
// assume 0<=a<p, p is coprime to m
{
    long b(m);
    a = SubMod(a, x%p, p);
    b *= MulMod(a, InvMod(m%p, p), p);
    x += b;
    m *= p;
}