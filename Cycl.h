// uses NTL
//   http://www.shoup.net/ntl

#ifndef __Cycl_h__
#define __Cycl_h__

#include<NTL/ZZ_pX.h>
using namespace NTL;

struct Cycl {// cyclotomic integer
    ZZ_pX a;// polynomial representation
    static long p;// degree(modulus)+1 (prime)
    inline static void init(long p1) { p=p1; }
    inline void operator=(long x) { a=x; }
    inline void operator*=(long x) { a*=x; }
    void set(const long*);
    void operator*=(const Cycl&);
    long IsRoot1();
};

inline const ZZ_pX& rep(const Cycl& x) { return x.a; }

void power(Cycl&, const Cycl&, const ZZ&);
void power(Cycl&, const Cycl&, long);
void jsum1(Cycl&, long, const long*, long);
void jsum2(Cycl&, long, const long*, long);
void gsum(Cycl&, long, const long*, long=Cycl::p);

#endif // __Cycl_h__