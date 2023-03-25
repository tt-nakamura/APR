#include<NTL/ZZ.h>
using namespace NTL;

long apr(const ZZ&);

main() {
    long n(5000),l(64),i(0),j(0),k(0),m;
    ZZ p;
    while(i+j<n) {
        RandomBits(p,l);
        if(p<2) continue;
        if(m = apr(p)) i++; else j++;
        if(m != ProbPrime(p)) k++;
    }
    std::cout << (k ? "failure" : "successful");
    std::cout << std::endl;
    std::cout << "primes: " << i << std::endl;
    std::cout << "composites: " << j << std::endl;
}