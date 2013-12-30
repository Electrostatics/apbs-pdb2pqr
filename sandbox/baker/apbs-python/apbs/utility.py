""" This contains some helper functions for parsing APBS input files """
import math
from operator import mul

def sieve(n):
    """ This implements the Sieve of Eratosthenes as described at http://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
    The input is an integer and the output is a list of primes up to that integer """
    if n == 1:
        return n
    sqrt_n = int(math.sqrt(n))
    isPrime = [True] * (n+1)
    for i in range(2, sqrt_n+1):
        if isPrime[i]:
            for j in range(i*i, n+1, i):
                isPrime[j] = False
    primes = []
    for i, value in enumerate(isPrime):
        if (i>1) and value:
            primes.append(i)
    return primes

def factors(n):
    """ Get all factors of n """
    primes = sieve((n//2)+1)
    primes.reverse()
    factors = []
    for prime in primes:
        mult = 0
        m = prime
        while True:
            if n%m == 0:
                m = m*prime
                factors.append(prime)
            else:
                break
        if mult > 0:
            factors.append((prime, mult))
    factors.sort()
    return factors

def product(prodlist):
    """ Return the product of a list """
    return reduce(mul,prodlist,1)

if __name__ == '__main__':
    inFactors = [5, 5, 2, 2, 2, 2, 2, 2, 2]
    inFactors.sort()
    n = product(inFactors)
    outFactors = factors(n)
    outFactors.sort()
    print "The factors of %d are: %s" % (n, outFactors)
    print "Checking the answer..."
    for x in inFactors:
        outFactors.remove(x)
    if len(outFactors) > 0:
        raise IndexError
    else:
        print "Everything is OK"