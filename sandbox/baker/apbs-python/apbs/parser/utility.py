""" This contains some helper functions for parsing APBS input files """
import math
import unittest

def sieve(num):
    """ This implements the Sieve of Eratosthenes as described at
    http://en.wikipedia.org/wiki/Sieve_of_Eratosthenes.  The input is an integer and the output is
    a list of primes up to that integer """
    if num == 1:
        return num
    sqrt_n = int(math.sqrt(num))
    is_prime = [True] * (num+1)
    for i in range(2, sqrt_n+1):
        if is_prime[i]:
            for j in range(i*i, num+1, i):
                is_prime[j] = False
    primes = []
    for i, value in enumerate(is_prime):
        if (i > 1) and value:
            primes.append(i)
    return primes

def factors(num):
    """ Get all factors of n """
    primes = sieve((num//2)+1)
    primes.reverse()
    num_factors = []
    for prime in primes:
        mult = 0
        mnum = prime
        while True:
            if num % mnum == 0:
                mnum = mnum*prime
                num_factors.append(prime)
            else:
                break
        if mult > 0:
            num_factors.append((prime, mult))
    num_factors.sort()
    return num_factors

def product(prodlist):
    """ Return the product of a list """
    prod = 1
    for term in prodlist:
        prod = prod * term
    return prod

class _TestUtilities(unittest.TestCase):
    """ Test the utilities """
    def test_factors(self):
        """ Test the factoring functions """
        in_factors = [5, 5, 2, 2, 2, 2, 2, 2, 2]
        in_factors.sort()
        num = product(in_factors)
        out_factors = factors(num)
        out_factors.sort()
        for fac in in_factors:
            out_factors.remove(fac)
        if len(out_factors) > 0:
            raise IndexError
