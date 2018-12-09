import math

# number theory

def gcd(*arg):
    if len(arg) == 2:
        a, b = arg
        if a < 0:
            a = -a
        if b < 0:
            b = -b
        
        if b == 0:
            return a
        return gcd(b, a%b)
    
    r = gcd(arg[0], arg[1])
    for i in range(2, len(arg)):
        r = gcd(r, arg[i])
    return r

_isPrime = []
_primes = []
    
def sieve(n):
    global _isPrime, _primes
    _isPrime = [True] * n
    _isPrime[0] = _isPrime[1] = False
    
    i = 2
    while i*i < n:
        if _isPrime[i]:
            j = i*i
            while j < n:
                _isPrime[j] = False
                j += i
        i += 1
        
    _primes = []
    for i in range(n):
        if _isPrime[i]:
            _primes.append(i)

def primeList():
    global _primes
    return _primes

def isPrime(n):
    global _isPrime
    if n >= 0 and n < len(_isPrime):
        return _isPrime[n]
    return None

def isPrimeList():
    global _isPrime
    return _isPrime
            
def modpow(b, e, m):
    r = 1
    while e > 0:
        if e % 2 == 1:
            r = r*b % m
        b = b*b % m
        e //= 2
    return r

# run Miller-Rabin primality test k times
def probablePrime(n, k):
    pass

# return a map from primes to exponents in the prime factorization of n
def factorize(n):
    pass

# return the exponents of the prime factorization of n (up to the sieve size)
def primeExponents(n):
    pass
    
# return a sorted list of the divisors of n
def divisors(n):
    pass

# return sum[p^e||n] e^k
def sigma(n, d):
    pass
    
def totient(n):
    pass

def mobius(n):
    r = 1
    factorization = factorize(n)
    for p, e in factorization:
        if e > 1:
            return 0
        r *= -1
    return r
    
    
# combinatorics

def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r

def choose(n, r):
    pass

def permute(n, r):
    r = 1
    for i in range(n-r, n+1):
        r *= i
    return r


# computer science

# note: nextPermutation and prevPermutation do not work with strings
def nextPermutation(v):
    i = len(v) - 1
    while i > 0 and v[i-1] >= v[i]:
        i = i - 1
    if i <= 0:
        return False

    j = len(v) - 1
    while v[j] <= v[i -1]:
        j = j - 1
    v[i-1],v[j] = v[j],v[i-1]

    v[i:] = v[len(v) - 1 : i-1 : -1]
    return True

def prevPermutation(v):
    i = len(v) - 1
    while i > 0 and v[i-1] <= v[i]:
        i = i - 1
    if i <= 0:
        return False

    j = len(v) - 1
    while v[j] >= v[i -1]:
        j = j - 1
    v[i-1],v[j] = v[j],v[i-1]

    v[i:] = v[len(v) - 1 : i-1 : -1]
    return True
    
    
