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
    
def sieve(n=10000000):
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
    if len(_primes) < 10:
        sieve()
    return _primes

def isPrime(n):
    global _isPrime
    if len(_isPrime) < 10:
        sieve()
    if n < 0:
        return False
    if n < len(_isPrime):
        return _isPrime[n]
    return sigma(n, 0) == 2

def isPrimeList():
    global _isPrime
    if len(_isPrime) < 10:
        sieve()
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
    primes = primeList()
    primev = isPrimeList()
    
    r = {}
    for p in primes:
        if p ** 2 > n:
            break
        if n < len(primev) and primev[n]:
            break
        while n % p == 0:
            r[p] = r.get(p, 0) + 1
            n //= p
    if n != 1:
        r[n] = r.get(n, 0) + 1
    return r

# return the exponents of the prime factorization of n (up to the sieve size)
def primeExponents(n):
    pass
    
# return a sorted list of the divisors of n
def divisors(n):
    pass

# return sum[d|n] d^k
def sigma(n, k):
    v = factorize(n)
    r = 1
    for p, e in v.items():
        if k != 0:
            r *= (p ** (k*(e+1)) - 1) // (p**k - 1)
        else:
            r *= e+1
    return r

# return sum[p^e||n] e^k
def omega(n, k):
    v = factorize(n)
    r = 0
    for p, e in v.items(): 
        r += e ** k
    return r
    
def totient(n):
    v = factorize(n)
    r = n
    for p, e in v.items():
        r *= p-1
        r //= p
    return r

def mobius(n):
    r = 1
    v = factorize(n)
    for p, e in v.items():
        if e > 1:
            return 0
        r *= -1
    return r
    
def radical(n):
    v = factorize(n)
    r = 1
    for p in v:
        r *= p
        
    return r

    
# combinatorics

def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r

def choose(n, r):
    if 2*r > n:
        return choose(n, n-r)
    
    res = 1
    
    for i in range(0, r):
        res *= n-i
        res //= i+1
    
    return res

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
    
    
