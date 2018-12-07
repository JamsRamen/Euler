import math

# number theory

def gcd(a, b):
    if a < 0:
        a = -a;
    if b < 0:
        b = -b;
    
    if b == 0:
        return a
    return gcd(b, a%b);

def sieve(n):
    isPrime = [True] * n
    isPrime[0] = isPrime[1] = False
    
    i = 2
    while i*i < n:
        if isPrime[i]:
            j = i*i
            while j < n:
                isPrime[j] = False
                j += i
        i += 1
        
    primes = []
    for i in range(n):
        if isPrime[i]:
            primes.append(i)
    
    return isPrime, primes

def modpow (b, e, m):
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

# return the exponents of the prime factorization of n
# return a map of primes to their exponents
def factorize(n):
    pass
    
# return a sorted list of the divisors of n
def divisors(n):
    pass

# return sum[p^e||n] e^k
def sigma(n, k):
    return 
    
def totient(n):
    pass

def mobius(n):
    pass

    
    
# combinatorics

def factorial(n):
    pass

def choose(n, r):
    pass

def permute(n, r):
    pass


# computer science

# note: nextPermutation does not work with strings
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
    ans = [chr(i) for i in v]
    v = "".join(ans)
    return True


