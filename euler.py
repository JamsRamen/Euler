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
            
            