def xor(a:bytes, b:bytes):
     return bytes([ai^bi for (ai,bi) in zip(a,b)])

def gcd(a:int, b:int):
    # Euclidean algorithm
    while b != 0:
        (a,b) = (b, a%b)
    return a

def gcdext(a:int, b:int):
    # Extended Euclidean algorithm
    (sold, snew) = (1,0)
    (told, tnew) = (0,1)
    while b != 0:
        q = a // b
        (a,b) = (b, a-q*b)
        (sold, snew) = (snew, sold - q*snew)
        (told, tnew) = (tnew, told - q*tnew)
    return (sold, told)

def modinv(a:int, p:int):
    s,t = gcdext(a,p)
    return s % p #take it mod p again just in case s is negative

def legendre(a:int, p:int):
    # Legendre symbol via Euler's Criterion
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def modsqrt(a:int, p:int):
    #Tonelli-Shanks algorithm
    #Step 1: Write p-1 = Q*2^M
    (Q,M) = (p-1,0)
    while Q%2 == 0:
        (Q,M) = (Q//2, M+1)
    
    #Step 2: Find a nonresidue z
    z = 2
    while pow(z, (p-1)//2, p) == 1:
        z += 1
	c = pow(z, Q, p)

    #Step 3: Iterate
    R = pow(a, (Q+1)//2, p)
    t = pow(a, Q, p)
    # R^2 = a^(Q+1) = a*a^Q, we eventually want a^Q=1modp so R^2=amodp

    if t == 0: #handle the case if a is a nonresidue
    	return 0

    while t != 1:
        # Find i so t^(2^i) = 1modp
        # t^2^(M-1) = a^(Q*2^(M-1)) = a^((p-1)/2) = 1 by Euler's Criterion
        ind = -1
        tt = t #create secondary variable for t that we can manipulate
        for i in range(1,M):
            tt = pow(tt,2,p)
            if tt == 1:
                ind = i
                break
            if ind == -1: #no i found, therefore a is a nonsqrt
                return 0

        b = pow(c, pow(2, M - 1 - i), p)
        M = i
        c = (b*b) % p
        t = (t*b*b) % p
        R = (R*b) % p
    return R

from functools import reduce
def chinrem(n, a):
    sum = 0
    prod = reduce(lambda a, b: a*b, n)
    for n_i, a_i in zip(n,a):
        p = prod/n_i
        sum += a_i * modinv(p, n_i) * p
    return sum % prod
