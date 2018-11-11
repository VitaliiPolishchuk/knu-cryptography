from long_ariphmetic import *
import random
from fractions import gcd
from collections import Counter

def bsgs(g, h, p):
    '''
    Solve for x in h = g^x mod p given a prime p.
    If p is not prime, you shouldn't use BSGS anyway.
    '''
    N = int(sqrt(sub(p, 1))) # phi(p) is p-1 if p is prime

    # Store hashmap of g^{1...m} (mod p). Baby step.
    tbl = {deg_mod(g, i, p):i for i in range(N)}

    # Precompute via Fermat's Little Theorem
    c = deg_mod(g, mul(N, sub(p, 2)), p)

    # Search for an equivalence in the table. Giant step.
    for j in range(N):
        y = mul_mod(h, deg_mod(c, j, p), p)
        if y in tbl: return add(mul(j, N), tbl[y])

    # Solution not found
    return None

# print(bsgs(5, 16, 59))

def miller_rabin(n, k = 10):

    # Implementation uses the Miller-Rabin Primality Test
    # The optimal number of rounds for this test is 40
    # See http://stackoverflow.com/questions/6325576/how-many-iterations-of-rabin-miller-should-i-use-for-cryptographic-safe-primes
    # for justification

    # If number is even, it's a composite number

    if eq(n,2):
        return True

    if eq(mod(n, 2), 0):
        return False

    r, s = 0, sub(n, 1)
    while eq(mod(s, 2),0):
        r = add(r, 1)
        s = div(s, 2)
    for _ in range(k):
        a = random.randrange(2, sub(n, 1))
        x = deg_mod(a, s, n)
        if x == 1 or x == sub(n, 1):
            continue
        for _ in range(sub(r, 1)):
            x = deg_mod(x, 2, n)
            if x == sub(n, 1):
                break
        else:
            return False
    return True

def mobius(n):
    pr = []
    sz = add(n, 1)
    lpd = [i for i in range(sz)]
    mi = [-1 for i in range(sz)]
    mi[1] = 1
 
    for i in range(2, int(add(div(n, 2),1))):
        cp = lpd[i]
        if eq(cp, i):
            pr.append(i)
        for p in pr:
            j = mul(i, p)
            if greater(j, n):
                break
            lpd[j] = p
 
            if less(p, cp):    
                mi[j] = -mi[i]
            else:
                mi[j] = 0
                break
    return mi[n]

def eulers_loop(number):
    coprime_ints = [1]
    for i in range(2, number):
        if eq(gcd(number, i), 1):
            coprime_ints.append(i)
    return len(set(coprime_ints))

def count_euler(number):
    if not isinstance(number, int):
        raise TypeError('Euler\'s totient finction could not be called for type: {0} of object {1}'.format(
            type(number), number
        ))
    elif less(number, 1):
        raise ValueError('Euler\'s function could be called for natural number')
    elif eq(number, 1):
        return 1
    elif miller_rabin(number):
        return sub(number, 1)
    for num in range(sub(number, 1)):
        return eulers_loop(number)

# print(miller_rabin(17,100))
# print(count_euler(46))
# print(mobius(4))

def jacobi(a, n):
  """Jacobi symbol"""

  # Based on the Handbook of Applied Cryptography (HAC), algorithm 2.149.

  # This function has been tested by comparison with a small
  # table printed in HAC, and by extensive use in calculating
  # modular square roots.

  assert greater(n, 3)
  assert eq(mod(n, 2), 1)
  a = mod(a, n)
  if eq(a, 0):
    return 0
  if eq(a, 1):
    return 1
  a1, e = a, 0
  while eq(mod(a1, 2), 0):
    a1, e = int(div(a1, 2)), add(e, 1)
  if eq(mod(e, 2), 0) or eq(mod(n, 8), 1) or eq(mod(n, 8), 7):
    s = 1
  else:
    s = -1
  if eq(a1, 1):
    return s
  if eq(mod(n, 4), 3) and eq(mod(a1, 4), 3):
    s = -s
  return s * jacobi(mod(n, a1), a1)

# print(jacobi(10, 77))

# first
def primesbelow(N):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    #""" Input N>=6, Returns a list of primes, 2 <= p < N """
    correction = greater(mod(N, 6), 1)
    N = {0:N, 1:sub(N,1), 2:add(N,4), 3:add(N,3), 4:add(N,2), 5:add(N,1)}[mod(N,6)]
    sieve = [True] * (int(div(N, 3)))
    sieve[0] = False
    for i in range(add(int(div(int(deg(N, .5)), 3)), 1)):
        if sieve[i]:
            k = (add(mul(3, i), 1)) | 1
            sieve[int(div(mul(k,k), 3))::mul(2, k)] = mul([False], int(add(div(sub(sub(div(N,6), int(div(mul(k,k),6))), 1), k), 1)))
            sieve[int(div(sub(add(mul(k,k), mul(4,k)), mul(mul(2,k),(mod(i,2)))), 3))::mul(2,k)] = mul([False], int(add(div(sub(sub(div(N, 6), int(div((  sub(add(mul(k,k), mul(4,k)), mul(mul(2,k),(mod(i,2))))  )  ,6))), 1), k), 1)))
            # sieve[int(div(mul(k,k), 3))::mul(2, k)] = [False] * (add(int(div(sub(sub(int(div(N,6)), int(div(mul(k,k), 6))), 1), k)), 1))
            # sieve[int(div(sub(add(mul(k,k), mul(4,k)), mul(mul(2,k),(mod(i,2)))), 3))::mul(2,k)] = [False] * ((N // 6 - (k*k + 4*k - 2*k*(i%2))//6 - 1) // k + 1)
    return [2, 3] + [(add(mul(3, i), 1)) | 1 for i in range(1, int(div(N,3)) - correction) if sieve[i]]

smallprimeset = set(primesbelow(100000))
_smallprimeset = 100000
def isprime(n, precision=7):
    # http://en.wikipedia.org/wiki/Miller-Rabin_primality_test#Algorithm_and_running_time
    if less(n, 1):
        raise ValueError("Out of bounds, first argument must be > 0")
    elif less(n, 4):
        return greater(n, 1)
    elif eq(mod(n,2), 0):
        return False
    elif less(n, _smallprimeset):
        return n in smallprimeset


    d = sub(n, 1)
    s = 0
    while eq(mod(d, 2), 0):
        d = int(div(d, 2))
        s = add(s, 1)

    for repeat in range(precision):
        a = random.randrange(2, sub(n, 2))
        x = deg_mod(a, d, n)

        if eq(x, 1) or eq(x, sub(n, 1)): continue

        for r in range(sub(s, 1)):
            x = deg_mod(x, 2, n)
            if eq(x, 1): return False
            if eq(x, sub(n, 1)): break
        else: return False

    return True

smallprimes = primesbelow(1000) # might seem low, but 1000*1000 = 1000000, so this will fully factor every composite < 1000000
def primefactors(n, sort=False):
    factors = []

    for checker in smallprimes:
        while eq(mod(n, checker), 0):
            factors.append(checker)
            n = int(div(n, checker))
        if greater(checker, n): break

    if less(n, 2): return factors

    while greater(n, 1):
        if miller_rabin(n):
            factors.append(n)
            break
        factor = pollard_brent(n) # trial division did not fully factor, switch to pollard-brent
        print(factor)
        factors.extend(primefactors(factor)) # recurse to factor the not necessarily prime factor returned by pollard-brent
        n = int(div(n, factor))

    if sort: factors.sort()

    return factors

def factorization(n):
    factors = {}
    for p1 in primefactors(n):
        try:
            factors[p1] = add(factors[p1], 1)
        except KeyError:
            factors[p1] = 1
    return factors.keys()

# print(factorization(3237))

def calculate_legendre(a, p):
    """
    Calculate the legendre symbol (a, p) with p is prime.
    The result is either -1, 0 or 1
    >>> calculate_legendre(3, 29)
    -1
    >>> calculate_legendre(111, 41) # Beispiel aus dem Skript, S. 114
    -1
    >>> calculate_legendre(113, 41) # Beispiel aus dem Skript, S. 114
    1
    >>> calculate_legendre(2, 31)
    1
    >>> calculate_legendre(5, 31)
    1
    # http://math.stackexchange.com/q/221223/6876
    >>> calculate_legendre(150, 1009)
    1
    # http://math.stackexchange.com/q/221223/6876
    >>> calculate_legendre(25, 1009)
    1
    # http://math.stackexchange.com/q/221223/6876
    >>> calculate_legendre(2, 1009)
    1
    # http://math.stackexchange.com/q/221223/6876
    >>> calculate_legendre(3, 1009)
    1
    """
    if greater(a, sub(p,1)) or less(a, 0):
        return calculate_legendre(mod(a, p), p)
    elif eq(a, 0) or eq(a, 1):
        return a
    elif eq(a, 3):
        if eq(mod(p, 8), 1) or eq(mod(p, 8), 7):
            return 1
        else:
            return -1
    elif eq(a, sub(p,1)):
        if eq(mod(p, 4), 1):
            return 1
        else:
            return -1
    elif not is_prime(a):
        factors = factorization(a)
        product = 1
        for pi in factors:
            product = mul(product, calculate_legendre(pi, p))
        return product
    else:
        if eq(mod(div(sub(p,1),2), 2), 0) or eq((mod(div(sub(a,1),2), 2)), 0):
            return calculate_legendre(p, a)
        else:
            return mul(-1, calculate_legendre(p, a))

def is_prime(a):
    """
    Check if `a` is a prime number.
    Parameters
    ----------
    a : int, a >= 2
    """
    return all(a % i for i in range(2, a))

# print(calculate_legendre(3, 19))

""" This is the main function used in the Cipolla-Lehmer algorithm for
calculating square roots mod a prime p """
def chipolla(c,b,p):
    t1 = mul_mod(b,b,p)
    t2 = mul_mod(4,c,p)
    t2 = sub_mod(p,t2,p)
    g = add_mod(t1,t2,p)
    e = div(sub(p, 1), 2)
    h = exp1(e,g,p)
    s = 1
    if (eq(h,0) or eq(h,1)):
        s = 0
    e = add(e,1)
    t1 = mul_mod(sub(p,1),b,p)
    t2 = mod(c,p)
    q = [t1,t2]
    a = [1,0]
    b = exp2(e,a,q,p)
    t = mul(s,b[1])
    return(t)

# print(chipolla(2575*2575, 100000000, 122342300429))

def el_Gamal():
    print("ElGamal based Elliptic Curve Cryptography\nby Dhruv Dixit:\t 15BCE1324\nVIT University, Chennai\n\nElliptic Curve General Form:\t y^2 mod n=(x^3  + a*x + b)mod n\nEnter 'n':")

    def polynomial(LHS,RHS,n):
        for i in range(0,n):
            LHS[0].append(i)
            RHS[0].append(i)
            LHS[1].append(mod(add(add(mul(mul(i,i),i), mul(a,i)), b),n))
            RHS[1].append(mul_mod(i,i,n))


    def points_generate(arr_x,arr_y,n):
        count=0
        for i in range(0,n):
            for j in range(0,n):
                if(LHS[1][i]==RHS[1][j]):
                    count+=1
                    arr_x.append(LHS[0][i])
                    arr_y.append(RHS[0][j])
        return count

    #main
    n=int(input())
    LHS=[[]]
    RHS=[[]]
    LHS.append([])
    RHS.append([])
    print("Enter value of 'a':")
    a=int(input())
    print("Enter value of 'b':")
    b=int(input())

    #Polynomial
    polynomial(LHS,RHS,n)

    arr_x=[]
    arr_y=[]
    #Generating base points
    count=points_generate(arr_x,arr_y,n)
        
    #Print Generated Points
    print("Generated points are:")
    for i in range(0,count):
        print(add(i,1)," (",arr_x[i],",",arr_y[i],")\n")



    #Calculation of Base Point
    bx=arr_x[0]
    by=arr_y[0]
    print("Base Point taken is:\t(",bx,",",by,")\n")


    print("Enter the random number 'd' i.e. Private key of Sender (d<n):")
    d=int(input())
    if(d>=n):
        print("'d' should be less than 'n'.")
    else:
        #Q i.e. sender's public key generation
        Qx=d*bx
        Qy=d*by
        print("Public key of sender is:\t(",Qx,",",Qy,")\n")

        #Encrytion
        print("Enter the random number 'k' (k<n):\n")
        k=int(input())
        if(greater(k, sub(n,1))):
            print("'k' should be less than 'n'")
        else:
            print("Enter the message to be sent:\n")
            M=int(input())

            #Cipher text 1 generation
            C1x=mul(k,bx)
            C1y=mul(k,by)
            print("Value of Cipher text 1 i.e. C1:\t(",C1x,",",C1y,")\n")

            #Cipher text 2 generation
            C2x=add(mul(k,Qx), M)
            C2y=add(mul(k,Qy),M)
            print("Value of Cipher text 2 i.e. C2:\t(",C2x,",",C2y,")\n")

        #Decryption
            Mx=sub(C2x,mul(d,C1x))
            My=sub(C2y,mul(d,C1y))
            print("The message recieved by reciever is:\t",Mx)  
el_Gamal()