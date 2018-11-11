import math 

def gauss(A):
    m = len(A)
    assert all([len(row) == m + 1 for row in A[1:]]), "Matrix rows have non-uniform length"
    n = m + 1
    
    for k in range(m):
        pivots = [abs(A[i][k]) for i in range(k, m)]
        i_max = pivots.index(max(pivots)) + k
        
        # Check for singular matrix
        assert A[i_max][k] != 0, "Matrix is singular!"
        
        # Swap rows
        A[k], A[i_max] = A[i_max], A[k]

        
        for i in range(k + 1, m):
            f = A[i][k] / A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[k][j] * f

            # Fill lower triangular matrix with zeros:
            A[i][k] = 0
    
    # Solve equation Ax=b for an upper triangular matrix A         
    x = []
    for i in range(m - 1, -1, -1):
        x.insert(0, A[i][m] / A[i][i])
        for k in range(i - 1, -1, -1):
            A[k][m] -= A[k][i] * x[0]
    return x

def add(a, b):
    return a + b

def sub(a, b):
    return a - b

def mul(a, b):
    return a * b

def div(a, b):
    return a / b

def deg(a, b):
    return a ** b

def mod(a, m):
    return a % m

def add_mod(a, b, m):
    return (a + b) % m

def sub_mod(a, b, m):
    return (a - b) % m

def mul_mod(a, b, m):
    return (a * b) % m

def div_mod(a, b, m):
    return (a / b) % m

def deg_mod(a, b, m):
    return (a ** b) % m

def eq(a, b):
    return a == b

def greater(a, b):
    return a > b

def less(a, b):
    return a < b

def sqrt(a):
    return round(math.sqrt(a), 0)

""" exp2 exponentiates in GF(p^2) by calculating (g(x))^e mod <q(x)> """
def exp2(e,g,q,n):
    t = [0,1]
    sq = g
    e1 = e
    while(e1!=0):
        if (e1%2)==1:
            t = mult2(sq,t,q,n)
            e1 = (e1-1)//2
        else:
            e1 = e1//2
        sq = mult2(sq,sq,q,n)
    return(t)

""" This calculates g^e (mod n) for integers e,g, and n """
def exp1(e,g,n):
    t = 1
    sq = g
    e1 = e
    while(e1!=0):
        if (e1%2)==1:
            t = (sq*t)%n
            e1 = (e1-1)//2
        else:
            e1 = e1//2
        sq = (sq*sq)%n
    return(t)

""" This multiplies two polynomials in GF(p^2)
that is it calculates c(x) = a(x)*b(x) <mod q(x)>
where a(x) = a[0]x + a[1] and b(x) = b[0]x + b[1] and
q(x) = x^2 + q[0]x + q[1]
"""
def mult2(a,b,q,n):
    t1 = (a[0]*b[1])%n
    t2 = (a[1]*b[0])%n
    t1 = (t1+t2)%n
    t2 = (a[1]*b[1])%n
    t3 = ((n-1)*q[0])%n
    t4 = ((n-1)*q[1])%n
    t5 = (a[0]*b[0])%n
    t3 = (t5*t3)%n
    t4 = (t5*t4)%n
    c = [(t1+t3)%n , (t2+t4)%n]
    return(c)

# print(add_mod(1234567890, 1, 6457357))
# print(mul(1234567890, 1))
# print(deg_mod(1234567890, 4, 14325))
# print(div(1234567890, 100))
# print(greater(1234567890, 1))
# print(less(1234567890, 1))
# print(sqrt(1234567890))