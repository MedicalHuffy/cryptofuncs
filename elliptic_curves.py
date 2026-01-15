class EllipticCurve(object):
    # Weierstrass Elliptic Curve
    def __init__(self, a:int, b:int, p:int):
        # y^2 = x^3 + ax + b mod p
        self.a = a
        self.b = b
        self.p = p
        self.inf = [0,float('inf')]
        
    def y(self, x:int):
        ysquared = x**3 + self.a*x + self.b
        return modsqrt(ysquared)

    def add(self, P:list, Q:list):
        if P == self.inf:
            return Q
        if Q == self.inf:
            return P
        Px, Py = P
        Qx, Qy = Q
        if Px == Qx and Py == -Qy:
            return self.inf
        if P == Q:
            l = (3*Px**2 + self.a) * modinv(2*Py, self.p)
        else:
            l = (Qy - Py) * modinv(Qx - Px, self.p)
        x = (l**2 - Px - Qx) % self.p
        y = (l*(Px - x) - Py) % self.p
        return (x,y)
        
    def mult(self, n:int, P:list):
        Q = P
        R = self.inf
        while n > 0:
            if n % 2:
                R = self.add(R, Q)
            Q = self.add(Q,Q)
            n //= 2
        return R
