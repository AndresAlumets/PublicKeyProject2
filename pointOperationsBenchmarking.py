#Library where I define various operations of points on a Weiserstrass curve
import sys
import math
import random
import time

operations= 0

class WeierstrassCurve:
    """Represents an elliptic curve in Weierstrass form: y^2 = x^3 + ax + b"""
    
    def __init__(self, a, b, p=None):
        """
        Initialize a Weierstrass curve.
        
        Args:
            a: coefficient a in y^2 = x^3 + ax + b
            b: coefficient b in y^2 = x^3 + ax + b
            p: prime modulus (optional, for finite field arithmetic)
        """
        self.a = a
        self.b = b
        self.p = p
    
    def is_on_curve(self, x, y):
        """Check if point (x, y) lies on the curve"""
        left = (y ** 2) % self.p
        right = (x ** 3 + self.a * x + self.b) % self.p
        return left == right
    
    def writeOut(self):
        return f"WeierstrassCurve(a={self.a}, b={self.b}, p={self.p})"


class AffinePoint:
    """Represents a point in affine coordinates on a Weierstrass curve"""
    
    def __init__(self, x, y, curve):
        """
        Initialize an affine point.
        
        Args:
            x: x-coordinate
            y: y-coordinate
            curve: instance of WeierstrassCurve
        """
        self.x = x
        self.y = y
        self.curve = curve

    
    def writeOut(self):
        return f"AffinePoint(x={self.x}, y={self.y}) on {self.curve.writeOut()}"
    
class ProjectivePoint:
    """Represents a point in projective coordinates on a Weierstrass curve"""
    
    def __init__(self, X, Y, Z, curve):
        """
        Initialize a projective point.
        
        Args:
            X: X-coordinate
            Y: Y-coordinate
            Z: Z-coordinate
            curve: instance of WeierstrassCurve
        """
        self.X = X
        self.Y = Y
        self.Z = Z
        self.curve = curve

        
    def writeOut(self):
        return f"ProjectivePoint(X={self.X}, Y={self.Y}, Z={self.Z}) on {self.curve.writeOut()}"
    
class JacobianPoint:
    """Represents a point in Jacobian coordinates on a Weierstrass curve"""
    
    def __init__(self, X, Y, Z, curve):
        """
        Initialize a Jacobian point.
        
        Args:
            X: X-coordinate
            Y: Y-coordinate
            Z: Z-coordinate
            curve: instance of WeierstrassCurve
        """
        self.X = X
        self.Y = Y
        self.Z = Z
        self.curve = curve

    def writeOut(self):
        return f"JacobianPoint(X={self.X}, Y={self.Y}, Z={self.Z}) on {self.curve.writeOut()}"

def to_projective(P):
    """Convert an affine point P to projective coordinates."""
    if P is None:
        return ProjectivePoint(0, 1, 0, None)  # Point at infinity
    return ProjectivePoint(P.x, P.y, 1, P.curve)

def to_jacobian(P):
    """Convert an affine point P to Jacobian coordinates."""
    if P is None:
        return JacobianPoint(1, 1, 0, None)  # Point at infinity
    return JacobianPoint(P.x, P.y, 1, P.curve)

def to_affine(P, type):
    if type == "projective":
        if P.Z == 0:
            return None  # Point at infinity
        x = (P.X * pow(P.Z, -1, P.curve.p)) % P.curve.p
        y = (P.Y * pow(P.Z, -1, P.curve.p)) % P.curve.p
        return AffinePoint(x, y, P.curve)
    elif type == "jacobian":
        if P.Z == 0:
            return None  # Point at infinity
        z_inv = pow(P.Z, -1, P.curve.p)
        z_inv2 = (z_inv ** 2) % P.curve.p
        z_inv3 = (z_inv2 * z_inv) % P.curve.p
        x = (P.X * z_inv2) % P.curve.p
        y = (P.Y * z_inv3) % P.curve.p
        return AffinePoint(x, y, P.curve)
    
def mul(a, b, p):
    global operations
    operations += 1
    return (a * b) % p

def add(a, b, p):
    global operations
    operations += 1
    return (a + b) % p

def sub(a, b, p):
    global operations
    operations += 1
    return (a - b) % p

def square(a, p):
    global operations
    operations += 1
    return (a * a) % p

def cube(a, p):
    global operations
    operations += 2
    return (a * a * a) % p

def inverse(a, p):
    global operations
    operations += math.ceil(math.log2(p))
    return pow(a, -1, p)

def double_affine_point(P, p):
    if P == None:
        return None    #Point at infinity
    if P.y == 0:
        return None    #Point at infinity
    # Point doubling
    m = mul(3, square(P.x, p), p)
    m = add(m, P.curve.a, p)
    m = mul(m, pow(mul(2, P.y, p), -1, p), p)

    x_r = sub(square(m, p), mul(2, P.x, p), p)
    y_r = sub(mul(m, sub(P.x, x_r, p), p), P.y, p)
    return AffinePoint(x_r, y_r, P.curve)


def add_affine_points(P, Q, p):
    """Add two points P and Q in affine coordinates on a Weierstrass curve."""
    if P== None:        #Lets define None as the point at infinity
        return Q
    if Q== None:
        return P
    if P.x == Q.x and P.y != Q.y:
        return None    #Point at infinity
    if P.x == Q.x and P.y == Q.y:
        return double_affine_point(P, p)
    # Point addition
    m = mul(sub(Q.y, P.y, p), pow(sub(Q.x, P.x, p), -1, p), p)
    x_r = sub(square(m, p), add(P.x, Q.x, p), p)
    y_r = sub(mul(m, sub(P.x, x_r, p), p), P.y, p)
    
    return AffinePoint(x_r, y_r, P.curve)

def double_projective_point(P, p):
    if P.Z == 0:
        return P  # Point at infinity
    if (2 * P.Y) % p == 0:
        return ProjectivePoint(0, 1, 0, P.curve)  # Point at infinity
    # Point doubling in projective coordinates
    A = add(mul(P.curve.a, square(P.Z, p), p), mul(3, square(P.X, p), p), p)
    B = mul(P.Y, P.Z, p)
    C = mul(mul(P.X, P.Y, p), B, p)
    D = sub(square(A, p), mul(8, C, p), p)
    X_r = mul(2, B, p)
    X_r = mul(X_r, D, p)
    Y_r = mul(A, sub(mul(4, C, p), D, p), p)
    Y_r = sub(Y_r, mul(8, square(mul(P.Y, B, p), p), p), p)
    Z_r = mul(8, cube(B, p), p)

    return ProjectivePoint(X_r, Y_r, Z_r, P.curve)

def add_projective_points(P, Q, p):
    if P.Z == 0:
        return Q
    if Q.Z == 0:
        return P
    if equals(P, Q, "projective"):
        return double_projective_point(P, p)
    if (P.X * Q.Z) % p == (Q.X * P.Z) % p and (P.Y * Q.Z) % p != (Q.Y * P.Z) % p:
        return ProjectivePoint(0, 1, 0, P.curve)  # Point at infinity
    # Point addition in projective coordinates
    A = sub(mul(Q.Y, P.Z, p), mul(P.Y, Q.Z, p), p)
    B = sub(mul(Q.X, P.Z, p), mul(P.X, Q.Z, p), p)
    C = square(A, p)
    C = mul(C, P.Z, p)
    C = mul(C, Q.Z, p)
    C = sub(C, cube(B, p), p)
    C = sub(C, mul(2, mul(square(B, p), mul(P.X, Q.Z, p), p), p), p)
    X_r = mul(B, C, p)
    Y_r = mul(A, sub(mul(square(B, p), mul(P.X, Q.Z, p), p), C, p), p)
    Y_r = sub(Y_r, mul(cube(B, p), mul(P.Y, Q.Z, p), p), p)
    Z_r = mul(cube(B, p), mul(P.Z, Q.Z, p), p)
    return ProjectivePoint(X_r, Y_r, Z_r, P.curve)

def double_jacobian_point(P, p):
    if P.Z == 0:
        return P  # Point at infinity
    if (2 * P.Y) % p == 0:
        return JacobianPoint(1, 1, 0, P.curve)  # Point at infinity
    # Point doubling in Jacobian coordinates
    #Taken from https://eprint.iacr.org/2014/1014.pdf
    L = mul(3, square(P.X, p), p)
    L = add(L, mul(P.curve.a, mul(square(P.Z, p), square(P.Z, p), p), p), p)
    Z_r = mul(2, mul(P.Y, P.Z, p), p)
    X_r = sub(square(L, p), mul(8, mul(P.X, square(P.Y, p), p), p), p)
    Y_r = sub(mul(L, sub(mul(4, mul(P.X, square(P.Y, p), p), p), X_r, p), p), mul(8, square(square(P.Y, p), p), p), p)
    return JacobianPoint(X_r, Y_r, Z_r, P.curve)

def add_jacobian_points(P, Q, p):
    if P.Z == 0:
        return Q
    if Q.Z == 0:
        return P
    if equals(P, Q, "jacobian"):
        return double_jacobian_point(JacobianPoint(P.X, P.Y, P.Z, P.curve), p)
    if (P.X * Q.Z ** 2) % p == (Q.X * P.Z ** 2) % p and (P.Y * Q.Z ** 3) % p != (Q.Y * P.Z ** 3) % p:
        return JacobianPoint(1, 1, 0, P.curve)  # Point at infinity
    # Point addition in Jacobian coordinates
    Z1Z1 = mul(P.Z, P.Z, p)
    Z2Z2 = mul(Q.Z, Q.Z, p)
    U1 = mul(P.X, Z2Z2, p)
    U2 = mul(Q.X, Z1Z1, p)
    Z1Z1Z1 = mul(P.Z, Z1Z1, p)
    Z2Z2Z2 = mul(Q.Z, Z2Z2, p)
    S1 = mul(P.Y, Z2Z2Z2, p)
    S2 = mul(Q.Y, Z1Z1Z1, p)
    H = sub(U2, U1, p)
    r = sub(S2, S1, p)
    H2 = square(H, p)
    H3 = mul(H, H2, p)
    U1H2 = mul(U1, H2, p)
    X_r = sub(sub(square(r, p), H3, p), mul(2, U1H2, p), p)
    Y_r = sub(mul(r, sub(U1H2, X_r, p), p), mul(S1, H3, p), p)
    Z_r = mul(H, mul(P.Z, Q.Z, p), p)
    return JacobianPoint(X_r, Y_r, Z_r, P.curve)

def add_points(P, Q, type, p):
    if type == "affine":
        return add_affine_points(P, Q, p)
    elif type == "projective":
        return add_projective_points(P, Q, p)
    elif type == "jacobian":
        return add_jacobian_points(P, Q, p)
    
 #Double point function
 # Input type of P can be multiple types
 # Returns point in same system as input   
def double_point(P, type, p):
    if type == "affine":
        return double_affine_point(P, p)
    elif type == "projective":
        return double_projective_point(P, p)
    elif type == "jacobian":
        return double_jacobian_point(P, p)
    
def scalar_multiply(P, k, type, p):
    R = infinity(type, P.curve)
    Q = P
    while k > 0:
        if k % 2 == 1:
            R = add_points(R, Q, type, p)
        Q = double_point(Q, type, p)
        k = k // 2
    return R

def infinity(type, curve):
    if type == "affine":
        return None
    elif type == "projective":
        return ProjectivePoint(0, 1, 0, curve)
    elif type == "jacobian":
        return JacobianPoint(1, 1, 0, curve)
    
def isInfinity(P, type):
    if type == "affine":
        return P == None
    elif type == "projective":
        return P.Z == 0
    elif type == "jacobian":
        return P.Z == 0
    
def equals(P, Q, type):
    if type == "affine":
        if P == None and Q == None:
            return True
        if P == None or Q == None:
            return False
        return P.x == Q.x and P.y == Q.y
    elif type == "projective":
        if P.Z == 0 and Q.Z == 0:
            return True
        if P.Z == 0 or Q.Z == 0:
            return False
        return (P.X * Q.Z) % P.curve.p == (Q.X * P.Z) % P.curve.p and (P.Y * Q.Z) % P.curve.p == (Q.Y * P.Z) % P.curve.p
    elif type == "jacobian":
        if P.Z == 0 and Q.Z == 0:
            return True
        if P.Z == 0 or Q.Z == 0:
            return False
        return (P.X * Q.Z ** 2) % P.curve.p == (Q.X * P.Z ** 2) % P.curve.p and (P.Y * Q.Z ** 3) % P.curve.p == (Q.Y * P.Z ** 3) % P.curve.p
    


def init_curve(p):
    return WeierstrassCurve(p-3, 18958286285566608000408668544493926415504680968679321075787234672564, p=p)

#point_type = sys.argv[1]
p = int(sys.argv[2])

#I just want to create a curve and a point on it for testing
#Because those are not given as arguments
# So it is all tuned for a specific curve

curve = init_curve(p)
Gx = 19277929113566293071110308034699488026831934219452440156649784352033
Gy = 19926808758034470970197974370888749184205991990603949537637343198772

runs = 50000

#run_tests(Gx, Gy, curve, p, point_type)
for point_type in ["affine", "projective", "jacobian"]:
    operations = 0
    print("Running tests for", point_type)
    fixedPoint = AffinePoint(Gx, Gy, curve)
    if point_type == "projective":
        fixedPoint = to_projective(fixedPoint)
    elif point_type == "jacobian":
        fixedPoint = to_jacobian(fixedPoint)
    P1 = scalar_multiply(fixedPoint, 123456789, point_type, p)
    P2 = scalar_multiply(fixedPoint, 987654321, point_type, p)
    time_start = time.time()
    for i in range(runs):
        double_point(P1, point_type, p)
    time_end = time.time()
    print("Operations done in doubling:", operations, "for", runs, "doublings with point type", point_type)
    print("Time taken for", runs, "doublings with point type", point_type, ":", time_end - time_start, "seconds")
    operations = 0
    time_start = time.time()
    for i in range(runs):
        add_points(P1, P2, point_type, p)
    time_end = time.time()
    print("Operations done in adding:", operations, "for", runs, "additions with point type", point_type)
    print("Time taken for", runs, "additions with point type", point_type, ":", time_end - time_start, "seconds")