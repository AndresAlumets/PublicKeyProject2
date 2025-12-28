#Library where I define various operations of points on a Weiserstrass curve
import sys
import math
import random

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
        
        if not curve.is_on_curve(x, y):
            raise ValueError("The point is not on the given curve.")
    
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

        if Z == 0:
            return  # Point at infinity, no need to check
        left = (Y ** 2)*Z % curve.p
        right = (X ** 3 + curve.a * X * Z ** 2 + curve.b * Z ** 3) % curve.p
        if left != right:
            raise ValueError("The point is not on the given curve.")
        
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
        if Z == 0:
            return  # Point at infinity, no need to check
        left = (Y ** 2) % curve.p
        right = (X ** 3 + curve.a * X * Z ** 4 + curve.b * Z ** 6) % curve.p
        if left != right:
            raise ValueError("The point is not on the given curve.")

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

def double_affine_point(P, p):
    if P == None:
        return None    #Point at infinity
    if P.y == 0:
        return None    #Point at infinity
    # Point doubling
    m= (3 * P.x ** 2 + P.curve.a) * pow(2 * P.y, -1, p) % p
    x_r = (m ** 2 - 2 * P.x) % p
    y_r = (m * (P.x - x_r) - P.y) % p
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
    m = (Q.y - P.y) * pow(Q.x - P.x, -1, p) % p
    
    x_r = (m ** 2 - P.x - Q.x) % p
    y_r = (m * (P.x - x_r) - P.y) % p
    
    return AffinePoint(x_r, y_r, P.curve)

def double_projective_point(P, p):
    if P.Z == 0:
        return P  # Point at infinity
    if (2 * P.Y) % p == 0:
        return ProjectivePoint(0, 1, 0, P.curve)  # Point at infinity
    # Point doubling in projective coordinates
    A = (P.curve.a * P.Z ** 2 + 3* P.X ** 2) % p
    B = (P.Y * P.Z) % p
    C = (P.X * P.Y * B) % p
    D = (A ** 2 - 8 * C) % p
    X_r = (2 * B * D) % p
    Y_r = (A * (4 * C - D) - 8 * (B ** 2) * P.Y ** 2) % p
    Z_r = (8 * B ** 3) % p
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
    A = (Q.Y * P.Z - P.Y * Q.Z) % p
    B = (Q.X * P.Z - P.X * Q.Z) % p
    C = (A ** 2 * P.Z * Q.Z - B ** 3 - 2* B**2* P.X*Q.Z) % p
    X_r = (B * C) % p
    Y_r = (A * (B ** 2 * P.X * Q.Z - C) - B ** 3 * P.Y * Q.Z) % p
    Z_r = (B ** 3 * P.Z * Q.Z) % p
    return ProjectivePoint(X_r, Y_r, Z_r, P.curve)

def double_jacobian_point(P, p):
    if P.Z == 0:
        return P  # Point at infinity
    if (2 * P.Y) % p == 0:
        return JacobianPoint(1, 1, 0, P.curve)  # Point at infinity
    # Point doubling in Jacobian coordinates
    #Taken from https://eprint.iacr.org/2014/1014.pdf
    L = (3 * P.X ** 2 + P.curve.a * P.Z ** 4) % p
    Z_r = (2 * P.Y * P.Z) % p
    X_r = (L**2 - 8 * P.X * P.Y ** 2) % p
    Y_r = (L * (4 * P.X * P.Y ** 2 - X_r) - 8 * P.Y ** 4) % p
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
    Z1Z1 = (P.Z ** 2) % p
    Z2Z2 = (Q.Z ** 2) % p
    U1 = (P.X * Z2Z2) % p
    U2 = (Q.X * Z1Z1) % p
    Z1Z1Z1 = (P.Z * Z1Z1) % p
    Z2Z2Z2 = (Q.Z * Z2Z2) % p
    S1 = (P.Y * Z2Z2Z2) % p
    S2 = (Q.Y * Z1Z1Z1) % p
    H = (U2 - U1) % p
    r = (S2 - S1) % p
    H2 = (H ** 2) % p
    H3 = (H * H2) % p
    U1H2 = (U1 * H2) % p
    X_r = (r ** 2 - H3 - 2 * U1H2) % p
    Y_r = (r * (U1H2 - X_r) - S1 * H3) % p
    Z_r = (H * P.Z * Q.Z) % p
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
    
def test1(fixedPoint, p, point_type):
    #Simple addition test
    P = fixedPoint
    Q = fixedPoint
    R = add_points(P, Q, point_type, p)
    expected_R = double_point(P, point_type, p)
    if equals(R, expected_R, point_type):
        print("Test1 Doubling : passes")
    else:
        print("Test1 Doubling : fails")

def test2(fixedPoint, p, point_type):
    P = fixedPoint
    Q = infinity(point_type, fixedPoint.curve)
    R = add_points(P, Q, point_type, p)
    if equals(R, P, point_type):
        print("Test2 Addition with infinity : passes")
    else:
        print("Test2 Addition with infinity : fails")

def test3(fixedPoint, p, point_type):
    P = infinity(point_type, fixedPoint.curve)
    Q = fixedPoint
    R = add_points(P, Q, point_type, p)
    if equals(R, Q, point_type):
        print("Test3 Addition with infinity : passes")
    else:
        print("Test3 Addition with infinity : fails")

def test4(fixedPoint, p, point_type):
    P = fixedPoint
    if point_type == "affine":
        Q = AffinePoint(P.x, (-P.y) % p, P.curve)
    elif point_type == "projective":
        Q = ProjectivePoint(P.X, (-P.Y) % p, P.Z, P.curve)
    elif point_type == "jacobian":
        Q = JacobianPoint(P.X, (-P.Y) % p, P.Z, P.curve)
    R = add_points(P, Q, point_type, p)
    if isInfinity(R, point_type):
        print("Test4 Addition with negative point : passes")
    else:
        print("Test4 Addition with negative point : fails")

def test5(fixedPoint, p, point_type):
    k1 = 4
    k2 = 5
    k_total = 10
    P = fixedPoint
    R1 = scalar_multiply(P, k1, point_type, p)
    R2 = scalar_multiply(P, k_total - k1, point_type, p)
    R_res1 = add_points(R1, R2, point_type, p)
    R3 = scalar_multiply(P, k2, point_type, p)
    R4 = scalar_multiply(P, k_total - k2, point_type, p)
    R_res2 = add_points(R3, R4, point_type, p)
    if equals(R_res1, R_res2, point_type):
        print("Test5 Distributivity : passes")
    else:
        print("Test5 Distributivity : fails")

def test6(fixedPoint, p, point_type):
    k1 = 4
    k2 = 5
    P = fixedPoint
    R1 = scalar_multiply(P, k1, point_type, p)
    R2 = scalar_multiply(P, k2, point_type, p)
    R_res1 = add_points(R1, R2, point_type, p)
    R_res2 = add_points(R2, R1, point_type, p)
    if equals(R_res1, R_res2, point_type):
        print("Test6 Commutativity : passes")
    else:
        print("Test6 Commutativity : fails")

def test7(fixedPoint, p, point_type):
    k1 = 4
    k2 = 5
    P = fixedPoint
    R1 = scalar_multiply(P, k1, point_type, p)
    R2 = scalar_multiply(P, k2, point_type, p)
    R3 = scalar_multiply(P, k1 + k2, point_type, p)
    R_res = add_points(R1, R2, point_type, p)
    if equals(R_res, R3, point_type):
        print("Test7 Addition of Scalars : passes")
    else:
        print("Test7 Addition of Scalars : fails")

def test8(fixedPoint, p, point_type):
    k1 = random.randint(1, p-1)
    k2 = random.randint(1, p-1)
    k3 = random.randint(1, p-1)
    R1 = scalar_multiply(fixedPoint, k1, point_type, p)
    R2 = scalar_multiply(R1, k2, point_type, p)
    R3 = scalar_multiply(fixedPoint, k3, point_type, p)
    R_res1 = add_points(R1, add_points(R2, R3, point_type, p), point_type, p)
    R_res2 = add_points(add_points(R1, R2, point_type, p), R3, point_type, p)
    if equals(R_res1, R_res2, point_type):
        print("Test8 Associativity : passes")
    else:
        print("Test8 Associativity : fails")


def run_tests(Gx, Gy, curve, p, point_type):
    fixedPoint = AffinePoint(Gx, Gy, curve)
    if point_type == "affine":
        print("Running tests for Affine Points")
    elif point_type == "projective":
        print("Running tests for Projective Points")
        fixedPoint = to_projective(fixedPoint)
    elif point_type == "jacobian":
        print("Running tests for Jacobian Points")
        fixedPoint = to_jacobian(fixedPoint)
    test1(fixedPoint, p, point_type)
    test2(fixedPoint, p, point_type)
    test3(fixedPoint, p, point_type)
    test4(fixedPoint, p, point_type)
    test5(fixedPoint, p, point_type)
    test6(fixedPoint, p, point_type)
    test7(fixedPoint, p, point_type)
    test8(fixedPoint, p, point_type)

def init_curve(p):
    return WeierstrassCurve(p-3, 18958286285566608000408668544493926415504680968679321075787234672564, p=p)

point_type = sys.argv[1]
p = int(sys.argv[2])

#I just want to create a curve and a point on it for testing
#Because those are not given as arguments
# So it is all tuned for a specific curve

curve = init_curve(p)
print(curve.writeOut())
Gx = 19277929113566293071110308034699488026831934219452440156649784352033
Gy = 19926808758034470970197974370888749184205991990603949537637343198772

run_tests(Gx, Gy, curve, p, point_type)
