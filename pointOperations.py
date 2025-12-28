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
    L = (P.Y - Q.Y) % p
    Z_r = ((P.X - Q.X) * P.Z) % p
    X_r = (L ** 2 -  (P.X + Q.X) * (P.X - Q.X) ** 2 ) % p
    Y_r = (L* (P.X * (P.X - Q.X) ** 2 - X_r) - P.Y * (P.X - Q.X) ** 3) % p
    return JacobianPoint(X_r, Y_r, Z_r, P.curve)

def add_points(P, Q, type, p):
    if type == "affine":
        return add_affine_points(P, Q, p)
    elif type == "projective":
        P_proj = to_projective(P)
        Q_proj = to_projective(Q)
        R_proj = add_projective_points(P_proj, Q_proj, p)
        return to_affine(R_proj, type)
    elif type == "jacobian":
        P_jac = to_jacobian(P)
        Q_jac = to_jacobian(Q)
        R_jac = add_jacobian_points(P_jac, Q_jac, p)
        return to_affine(R_jac, type)
    
def double_point(P, type, p):
    if type == "affine":
        return double_affine_point(P, p)
    elif type == "projective":
        P_proj = to_projective(P)
        R_proj = double_projective_point(P_proj, p)
        return to_affine(R_proj, type)
    elif type == "jacobian":
        P_jac = to_jacobian(P)
        R_jac = double_jacobian_point(P_jac, p)
        return to_affine(R_jac, type)
    
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
    
def test1(Gx, Gy, curve, p, point_type):
    #Simple addition test
    P = AffinePoint(Gx, Gy, curve)
    Q = AffinePoint(Gx, Gy, curve)
    R = add_points(P, Q, point_type, p)
    expected_R = double_point(P, point_type, p)
    if equals(R, expected_R, "affine"):
        print("Test1 Doubling : passes")
    else:
        print("Test1 Doubling : fails")

def test2(Gx, Gy, curve, p, point_type):
    P = AffinePoint(Gx, Gy, curve)
    Q = None
    R = add_points(P, Q, point_type, p)
    if equals(R, P, "affine"):
        print("Test2 Addition with infinity : passes")
    else:
        print("Test2 Addition with infinity : fails")

def test3(Gx, Gy, curve, p, point_type):
    P = None
    Q = AffinePoint(Gx, Gy, curve)
    R = add_points(P, Q, point_type, p)
    if equals(R, Q, "affine"):
        print("Test3 Addition with infinity : passes")
    else:
        print("Test3 Addition with infinity : fails")

def test4(Gx, Gy, curve, p, point_type):
    P = AffinePoint(Gx, Gy, curve)
    Q = AffinePoint(Gx, (-Gy) % p, curve)
    R = add_points(P, Q, point_type, p)
    if R == None:
        print("Test4 Addition with negative point : passes")
    else:
        print("Test4 Addition with negative point : fails")

def test5(Gx, Gy, curve, p, point_type):
    k1 = 4
    k2 = 5
    k_total = 10
    P = AffinePoint(Gx, Gy, curve)
    R1 = P
    for _ in range(k1 - 1): #Would be faster to do double and add but this is just for testing
        R1 = add_points(R1, P, point_type, p)
    R2 = P
    for _ in range(k_total - k1 -1):
        R2 = add_points(R2, P, point_type, p)
    R_res1 = add_points(R1, R2, point_type, p)
    R3 = P
    for _ in range(k2 - 1):
        R3 = add_points(R3, P, point_type, p)
    R4 = P
    for _ in range(k_total - k2 -1):
        R4 = add_points(R4, P, point_type, p)
    R_res2 = add_points(R3, R4, point_type, p)
    if equals(R_res1, R_res2, "affine"):
        print("Test5 Distributivity : passes")
    else:
        print("Test5 Distributivity : fails")

def test6(Gx, Gy, curve, p, point_type):
    k1 = 4
    k2 = 5
    P = AffinePoint(Gx, Gy, curve)
    R1 = P
    for _ in range(k1 - 1):
        R1 = add_points(R1, P, point_type, p)
    R2 = P
    for _ in range(k2 -1):
        R2 = add_points(R2, P, point_type, p)
    R_res1 = add_points(R1, R2, point_type, p)
    R_res2 = add_points(R2, R1, point_type, p)
    if equals(R_res1, R_res2, "affine"):
        print("Test6 Commutativity : passes")
    else:
        print("Test6 Commutativity : fails")

def run_tests(Gx, Gy, curve, p, point_type):
    test1(Gx, Gy, curve, p, point_type)
    test2(Gx, Gy, curve, p, point_type)
    test3(Gx, Gy, curve, p, point_type)
    test4(Gx, Gy, curve, p, point_type)
    test5(Gx, Gy, curve, p, point_type)
    test6(Gx, Gy, curve, p, point_type)

point_type = sys.argv[1]
p = int(sys.argv[2])

#I just want to create a curve and a point on it for testing
#Because those are not given as arguments
# So it is all tuned for a specific curve

a= p-3
b= 18958286285566608000408668544493926415504680968679321075787234672564
while (4 * a**3 + 27 * b**2) % p == 0:
    b+=1

curve = WeierstrassCurve(a, b, p=p)
print(curve.writeOut())
Gx = 19277929113566293071110308034699488026831934219452440156649784352033
Gy = 19926808758034470970197974370888749184205991990603949537637343198772
left = (Gy ** 2) % p
right = (Gx ** 3 + a * Gx + b) % p


run_tests(Gx, Gy, curve, p, point_type)
