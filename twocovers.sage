"""
This code is a Sage port of some of the Magma code in "twocovers.m". Useful references:
[1] E.V. Flynn, D. Testa and R. van Luijk. ``Two-coverings of Jacobians of curves of genus two''.
    Proc. London Math. Soc. (3) 104 (2012), 387—429.
[2] E.V. Flynn, B. Poonen and E. Schaefer. ``Cycles of Quadratic Polynomials and Rational Points
    on a Genus 2 Curve''. Duke Math. J. 90 (1997), 435—463.
[3] Cassels, J. W. S., ``The Mordell-Weil group of curves of genus 2'', in: M. Artin, J. Tate (eds.),
    Arithmetic and Geometry I, Birkhauser, Boston, (1983), 27–60.
[4] E.V. Flynn. ``The Group Law on the Jacobian of a Curve of Genus 2''. J. reine angew. Math. 439 (1993), 45—69.

See http://magma.maths.usyd.edu.au/magma/handbook/text/1530 for details of how points on the Jacobian are implemented in Magma.

Table of contents:
   * Section 1: The Cassels map in genus 2: not implemented in Sage yet
   * Section 2: Desingularized twisted Kummer surfaces in P^5
   * Section 3: Maps on Kummer surfaces and their twists
   * Section 4: Embedding the genus 5 curve
   * Section 5: Searching for points on the curve
   * Appendix: Test code
"""

load("twocovers-data.sage") # contains some precomputed data
from sage.schemes.projective.projective_rational_point import sieve

"""
*************************************************************
  Section 2: Desingularized twisted Kummer surfaces in P^5
*************************************************************
"""

def kummer_matrices(f):
    """See Flynn-Testa-van Luijk [1], Remark 3.7. Must have deg(f) = 6."""
    assert f.degree() == 6
    f0, f1, f2, f3, f4, f5, f6 = list(f)
    R = matrix(f.base_ring(), 
        [[0, 0, 0, 0, 0, -f0/f6],
         [1, 0, 0, 0, 0, -f1/f6],
         [0, 1, 0, 0, 0, -f2/f6],
         [0, 0, 1, 0, 0, -f3/f6],
         [0, 0, 0, 1, 0, -f4/f6],
         [0, 0, 0, 0, 1, -f5/f6]])
    T = matrix(f.base_ring(), 6, 6,
        [[f1, f2, f3, f4, f5, f6],
         [f2, f3, f4, f5, f6, 0],
         [f3, f4, f5, f6, 0,  0],
         [f4, f5, f6, 0,  0,  0],
         [f5, f6, 0,  0,  0,  0],
         [f6, 0,  0,  0,  0,  0]])
    return R, T

def kummer_quadratic_form(f, delta, j):
    """Quadratic form Q_j^{(delta}) as defined in [1], section 4."""
    R, T = kummer_matrices(f)
    A.<t> = f.base_ring()[]
    L.<t> = A.quotient(f)
    try:
        d = list(L(A(list(delta))))
    except TypeError:
        d = [delta, 0, 0, 0, 0, 0]
    M = f.leading_coefficient() * sum([d[i] * R^(i+j) * T for i in range(6)])
    return QuadraticForm(M).polynomial()

def twisted_kummer(f, delta):
    """Construct the desingularized twisted Kummer surface in P^5"""
    Q0, Q1, Q2 = (kummer_quadratic_form(f, delta, j) for j in [0, 1, 2])
    P5.<v1, v2, v3, v4, v5, v6> = ProjectiveSpace(f.base_ring(), 5)
    return P5.subscheme([Q(*P5.gens()) for Q in [Q0, Q1, Q2]])

"""
*************************************************************
  Section 3: Maps on Kummer surfaces and their twists
*************************************************************
"""

# Source: https://people.maths.ox.ac.uk/flynn/genus2/kummer/duplication
def kummer_duplication(f):
    """Duplication map on (singular) Kummer surface in P^3, as computed by Flynn (cf. [4])"""
    assert f.degree() == 6
    f0, f1, f2, f3, f4, f5, f6 = list(f)
    K = KummerSurface(Jacobian(HyperellipticCurve(f)))
    P3 = K.ambient_space()
    k1, k2, k3, k4 = P3.gens()
    # For convenience, the polynomials are stored separately in "twocovers-data.sage"
    delta = [KUMMER_DUPLICATION_POLY[i](*list(f) + [k1, k2, k3, k4]) for i in range(4)]
    return K.Hom(K)(delta)

# Source: [1], Remark 3.8
def kummer_blowup_map(f):
    """The map from the desingularized Kummer in P^5 to the singular Kummer surface in P^3"""
    assert f.degree() == 6
    f0, f1, f2, f3, f4, f5, f6 = list(f)
    Y1 = twisted_kummer(f, 1)
    b1, b2, b3, b4, b5, b6 = Y1.ambient_space().gens()
    X1 = KummerSurface(Jacobian(HyperellipticCurve(f)))
    k1 = b1*b3 - b2^2
    k2 = b1*b4 - b2*b3
    k3 = b2*b4 - b3^2
    k4 = f0*b1^2 + f1*b1*b2 + f2*b2^2 + f3*b2*b3 + f4*b3^2 + f5*b3*b4 + f6*b4^2
    return Y1.Hom(X1)([k1, k2, k3, k4])

def kummer_smooth_duplication(f):
    """Combine the blowup and duplication maps"""
    return kummer_duplication(f) * kummer_blowup_map(f)

# Compute divisor on Jacobian given by applying the twisted duplication map to the
# point P on twisted_kummer(f, delta).
def twocover_eval(f, delta, P):
    A.<x> = (f.base_ring())[]
    L.<T> = A.quotient(f)
    delta = L(list(delta))
    g = [sum([list(f)[j+i] * T^j for j in range(6 - i + 1)]) for i in range(1, 7)]
    xi = sum([list(P)[i] * g[i] for i in range(6)])
    H = A(list(delta * xi^2))
    return H

"""
*************************************************************
  Section 4: Embedding the genus 5 curve
*************************************************************
"""

def base_kummer_locus(f, root=None):
    """Compute the image of the map from C: y^2 = f(x) to the (singular) Kummer surface of Jac(C) under
    the Abel-Jacobi embedding of the point (0, 0), assumed to lie on C."""
    assert f.degree() == 6
    if root is None:
        roots = f.roots()
        if len(roots) == 0:
            raise ValueError("Polynomial must have a root over the base field")
        else:
            root = roots[0][0]
    else:
        assert f(root) == 0
    f0, f1, f2, f3, f4, f5, f6 = list(f)
    K = KummerSurface(Jacobian(HyperellipticCurve(f)))
    P3 = K.ambient_space()
    k1, k2, k3, k4 = P3.gens()
    return P3.subscheme(K.defining_polynomials() + (root^2 * k1 - root * k2 + k3,)).reduce()

def genus5_in_P5(f, root=None, delta=1):
    """The genus 5 curve embedded in the desingularized Kummer surface Y1 in P^5"""
    assert f.degree() == 6
    if root is None:
        roots = f.roots()
        if len(roots) == 0:
            raise ValueError("Polynomial must have a root over the base field")
        else:
            root = roots[0][0]
    else:
        assert f(root) == 0
    g1, g2, g3, g4, g5, g6 = list(f // (f.parent().gens()[0] - root))
    Y = twisted_kummer(f, delta)
    P5 = Y.ambient_space()
    b1, b2, b3, b4, b5, b6 = P5.gens()
    Z = P5.subscheme(Y.defining_polynomials() + (g1*b1 + g2*b2 + g3*b3 + g4*b4 + g5*b5 + g6*b6,))
    return Z

def genus5_canonical(f, root, delta):
    """The genus 5 curve, canonically embedded in P^4 (as a hyperplane section of Y_delta in P^5)"""
    Z = genus5_in_P5(f, root, delta)
    k = f.base_ring()
    P4.<v1, v2, v3, v4, v5> = ProjectiveSpace(k, 4)
    lin_rel = [p for p in Z.defining_polynomials() if p.degree() == 1][0]
    # We choose to eliminate v6 because, by construction, its coefficient is nonzero.
    v6 = -lin_rel(v1, v2, v3, v4, v5, 0) / lin_rel(0, 0, 0, 0, 0, 1)
    v = [v1, v2, v3, v4, v5, v6]
    can = P4.Hom(Z.ambient_space())(v)
    Z_can = P4.subscheme([poly(v) for poly in Z.defining_polynomials() if poly.degree() > 1])
    return Z_can, can

def twisted_duplication_map(f, root, delta):
    """Twisted duplication map on the genus 5 canonical curve Z_delta in P^4 = P(ker(eval @ root)) subset P^5"""
    k = f.base_ring()
    Z, can = genus5_canonical(f, root, delta)
    v = can.defining_polynomials()
    P1.<x, z> = ProjectiveSpace(k, 1)
    Q3 = kummer_quadratic_form(f, delta, 3)(v)
    Q4 = kummer_quadratic_form(f, delta, 4)(v)
    f5 = list(f)[5]
    f6 = f.leading_coefficient()
    return Z.Hom(P1)([-(f5 + f6*root)*Q3 - f6*Q4, f6*Q3])

"""
*************************************************************
  Section 5: Searching for points on the curve
*************************************************************
"""
def poly_to_x_coord(H, root):
    A.<x> = (H.base_ring())[]
    factor = H(x) // (H(x).leading_coefficient() * (x - root))
    if factor == 1:
        return Infinity
    else:
        return -factor.constant_coefficient()

def point_search_on_twist(f, root, delta, bound):
    Z = genus5_in_P5(f, root, delta)
    twist_pts = sieve(Z, bound) # warning: this seems to be really slow
    x_coords = [poly_to_x_coord(twocover_eval(f, delta, P), root) for P in twist_pts]
    finite_pts1 = {(x, sqrt(f(x)), 1) for x in x_coords if x != Infinity}
    finite_pts2 = {(x, -sqrt(f(x)), 1) for x in x_coords if x != Infinity}
    # can't give coordinates for these because Sage doesn't embed hyperelliptic curves into the right weighted projective space
    infinite_pts = {Infinity, -Infinity}
    return finite_pts1 | finite_pts2 | infinite_pts

def twocover_point_search(C, w, bound, twists):
    f, h = C.hyperelliptic_polynomials()
    assert h == 0
    root = list(w)[0]
    pt_coords = set.union(*[point_search_on_twist(f, root, delta, bound) for delta in twists])
    finite_pts = {C(P) for P in pt_coords if P not in {Infinity, -Infinity}}
    infinite_pts = {P for P in pt_coords if P in {Infinity, -Infinity}}
    return finite_pts | infinite_pts

"""
*************************************************************
  Appendix: test code
*************************************************************
"""
A.<x> = QQ[]
# Genus 2 curve 533963.a.533963.1
f = x^6 + 2*x^5 + 7*x^4 + 4*x^3 - 13*x^2 - 10*x + 1
root = -1
L.<T> = A.quotient(f)
# twist parameters, computed in Magma
deltas = [1 * T^0,
    T^2 - 2,
    T^5 + T^4 + 6*T^3 - T^2 - 10*T + 1,
    -T^5 - T^4 - 7*T^3 + T^2 + 13*T + 1,
    -5*T^5 - 6*T^4 - 23*T^3 + 12*T^2 + 46*T - 4,
    4*T^5 + 4*T^4 + 23*T^3 - 10*T^2 - 45*T + 4,
    T^5 + T^4 + 6*T^3 - 2*T^2 - 12*T,
    -T^5 - 5*T^3 + 9*T - 1]

B.<alpha,g1,g2,g3,g4,g5,g6> = QQ[]
C.<X> = FractionField(B)[]
f_gen = (X - alpha)*(g1 + g2*X + g3*X^2 + g4*X^3 + g5*X^4 + g6*X^5)

