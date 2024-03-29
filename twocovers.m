/*
This code makes explicit the construction of two-coverings with rational points in [1]. Useful references:
[1] E.V. Flynn, D. Testa and R. van Luijk. ``Two-coverings of Jacobians of curves of genus two''.
    Proc. London Math. Soc. (3) 104 (2012), 387—429.
[2] E.V. Flynn, B. Poonen and E. Schaefer. ``Cycles of Quadratic Polynomials and Rational Points
    on a Genus 2 Curve''. Duke Math. J. 90 (1997), 435—463.
[3] Cassels, J. W. S., ``The Mordell-Weil group of curves of genus 2'', in: M. Artin, J. Tate (eds.),
    Arithmetic and Geometry I, Birkhauser, Boston, (1983), 27–60.
[4] E.V. Flynn. ``The Group Law on the Jacobian of a Curve of Genus 2''. J. reine angew. Math. 439 (1993), 45—69.

See http://magma.maths.usyd.edu.au/magma/handbook/text/1530 for details of how points on the Jacobian are implemented in Magma.

Table of contents:
   * Section 1: The Cassels map in genus 2
   * Section 2: Desingularized twisted Kummer surfaces in P^5
   * Section 3: Embedding the genus 5 curve
   * Section 4: Genus one quotients
   * Section 5: Searching for points on the curve
   * Section 6: Applying elliptic Chabauty to twists
   * Section 7: Numerical computation of degree of isogeny
*/

// Patch for a couple bugs in elliptic Chabauty.
// Fix due to Michael Stoll; should be included in future version of Magma.
AttachSpec("patch/patch.spec");

/*************************************************************
  Section 1: The Cassels map in genus 2
 *************************************************************/

// Polynomial computations to compute Cassels map
function cassels_map_poly(a, f)
    if Degree(a) eq 0 then
        return 1;
    elif GCD(a, f) eq 1 then
        return ((-1)^(Degree(a)) * a) mod f;
    elif Degree(a) eq 1 then
        F0 := f div a;
        // Solve g = -a (mod F0), g = F0 (mod a)
        g := CRT([-a, F0], [F0, a]);
        return g mod f;
    elif Degree(a) eq 2 then
        if IsIrreducible(a) then
            F0 := f div a;
            // Solve g = a (mod F0), g = -F0 (mod a)
            g := CRT([a, -F0], [F0, a]);
            return g mod f;
        else
            factors := Factorization(a);
            if IsSquare(a) then
                return 1;
            else
                // Call recursively on the factors
                return $$(factors[1][1], f) * $$(factors[2][1], f);
            end if;
        end if;
    end if;
end function;

// Let C: y^2 = f(x) be a genus 2 curve over a number field k, let J = Jac(C), and let L = k[T]/(f).
// Input: a point P in J(k).
// Output: a polynomial g such that the class of g(T) in L*/(L*)^2 k* (or L*/(L*)^2 if deg(f) = 5) is
// the value of the Cassels map at P.
function cassels_map(P)
    if Type(P) ne JacHypPt then
        error "Input to cassels_map must be point of Jacobian";
    end if;
    a := P[1];
    J := Parent(P);
    f, h := HyperellipticPolynomials(Curve(J));
    if h ne 0 then
        error "cassels_map only implemented for genus 2 curves in form y^2 = f(x)";
    end if;
    assert Dimension(J) eq 2;
    assert P[3] le 2;
    assert Degree(a) le P[3];
    return cassels_map_poly(a, f);
end function;

// Given an element delta of L = K[X]/(f), simplify f as much as possible mod squares
// Currently only simplifies the linear factors, and only over Q
function reduce_twist_param(f, delta)
    assert IsSquarefree(f);
    assert BaseRing(f/1) eq Rationals();
    R := Parent(f/1);
    if IsCoercible(Rationals(), delta) then
        delta := R!delta;
    else
        delta := R!Eltseq(delta);
    end if;
    fact := Factorization(f);
    factors := [p[1] : p in fact];
    ds := [];
    for i in [1 .. #factors] do
        if Degree(factors[i]) eq 1 then
            d := Evaluate(delta, -Roots(factors[i])[1][1]);
            num_sqfree := SquarefreeFactorization(Numerator(d));
            denom_sqfree := SquarefreeFactorization(Denominator(d));
            ds[i] := R!(num_sqfree / denom_sqfree);
        else
            ds[i] := delta mod factors[i];
        end if;
    end for;
    // Solve poly = ds[i] (mod factors[i])
    poly := CRT(ds, factors);
    return poly mod f;
end function;

// Given degree 5 or 6 polynomial f over Q, output the 2-cover twist parameters associated to C: y^2 = f(x).
// More precisely: let L = k[X]/(f), and let J be the Jacobian of C. We compute generators in L* for the image
// of the Cassels map from J(k)/2J(k) to L*/(L*)^2 k* (or L*/(L*)^2 if deg(f) = 5).
function twist_param_generators(f)
    k := CoefficientRing(Parent(f));
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    G, i := MordellWeilGroup(J); // this is only currently implemented over Q
    // Get the free or even-order generators of J(k)
    gen_list := [i(g) : g in Generators(G) | IsEven(Order(g))];

    R<t> := PolynomialRing(k);
    L<X> := quo< R | f>;
    twist_params := {Evaluate(cassels_map(P), X) : P in gen_list};
    return twist_params;
end function;

// Multiplies out all subsets of a given set
function products_of_subsets(set)
    type := Parent(Rep(set));
    return Include({&*S : S in Subsets(set) | not IsEmpty(S)}, type!1);
end function;

/*************************************************************
  Section 2: Desingularized twisted Kummer surfaces in P^5
 *************************************************************/

// See reference [1], Remark 3.7. Must have deg(f) = 6.
function kummer_matrices(f)
    assert Degree(f) eq 6;
    f0, f1, f2, f3, f4, f5, f6 := Explode(Coefficients(f));
    R := Matrix(CoefficientRing(f), 6, 6,
        [0, 0, 0, 0, 0, -f0/f6,
         1, 0, 0, 0, 0, -f1/f6,
         0, 1, 0, 0, 0, -f2/f6,
         0, 0, 1, 0, 0, -f3/f6,
         0, 0, 0, 1, 0, -f4/f6,
         0, 0, 0, 0, 1, -f5/f6]);
    T := Matrix(CoefficientRing(f), 6, 6,
        [f1, f2, f3, f4, f5, f6,
         f2, f3, f4, f5, f6, 0,
         f3, f4, f5, f6, 0,  0,
         f4, f5, f6, 0,  0,  0,
         f5, f6, 0,  0,  0,  0,
         f6, 0,  0,  0,  0,  0]);
    return R, T;
end function;

// Quadratic form Q_j^{(\delta}) as defined in [1], section 4.
function kummer_quadratic_form(f, delta, j)
    R, T := kummer_matrices(f);
    A := PolynomialRing(CoefficientRing(f));
    L<t> := quo< A | Evaluate(f, A.1) >;
    d := (Coefficients(L!Eltseq(delta)) cat [0 : i in [1 .. 6]])[1 .. 6];
    M := LeadingCoefficient(f) * &+[d[i+1] * R^(i + j) * T : i in [0 .. 5]];
    return QuadraticForm(M);
end function;

// Construct the desingularized twisted Kummer surface in P^5
function twisted_kummer(f, delta)
    Q0 := kummer_quadratic_form(f, delta, 0);
    Q1 := kummer_quadratic_form(f, delta, 1);
    Q2 := kummer_quadratic_form(f, delta, 2);
    P5<v1, v2, v3, v4, v5, v6> := ProjectiveSpace(CoefficientRing(f), 5);
    return Scheme(P5, [Evaluate(Q, [P5.i : i in [1 .. 6]]) : Q in [* Q0, Q1, Q2 *]]);
end function;

/*************************************************************
  Section 3: Embedding the genus 5 curve
 *************************************************************/

// The genus 5 curve embedded in the desingularized Kummer surface Y_delta in P^5
function genus5_in_P5(f, root, delta)
    assert Degree(f) eq 6 and Evaluate(f, root) eq 0;
    g1, g2, g3, g4, g5, g6 := Explode(Coefficients(f div (Parent(f).1 - root)));
    Y<v1, v2, v3, v4, v5, v6> := twisted_kummer(f, delta);
    Z := Curve(Y, [g1*v1 + g2*v2 + g3*v3 + g4*v4 + g5*v5 + g6*v6]);
    return Z;
end function;

// The genus 5 curve, canonically embedded in P^4 (as a hyperplane section of Y_delta in P^5)
function genus5_canonical(f, root, delta)
    Z := genus5_in_P5(f, root, delta);
    k := BaseRing(f);
    P4<v1, v2, v3, v4, v5> := ProjectiveSpace(k, 4);
    lin_rel := [p : p in DefiningEquations(Z) | Degree(p) eq 1][1];
    // We choose to eliminate v6 because, by construction, its coefficient is nonzero.
    v6 := -Evaluate(lin_rel, [v1, v2, v3, v4, v5, 0]) / k!Coefficient(lin_rel, Z.6, 1);
    v := [v1, v2, v3, v4, v5, v6];
    can := map<P4 -> Ambient(Z) | v>;
    Z_can := Curve(P4, [Evaluate(p, v) : p in DefiningEquations(Z) | Degree(p) ne 1]);
    return Z_can, can;
end function;

// Twisted duplication map on the genus 5 canonical curve Z_delta in P^4 = P(ker(eval @ root)) \subset P^5
function twisted_duplication_map(f, root, delta)
    k := BaseRing(f);
    Z<v1, v2, v3, v4, v5>, can := genus5_canonical(f, root, delta);
    v := DefiningPolynomials(can);
    Q3 := Evaluate(kummer_quadratic_form(f, delta, 3), v);
    Q4 := Evaluate(kummer_quadratic_form(f, delta, 4), v);
    f5 := Coefficient(f, 5);
    f6 := Coefficient(f, 6);
    P1<x, z> := ProjectiveSpace(k, 1);
    return map<Z -> P1 | [-(f5 + f6*root)*Q3 - f6*Q4, f6*Q3]>;
end function;

/*************************************************************
  Section 4: Genus one quotients
 *************************************************************/

// Computes delta(b_1)*delta(b_2)*delta(b_3)*delta(b_4) where b_1, ..., b_4 are the roots of the quartic h
function elliptic_twist_param(h, delta)
    assert Degree(h) eq 4;
    K := BaseRing(h);
    d := Parent(h)!Eltseq(delta) mod h;
    A<b1, b2, b3, b4> := PolynomialRing(K, 4);
    poly := Evaluate(d, b1) * Evaluate(d, b2) * Evaluate(d, b3) * Evaluate(d, b4);
    _, sym := IsSymmetric(poly);
    h0, h1, h2, h3, h4 := Explode(Eltseq(h));
    return Evaluate(sym, [-h3/h4, h2/h4, -h1/h4, h0/h4]);
end function;

// The genus 1 hyperelliptic curve that the twisted duplication map factors through
function genus1_quotient(h, root, delta)
    twist := elliptic_twist_param(h, delta);
    return HyperellipticCurve(Evaluate(h, root)/twist * h);
end function;

// The quadric form Y on L = K[X]/(f) defined for each s in L by:
// Y(s) = s(b1)*s(b2)*s(b3)*s(b4), where {b1,b2,b3,b4} are the roots of h.
function genus1_quadric(h, root1, root2)
    assert Degree(h) eq 4;
    K := BaseRing(h);
    x := Parent(h).1;
    f := h * (x - root1) * (x - root2);
    _, T := kummer_matrices(f);
    g := [Parent(h)!Eltseq(row) : row in Rows(T)];
    R<[b]> := PolynomialRing(K, 4);
    S<[vv]> := PolynomialRing(R, 6);
    ev_prod := &*[&+[Evaluate(g[i], b[j])*vv[i] : i in [1 .. 6]] : j in [1 .. 4]];
    coeffs_poly, monomials := CoefficientsAndMonomials(ev_prod);
    h0, h1, h2, h3, h4 := Explode(Eltseq(h));
    inputs := [-h3/h4, h2/h4, -h1/h4, h0/h4];
    coeffs := [];
    for poly in coeffs_poly do
        bool, sym := IsSymmetric(poly);
        assert bool;
        Append(~coeffs, Evaluate(sym, inputs));
    end for;
    B := PolynomialRing(K, 6);
    return &+[coeffs[i] * Monomial(B, Exponents(monomials[i])) : i in [1 .. #monomials]];
end function;

// The map from Z_delta to the genus 1 curve factoring the twisted duplication map
function genus1_map(f, root1, root2, delta)
    assert Degree(f) eq 6;
    assert Evaluate(f, root1) eq 0 and Evaluate(f, root2) eq 0;

    K := Parent(root2/1);
    A := ChangeRing(Parent(f), K);
    h := A!f div ((A.1 - root1)*(A.1 - root2));
    D := genus1_quotient(h, root1, delta);

    Z, can := genus5_canonical(f, root1, delta);
    ZK<v1,v2,v3,v4,v5> := ChangeRing(Z, K);
    v := [CoordinateRing(Ambient(ZK))!poly : poly in DefiningPolynomials(can)];
    Q3 := Evaluate(kummer_quadratic_form(f, delta, 3), v);
    Q4 := Evaluate(kummer_quadratic_form(f, delta, 4), v);
    f5 := Coefficient(f, 5);
    f6 := Coefficient(f, 6);
    xx := -(f5 + f6*root1)*Q3 - f6*Q4;
    yy := f6^3 * Evaluate(genus1_quadric(h, root1, root2), v);
    zz := f6*Q3;
    return map<ZK -> D | [xx, yy, zz]>;
end function;

/*************************************************************
  Section 5: Searching for points on the curve
 *************************************************************/

function lift_base_point(f, P)
    if P[2] eq 0 then
        square, root := IsSquare(LeadingCoefficient(f));
        if square then
            return {[1, root, 0], [1, -root, 0]};
        else
            return {};
        end if;
    else
        x := P[1]/P[2];
        square, root := IsSquare(Evaluate(f, x));
        if square then
            return {[x, root, 1], [x, -root, 1]};
        else
            return {};
        end if;
    end if;
end function;

function point_search_on_twist(f, root, delta : Bound := 10000)
    phi := twisted_duplication_map(f, root, delta);
    Z := Domain(phi);
    P1 := Codomain(phi);
    twist_pts := PointSearch(Z, Bound);
    base_pts := [phi(P) : P in twist_pts];
    return &join[lift_base_point(f, P) : P in base_pts];
end function;

function even_weierstrass_model(C)
    C_simp, phi := SimplifiedModel(C);
    f := HyperellipticPolynomials(C_simp);
    if IsEven(Degree(f)) then
        return C_simp, phi;
    else
        i := 0;
        while Evaluate(f, i) eq 0 do
            i := i + 1;
        end while;
        x := Parent(f).1;
        f_good := Numerator(x^6 * Evaluate(f, 1/x + i));
        C_good := HyperellipticCurve(f_good);
        coord_change := map< C_simp -> C_good | [C_simp.3, C_simp.2, C_simp.1 - i*C_simp.3] >;
        assert IsIsomorphism(coord_change);
        phi_good := Expand(phi * coord_change);
        return C_good, phi_good;
    end if;
end function;

function twocover_point_search(C, w : Bound := 10000)
    C_good, phi := even_weierstrass_model(C);
    f := HyperellipticPolynomials(C_good);
    root := Eltseq(phi(w))[1];
    twists := products_of_subsets(twist_param_generators(f));
    pt_coords := &join[point_search_on_twist(f, root, delta : Bound := Bound) : delta in twists];
    return {Inverse(phi)(C_good!P) : P in pt_coords | P in C_good};
end function;

/*************************************************************
  Section 6: Applying elliptic Chabauty to twists
 *************************************************************/

// Determine if univariate polynomial has real root. Errs on the side of false positives.
// TODO: make this numerically stable
function has_real_root(f)
    has_root, root := HasRoot(f);
    if has_root then
        return has_root, root;
    end if;
    complex_roots, errors := RootsNonExact(f);
    for i in [1 .. #complex_roots] do
        z := complex_roots[i];
        err := Abs(errors[i]);
        if Abs(Im(z)) le err then
            // Complex root can't be bounded away from the real line, assume it's real
            return true, Re(z);
        end if;
    end for;
    return false, _;
end function;

// Test if a plane curve has real points
function curve_has_real_points(C : precision := 200)
    RR := RealField(precision);
    F := DefiningPolynomial(C);
    // Check point at infinity
    if Evaluate(F, [1, 0, 0]) eq 0 then
        return true, [1, 0, 0];
    end if;
    // Check line at infinity
    A<T> := PolynomialRing(RR);
    f_inf := Evaluate(F, [T, 1, 0]);
    has_root, root := has_real_root(f_inf);
    if has_root then
        return true, [root, 1, 0];
    end if;
    // Check affine plane
    A<x, y> := PolynomialRing(Rationals(), 2);
    f := Evaluate(F, [x, y, 1]);
    D := Evaluate(Discriminant(f, y), [T, 0]);
    D_roots := [root[1] : root in Roots(D)];
    Sort(~D_roots);
    r := #D_roots;
    if r eq 0 then
        test_x_coords := [0];
    else
        test_x_coords := [D_roots[1] - 1] cat [(D_roots[i] + D_roots[i+1])/2 : i in [1 .. r-1]] cat [D_roots[r] + 1];
    end if;
    for b in test_x_coords do
        has_root, root := has_real_root(Evaluate(f, [b, T]));
        if has_root then
            return true, [b, root, 1];
        end if;
    end for;
    return false, _;
end function;

// Test if the twist is locally solvable. If not, return a prime that fails.
function is_genus5_locally_solvable(Z, bad_primes : Skip := [], CheckReal := true)
    // Need to check primes for which the Weil lower bound p+1-10*sqrt(p) is non-positive
    prime_threshold := 98;
    small_primes := PrimesUpTo(prime_threshold);
    // Also need to check primes of bad reduction
    large_bad_primes := [p : p in bad_primes | p gt prime_threshold];
    // Create a singular plane model for Z by projecting away from a point twice
    Z_planar := Curve(Projection(Projection(Z)));
    if Genus(Z_planar) ne 5 then
        // in the rare event we get unlucky with the first projection
        Z_space := Projection(Z, Ambient(Z)![0, 1, 0, 0, 0]);
        Z_planar := Curve(Projection(Z_space, Ambient(Z_space)![0, 1, 1, 0]));
        // check just in case we get unlucky twice
        assert Genus(Z_planar) eq 5;
    end if;
    for p in small_primes cat large_bad_primes do
        if p in Skip then
            continue;
        end if;
        if not IsLocallySolvable(Z_planar, p : Smooth, AssumeIrreducible) then
            return false, p;
        end if;
    end for;
    if CheckReal then
        if curve_has_real_points(Z_planar) then
            return true;
        else
            return false, Infinity();
        end if;
    else
        return true;
    end if;
end function;

// Test if the twist is locally solvable at 2.
function is_genus5_locally_solvable_at_2(Z)
    // Create a singular plane model for Z by projecting away from a point twice
    Z_planar := Curve(Projection(Projection(Z)));
    if Genus(Z_planar) ne 5 then
        // in the rare event we get unlucky with the first projection
        Z_space := Projection(Z, Ambient(Z)![0, 1, 0, 0, 0]);
        Z_planar := Curve(Projection(Z_space, Ambient(Z_space)![0, 1, 1, 0]));
        // check just in case we get unlucky twice
        assert Genus(Z_planar) eq 5;
    end if;
    return IsLocallySolvable(Z_planar, 2 : Smooth, AssumeIrreducible);
end function;

// Returns the fields and class groups needed for computing the 2-Selmer group, conditional on GRH.
// Assumes E is given in the form y^2 = p(x).
function two_selmer_class_groups(E)
    p := HyperellipticPolynomials(E);
    A := AbsoluteAlgebra(quo< Parent(p) | p >);
    Cl := [ClassGroup(F : Proof := "GRH") : F in Components(A)];
    cls := [<A[i], Cl[i]> : i in [1 .. #Cl]];
    return [c : c in cls | Degree(c[1]) gt 1];
end function;

function twist_ell_curve(f, root, g, delta)
    if Degree(g) eq 1 then
        K := RationalField();
        w := Roots(g)[1][1];
    else
        K<w> := NumberField(g);
    end if;
    A := ChangeRing(Parent(f), K);
    h := A!f div ((A.1 - root) * (A.1 - w));
    D := genus1_quotient(h, root, delta);
    // Minimise and reduce D
    D_min := Reduce(Minimise(GenusOneModel(D)));
    E := Jacobian(D_min);
    return E;
end function;

/* Encodes a polynomial in several variables over a number field in a JSON-compatible format */
function pack_polynomial(poly)
    coeffs, monomials := CoefficientsAndMonomials(poly);
    coeff_data := [Eltseq(c) : c in coeffs];
    exponents := [Exponents(m) : m in monomials];
    return <coeff_data, exponents>;
end function;

function pack_map(fn)
    return [pack_polynomial(poly) : poly in DefiningEquations(fn)];
end function;

/* Reconstructs a polynomial from JSON-compatiable data */
function unpack_polynomial(R, coeff_data, exponents)
    K := BaseRing(R);
    monomials := [Monomial(R, [Integers()!e : e in exps]) : exps in exponents];
    coeffs := [K!c : c in coeff_data];
    return Polynomial(coeffs, monomials);
end function;

function unpack_map(R, eqn_data)
    return [unpack_polynomial(R, eqn[1], eqn[2]) : eqn in eqn_data];
end function;

function twist_chabauty_map(f, root, g, delta, base_pt, E_aInv)
    if Degree(g) eq 1 then
        K := RationalField();
        w := Roots(g)[1][1];
    else
        K<w> := NumberField(g);
    end if;
    phi := genus1_map(f, root / 1, w, delta);
    ZK := Domain(phi);
    D := Codomain(phi);
    // Minimise and reduce D
    D_min := Reduce(Minimise(GenusOneModel(D)));
    D_min_hyp := HyperellipticCurve(D_min);
    is_iso, min_map := IsIsomorphic(D, D_min_hyp);
    assert is_iso;
    E := EllipticCurve([K!a : a in E_aInv]);
    E1, j := EllipticCurve(D_min_hyp, min_map(phi(ZK!base_pt)));
    E1_to_D := Inverse(j) * Inverse(min_map);
    // This model of the Jacobian isn't optimal, so identify it with E for faster Mordell-Weil computations.
    is_iso, E_to_E1 := IsIsomorphic(E, E1);
    assert is_iso;
    // Construct the map to P^1 for elliptic Chabauty
    Ecov := Expand(E_to_E1 * E1_to_D * map< D -> ProjectiveSpace(Rationals(), 1) | [D.1, D.3] >);
    return Ecov;
end function;

function twist_chabauty(f, root, g, delta, E_aInv, Ecov_data, MW_orders, MW_gens)
    if Degree(g) eq 1 then
        K := RationalField();
        w := Roots(g)[1][1];
    else
        K<w> := NumberField(g);
    end if;
    E := EllipticCurve([K!a : a in E_aInv]);
    P1 := ProjectiveSpace(Rationals(), 1);
    Ecov := map< E -> P1 | unpack_map(CoordinateRing(Ambient(E)), Ecov_data) >;

    A := AbelianGroup(MW_orders);
    gens := [E![K!c : c in coords] : coords in MW_gens];
    mw := map< A -> E(K) | a :-> &+[E(K) | e[i] * gens[i] : i in [1 .. #gens]] where e is Eltseq(a) >;
    mw_map := map< A -> E | a :-> mw(a) >;
    // Run elliptic Chabauty and record the set of points and the index parameter R
    V, R := Chabauty(mw_map, Ecov);
    x_coords := [Ecov(mw_map(pt)) : pt in V];

    dup := twisted_duplication_map(f, root, delta);
    preimages := [Pullback(dup, Codomain(dup)!Eltseq(pt)) : pt in x_coords];
    rational_pts := &join[{pt : pt in Points(S)} : S in preimages];
    return x_coords, rational_pts;
end function;

/*************************************************************
  Section 7: Numerical computation of degree of isogeny
 *************************************************************/

/*
   Given a list of six distinct roots (e.g. [0, 1, -1, 2, -2, 3]), analytically compute an isogeny between two abelian varieties:
   J = Jac(Z), where Z is a genus 5 curve associated by two-cover descent to the genus 2 curve C: y^2 = f(x), where the roots of f are the input;
   A = E_1 x E_2 x E_3 x E_4 x E_5, where the E_i are the elliptic curves constructed as the codomain of genus1_map.
*/
function analytic_isogeny_computation(roots : precision := 30)
    assert #roots eq 6 and #Seqset(roots) eq 6;
    precision := Max(precision, 30); // precision < 30 doesn't work
    root := roots[1];
    other_roots := roots[2 .. 6];
    QQ := Rationals();
    CC := ComplexField(precision);
    Rt<t> := PolynomialRing(QQ);
    f := &*[(t - r) : r in roots];

    // Construct an equation of a planar model of the genus 5 curve
    Z := genus5_canonical(f, root, 1);
    bool, pr := Genus5PlaneCurveModel(Z);
    if bool then
        C := Codomain(pr);
    else
        C, pr1 := Projection(Z);
        C, pr2 := Projection(C);
        pr := Expand(pr1 * pr2);
    end if;
    Qxy<x, y> := PolynomialRing(QQ, 2);
    poly_S := Evaluate(DefiningEquation(C), [x, y, 1]);

    // Construct the genus 1 curves algebraically
    Ds := [genus1_quotient(f div ((t - root)*(t - w)), root, 1) : w in other_roots];

    // Construct the Riemann surface associated to Z
    S := RiemannSurface(poly_S : Precision := precision);
    // Construct the elliptic curves as Riemann surfaces
    Es := [RiemannSurface(HyperellipticPolynomials(D), 2 : Precision := precision) : D in Ds];

    // Small period matrix of the analytic Jacobian of Z
    P_J := SmallPeriodMatrix(S);
    // Small period matrix of the product of the elliptic curves
    P_A := DiagonalMatrix(CC, [SmallPeriodMatrix(E)[1][1] : E in Es]);

    // Find an isogeny between P_J and P_A
    bool, M := IsIsogenousPeriodMatrices(P_J, P_A);
    // The determinant of M should be 32
    return bool, M, Determinant(M);
end function;

