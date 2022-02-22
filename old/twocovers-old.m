/*
This file contains some functions that were used in an older iteration of the code, but are no longer used.
*/

/*************************************************************
  Section 6: Applying elliptic Chabauty to twists
 *************************************************************/

// Run elliptic Chabauty to determine points of Z_delta.
/*
   The first return value is a set of points, the second return value a boolean.
   The boolean indicates whether the set of points is proven to be complete.
   There are three reasons why the boolean could be false:
   (1) The curve is locally solvable but we can't find any rational points.
   (2) The rank of the elliptic curve is greater than [Q(root2) : Q].
   (3) We can't provably compute the Mordell-Weil group of an elliptic curve.
   In case (1), we give up and return the empty set.
   In case (2), we give up and return whatever points we can find with a search.
   In case (3), we return the set of points that maps to the subgroup that we're able to find.
   The third return value is a list of number fields and conditionally computed class groups.
*/
function chabauty_on_twist(f, root1, root2, delta : Effort := 1, Bound := 10000, PrimeBound := -1)
    Z := genus5_canonical(f, root1, delta);
    // Search for points up to the given height bound on the genus 5 curve
    pts_Z := PointSearch(Z, Bound);
    if IsEmpty(pts_Z) then // check if Z is locally solvable (looking only at primes under the bound)
        fact := Factorization(Conductor(HyperellipticCurve(f)));
        bad_primes_under_bound := [p[1] : p in fact | PrimeBound eq -1 or p[1] le PrimeBound];
        if not is_genus5_locally_solvable(Z, bad_primes_under_bound) then
            return {}, 0, true, []; // Z is not locally solvable, so Z(Q) is empty
        end if;
        // Z is locally solvable at finite primes but appears to have no points, so we give up.
        return {}, 0, false, [];
        // TODO: try Mordell-Weil sieve before giving up
    end if;
    // Since Z has a rational point, we use elliptic Chabauty
    phi := genus1_map(f, root1, root2, delta);
    ZK := Domain(phi);
    D := Codomain(phi);
    K := BaseRing(D);
    // Choose a point of minimum height
    min_ht := Minimum([HeightOnAmbient(pt) : pt in pts_Z]);
    base_pt := [ZK!Eltseq(pt) : pt in pts_Z | HeightOnAmbient(pt) eq min_ht][1];
    // Minimise and reduce D
    D_min := Reduce(Minimise(GenusOneModel(D)));
    D_min_hyp := HyperellipticCurve(D_min);
    is_iso, min_map := IsIsomorphic(D, D_min_hyp);
    assert is_iso;
    E := Jacobian(D_min);
    E1, j := EllipticCurve(D_min_hyp, min_map(phi(base_pt)));
    E1_to_D := Inverse(j) * Inverse(min_map);
    // This model of the Jacobian isn't optimal, so identify it with E for faster Mordell-Weil computations.
    is_iso, E_to_E1 := IsIsomorphic(E, E1);
    assert is_iso;
    // Construct the map to P^1 for elliptic Chabauty
    Ecov := E_to_E1 * E1_to_D * map< D -> ProjectiveSpace(Rationals(), 1) | [D.1, D.3] >;
    // Compute the Mordell-Weil group, and record whether it was provably computed
    A, mw, bool_rank, bool_index := MordellWeilGroup(E : Effort := Effort);
    mw_map := map< A -> E | a :-> mw(a) >;
    if TorsionFreeRank(A) eq 0 and Degree(MinimalPolynomial(root2)) eq 1 then
        V := Setseq(Set(A));
    elif TorsionFreeRank(A) lt Degree(MinimalPolynomial(root2)) then
        // Run elliptic Chabauty and record the set of points and the index parameter R
        V, R := Chabauty(mw_map, Ecov);
    else // Rank is too high, give up and return any points found
        return pts_Z, A, false, [];
    end if;
    VD := [E1_to_D(E_to_E1(mw_map(pt))) : pt in V];
    preimages := [Pullback(phi, pt) : pt in VD];
    rational_pts := &join[{Z!Eltseq(pt) : pt in Points(S) | IsCoercible(Z, Eltseq(pt))} : S in preimages];
    return rational_pts, A, (bool_rank and bool_index), two_selmer_class_groups(E);
end function;

// Input: a degree 6 polynomial f over Q, a rational root root1, an irrational root root2, and a height bound.
// Output: a list of possible x-coordinates of rational points on the genus 2 curve y^2 = f(x).
// The list might include some extraneous points. If the boolean second output is true, the list is complete.
// The third output is a tuple <points found on Z_delta, the twist parameter delta, whether the list is complete>.
function twocover_chabauty(f, root1, root2 : Bound := 10000)
    QQ := Rationals();
    assert Degree(f) eq 6 and BaseRing(f) eq QQ;
    assert IsCoercible(QQ, root1) and Evaluate(f, root1) eq 0;
    assert Evaluate(f, root2) eq 0 and Degree(MinimalPolynomial(root2)) ge 2;
    twists := products_of_subsets(twist_param_generators(f));
    Z_pts := [* *];
    for delta in twists do
        pts, _, bool := chabauty_on_twist(f, root1, root2, delta : Bound := Bound);
        Append(~Z_pts, <pts, delta, bool>);
    end for;
    verified := &and[t[3] : t in Z_pts];
    duplication_maps := [twisted_duplication_map(f, root1, delta) : delta in twists];
    P1 := ProjectiveSpace(QQ, 1);
    x_coords := &join[ {P1!Eltseq(pi(Eltseq(P))) where pi is duplication_maps[i] : P in Z_pts[i][1]} : i in [1 .. #twists] ];
    return x_coords, verified, Z_pts;
end function;

// Provably computes rational points on a curve using two-cover descent and elliptic Chabauty.
// Writes the output to the provided file name.
procedure descent_procedure(C : SearchBound := 10000, AssumeGRH := true, OutputFile := "", Effort := 1)
    t_total := Time();
    if #OutputFile eq 0 then
        OutputFile := Sprintf("twocover-descent-output-%o.txt", Realtime());
    end if;
    fprintf OutputFile, "Curve:\n%o\n\n", C;
    C_even, phi := even_weierstrass_model(C);
    f := HyperellipticPolynomials(C_even);
    fprintf OutputFile, "Model:\ny^2 = f(x) := %o\nChange of coordinates:\n%o\n\n", f, phi;
    root := Roots(f)[1][1];
    twists := Setseq(products_of_subsets(twist_param_generators(f)));
    fprintf OutputFile, "Twists:\n%o\n\n", twists;
    Zs := [genus5_canonical(f, root, delta) : delta in twists];
    pts_search := [* PointSearch(Z, SearchBound) : Z in Zs *];
    fprintf OutputFile, "Results of initial point search:\n%o\n\n", pts_search;
    P1 := ProjectiveSpace(RationalField(), 1);
    class_group_fields := [];
    // We assume GRH to start, then (if AssumeGRH is false) make the computation unconditional at the end.
    SetClassGroupBounds("GRH");

    // We'll try using each factor of the quintic f/(x - root)
    g_quintic := f div (Parent(f).1 - root);
    gs := [fact[1] : fact in Factorization(g_quintic)];
    for delta in twists do
        t0 := Time();
        for g in gs do
            if Degree(g) eq 1 then
                K := RationalField();
                w := Roots(g)[1][1];
            else
                K<w> := NumberField(g);
            end if;
            pts, A, verified, cls := chabauty_on_twist(f, root, w, delta : Effort := Effort);
            // Record new fields and class numbers to be verified later
            for c in cls do
                if not &or[IsIsomorphic(cgf[1], c[1]) : cgf in class_group_fields] then
                    Append(~class_group_fields, c);
                end if;
            end for;
            if verified then
                if IsEmpty(cls) then
                    verified_str := "yes (unconditionally)";
                else
                    verified_str := "yes (assuming GRH)";
                end if;
            else
                verified_str := "no";
            end if;
            pi := twisted_duplication_map(f, root, delta);
            x_coords := "";
            for P in pts do
                x := P1!Eltseq(pi(Eltseq(P)));
                x_coords cat:= Sprintf("  pi(%o) = %o\n", P, x);
            end for;
            if Type(A) eq GrpAb then
                mw_string := Sprintf("%o", A);
            else
                mw_string := "[not computed]";
            end if;
            fprintf OutputFile,
                "delta = %o\n" cat
                "g = %o\n" cat
                "Found rational points:\n%o\n" cat
                "Corresponding x-coordinates:\n%o" cat
                "List of points provably complete? %o\n" cat
                (verified select "" else "Subgroup of ") cat
                "Mordell-Weil group computed:\n%o\n" cat
                "Time: %o\n\n",
                delta, g, pts, x_coords, verified_str, mw_string, Time(t0);
            if verified then
                break;
            end if;
        end for;
    end for;
    for cgf in class_group_fields do
        F := cgf[1];
        Cl := cgf[2];
        fprintf OutputFile, "To remove dependence on GRH, need to verify that the number field with defining polynomial\n%o\n" cat
            "has class number %o.\n", Parent(f)!DefiningPolynomial(F), #Cl;
        if not AssumeGRH then
            t0 := Time();
            fprintf OutputFile, "Checking primes up to Minkowski bound %o...\n", MinkowskiBound(F);
            Cl_actual := ClassGroup(F : Proof := "Full");
            assert #Cl eq #Cl_actual;
            fprintf OutputFile, "Verified in %o seconds.\n", Time(t0);
        end if;
    end for;
    if not AssumeGRH then
        fprintf OutputFile, "\nThe results depending on GRH have now been unconditionally proven.\n";
    end if;
    fprintf OutputFile, "\n";
    fprintf OutputFile, "Done.\nTotal time: %o\n", Time(t_total);
end procedure;

