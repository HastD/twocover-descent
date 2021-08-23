load "twocovers.m";

t_total := Time();

C := HyperellipticCurve(); // TODO
SearchBound := 10000;
Effort := 1;
OutputFile := LABEL * "-results.txt";
SetClassGroupBounds("GRH");

fprintf OutputFile, "Curve:\n%o\n\n", C;
C_even, phi := even_weierstrass_model(C);
f := HyperellipticPolynomials(C_even);
fprintf OutputFile, "Model:\ny^2 = f(x) := %o\nChange of coordinates:\n%o\n\n", f, phi;
root := Roots(f)[1][1];
twists := Setseq(products_of_subsets(twist_param_generators(f)));
fprintf OutputFile, "Twists:\n%o\n\n", twists;
Zs := [genus5_canonical(f, root, delta) : delta in twists];
pts_search := [* PointSearch(Z, 10000) : Z in Zs *];
fprintf OutputFile, "Results of initial point search:\n%o\n\n", pts_search;
P1 := ProjectiveSpace(RationalField(), 1);
class_group_fields := [];

// We'll try using each factor of the quintic f/(x - root)
g_quintic := f div (Parent(f).1 - root);
gs := [fact[1] : fact in Factorization(g_quintic)];
verified_twist_count := 0;
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
            verified_twist_count :+= 1;
            break;
        end if;
    end for;
end for;

success := (verified_twist_count eq #twists);

for cgf in class_group_fields do
    F := cgf[1];
    Cl := cgf[2];
    fprintf OutputFile, "To remove dependence on GRH, need to verify that the number field with defining polynomial\n%o\n" cat
        "has class number %o.\n", Parent(f)!DefiningPolynomial(F), #Cl;
end for;
fprintf OutputFile, "\n";
if success then
    fprintf OutputFile, "SUCCESS\n";
else
    fprintf OutputFile, "FAILURE\n";
end if;
fprintf OutputFile, "Done.\nTotal time: %o\n", Time(t_total);

quit;

