/*************************************************************
  Test code for twocovers.m
 *************************************************************/

load "twocovers.m";

SetClassGroupBounds("GRH");
SetVerbose("MordellWeilGroup", 1);
SetVerbose("MWSha", 1);
SetVerbose("Selmer", 1);
SetVerbose("TwoDescent", 1);
SetVerbose("CasselsTate", 1);
SetVerbose("EllChab", 1);

f := HyperellipticPolynomials(even_weierstrass_model(C));
root1 := Roots(f)[1][1]; root1;
g := f div (x - root1);   
K<w> := NumberField(g);
root2 := w;
bound := 10000;

twists := Setseq(products_of_subsets(twist_param_generators(f)));

Zs := [genus5_canonical(f, root1, delta) : delta in twists];
[* PointSearch(Z, 10000) : Z in Zs *];

phis := [genus1_map(f, root1, root2, delta) : delta in twists];
Ds := [Reduce(Minimise(GenusOneModel(Codomain(phi)))) : phi in phis];
Es := [Jacobian(D) : D in Ds];

// Genus 2 curve 533963.a.533963.1
R<x> := PolynomialRing(Rationals());
f := x^6 + 2*x^5 + 7*x^4 + 4*x^3 - 13*x^2 - 10*x + 1;
root := -1;
g := f div (x - root);
K<w> := NumberField(g);
h := ChangeRing(g, K) div (ChangeRing(x, K) - K.1);

//B<xx, zz> := PolynomialRing(K, 2);
//h_homog := &+[Eltseq(h)[i+1] * xx^i * zz^(4-i) : i in [0 .. 4]];


//root_rand := Random(-10^4, 10^4);
//g_rand := A![Random(-10^5, 10^5) : i in [0 .. 5]];
//f_rand := (x - root_rand) * g_rand;

// Genus 2 curve 6982.a.13964.1
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(R![0, 2, 0, -3, 0, 1], R![1, 0, 0, 1]);

// Genus 2 curve 96431.a.96431.1
C := HyperellipticCurve(R![0, -1, 3, 0, -3, -1], R![1, 1, 1]);
f := HyperellipticPolynomials(even_weierstrass_model(C));

// Genus 2 curve 6443.a.6443.1
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 1, 1, -2, -1, 1], R![1]);
f := HyperellipticPolynomials(even_weierstrass_model(C));

// Genus 2 curve 141991.b.141991.1
R<x> := PolynomialRing(Rationals()); C := HyperellipticCurve(R![0, 0, 1, -2, -2, 1], R![1, 1, 1]);

SetVerbose("LocSol", 1);
SetVerbose("ClassGroup", 2);

// Genus 2 curve with a rational Weierstrass point and at least 303 rational points
load "twocovers.m";
SetClassGroupBounds("GRH");
R<x> := PolynomialRing(Rationals());
f2 := 98017920*x^5 - 3192575*x^4 - 274306650*x^3 + 256343425*x^2 - 76075320*x + 2740^2;
X2 := HyperellipticCurve(f2);
C := even_weierstrass_model(X2);
f := HyperellipticPolynomials(C);
pts := Points(C : Bound:=10000);
pts_J := [P - pts[1] : P in pts];
// time bas := ReducedBasis(pts_J); 
L<X> := quo< R | f>;
// gens := [Evaluate(cassels_map(P), X) : P in bas];
gens := [ // precomputed
    7507600*X^5 - 76075320*X^4 + 256343425*X^3 - 274306650*X^2 - 3192576*X + 98017920,
    7507600*X^5 - 76075320*X^4 + 256343425*X^3 - 274306649*X^2 - 3192576*X + 98017920,
    -7507600*X^5 + 76075320*X^4 - 256343425*X^3 + 274306651*X^2 + 3192576*X - 98017920,
    22522800*X^5 - 228225960*X^4 + 769030275*X^3 - 822919949*X^2 - 9577728*X + 294053760,
    30030400/3*X^5 - 101433760*X^4 + 1025373700/3*X^3 - 365742199*X^2 - 4256768*X + 130690560,
    -5630700*X^5 + 57056490*X^4 - 769030275/4*X^3 + 411459977/2*X^2 + 2394432*X - 73513440,
    1876900*X^5 - 19018830*X^4 + 256343425/4*X^3 - 137153323/2*X^2 - 798144*X + 24504480,
    37538000*X^5 - 380376600*X^4 + 1281717125*X^3 - 1371533249*X^2 - 15962880*X + 490089600,
    9384500*X^5 - 95094150*X^4 + 1281717125/4*X^3 - 685766623/2*X^2 - 3990720*X + 122522400,
    -10510640*X^5 + 106505448*X^4 - 358880795*X^3 + 384029311*X^2 + 22348032/5*X - 137225088,
    60060800/7*X^5 - 608602560/7*X^4 + 2050747400/7*X^3 - 2194453193/7*X^2 - 25540608/7*X + 112020480,
    6569150*X^5 - 66565905*X^4 + 1794403975/8*X^3 - 960073271/4*X^2 - 2793504*X + 85765680,
    X^2 - 453/217*X + 176/217,
    82583600*X^5 - 836828520*X^4 + 2819777675*X^3 - 3017373149*X^2 - 35118336*X + 1078197120,
    X^2 + 5/7*X - 33/2,
    24024320*X^5 - 243441024*X^4 + 820298960*X^3 - 877781279*X^2 - 51081216/5*X + 313657344,
    X^2 + 241/25*X - 34/25,
    210212800/19*X^5 - 2130108960/19*X^4 + 7177615900/19*X^3 - 7680586181/19*X^2 - 89392128/19*X + 2744501760/19
];
g := f div x;   
K<w> := NumberField(g);
bad_primes := [2, 3, 5, 7, 11, 13, 17, 18517, 89241827, 322251693061, 15307704844889];
small_primes := PrimesUpTo(98);
large_bad_primes := [p : p in bad_primes | p gt 98];

