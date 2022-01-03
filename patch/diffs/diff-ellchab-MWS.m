*** /opt/magma/magma2.26-10/package/Geometry/Algaem/ellchab-MWS.m	2021-03-30 01:23:06.000000000 +0200
--- ellchab-MWS.m	2022-01-02 17:29:44.416628555 +0100
***************
*** 32,36 ****
    nEp:=[#e: e in Ep];
    vprint EllChab,2: "Group orders:",nEp;
!   idx:=[i: i in [1..#Ep] | Max(PrimeBasis(nEp[i])) le SmoothBound];
    if #idx eq 0 then
      vprint EllChab,2: "-"^40;
--- 32,38 ----
    nEp:=[#e: e in Ep];
    vprint EllChab,2: "Group orders:",nEp;
!   // Changed by Michael Stoll, 2022-01-02, to allow nEP[i] = 1.
!   idx:=[i: i in [1..#Ep] | forall{p : p in PrimeBasis(nEp[i]) | p le SmoothBound}];
! //   idx:=[i: i in [1..#Ep] | Max(PrimeBasis(nEp[i])) le SmoothBound];
    if #idx eq 0 then
      vprint EllChab,2: "-"^40;
***************
*** 317,328 ****
     InitialPrimes: Information at all good primes below this bound will be used initially.}

    if ISA(Type(IndexBound),RngIntElt) and IndexBound ne -1 then
       IndexBound:=PrimeBasis(IndexBound);
    end if;

-   // MW: should have a bunch of requires, ... ?
    if IsFinite(Domain(MWmap)) then return finite_chabauty(MWmap,Ecov); end if;

-   E:=Codomain(MWmap);
    BadPrimes:=Seqset(PrimeBasis(&*[Denominator(a): a in aInvariants(E)]) cat
               PrimeBasis(Discriminant(IntegerRing(BaseRing(E)))) cat
--- 319,344 ----
     InitialPrimes: Information at all good primes below this bound will be used initially.}

+   // Added by Michael Stoll, 2022-01-02:
+   // Some argument checking, plus replacing point set by elliptic curve if necessary.
+   // This should allow for using the map returned by MordellWeilGroup directly as MWmap.
+   E := Codomain(MWmap);
+   if Type(E) eq SetPtEll then
+     // Change codomain of MWmap (and E) to be the curve, not the point set.
+     E := Curve(E);
+     MWmap := map<Domain(MWmap) -> E | x :-> MWmap(x)>;
+   end if;
+   require Domain(Ecov) eq E: "Cover Ecov must have domain the codomain of MWmap";
+   require ISA(Type(Codomain(Ecov)),Prj) and Dimension(Codomain(Ecov)) eq 1:
+     "Cover Ecov should be to a projective line";
+   require BaseRing(Codomain(Ecov)) eq Rationals():
+     "The cover Ecov should map to a projective line defined over the rationals.";
+   //
+
    if ISA(Type(IndexBound),RngIntElt) and IndexBound ne -1 then
       IndexBound:=PrimeBasis(IndexBound);
    end if;

    if IsFinite(Domain(MWmap)) then return finite_chabauty(MWmap,Ecov); end if;

    BadPrimes:=Seqset(PrimeBasis(&*[Denominator(a): a in aInvariants(E)]) cat
               PrimeBasis(Discriminant(IntegerRing(BaseRing(E)))) cat
