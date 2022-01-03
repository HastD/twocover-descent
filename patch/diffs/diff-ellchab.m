*** /opt/magma/magma2.26-10/package/Geometry/Algaem/ellchab.m	2021-03-30 01:23:06.000000000 +0200
--- ellchab.m	2022-01-02 17:17:28.556617455 +0100
***************
*** 247,253 ****

    for p in Prs do
!     Kp,toKp:=MyCompletion(p);
!     error if Valuation(Kp!Minimum(p))/(Minimum(p)-1) ge 1,
        "Prime has too high ramification compared to residue characteristic.";

      // ensure BaseRing(Kp) is pAdicField(Minimum(p))
--- 247,265 ----

    for p in Prs do
!     // Changed the following to ensure that the increased precision
!     // is actually used in the computation.
!     // Michael Stoll, 2022-01-02
!
!     // Determine (absolute) ramification index of p.
!     e_p := RamificationIndex(p, Minimum(p));
!     // Check that ramification is not too high.
!     error if e_p ge Minimum(p) - 1,
        "Prime has too high ramification compared to residue characteristic.";
+     // Define starting precision.
+     prc := 2*(pAdicPrec+1)*e_p;
+     repeat
+       // Compute with current precision.
+       // If not successful, double precision and repeat.
+       Kp,toKp:=MyCompletion(p : Precision := prc);

        // ensure BaseRing(Kp) is pAdicField(Minimum(p))
***************
*** 256,275 ****
      if Degree(Kp) eq 1 then
        Kp := BaseRing(Kp);
      end if;

      pi:=UniformizingElement(Kp);
-     prc:=2*(pAdicPrec+1)*AbsoluteRamificationDegree(Kp);
-     oldprec:=Kp`DefaultPrecision;
-     repeat
-       Kp`DefaultPrecision:=prc;
        EKp:=PointSet(E,map<Domain(toKp)->Codomain(toKp)|a:->toKp(a)>);
-       error if prc gt Kp`DefaultPrecision, "Insufficient precision available";
        Gp:=[EKp![toKp(c)+O(pi^prc):c in Eltseq(MWmap(g))]:
                     g in OrderedGenerators(Domain(MWmap))];
!       L:=[p where p:=&+[v[i]*Gp[i]:i in [1..#Gp]]: v in KerPrsVec];
        zBp:=[-p[1]/p[2]:v in KerPrsVec | RelativePrecision(p[2]) gt 0 where p:=&+[v[i]*Gp[i]:i in [1..#Gp]]];
        prc:=prc*2;
      until #zBp eq #KerPrsVec and MinPrec(zBp) ge (pAdicPrec+1)*RamificationDegree(Kp);
-     Kp`DefaultPrecision:=oldprec;

      zBmat:=zBmat cat [[Integers()!(Eltseq(a)[i]/Minimum(p)): a in zBp]:
--- 268,285 ----
        if Degree(Kp) eq 1 then
          Kp := BaseRing(Kp);
+         // change toKp accordingly
+         toKp := map<Domain(toKp) -> Kp | a :-> Kp!toKp(a)>;
        end if;

+       vprintf EllChab, 3: "==> p-adic precision %o\n", prc;
        pi:=UniformizingElement(Kp);
        EKp:=PointSet(E,map<Domain(toKp)->Codomain(toKp)|a:->toKp(a)>);
        Gp:=[EKp![toKp(c)+O(pi^prc):c in Eltseq(MWmap(g))]:
                     g in OrderedGenerators(Domain(MWmap))];
!       L:=[&+[v[i]*Gp[i]:i in [1..#Gp]]: v in KerPrsVec];
        zBp:=[-p[1]/p[2]:v in KerPrsVec | RelativePrecision(p[2]) gt 0 where p:=&+[v[i]*Gp[i]:i in [1..#Gp]]];
+       vprintf EllChab, 3: "==> zBp = %o\n", zBp;
        prc:=prc*2;
      until #zBp eq #KerPrsVec and MinPrec(zBp) ge (pAdicPrec+1)*RamificationDegree(Kp);

      zBmat:=zBmat cat [[Integers()!(Eltseq(a)[i]/Minimum(p)): a in zBp]:
