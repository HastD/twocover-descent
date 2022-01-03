*** /opt/magma/magma2.26-10/package/Geometry/CrvG2/chabauty-MWS.m	2021-05-17 09:58:53.000000000 +0200
--- chabauty-MWS.m	2022-01-02 17:28:20.888627295 +0100
***************
*** 93,96 ****
--- 93,102 ----
      -- Some slight editing of verbose messages.
  
+   January 2022, Michael Stoll:
+     -- prompted by a bug report by Daniel Hast that indicated that
+        DiscreteLogMapSmooth() gives an error when the group is trivial,
+        fixed this intrinsic so that it works correctly also in the
+        trivial case.
+ 
   ***********************************************************************/
  
***************
*** 150,154 ****
--- 156,168 ----
   (This intrinsic is temporary and will be removed in a later release)}
  
+   require Domain(m) eq G: "the second argument must be a map with domain the first argument";
    n := #Invariants(G);
+   // added 2022-01-02, Michael Stoll
+   if n eq 0 then
+     // Group is trivial: return unique map from X to G.
+     // Without catching this special case, the assert below would throw an error.
+     return map<Codomain(m) -> G | x :-> G!0>;
+   end if;
+   //
    f := Factorization(#G);
    cofs := [&*[Integers()|f[i,1]^f[i,2] : i in [1..#f] | i ne j] : j in [1..#f]];
