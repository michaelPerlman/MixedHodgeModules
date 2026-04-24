--Copyright 2026 by Andras Lorincz and Michael Perlman



--test functions--

testHodgeIdealforWeightedHomog = method(Options => {UseGenLevel => True});

--compares hodgeIdeal to hodgeIdealWeightedHomogIsolated for all
--interesting alpha for p steps (not just p).
--Use "UseGenLevel => False" to turn off generation level use in hodgeIdeal

testHodgeIdealforWeightedHomog(RingElement,ZZ,List) := options -> (f,p,w) -> (

--w is the weight of f such that w*(exponent f)=1
--p is the max index you want to check

--output: triples {a,i,eqCheck} where eqCheck is falset
--prints all ideals and truth values along the way

   alphaToCheck := keys hodgeOnV(f,0);--just uses p=0 to find "interesting alphas". Change?
   badPairs := for i from 0 to p list (
       for a in alphaToCheck list (
	   I := hodgeIdeal(f,a,i, UseGenLevel => options.UseGenLevel);
	   print I;
	   J := hodgeIdealWeightedHomogIsolated(f,a,i,w);
	   print J;
	   eqCheck := I == J;
	   {a,i,eqCheck}));
   badPairs = flatten badPairs;
  select(badPairs, ai -> (ai_2) == false)
   )
   

------------------------------------------------------------------------

testBrieskornPham = method(Options => {UseGenLevel => True});

--compares hodgeIdeal to hodgeIdealBrieskornPham for all
--interesting alpha for p steps (not just p).
--Use "UseGenLevel => False" to turn off generation level us in hodgeIdeal

testBrieskornPham(List,ZZ) := options -> (L,p) -> (

--L is a list of exponents
--p is the max index you want to check

--output: triples {a,i,eqCheck}  where eqCheck is false
--prints all ideals and truth values along the way

   n := #L;
   z := local z;
   S := QQ[z_1..z_n];
   f := sum apply(toList(0..n-1), i-> S_i^(L_i));

   alphaToCheck := keys hodgeOnV(f,0);--just uses p=0 to find "interesting alphas". Change?
   badPairs := for i from 0 to p list (
       for a in alphaToCheck list (
	   I := hodgeIdeal(f,a,i, UseGenLevel => options.UseGenLevel);
	   print I;
	   J := hodgeIdealBrieskornPham(L,a,i);
	   print J;
	   R := ring J;
           F := map(R,S,gens R);
           JS := preimage(F,J);
	   eqCheck := I == JS;
	   {a,i,eqCheck}));
   badPairs = flatten badPairs;
   select(badPairs, ai -> (ai_2) == false)
   )


--------------------------------------------------------------------

testWeightedHodge = (f,pMax,mMax) -> (
--uses weightCheck to check weights of elements of weighted Hodge ideals
--checks for all interesting alpha and all levels \leq pMax,mMax
--returns a list of all the failures
--prints ALL ideals along the way

  alphaToCheck := keys hodgeOnV(f,0);
  bad := {};
  for wt from 0 to mMax do (
    for p from 0 to pMax do (
      for a in alphaToCheck do (
        I := weightedHodgeIdeal(f,a,p,wt);
        L := flatten entries gens I;
	print L;
        failures := select(L, g -> not weightCheck(f,g,a+p,wt));
        if #failures > 0 then bad = append(bad, {a,p,wt,failures});
      );
    );
  );
  bad
);


------------------------------------------------------------------

compareHodgeIequalHigherMultI = method(Options => {UseGenLevel => True});

--compares hodgeIdeal_k to higher multiplier ideal_k for all
-- 0\leq k\leq p and all alpha in keys hodgeOnV
--Use "UseGenLevel => False" to turn off generation level use in hodgeIdeal

--Important: these should only be equal for "nice" classes of f
-- eg semi-invariant functions

--see cusp example below and compare to Davis--Yang page 4
--compare to Zhang Conjecture E

compareHodgeIequalHigherMultI(RingElement,ZZ) := options -> (f,p) -> (

--p is the max index you want to check

--output: triples {a,i,eqCheck} where eqCheck is false
--prints all ideals and truth values along the way

   alphaToCheck := keys hodgeOnV(f,0);--just uses p=0 
   neqPairs := for i from 0 to p list (
       for a in alphaToCheck list (
	   I := hodgeIdeal(f,a,i, UseGenLevel => options.UseGenLevel);
	   print I;
	   J := higherMultiplierIdeal(f,a,i);
	   print J;
	   eqCheck := I == J;
	   {a,i,eqCheck}));
   neqPairs = flatten neqPairs;
   select(neqPairs, ai -> (ai_2) == false)
   )



------------------------------------------------------------------

testHodgeIhigherMultImodf = method(Options => {UseGenLevel => True});

-- compares hodgeIdeal_k+(f) to higher multiplier ideal_k + (f) for all
-- 0\leq k\leq p and all alpha in keys hodgeOnV
-- should be equal by MP20 Theorem A
-- checks equality mod (f): I + (f) == J + (f)
-- prints minimal generators of the sums along the way

--Use "UseGenLevel => False" to turn off generation level use in hodgeIdeal

testHodgeIhigherMultImodf(RingElement,ZZ) := options -> (f,p) -> (

   R := ring f;
   fIdeal := ideal(f);

   alphaToCheck := keys hodgeOnV(f,0); -- just uses p=0

   neqTriples := for i from 0 to p list (
       for a in alphaToCheck list (

           I := hodgeIdeal(f,a,i, UseGenLevel => options.UseGenLevel);
           J := higherMultiplierIdeal(f,a,i);

           Iplus := ideal mingens (I + fIdeal);
           Jplus := ideal mingens (J + fIdeal);

           print Iplus;
           print Jplus;

           eqCheck := Iplus == Jplus;

           {a,i,eqCheck}
       )
   );

   neqTriples = flatten neqTriples;
   select(neqTriples, ai -> (ai_2) == false)
)


