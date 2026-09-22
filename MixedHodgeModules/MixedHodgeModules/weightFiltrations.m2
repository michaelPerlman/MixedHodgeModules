--Copyright 2026 by Andras Lorincz and Michael Perlman
--
-- Weighted-Hodge-ideal code and related downstream invariants.  The V-filtration
-- helpers and the core weightHodgeOnV computations
-- live in vFiltrations.m2.



---------------------------------------------------------------
---------------------------------------------------------------

--ref: Theorem A of Olano WEIGHTED HODGE IDEALS OF REDUCED DIVISORS

weightedHodgeIdeal = method();

weightedHodgeIdeal(RingElement, ZZ, ZZ, ZZ) :=
weightedHodgeIdeal(RingElement, QQ, ZZ, ZZ) := (f,alpha,p,m) -> (
--alpha is a rational number in (0,1]
--p and m and non-negative integers
--calculates the weighted Hodge ideal I^m_p(alpha*f)

     checkAlpha alpha;
     checkP p;
     if m < 0 then error "expected weight m to be a non-negative integer";
     R := ring f;
     if not R.cache#?WeightedHodgeIdealCache then R.cache#WeightedHodgeIdealCache = new MutableHashTable;
     tbl := R.cache#WeightedHodgeIdealCache;
     key := (f, sub(alpha,QQ), p, m);
     if tbl#?key then return tbl#key;
     WF := weightHodgeOnV(f,alpha,p,m);
     I := ideal mingens ideal malgrangeEval(WF, f, alpha, p);

     tbl#key = I;
     I
     )


---------------------------------------------------------------
---------------------------------------------------------------

--ref: Theorem A of Olano WEIGHTED MULTIPLIER IDEALS OF REDUCED DIVISORS

adjointIdeal = method();

adjointIdeal(RingElement) := f -> (
--gives the adjoint ideal of f
--see Lazarsfeld Positivity II, Section 9.3E
--the divisor f has rat sing iff = ideal(1_R)

   weightedHodgeIdeal(f,1,0,1)
   )


---------------------------------------------------------------
---------------------------------------------------------------



bMinMax= method();
--calculates k1 and k2 from Lorincz--Yang paper

--used in: weightCheck, weightLevel


bMinMax(RingElement, ZZ) :=
bMinMax(RingElement, QQ) := (bfs,alpha) -> (-----bfs is the b-function of f
    Roots := bFunctionRoots(bfs);
    alphaInts := sort select(Roots, i-> denominator(sub(i+alpha,QQ)) == 1);
    sizeAI := #alphaInts;
    k1 := 0;
    k2 := 0;
    if sizeAI > 0 then (
        if alphaInts_0 <= -alpha then (
            k1 = -(alphaInts_0+alpha));
        if alphaInts_(sizeAI-1 ) >= -alpha then (
            k2 = alphaInts_(sizeAI-1) + alpha + 1));
    {k1, k2}
    )



----------------------------------------------------------------
---------------------------------------------------------------


weightCheck=method();
--checks if g/f^(alpha) is in weight w piece of S_f*f^(-alpha)
--if beta\in (0,1] and k\in ZZ  with alpha = k + beta then
--checks if (g/f^k)*f^(-beta) is in weight w piece of S_f*f^(-beta)

--uses: bMinMax

weightCheck(RingElement,RingElement,ZZ,ZZ) :=
weightCheck(RingElement,RingElement,QQ,ZZ) := (f,g,alpha,w)-> (
  if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
  if w < 0 then error "expected weight w to be a non-negative integer";
  b := cachedGlobalBFunction f;
  bMM := bMinMax(b,alpha);
  k1 := floor bMM_0;
  k2 := floor bMM_1;
  If := cachedAnnFs f;
  varDsFs := flatten entries vars ring If;
  s := varDsFs_(#varDsFs-1);
  newf := substitute(f, ring If);
  newg := substitute(g, ring If);
  use ring If;
  J := If + ideal((s + alpha + k1)^(w + 1), newf^(k1 + k2));
  result := isMember((((s + alpha + k1)^w) * newg * newf^(k1))_(ring If),J);
  use ring f;--restore caller's ring after `use ring If`
  result
  )


----------------------------------------------------------------
---------------------------------------------------------------



weightLevel= method();
--finds minimal weight level of g/f^(alpha) in S_f*f^(-alpha)
--if beta\in (0,1] and k\in ZZ  with alpha = k + beta then
--finds minimal weight level of (g/f^k)*f^(-beta) in S_f*f^(-beta)


--uses: bMinMax
--used by: weightLength(ByWeightLevel)

weightLevel(RingElement, RingElement, ZZ) :=
weightLevel(RingElement, RingElement, QQ) := (f,g,alpha) -> (
    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    bfs := cachedGlobalBFunction f;
    weightLevel(f,g,bfs,alpha)
    )


weightLevel(RingElement, RingElement, RingElement, ZZ) :=
weightLevel(RingElement, RingElement, RingElement, QQ) := (f,g,bfs,alpha) -> (
    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    bMM := bMinMax(bfs,alpha);
    k1 := floor bMM_0;
    --Test from weight zero: polynomial numerators can lower the weight level.
    w := -1;
    Roots := bFunctionRoots bfs;
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i + alpha,QQ)) == 1);
    if (#AlphaInts == 0) or (AlphaInts_0 < -alpha) then w = 0 else (
    	k2 := floor bMM_1;
    	If := cachedAnnFs f;-----hard part
    	varDsFs := flatten entries vars ring If;
    	s := varDsFs_(#varDsFs-1);
    	newf := substitute(f, ring If);
    	newg := substitute(g, ring If);
    	use ring If;
    	mem := false;
    	while mem == false do (
            w = w + 1;
            J := If + ideal((s + alpha + k1)^(w + 1), newf^(k1 + k2));
	    h := (((s + alpha + k1)^w) * newg * newf^(k1))_(ring If);
            mem = ((h%J) == 0_(ring If))));
    	use ring f;--restore caller's ring after `use ring If`
    	w
    	)



----------------------------------------------------------------
----------------------------------------------------------------


weightLength = method(Options => {LengthStrategy => ByNuAlpha, NuMethod => ByAnnFs})
--finds maximal weight w such that Gr^W_w(S_f*f^(-alpha)) is nonzero
--may assume alpha\in (0,1]

--uses weightLevel (if ByWeightLevel) or nuAlpha (if ByAnnFs)
--default: ByAnnFs


weightLength(RingElement, ZZ) :=
weightLength(RingElement, QQ) := options -> (f,alpha) -> (
    checkAlpha alpha;
    if options.LengthStrategy == ByWeightLevel then (
        weightLengthByLevel(f,alpha))
    else weightLengthByNu(f,alpha, NuMethod => options.NuMethod)
    )



--strategies


weightLengthByLevel= method();
--formerly called weightLength

--uses: weightLevel

weightLengthByLevel(RingElement, ZZ) :=
weightLengthByLevel(RingElement, QQ) := (f,alpha) -> (
    bfs := cachedGlobalBFunction f;
    w := 0;
    Roots := bFunctionRoots(bfs);
    alphaInts := sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
    if #alphaInts >0 then (
	 k0 := alphaInts_0;
	 w = weightLevel(f,1_(ring f),bfs,-k0));
     w
     )


 ----------------------------------------------------------------



weightLengthByNu=method(Options => {NuMethod => ByAnnFs});
--formerly called slowWeightLength
--user may select NuMethod using syntax above
--default is ByAnnFs

--uses: nuAlpha

weightLengthByNu(RingElement, ZZ) :=
weightLengthByNu(RingElement, QQ) := options -> (f,alpha) -> (
     bfs := cachedGlobalBFunction f;
     w := 0;
     Roots := bFunctionRoots(bfs);
     alphaInts :=sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
     if #alphaInts > 0 then (
	 k0 := alphaInts_0;
	 w = nuAlpha(f,1_(ring f), bfs, -k0, NuMethod => options.NuMethod);
	 );
     w
     )



----------------------------------------------------------------
----------------------------------------------------------------

localCohomFW = method();

localCohomFW(Ideal, ZZ, ZZ, ZZ) := (I,q,p,m) -> (

    S := ring I;
    n := numgens S;

    if m<n then error "the weight should be at least numgens(ring I)";

    gensI := flatten entries gens I;
    fJs := subsets(gensI,q);

    WFCqSmds := apply(fJs, J -> module(weightedHodgeIdeal(product(J),1,p,m-n)));
    WFCq := directSum WFCqSmds;

    gensIpMat := matrix{ apply(gensI, g -> g^(p+1))};

    dq := matrix entries transpose koszul(q+1, gensIpMat);
    dqm1 := matrix entries transpose koszul(q, gensIpMat);

    if ((rank target dq) == 0) then dq = matrix{{0_S}};

    Dq := map(target dq, super WFCq, dq);
    Dqm1 := map(super WFCq, , dqm1);

    Kq := intersect(ker(Dq), WFCq);
    Imq := intersect(image(Dqm1), WFCq);
    Hq := Kq/Imq;

    prune Hq
    )



end
