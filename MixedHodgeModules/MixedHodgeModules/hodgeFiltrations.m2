-- Copyright 2026 by Andras Lorincz and Michael Perlman
--
-- Hodge ideals and related downstream invariants.  The V-filtration helpers
-- and the core hodgeOnV / HRH computations live in vFiltrations.m2.



---------------------------------------------------------------
---------------------------------------------------------------



--if you want to test this function without generation level tricks,
--use strategy "UseGenLevel => False"

hodgeIdeal = method(Options => {UseGenLevel => True})


hodgeIdeal(RingElement,ZZ,ZZ) :=
hodgeIdeal(RingElement,QQ,ZZ) := options -> (f,alpha,p) -> (
-- f is a polynomial in polynomial ring R
-- alpha is a rational number in (0,1]
-- p is a non-negative integer
-- outputs the Hodge ideal I_p(alpha f)

    checkAlpha alpha;
    checkP p;
    R := ring f;
    if not R.cache#?HodgeIdealCache then R.cache#HodgeIdealCache = new MutableHashTable;
    tbl := R.cache#HodgeIdealCache;
    key := (f, sub(alpha,QQ), p, options.UseGenLevel);
    if tbl#?key then return tbl#key;
    n := numgens R;
    I := ideal(1_R);
    FpValpha := {};

    if options.UseGenLevel == False then (
       FpValpha = hodgeOnV(f,alpha,p);
       I = ideal mingens ideal malgrangeEval(FpValpha, f, alpha, p);
       )
       else (
             if (alpha == 1_QQ or alpha == 1_ZZ) and (p > max{n-2,0}) then (--MP Theorem B
                 I = generateNext(f,1,p-1))

               else if p > max{n-1,0} then (--MP20 Theorem E
                       I = generateNext(f,alpha,p-1))

               else (--do Blanco

                    FpValpha = hodgeOnV(f,alpha,p);
	            I = ideal mingens ideal malgrangeEval(FpValpha, f, alpha, p);
	            );
		);

    tbl#key = I;
    I
    )




---------------------------------------------------------------
---------------------------------------------------------------


generateNext = method()

generateNext(RingElement, ZZ, ZZ) :=
generateNext(RingElement, QQ, ZZ) := (f,alpha,p) -> (
--calculates the image of F_1(D)*F_p(R_f*f^(-alpha)) in F_(p+1)(R_f*f^(-alpha)) as an ideal in R
--the output will be a subset of I_(p+1)(alpha*f)
--alpha should be in (0,1]

    checkAlpha alpha;
    checkP p;
    R := ring f;
    n := numgens R;

    Ip := flatten entries gens hodgeIdeal(f,alpha,p);

    fgens := apply(Ip, m -> f*m);
    partialGens := flatten for m in Ip list (
	    for j from 0 to n-1 list f*diff(R_j,m)-(p+alpha)*m*diff(R_j,f));
    F1DFp := fgens|partialGens;

    J := ideal mingens ideal(F1DFp);
    J
    )




---------------------------------------------------------------
---------------------------------------------------------------

doesGenerateNext = method();

doesGenerateNext(RingElement, ZZ, ZZ) :=
doesGenerateNext(RingElement, QQ, ZZ) := (f,alpha,p) -> (
--tests if F_p(R_f*f^{-alpha}) generates F_(p+1)(R_f*f^{-alpha})

    checkAlpha alpha;
    checkP p;
    Ip1 :=  hodgeIdeal(f,alpha,p+1);
    J := generateNext(f,alpha,p);
    J == Ip1
    )


---------------------------------------------------------------
---------------------------------------------------------------

generationLevel = method();
--finds generation level on Hodge filtration of R_f*f^{-alpha}
--f should belong to a ring with at least two variables
--if alpha omitted, sets alpha=1

generationLevel(RingElement) := f -> (
    generationLevel(f,1)
    )


generationLevel(RingElement,ZZ) :=
generationLevel(RingElement,QQ) := (f,alpha) -> (

    checkAlpha alpha;
    R := ring f;
    n := numgens R;

    maxg := n-1;--MP20 Theorem E
    if (alpha == 1_QQ or alpha == 1_ZZ) then maxg = n-2;--MP Theorem B

    gend := true;
    p := maxg;

    while gend and p>=1 do (
	p = p-1;
	gend = doesGenerateNext(f,alpha,p);
	);

    if gend then p=0 else p=p+1;
    p
    )





---------------------------------------------------------------
---------------------------------------------------------------
-- weighted homogeneous isolated code

--helper function
weightGeqAlpha = (R,w,alpha) -> (
--R is poly ring with n vars
--w is a list, a weight vector with n many rational entries
--alpha is a rational number in (0,1]
--outputs R\geq Alpha as on page 2 of Zhang, as an ideal of R

--Acknowledgment: Timothy Duff told me how to code this function

   if #w != numgens R then
       error("weightGeqAlpha: expected length of w to equal numgens R = ", toString numgens R);
   if not all(w, wi -> sub(wi, QQ) > 0) then
       error "weightGeqAlpha: expected all weights to be positive rationals";

   n := numgens R;
   sumW := sum(w);

   ------------------------------------------------------------------------
   -- Step 1. find bounds for variables on coordinate axes
   ------------------------------------------------------------------------
   simplexVerts := apply( toList(0..n-1), i -> max(0,ceiling((-sumW+alpha)/w_i)+1));--+1 needed?

   ------------------------------------------------------------------------
   -- Step 2. find lattice points of polyhedron with faces: coordinate planes,
   -- bounds from Step 1, and plane pho(g)\geq alpha
   ------------------------------------------------------------------------

   zeroList := apply(toList (0..n-1), i -> 0);
   vMatrix := transpose matrix {append(simplexVerts|zeroList, sumW-alpha)};
   negWvector := transpose matrix {apply(w, i-> -i)};
   Imatrix := transpose ((id_(QQ^n))|(-(id_(QQ^n)))|negWvector);
   simplexLP := latticePoints polyhedronFromHData(Imatrix,vMatrix);
   simplexLPLists := apply(simplexLP, i -> flatten entries transpose i);

   ------------------------------------------------------------------------
   -- Step 3. create monomial ideal in R generated by m with rho(m)\geq alpha
   ------------------------------------------------------------------------

   ideal mingens ideal apply(simplexLPLists, i -> R_i)
   )


---------------------------------------------------------------

--private recursive helper: assumes inputs already validated.  Cached on
--R.cache#HodgeIdealWeightedHomogIsolatedCache keyed by (f, alpha, p, w), so
--recursive calls and surrounding loops (e.g. deRhamInterval) reuse work.
hodgeIdealWeightedHomogIsolatedHelper = (f,alpha,p,w) -> (

  R := ring f;
  if not R.cache#?HodgeIdealWeightedHomogIsolatedCache then
      R.cache#HodgeIdealWeightedHomogIsolatedCache = new MutableHashTable;
  tbl := R.cache#HodgeIdealWeightedHomogIsolatedCache;
  key := (f, sub(alpha, QQ), p, w);
  if tbl#?key then return tbl#key;

  hodgeI := ideal(1_R);

  if p == 0 then hodgeI = weightGeqAlpha(R, w, alpha)

      else  (
	  RgeqPplusAlpha := weightGeqAlpha(R,w,p+alpha);--the first summand of Zheng's formula
	  if RgeqPplusAlpha == ideal(1_R) then hodgeI = ideal(1_R)

                 else (

	               hodgeMinusOne := hodgeIdealWeightedHomogIsolatedHelper(f,alpha,p-1,w);
	               hodgeMinusOneGens := flatten entries gens gb hodgeMinusOne;
	               n := numgens R;
	               diffSummand := for i from 0 to n-1 list (
		       for a in hodgeMinusOneGens list f*diff(R_i,a)-(alpha+p-1)*a*diff(R_i,f));--the second summand
	               hodgeI = ideal mingens (RgeqPplusAlpha + ideal(flatten diffSummand) + f*hodgeMinusOne);

                 );
           );

  tbl#key = hodgeI;
  hodgeI
  )


hodgeIdealWeightedHomogIsolated = method();

--ref: Corollary B of Zhang "HODGE FILTRATION AND HODGE IDEALS FOR Q-DIVISORS WITH WEIGHTED HOMOGENEOUS ISOLATED SINGULARITIES"

hodgeIdealWeightedHomogIsolated(RingElement, ZZ, ZZ, List) :=
hodgeIdealWeightedHomogIsolated(RingElement, QQ, ZZ, List) := (f,alpha,p,w) -> (
-- f is a polynomial in polynomial ring R
-- f is weighted homog of weight w with isolated sing
-- w*(exponent f) = 1
-- alpha is a rational number in (0,1]
-- p is a non-negative integer
-- outputs the Hodge ideal I_p(alpha f)

  checkAlpha alpha;
  checkP p;
  R := ring f;
  if #w != numgens R then
      error "expected length of w to equal numgens(ring f)";
  if not all(exponents f, e -> sum apply(#e, i -> e_i * w_i) == 1) then
      error "expected f to be weighted homogeneous of weight 1 with respect to w";

  hodgeIdealWeightedHomogIsolatedHelper(f,alpha,p,w)
  )


---------------------------------------------------------------
---------------------------------------------------------------

hodgeIdealBrieskornPham = method();

hodgeIdealBrieskornPham(List, ZZ, ZZ) :=
hodgeIdealBrieskornPham(List, QQ, ZZ) := (L,alpha,p) -> (
--calculates Hodge ideals for x_1^(b_1)+\cdots+x_n^(b_n) in an internally
--created polynomial ring QQ[x_1..x_n].  Use (Ring, List, ZZ/QQ, ZZ) variant
--to compute in a caller-supplied ring.
--L is the list {b_1,\cdots,b_n}
--the b_i's should all be positive
--p is a non-negative integer

    if #L == 0 then error "expected L to be non-empty";
    if not all(L, b -> instance(b,ZZ) and b >= 1) then
        error "expected entries of L to be positive integers";
    checkAlpha alpha;
    checkP p;
    n := #L;
    x := local x;
    R := QQ[x_1..x_n];
    hodgeIdealBrieskornPham(R, L, alpha, p)
    )

hodgeIdealBrieskornPham(Ring, List, ZZ, ZZ) :=
hodgeIdealBrieskornPham(Ring, List, QQ, ZZ) := (R,L,alpha,p) -> (
--Variant that computes the Brieskorn-Pham Hodge ideal in a caller-supplied
--polynomial ring R = k[x_1,...,x_n].  Useful when the caller already has
--names for the variables and wants the output ideal to live in their ring.

    if #L == 0 then error "expected L to be non-empty";
    if not all(L, b -> instance(b,ZZ) and b >= 1) then
        error "expected entries of L to be positive integers";
    if #L != numgens R then
        error("expected length of L to equal numgens R = ", toString numgens R);
    checkAlpha alpha;
    checkP p;
    n := #L;
    f := sum apply(toList(0..n-1), i -> R_i^(L_i));
    w := apply(L, i -> 1/i);
    hodgeIdealWeightedHomogIsolated(f,alpha,p,w)
    )



---------------------------------------------------------------
---------------------------------------------------------------

higherMultiplierIdeal = method();

higherMultiplierIdeal(RingElement, ZZ, ZZ) :=
higherMultiplierIdeal(RingElement, QQ, ZZ) := (f,alpha,p) -> (
    checkAlpha alpha;
    checkP p;
    FpV := hodgeOnV(f,alpha,p);
    R := ring f;
    RDt := ring (FpV_0);
    Dt := last gens RDt;

    -- coefficient matrix: rows correspond to Dt^p, Dt^(p-1), ..., Dt^0
    DtPowersMat := (coefficients(matrix{FpV}, Variables => {Dt}))_1;

    -- take the top row (Dt^p row), as a 1 x (#cols) matrix
    topRowMat := DtPowersMat^{0};

    -- entries -> R
    topRowEntriesR := apply(flatten entries topRowMat, g -> sub(g, R));

    I := ideal mingens ideal topRowEntriesR;
    I
);



--------------------------------------------------------------
-------------------------------------------------------------

hodgeCheck = method();

hodgeCheck(RingElement,RingElement,ZZ,ZZ) :=
hodgeCheck(RingElement,RingElement,QQ,ZZ) := (f,g,alpha,p) -> (
--checks if g/f^(alpha) is in F_p(S_f*f^(-alpha)) for a rational number alpha>0
--if beta\in (0,1] and k\in ZZ  with alpha = k + beta then
--checks if (g/f^k)*f^(-beta) is in F_p(S_f*f^(-beta))

    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    checkP p;
    betaQQ := if (denominator(sub(alpha,QQ)) == 1) then 1_QQ else alpha-floor(alpha);
    poleOrder := if (denominator(sub(alpha,QQ)) == 1) then floor(alpha)-1 else floor(alpha);

    isMemberFp := true;

    if (p+1-poleOrder)<0 then isMemberFp = false

      else(
	  IpBeta := hodgeIdeal(f,betaQQ,p);
	  gTwist := g*f^(p+1-poleOrder);
	  isMemberFp = isMember(gTwist, IpBeta);
	  );
    isMemberFp
    )


--------------------------------------------------------------
-------------------------------------------------------------

hodgeLevel = method();

hodgeLevel(RingElement,RingElement,ZZ) :=
hodgeLevel(RingElement,RingElement,QQ) := (f,g,alpha) -> (
--finds the minimal p for which g/f^(alpha) is in F_p(S_f*f^(-alpha)) for alpha>0
--If beta\in (0,1] and k\in ZZ  with alpha = k + beta then finds the minimal p
--for which (g/f^k)*f^(-beta) is in F_p(S_f*f^(-beta))

   if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
   pNow := 0;
   isMemberFp := hodgeCheck(f,g,alpha,pNow);

   while not isMemberFp do (
       pNow = pNow+1;
       isMemberFp = hodgeCheck(f,g,alpha,pNow);
       );

   pNow
   )


---------------------------------------------------------------
---------------------------------------------------------------
-- Hodge ideals for the determinant hypersurface.
--
-- For each n, set S_n = (ZZ/32003)[x_(1,1)..x_(n,n)] and f_n = det of
-- the generic n x n matrix.  hodgeIdealDet(n, p) returns the p-th
-- Hodge ideal I_p(f_n) as an ideal in S_n.
---------------------------------------------------------------
---------------------------------------------------------------

--cache for hodgeIdealDet's session-global caches.  Cache
--symbols (DetGenericMatrixRingCache, ILambdaDetCache, HodgeIdealDetCache)
--are protected in MixedHodgeModules.m2 and used as keys on this ring.
detSessionRing = QQ[];

detEnsureCache = sym -> (
    if not detSessionRing.cache#?sym then
	detSessionRing.cache#sym = new MutableHashTable;
    detSessionRing.cache#sym
    )

--Returns Skk = (ZZ/32003)[x_(1,1)..x_(n,n)], cached per n.
--`local x` keeps the user's global symbol x untouched.
detGenericMatrixRing = n -> (detEnsureCache DetGenericMatrixRingCache)#n ??= (
    x := local x;
    (ZZ/32003)[x_(1,1)..x_(n,n)]
    );

detPadZeros = (L, n) -> L | apply(toList(1..n - #L), i -> 0);

--Partitions whose top p parts are equal (used for symbolic powers of (p x p minors)^d).
detSymbolicPowerParts = (p, d, n) -> (
    if d <= 0 then return {};
    partsList := {};
    for a from ceiling(d/(n - p + 1)) to d do (
	tailSum := d - a;
	tails := if n == p then (
	    if tailSum == 0 then {{}} else {}
	    )
	    else select(partitions(tailSum, n - p), mu ->
		#mu <= n - p and (#mu == 0 or mu_0 <= a));
	partsList = partsList | apply(tails, mu ->
	    join(apply(toList(1..p), i -> a), detPadZeros(toList mu, n - p)));
	);
    partsList
    );

--ILambda for partition lam (over ZZ/32003), cached by (lam, n).
ILambdaDet = (lam, n) -> (detEnsureCache ILambdaDetCache)#(lam, n) ??= (
    Skk := detGenericMatrixRing n;
    kk := coefficientRing Skk;
    r := local r;
    R := schurRing(r, n);
    conjlam := toList conjugate(new Partition from lam);
    e := dim r_lam;
    M := genericMatrix(Skk, n, n);
    ideal for i from 0 to e*e - 1 list (
	A := random(kk^n, kk^n);
	B := random(kk^n, kk^n);
	N := A * M * B;
	product for j from 0 to #conjlam - 1 list
	    det(N_{0..conjlam_j - 1}^{0..conjlam_j - 1})
	)
    );

--symbolic power (p x p minors)^(d) over ZZ/32003.
detSymbolicPower = (p, d, n) -> (
    Skk := detGenericMatrixRing n;
    if d <= 0 then return ideal(1_Skk);
    M := genericMatrix(Skk, n, n);
    if d == 1 then return minors(p, M);
    if p == 1 then return (minors(1, M))^d;
    sum apply(detSymbolicPowerParts(p, d, n), lam -> ILambdaDet(lam, n))
    );

hodgeIdealDet = method();
--uses formula from Perlman-Raicu
hodgeIdealDet(ZZ, ZZ) := (n, p) -> (detEnsureCache HodgeIdealDetCache)#(n, p) ??= (
    powers := apply(toList(1..n - 1), q ->
	{q, (n - q)*(p - 1) - binomial(n - q, 2)});
    ideals := apply(powers, g -> detSymbolicPower(g_0, g_1, n));
    if #ideals == 1 then ideals_0 else intersect ideals
    );

