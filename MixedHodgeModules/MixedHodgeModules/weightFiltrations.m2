--Copyright 2026 by Andras Lorincz and Michael Perlman



---------------------------------------------------------------
---------------------------------------------------------------
--weight filtration code

weightHodgeOnV = method(Options => {UseBasis => dtBasis})

--use "UseBasis => sBasis" to get s-basis (will have denominators)

weightHodgeOnV(RingElement, ZZ, ZZ, ZZ) :=
weightHodgeOnV(RingElement, QQ, ZZ, ZZ) := options -> (f,alpha,p,m) -> (
--alpha is a rational number in (0,1]
--p and m and non-negative integers
--calculates K^{\alpha}_m, the lift of ker(N^m)
--see Olano weighted Hodge ideals, proof of Theorem A

     checkAlpha alpha;
     checkP p;
     if m < 0 then error "expected weight m to be a non-negative integer";
     if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    R := ring f;
    if not R.cache#?WeightHodgeOnVCache then R.cache#WeightHodgeOnVCache = new MutableHashTable;
    tbl := R.cache#WeightHodgeOnVCache;
    key := (f, sub(alpha,QQ), p, m, options.UseBasis);
    if tbl#?key then return tbl#key;

    ------------------------------------------------------------------------
    -- Steps 1-3. Shared V-filtration setup: Ds, ss, DsF, rhoFp, Jp, M
    ------------------------------------------------------------------------

    negAlphaQQ := -sub(alpha,QQ);
    (Ds, ss, DsF, rhoFp, Jp, M) := prepareVfilt(f, p);

    ------------------------------------------------------------------------
    -- Step 4. Compute kernels
    ------------------------------------------------------------------------

    sLamMapsLess := {};
    sLamMapsAlpha := {};

    for i in rhoFp do (
        lam := sub(i_0, Ds);
        mult := i_1;
        a := (ss - lam)^mult;
	if i_0 < negAlphaQQ-p then sLamMapsLess = append(sLamMapsLess, a)
	else if i_0 == negAlphaQQ-p then sLamMapsAlpha = append(sLamMapsAlpha, a);
	);
    -- entries store (ss - lam)^mult; kernels via colon ideal

    Galpha := ideal flatten apply( sLamMapsLess, a -> flatten entries gens (Jp : ideal a));--hard step
    Klams := ideal flatten apply(sLamMapsAlpha, a -> flatten entries gens (Jp : ideal a));--hard step

    preVgAlpha := Galpha;--this corresponds to V^(>alpha)
    preVAlpha :=  (Galpha + Klams);--this corresponds to V^alpha

    ------------------------------------------------------------------------
    -- Step 5. Find weighted piece in preVAlpha
    ------------------------------------------------------------------------
    
    preVgAlpha = sub(preVgAlpha, Ds);

    K := ideal(1_(Ds));
    
    if preVgAlpha == ideal(0_(Ds)) then K = intersect(preVAlpha, (Jp : ideal((ss+p+alpha)^m)))--handles case when sLamMapsLess empty
       else (
	     K = intersect(preVAlpha, (preVgAlpha : ideal((ss+p+alpha)^m)));--find those in preVAlpha that are sent to preVgAlpha
                                                       
            );
    
    ------------------------------------------------------------------------
    -- Step 6. Eliminate, truncate, and change coordinates (Blanco 14 - 21)
    --this is faster than the corresponding step in hodgeOnV, as the modules are smaller
    ------------------------------------------------------------------------

    Halpha := DsToRs(K);
    Rs := ring Halpha_0;
    WFVpM := select(Halpha, g -> degree(Rs_0, g) <= p);--truncate
    -- this is W_mF_pV^{-alpha-p}(M)
    
    WFVpBf :=  WFVpM;
  
    if options.UseBasis == sBasis then WFVpBf = fromMRsToBf(WFVpBf, f, Rs,p);--twist to get B_f
    if options.UseBasis == dtBasis then WFVpBf = convertMStoDtBasisBf(WFVpBf, f, Rs, p);

    tbl#key = WFVpBf;
    WFVpBf
    )



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

--the following could be added as a new strategy for IHmodule in BernsteinSato
IHmoduleAdjoint = method();

IHmoduleAdjoint(RingElement) := f -> (
--gives presentation of IH using F_0(IC_f)
--useful for calculating b-functions on singular ambient varieties as in the recent
--work of Dirks (and my ongoing collaboration with Dirks).
--Different from current CompleteIntersection strategy, which uses Jacobian ideal

    R := ring f;
    gensIH := flatten entries gens adjointIdeal(f);
    r := #gensIH;
    D := makeWeylAlgebra R;
    fD :=sub(f,D);
    gensIHD := apply(gensIH, i-> sub(i,D));
    I := rationalFunctionAnnihilator(fD)+ideal(fD);--the submodule of R_f/R generated by [1/f]
    g := map(D^1/I, D^r, matrix{gensIHD});
    K := kernel g;
    D^r/K
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
  isMember((((s + alpha + k1)^w) * newg * newf^(k1))_(ring If),J)
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
    w := -1;
    Roots := bFunctionRoots bfs;
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i + alpha,QQ)) == 1);
    if (#AlphaInts == 0) or (AlphaInts_0 < -alpha) then w = 0 else (
	if k1 == 0 then (----lower bound in equation 4 (by Corollary 6.1)
	    B := factorBFunction bfs;
	    positionAlphaInts:=select(toList(0..(#B - 1)), i-> denominator( substitute(B#i#0 - (ring bfs)_0 - alpha,QQ)) == 1);
	    mults := apply(positionAlphaInts, i-> B#i#1);
	    w = max mults-1;
	   -- if ((alpha > 1) and (denominator(sub(alpha,QQ)) == 1) and (sub(g,ring f) == 1_(ring f))) then w = w+1);----Prop 1.6
	   );
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


----------------------------------------------------------------
----------------------------------------------------------------

monodromyWeightHodgeOnV = method(Options => {UseBasis => dtBasis})

--use "UseBasis => sBasis" to get s-basis (will have denominators)

monodromyWeightHodgeOnV(RingElement, ZZ, ZZ, ZZ) :=
monodromyWeightHodgeOnV(RingElement, QQ, ZZ, ZZ) := options -> (f,alpha,p,ell) -> (
--alpha is a rational number in (0,1]
--p is a non-negative integer
--ell is the monodromy weight index, centered at 0
--calculates W(N)_ell F_p Gr^alpha_V(B_f), lifted to F_p V^alpha(B_f)

    checkAlpha alpha;
    checkP p;
    if ell < 0 then error "expected input weight to be a non-negative integer";
    if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    ------------------------------------------------------------------------
    -- Steps 1-2. Shared V-filtration setup: Ds, ss, DsF, rhoFp, Jp, M
    ------------------------------------------------------------------------

    negAlphaQQ := -sub(alpha,QQ);
    (Ds, ss, DsF, rhoFp, Jp, M) := prepareVfilt(f, p);

    ------------------------------------------------------------------------
    -- Step 3. Compute eigenspace kernels for V^>alpha and V^alpha
    ------------------------------------------------------------------------

    sLamMapsLess := {};
    sLamMapsAlpha := {};

    for i in rhoFp do (
        lam := sub(i_0, Ds);
        mult := i_1;
        a := (ss - lam)^mult;

        if i_0 < negAlphaQQ-p then (
            sLamMapsLess = append(sLamMapsLess, a)
        )
        else if i_0 == negAlphaQQ-p then (
            sLamMapsAlpha = append(sLamMapsAlpha, a)
        );
    );
    -- entries store (ss - lam)^mult; kernels via colon ideal

    Galpha := ideal flatten apply(sLamMapsLess, a -> flatten entries gens (Jp : ideal a));
    Klams := ideal flatten apply(sLamMapsAlpha, a -> flatten entries gens (Jp : ideal a));

    preVgAlpha := Galpha;              -- corresponds to V^(>alpha)
    preVAlpha := Galpha + Klams;       -- corresponds to V^alpha

    ------------------------------------------------------------------------
    -- Step 4. Monodromy operator N = s + alpha
    -- In the shifted M-coordinate, this is ss + p + alpha
    ------------------------------------------------------------------------

    Nop := ss + p + alpha;

    ------------------------------------------------------------------------
    -- Step 5. Nilpotence index at alpha
    ------------------------------------------------------------------------

    nilIndex := 0;

    for i in rhoFp do (
        if i_0 == negAlphaQQ-p then nilIndex = i_1;
    );

    ------------------------------------------------------------------------
    -- Step 6. Helper functions for lifted kernels and images
    --
    -- kerLift(r) = lift of ker(N^r) inside V^alpha/V^>alpha
    -- imLift(b)  = lift of im(N^b) inside V^alpha/V^>alpha
    ------------------------------------------------------------------------

    preVgAlpha = sub(preVgAlpha, Ds);
    preVAlpha = sub(preVAlpha, Ds);

    kerLift := r -> (
        if r <= 0 then ideal(0_Ds)
        else intersect(preVAlpha, ((preVgAlpha + Jp) : ideal(Nop^r)))
    );

    imLift := b -> (
        if b == 0 then preVAlpha
        else (
            ideal flatten apply(flatten entries gens preVAlpha, g -> Nop^b*g)
        ) + preVgAlpha + Jp
    );

    ------------------------------------------------------------------------
    -- Step 7. Compute W(N)_ell, centered at 0:
    --
    -- W(N)_ell = sum_{b >= 0, ell+b >= 0}
    --            ker(N^(ell+b+1)) cap im(N^b)
    ------------------------------------------------------------------------

    Wlift := preVgAlpha + Jp;

    lowerB := max(0,-ell);
    upperB := nilIndex;

    if lowerB <= upperB then (
        for b from lowerB to upperB do (
            r := ell + b + 1;
            summand := intersect(kerLift(r), imLift(b));
            Wlift = Wlift + summand;
        );
    );

    K := intersect(preVAlpha, Wlift);

    ------------------------------------------------------------------------
    -- Step 8. Eliminate, truncate, and change coordinates
    ------------------------------------------------------------------------

    Halpha := DsToRs(K);
    Rs := ring Halpha_0;
    WFVpM := select(Halpha, g -> degree(Rs_0, g) <= p);
    -- this is W(N)_ell F_p Gr^alpha_V in the M-coordinate, lifted to V^alpha

    WFVpBf := WFVpM;

    if options.UseBasis == sBasis then WFVpBf = fromMRsToBf(WFVpBf, f, Rs,p);
    if options.UseBasis == dtBasis then WFVpBf = convertMStoDtBasisBf(WFVpBf, f, Rs, p);

    WFVpBf
    )



end

