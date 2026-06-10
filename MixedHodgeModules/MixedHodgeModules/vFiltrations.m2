-- Copyright 2026 by Andras Lorincz and Michael Perlman
--
-- V-filtration code.  This file collects the V-filtration computation kernel
-- (hodgeOnV, weightHodgeOnV, monodromyWeightHodgeOnV, HRHCheck, HRHLevel) 
--
--
-- Sections:
--   1. Input validation         checkAlpha, checkP
--   2. Caches                   cachedAnnFs, cachedGlobalBFunction, cachedRhoFp, canonicalAlpha
--   3. V-filtration setup       prepareVfilt
--   4. Elimination primitives   DsToRs, cachedJ0
--   5. truncation helpers       vSlicePoly, truncateBySDeg
--   6. Basis conversion         monomToShiftedQCoeffs, stdQpolysDesc,
--                               fromMRsToBf, convertMStoDtBasisBf,
--                               padDtPowersToP, malgrangeEval, toBfBasis
--   7. Main computations        hodgeOnV (both methods), weightHodgeOnV,
--                               monodromyWeightHodgeOnV
--   8. HRH                      HRHCheck, HRHLevel



---------------------------------------------------------------
---------------------------------------------------------------
--1. Input validation


checkAlpha = alpha -> (
    if sub(alpha,QQ) <= 0 or sub(alpha,QQ) > 1 then
        error "expected alpha to be a rational number in (0,1]";
    )

checkP = p -> (
    if p < 0 then error "expected p to be a non-negative integer";
    )


---------------------------------------------------------------
---------------------------------------------------------------
--2. Caches

cachedAnnFs = f -> (
    R := ring f;
    if not R.cache#?AnnFsCache then R.cache#AnnFsCache = new MutableHashTable;
    tbl := R.cache#AnnFsCache;
    if tbl#?f then tbl#f
    else tbl#f = AnnFs f
    )

cachedGlobalBFunction = f -> (
    R := ring f;
    if not R.cache#?GlobalBCache then R.cache#GlobalBCache = new MutableHashTable;
    tbl := R.cache#GlobalBCache;
    if tbl#?f then tbl#f
    else tbl#f = globalBFunction f
    )

--cache the sorted {root, multiplicity} list from the (p+1)-st generalized
--b-function of f.  Keyed on (f, p).  Roots are coerced to QQ.
cachedRhoFp = (f, p) -> (
    R := ring f;
    if not R.cache#?RhoFpCache then R.cache#RhoFpCache = new MutableHashTable;
    tbl := R.cache#RhoFpCache;
    key := (f, p);
    if tbl#?key then tbl#key
    else (
        bs := generalB({f}, 1_(ring f), Exponent => p+1);
        bsFac := factorBFunction bs;
        tbl#key = sort apply(toList(0..(#bsFac-1)), i ->
            {sub((ring bs)_0 - (bsFac#i#0), QQ), bsFac#i#1})
        )
    )

--internal: given alpha in (0,1] and p, return the upper endpoint of alpha's
--equivalence class under semi-continuity of F_p(V^{alpha}(B_f)) in alpha.  Two values
--of alpha that lie in the same interval (alpha_{i+1}, alpha_i] produce identical
--hodgeOnV(f, alpha, p) output, so they share the same representative.
--Breakpoints {alpha_1, ..., alpha_{k-1}} \cup {1} come from roots of the (p+1)-st
--generalized b-function via -root - p; rhoFp is already cached.
canonicalAlpha = (f, alpha, p) -> (
    alphaQQ := sub(alpha, QQ);
    rhoFp := cachedRhoFp(f, p);
    breakpoints := unique select(apply(rhoFp, i -> -i_0 - p), b -> b > 0 and b <= 1);
    if not member(1_QQ, breakpoints) then breakpoints = append(breakpoints, 1_QQ);
    breakpoints = sort breakpoints;
    canonical := 1_QQ;
    for b in breakpoints do (
        if b >= alphaQQ then (canonical = b; break);
        );
    canonical
    )


---------------------------------------------------------------
---------------------------------------------------------------
--3. V-filtration setup

--shared V-filtration setup used by hodgeOnV, weightHodgeOnV,
--monodromyWeightHodgeOnV, and HRHCheck.
--Returns the common quantities derived from f and p:
--  Ds    : Weyl algebra D[s] containing AnnFs(f)
--  ss    : the s-variable, Ds_(numgens Ds - 1)
--  DsF   : sub(f, Ds)
--  rhoFp : sorted list of {root, multiplicity} pairs from the (p+1)-st
--          generalized b-function (roots coerced to QQ; sorted increasing)
--  Jp    : ideal(gens gb AnnFsf) + ideal(DsF^(p+1))
--  M     : Ds^1 / Jp

prepareVfilt = (f, p) -> (
    AnnFsf := cachedAnnFs f;
    Ds := ring AnnFsf;
    Gp := flatten entries gens gb AnnFsf;
    ss := Ds_(numgens Ds - 1);
    DsF := substitute(f, Ds);
    rhoFp := cachedRhoFp(f, p);
    Jp := ideal(Gp) + ideal(DsF^(p+1));
    M := Ds^1 / Jp;
    (Ds, ss, DsF, rhoFp, Jp, M)
    )


---------------------------------------------------------------
---------------------------------------------------------------
--4. Elimination 

DsToRs = I -> (
--I is an ideal in Ds=R[dR_0..dR_n,s]--in this order
--eliminates dR_0..dR_n to give ideal in Rs
--outputs a list of elements of Grobner basis

    Ds := ring I;
    numGensDs := numgens Ds;
    ss := Ds_(numGensDs-1);
    partialVars := apply(toList(numGensDs//2..numGensDs-2), i -> Ds_i);
    elim := unique flatten entries gens eliminateWA(I,partialVars);--hard(est) step
    Rs := QQ[prepend(ss,apply(toList(0..numGensDs//2-1), i -> Ds_i)), MonomialOrder=>Lex];
    elimRs := apply(elim, i -> substitute(i,Rs));
    elimRsIdeal := ideal(elimRs);

    flatten entries gens gb elimRsIdeal
    )


--cache (Rs, J0) where J0 = Jp \cap R[s], keyed by (f, p).
--The one expensive Weyl-algebra elimination lives here, and the result
--is reused across all alpha by hodgeOnV / weightHodgeOnV.

cachedJ0 = (f, p) -> (
    R := ring f;
    if not R.cache#?J0Cache then R.cache#J0Cache = new MutableHashTable;
    tbl := R.cache#J0Cache;
    key := (f, p);
    if tbl#?key then tbl#key
    else (
        (Ds, ss, DsF, rhoFp, Jp, M) := prepareVfilt(f, p);
        H := DsToRs Jp;
        Rs := ring H_0;
        tbl#key = (Rs, ideal H)
        )
    )


---------------------------------------------------------------
---------------------------------------------------------------
--5. Polynomial / truncation helpers

--build a(s) = \prod (s - lambda)^mult in Rs, where Rs_0 = s
--and roots is a list of {lambda, mult} pairs (lambda in QQ).
--empty product returns 1_Rs.

vSlicePoly = (roots, Rs) -> (
    if #roots == 0 then 1_Rs
    else product apply(roots, i -> (Rs_0 - sub(i_0, Rs))^(i_1))
    )


--keep polys of s-degree <= p (where Rs_0 = s).

truncateBySDeg = (L, Rs, p) -> (
    select(L, g -> degree(Rs_0, g) <= p)
    )


---------------------------------------------------------------
---------------------------------------------------------------
--6. Basis conversion

-- Given polys in Rs truncated deg_s <= p, return the coefficients in the
-- shifted basis { t^{-p}, (t dt) t^{-p}, ..., (t^p dt^p) t^{-p} }
-- i.e. coefficients in the "shifted Q-basis" corresponding to (s+p)...(s+1).
monomToShiftedQCoeffs = (Htrunc, Rs, p) -> (
    HtruncMat := matrix{Htrunc};
    coeffsS := (coefficients(HtruncMat, Variables => {Rs_0}))_1;

    sMinuses := apply(toList(0..(p-1)), i -> Rs_0 + p - i); -- (s+p)...(s+1)
    shiftedQ := matrix{apply(toList(0..p), i ->
        ((-1)^i * (product drop(sMinuses, -p+i)))_(Rs)
    )};

    ch := (coefficients(shiftedQ, Variables => {Rs_0}))_1;
    inverse(ch) * coeffsS
)


--std Q polynomials for Malgrange isomorphism
stdQpolysDesc = (Rs, p) -> (
    -- returns [Q_p, ..., Q_0] where Q_i(s)=(-1)^i (s+i)...(s+1)
    Qlist := apply(toList(0..p), i -> (
        if i==0 then 1_Rs
        else ((-1)^i * product apply(toList(1..i), j -> Rs_0 + j))_(Rs)
    ));
    reverse Qlist
);


fromMRsToBf = (MRsList, f, Rs, p) -> (
    -- Input: list in Rs=QQ[s,x...] deg_s<=p coming from DsToRs truncation
    --maps M lattice to Bf lattice
    -- changes s representation of V-filt on M to that of B_f
    -- Output: list in Frac(R)[s] giving elements of F_p V^alpha(B_f)=R_f[s] f^s.

    dtCoeffs := monomToShiftedQCoeffs(MRsList, Rs, p);   -- (p+1)x(#gens) matrix

    -- Build the standard Q basis polynomials [Q_p,...,Q_0] as a row matrix
    Qdesc := stdQpolysDesc(Rs, p);
    Qrow  := matrix{Qdesc};  -- 1 x (p+1)

    -- Multiply: (1 x (p+1)) * ((p+1) x k) = (1 x k) polynomials in Rs
    polys := Qrow * dtCoeffs;

    --coerce to frac coefficient ring:
    Frs := frac Rs;
    fpow := (sub(f, Frs))^p;
    flatten entries (sub(polys, Frs) / fpow)
)


--this is for converting directly to dt basis
--helper function
--changing basis as in Blanco steps 17 - 19
convertMStoDtBasisBf = (Htrunc, f, Rs, p) -> (
--Htrunc is a list of elements of Rs of degree <=p in s
--f is an element of R
--p is a non-negative integer
--convention: s= -dt*t
--output: a list of elements of R[dt]

    HtruncMat := matrix{Htrunc};
    coeffsS := (coefficients(HtruncMat, Variables => {Rs_0}))_1;

    sMinuses := apply(toList(0..(p-1)), i -> Rs_0 + p - i);--{(s+p),(s+p-1),..,(s+1)}

    tDtBasisInS := matrix{apply(toList(0..p), i ->
        ((-1)^i * (product drop(sMinuses, -p+i)))_(Rs)
    )};

    coeffsdtmtmInv := inverse (coefficients(tDtBasisInS, Variables => {Rs_0}))_1;
    dtCoeffs := coeffsdtmtmInv * coeffsS;
    --we use the identity t^m*dt^m= (-1)^m * (s+m)*(s+m-1)*..*(s+1)
    --since we divide by f^p below (see fRPowers below), we use the identity
    --(t^m*dt^m)*t^(-p)= (-1)^m * t^(-p)* (s+p)*(s+p-1)*..*(s+p-m+1)

    dt := local dt;
    gensR := drop(gens Rs,1);--get rid of s
    RDt := QQ[append(gensR, dt)];
    fRDt := sub(f, RDt);

    fRPowers := matrix{apply(toList(0..p), i ->
        ((fRDt^i)/(fRDt^p)) * (dt^i)--multiply by appropriate powers of f (Blanco 18)
    )};

    dtCoeffsDt := sub(dtCoeffs, RDt);
    flatten entries sub(fRPowers * dtCoeffsDt, RDt)
);


--helper function for malgrangeEval to pad zeroes of the correct degree in DtPowers, if necessary
padDtPowersToP = (DtPowers, p) -> (
    RDt := ring DtPowers;

    -- current row degrees as integers (assumes they look like {k})
    currDegs := apply(degrees target DtPowers, d -> d#0);
    b := max currDegs;
    a := min currDegs;

    -- degrees missing above (top padding): p, p-1, ..., b+1
    topMissing := if p > b then toList(reverse (b+1..p)) else {};

    -- degrees missing below (bottom padding): a-1, a-2, ..., 0
    botMissing := if a > 0 then toList(reverse (0..a-1)) else {};

    -- build graded zero blocks with the right target shifts
    PadTopT := if #topMissing > 0
        then directSum apply(topMissing, k -> RDt^{1:-{k}})
        else null;

    PadBotT := if #botMissing > 0
        then directSum apply(botMissing, k -> RDt^{1:-{k}})
        else null;

    Ztop := if PadTopT === null then null else map(PadTopT, source DtPowers, 0_RDt);
    Zbot := if PadBotT === null then null else map(PadBotT, source DtPowers, 0_RDt);

    -- concatenate vertically in the correct order
    if Ztop === null and Zbot === null then DtPowers
    else if Ztop =!= null and Zbot === null then Ztop || DtPowers
    else if Ztop === null and Zbot =!= null then DtPowers || Zbot
    else Ztop || DtPowers || Zbot
);


--helper function to turn V-filtration in terms of dt to V-filtration in terms of s
--then substitute a rational number for s
malgrangeEval = (L,f,alpha,p) -> (
--p is a non-negative integer
--L is a list of elements of degree at most p in RDt (such as basis of F_pV^{alpha}(B_f))
--alpha is in (0,1]
--output is a list of elements of R

    R := ring f;
    alphaQQ := substitute(alpha,QQ);
    RDt := ring L_0;
    fRDt := sub(f,RDt);
    dt := last gens RDt;

    DtPowers :=  (coefficients(matrix{L}, Variables=> {dt}))_1;----higher power in the top row
    DtPowersPadded := padDtPowersToP(DtPowers, p);

    alphaPlusJs:=apply(toList(0..(p-1)), j-> alphaQQ+j);
    QAlphas := append(reverse apply(toList(0..(p-1)), i-> (product drop(alphaPlusJs, -p+1+i))),1);
    -- QAlphas = {(alpha)_p, (alpha)_(p-1), ..., (alpha)_0}
    -- where (alpha)_i = alpha*(alpha+1)*...*(alpha+i-1) is the rising factorial.

    fAlphaVector := matrix{apply(toList(0..p), i-> ((fRDt)^i)*(QAlphas_i))};
    newGens := fAlphaVector*DtPowersPadded;
    newGensR := sub(newGens, R);

    newGensR
    )


--change basis from M-coordinates in Rs to either sBasis (R_f[s] f^s)
--or dtBasis (R[dt]).  Dispatches to fromMRsToBf or convertMStoDtBasisBf.

toBfBasis = (L, f, Rs, p, useBasis) -> (
    if useBasis === sBasis then fromMRsToBf(L, f, Rs, p)
    else if useBasis === dtBasis then convertMStoDtBasisBf(L, f, Rs, p)
    else error "invalid UseBasis"
    )


---------------------------------------------------------------
---------------------------------------------------------------
--7. Main V-filtration computations

hodgeOnV = method(Options => {UseBasis => dtBasis});

--use "UseBasis => sBasis" to get s-basis (will have denominators) in R_f[s]f^s

hodgeOnV(RingElement,ZZ,ZZ) :=
hodgeOnV(RingElement,QQ,ZZ) := options -> (f,alpha,p) -> (
-- f is a polynomial in polynomial ring R
-- alpha is a rational number in (0,1]
-- p is a non-negative integer
-- outputs an R-basis for F_p(V^{alpha}i_{+}R) in terms of dt

--this differs from Blanco in two ways:
--(1) uses (p+1)-st generalized b-function instead of shifts of b-function,
--(2) only calculates for one alpha instead of all

--by the identity (Jp : a(s)) \cap R[s] = (Jp \cap R[s]) : a(s),
--the elimination W_alpha \cap R[s] reduces to a colon of J0 = Jp \cap R[s]
--by a(s) = \prod_{lam <= -alpha-p} (s - lam)^mult; J0 is cached on (f, p)

    checkAlpha alpha;
    checkP p;
    if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    R := ring f;
    if not R.cache#?HodgeOnVCache then R.cache#HodgeOnVCache = new MutableHashTable;
    tbl := R.cache#HodgeOnVCache;
    cAlpha := canonicalAlpha(f, alpha, p);
    key := (f, cAlpha, p, options.UseBasis);
    if tbl#?key then return tbl#key;
    -- if hodgeOnV(f, p) HashTable was cached, reuse its entry for cAlpha
    hashKey := (f, p, options.UseBasis);
    if tbl#?hashKey and (tbl#hashKey)#?cAlpha then (
        result := (tbl#hashKey)#cAlpha;
        tbl#key = result;
        return result;
        );

    ------------------------------------------------------------------------
    -- Steps 1-3. rhoFp via cachedRhoFp; J0 = Jp \cap R[s] via cachedJ0
    ------------------------------------------------------------------------

    negAlphaQQ := -substitute(alpha,QQ);
    rhoFp := cachedRhoFp(f, p);
    rhoFpAlpha := select(rhoFp, i -> ((i_0) <= negAlphaQQ-p));
    (Rs, J0) := cachedJ0(f, p);

    ------------------------------------------------------------------------
    -- Step 4. W_alpha \cap R[s] = J0 : a(s)  for a(s) = \prod (s - lam)^mult
    ------------------------------------------------------------------------

    aPoly := vSlicePoly(rhoFpAlpha, Rs);
    Halpha := flatten entries gens gb (J0 : ideal aPoly);

    ------------------------------------------------------------------------
    -- Step 5. Truncate by s-degree p and change basis (Blanco 14 - 21)
    ------------------------------------------------------------------------

    VFiltpM := truncateBySDeg(Halpha, Rs, p);
    --this is the R-basis of F_pV^{-alpha-p}M
    VFiltpBf := toBfBasis(VFiltpM, f, Rs, p, options.UseBasis);

    tbl#key = VFiltpBf;
    use R;--restore caller's ring after intermediate Rs / RDt creation
    VFiltpBf
)


---------------------------------------------------------------

hodgeOnV(RingElement,ZZ) := options -> (f,p) -> (
-- f is a polynomial in polynomial ring R
-- p is a non-negative integer
-- outputs a hash table with keys alpha,
-- where alpha\in (0,1] and -alpha-p is a root of
-- the (p+1)-st generalized b-function of f, and values
-- an R-basis for F_p(V^{alpha}i_{+}R)

--This is closer to Blanco because it does all alpha\in (0,1]
--but still uses (p+1)-st generalized b-function

--J0 = Jp \cap R[s] is cached on (f, p); for each alpha-root we extend the
--polynomial a(s) by one factor (s - lam)^mult and take J0 : a(s) in R[s]

    checkP p;
    if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    R := ring f;
    if not R.cache#?HodgeOnVCache then R.cache#HodgeOnVCache = new MutableHashTable;
    tbl := R.cache#HodgeOnVCache;
    key := (f, p, options.UseBasis);
    if tbl#?key then return tbl#key;

    ------------------------------------------------------------------------
    -- Steps 1-3. rhoFp via cachedRhoFp; J0 = Jp \cap R[s] via cachedJ0
    ------------------------------------------------------------------------

    rhoFp := cachedRhoFp(f, p);
    (Rs, J0) := cachedJ0(f, p);

    ------------------------------------------------------------------------
    -- Step 4. Partition roots: lam < -p-1 always in a(s);
    -- lam in [-p-1, -p) added one by one, yielding F_pV^alpha for each alpha
    ------------------------------------------------------------------------

    rhoLess := select(rhoFp, i -> i_0 < -p-1);
    rhoAlpha := select(rhoFp, i -> i_0 >= -p-1 and i_0 < -p);

    ------------------------------------------------------------------------
    -- Step 5. For each alpha-root, extend a(s) by (s - lam)^mult and compute
    -- Halpha = J0 : a(s); then truncate by s-degree p (Blanco 14 - 21)
    ------------------------------------------------------------------------

    aPoly := vSlicePoly(rhoLess, Rs);
    pairsAlphaB := for i in rhoAlpha list (
        aPoly = aPoly * (Rs_0 - sub(i_0, Rs))^(i_1);
        Halpha := flatten entries gens gb (J0 : ideal aPoly);
        {i_0, truncateBySDeg(Halpha, Rs, p)}
        );

    ------------------------------------------------------------------------
    -- Step 6. Change basis
    ------------------------------------------------------------------------

    VFiltpPairs := apply(pairsAlphaB, iH -> {-iH_0-p, toBfBasis(iH_1, f, Rs, p, options.UseBasis)});

    result := new HashTable from VFiltpPairs;
    tbl#key = result;
    use R;--restore caller's ring after intermediate Rs / RDt creation
    result
)


---------------------------------------------------------------

weightHodgeOnV = method(Options => {UseBasis => dtBasis})

--use "UseBasis => sBasis" to get s-basis (will have denominators)

weightHodgeOnV(RingElement, ZZ, ZZ, ZZ) :=
weightHodgeOnV(RingElement, QQ, ZZ, ZZ) := options -> (f,alpha,p,m) -> (
--alpha is a rational number in (0,1]
--p and m are non-negative integers
--calculates K^{\alpha}_m, the lift of ker(N^m)
--see Olano weighted Hodge ideals, proof of Theorem A
--
--K \cap R[s] = (J0 : a_<=(s)) \cap (J0 : a_<(s) * Nop^m)
--where a_<(s)  = \prod_{lam < -alpha-p} (s-lam)^mult,
--      a_<=(s) = a_<(s) * (s - lam_alpha)^mult_alpha   (boundary factor, if any),
--      Nop = s - lam_alpha = s + p + alpha

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
    -- Steps 1-3. rhoFp via cachedRhoFp; J0 = Jp \cap R[s] via cachedJ0
    ------------------------------------------------------------------------

    negAlphaQQ := -sub(alpha,QQ);
    rhoFp := cachedRhoFp(f, p);
    rhoLess := select(rhoFp, i -> i_0 < negAlphaQQ - p);
    rhoBoundary := select(rhoFp, i -> i_0 == negAlphaQQ - p);

    (Rs, J0) := cachedJ0(f, p);

    ------------------------------------------------------------------------
    -- Step 4. Build a_<(s), a_<=(s), and Nop = s - lam_alpha = s + p + alpha
    ------------------------------------------------------------------------

    aLessPoly := vSlicePoly(rhoLess, Rs);
    aLessEqPoly := aLessPoly * vSlicePoly(rhoBoundary, Rs);
    Nop := Rs_0 + p + alpha;

    ------------------------------------------------------------------------
    -- Step 5. K \cap R[s] = (J0 : a_<=) \cap (J0 : a_< * Nop^m)
    ------------------------------------------------------------------------

    KFinal := intersect(J0 : ideal aLessEqPoly, J0 : ideal (aLessPoly * Nop^m));
    Halpha := flatten entries gens gb KFinal;

    ------------------------------------------------------------------------
    -- Step 6. Truncate by s-degree p and change basis (Blanco 14 - 21)
    -- this is W_m F_p V^{-alpha-p}(M)
    ------------------------------------------------------------------------

    WFVpM := truncateBySDeg(Halpha, Rs, p);
    WFVpBf := toBfBasis(WFVpM, f, Rs, p, options.UseBasis);

    tbl#key = WFVpBf;
    use R;--restore caller's ring after intermediate Rs / RDt creation
    WFVpBf
    )


---------------------------------------------------------------

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

    use ring f;--restore caller's ring after intermediate Rs / RDt / Ds creation
    WFVpBf
    )


---------------------------------------------------------------
---------------------------------------------------------------
--8. HRH checks


HRHCheck = method();

HRHCheck(RingElement, ZZ) := (f,p) -> (
--checks if F_p(Gr^0_V(B_f)) = 0.  By [DOR, Thm H] and the can/var reduction,
--this holds iff N = s+p+1 annihilates F_p Gr^1_V(B_f)

    R := ring f;
    rhoFp := cachedRhoFp(f, p);
    (Rs, J0) := cachedJ0(f, p);

    -- roots <  -p-1  ->  V^{>1};   roots <= -p-1  ->  V^1
    aLess := vSlicePoly(select(rhoFp, i -> i_0 <  -p-1), Rs);
    aLeq  := vSlicePoly(select(rhoFp, i -> i_0 <= -p-1), Rs);

    Igt1 := J0 : ideal aLess;                                    -- lift of V^{>1}
    kerN := Igt1 : ideal(Rs_0 + p + 1);                          -- lift of ker(N) on Gr^1_V
    FpV1 := truncateBySDeg(                                       -- R-gens of F_p V^1
        flatten entries gens gb (J0 : ideal aLeq), Rs, p);

    use R;--restore caller's ring after intermediate Rs creation
    if #FpV1 == 0 then return true;                              -- F_p Gr^1_V = 0
    isSubset(ideal FpV1, kerN)                                   -- F_p V^1 subset ker N
    )


---------------------------------------------------------------

HRHLevel = method();

HRHLevel(RingElement) := f -> (

    n := numgens ring f;
    p := 0;
    HRHtrue := HRHCheck(f,p);

    while HRHtrue and p<(n-1)/2 do (
	p = p+1;
	HRHtrue = HRHCheck(f,p);
	);

    lev := 0;
    if HRHtrue then lev = "rational homology manifold"
      else lev = p-1;

    lev
    )


---------------------------------------------------------------
---------------------------------------------------------------
--cache management


--clearMHMCache R: reset all MixedHodgeModules caches stored on R, plus
--the session-global determinantal caches.  Useful for reclaiming memory
--between large sessions.

clearMHMCache = method()

clearMHMCache(Ring) := R -> (
    for sym in {AnnFsCache, GlobalBCache, RhoFpCache, J0Cache,
		HodgeOnVCache, WeightHodgeOnVCache, HodgeIdealCache,
		WeightedHodgeIdealCache, HodgeIdealWeightedHomogIsolatedCache,
		WeylAlgebraCache, PolyAnnCache} do
	if R.cache#?sym then remove(R.cache, sym);
    DetGenericMatrixRingCache = new MutableHashTable;
    ILambdaDetCache = new MutableHashTable;
    HodgeIdealDetCache = new MutableHashTable;
    )
