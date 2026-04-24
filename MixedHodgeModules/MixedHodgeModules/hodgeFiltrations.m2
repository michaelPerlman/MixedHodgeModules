-- Copyright 2026 by Andras Lorincz and Michael Perlman



---------------------------------------------------------------
---------------------------------------------------------------
--helper functions


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


---------------------------------------------------------------

-- Given polys in Rs truncated deg_s <= p, return the coefficients in the
-- shifted basis { t^{-p}, (t dt) t^{-p}, ..., (t^p dt^p) t^{-p} }
-- i.e. coefficients in the "shifted Q-basis" corresponding to (s+p)...(s+1).
-- need this to go from Mp to Bf
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

---------------------------------------------------------------

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

---------------------------------------------------------------

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

---------------------------------------------------------------

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


---------------------------------------------------------------


--this will be standard change of basis from s-basis to dt-basis given by Malgrange isom
-*
sBasisToDtBasisStandard = (sPolyList, f, Rs, p) -> (
    -- Input: list in Frac(R)[s] deg_s<=p that already represents elements of B_f
    -- Output: list in Frac(R)[dt] representing the same elements via Malgrange.

    sPolyMat := matrix{sPolyList};
    monomCoeffs := (coefficients(sPolyMat, Variables => {Rs_0}))_1;

    -- express in standard Q-basis
    Qdesc := stdQpolysDesc(Rs, p);
    QpolysDesc := matrix{Qdesc};                 -- [Q_p,...,Q_0]
    ch := (coefficients(QpolysDesc, Variables => {Rs_0}))_1;
    qCoeffs := (inverse ch) * monomCoeffs;       -- coefficients of Q_p,...,Q_0

    dt := local dt;
    gensR := drop(gens Rs, 1);
    RDt := (coefficientRing Rs)(monoid[append(gensR, dt)]); -- keep same coeff ring
    fRDt := sub(f, RDt);

    -- Q_p maps to f^p dt^p, ..., Q_0 maps to 1
    dtMonomsDesc := matrix{apply(toList(0..p), k -> (fRDt^(p-k)) * (dt^(p-k)))};

    flatten entries (dtMonomsDesc * sub(qCoeffs, RDt))
)
*-

---------------------------------------------------------------
---------------------------------------------------------------


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

    if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    ------------------------------------------------------------------------
    -- Step 1. Basic setup and annihilator data for f^s
    ------------------------------------------------------------------------

    negAlphaQQ := -substitute(alpha,QQ);
    AnnFsf := AnnFs f;
    Ds := ring AnnFsf;
    G := gens gb AnnFsf;
    numGensDs := numgens Ds;
    ss := Ds_(numGensDs-1);
    DsF := substitute(f,Ds);

    ------------------------------------------------------------------------
    -- Step 2. Compute the (p+1)-st generalized b-function and
    -- keep only the roots relevant for V^{alpha}
    ------------------------------------------------------------------------

    bs := generalB({f},1_(ring f), Exponent => p+1);
    bsFac := factorBFunction bs;
    rhoFp := apply(toList(0..(#bsFac-1)), i -> {sub((ring bs)_0-(bsFac#i#0),QQ),bsFac#i#1});--{(root,mult)}
     --the roots are rational numbers
    rhoFpAlpha := select(rhoFp, i -> ((i_0) <= negAlphaQQ-p));

    ------------------------------------------------------------------------
    -- Step 3. Build J_p and the module M = D/J_p
    --(Blanco step 5 (with a different shift on Jp, consistent with b^(p+1)_f))
    ------------------------------------------------------------------------

    Gp := flatten entries G;
    Jp := ideal(Gp) + ideal(DsF^(p+1));
    M := Ds^1 / Jp;
   
    ------------------------------------------------------------------------
    -- Step 4. Compute kernels and form W_alpha (Blanco 6 - 11)
    ------------------------------------------------------------------------

    
     sLamMaps := apply(rhoFpAlpha, i -> (
        lam := sub(i_0, Ds);
        mult := i_1;
        a := (ss - lam)^mult;
        map(M, M, matrix{{a}})
    ));
   
    
    Klams := apply(sLamMaps, i -> flatten entries gens kernel i);--hard step
    Walpha := trim ideal(flatten Klams);

    ------------------------------------------------------------------------
    -- Step 5. Eliminate, truncate, and (possibly) change basis (Blanco 14 - 21)
    ------------------------------------------------------------------------

    Halpha := DsToRs(Walpha);
    Rs := ring Halpha_0;
    VFiltpM := select(Halpha, g -> degree(Rs_0, g) <= p);--truncate
    --this is the R-basis of F_pV^{-alpha-p}M
    VFiltpBf := VFiltpM;

    if options.UseBasis == sBasis then VFiltpBf = fromMRsToBf(VFiltpBf, f, Rs, p);--twist to get Bf
    
    if options.UseBasis == dtBasis then VFiltpBf = convertMStoDtBasisBf(VFiltpBf, f, Rs, p);
   
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

    if (options.UseBasis != dtBasis) and  (options.UseBasis != sBasis) then error "invalid UseBasis";

    ------------------------------------------------------------------------
    -- Step 1. Basic setup and annihilator data for f^s
    ------------------------------------------------------------------------
    
    AnnFsf := AnnFs f;
    Ds := ring AnnFsf;
    G := gens gb AnnFsf;
    numGensDs := numgens Ds;
    ss := Ds_(numGensDs-1);
    DsF := substitute(f,Ds);

    ------------------------------------------------------------------------
    -- Step 2. Compute the (p+1)-st generalized b-function 
    ------------------------------------------------------------------------

    bs := generalB({f},1_(ring f), Exponent => p+1);
    bsFac := factorBFunction bs;
    rhoFp := sort apply(toList(0..(#bsFac-1)), i -> {sub((ring bs)_0-(bsFac#i#0),QQ),bsFac#i#1});--{(root,mult)}
    --important: increasing order
    --the roots are rational numbers
    
    ------------------------------------------------------------------------
    -- Step 3. Build J_p and the module M = D/J_p
    --(Blanco step 5 (with a different shift on Jp, consistent with b^(p+1)_f))
    ------------------------------------------------------------------------

    Gp := flatten entries G;
    Jp := ideal(Gp) + ideal(DsF^(p+1));
    M := Ds^1 / Jp;

    ------------------------------------------------------------------------
    -- Step 4. Compute kernels (Blanco 6 - 11)
    ------------------------------------------------------------------------

    sLamMapsLessNegOne := {};--will be G_{alpha'} in Blanco
    sLamMapsAlpha := {};--will be eigenspaces for alpha\in (0,1]

    for i in rhoFp do (
        lam := sub(i_0, Ds);
        mult := i_1;
        a := (ss - lam)^mult;
	if i_0 < -p-1 then sLamMapsLessNegOne = append(sLamMapsLessNegOne,  {i_0,map(M, M, matrix{{a}})})
	else if i_0 < -p then sLamMapsAlpha = append(sLamMapsAlpha, {i_0,map(M, M, matrix{{a}})});
	);

    Galpha := ideal flatten apply( sLamMapsLessNegOne , i -> flatten entries gens kernel (i_1));--hard step
    Klams := apply(sLamMapsAlpha, i -> {i_0, flatten entries gens kernel (i_1)});--hard step
   
    ------------------------------------------------------------------------
    -- Step 5. Eliminate differential variables, truncate, and (possibly) change basis (Blanco 14 - 21)
    ------------------------------------------------------------------------

    pairsAlphaB := for i in Klams list (
	 Galpha = Galpha + ideal(i_1);
	 Halpha := DsToRs(Galpha);
	 {i_0,Halpha});

    VFiltpPairs := apply(pairsAlphaB, iH -> {-iH_0-p, select(iH_1, g -> degree( (ring g)_0, g) <= p)});

    if options.UseBasis == sBasis then VFiltpPairs = apply(VFiltpPairs, iV -> {iV_0, fromMRsToBf(iV_1, f, ring (iV_1)_0, p)});

    if options.UseBasis == dtBasis then VFiltpPairs = apply(VFiltpPairs, iV -> {iV_0, convertMStoDtBasisBf(iV_1, f, ring (iV_1)_0, p)});

    new HashTable from VFiltpPairs
    
)



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

    R := ring f;
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

  R := ring f;
  hodgeI := ideal(1_R);
  
  if p == 0 then hodgeI = weightGeqAlpha(R, w, alpha)
  
      else  (
	  RgeqPplusAlpha := weightGeqAlpha(R,w,p+alpha);--the first summand of Zheng's formula
	  if RgeqPplusAlpha == ideal(1_R) then hodgeI = ideal(1_R)
	  
                 else (
		       
	               hodgeMinusOne := hodgeIdealWeightedHomogIsolated(f,alpha,p-1,w);
	               hodgeMinusOneGens := flatten entries gens gb hodgeMinusOne;
	               n := numgens R;
	               diffSummand := for i from 0 to n-1 list (
		       for a in hodgeMinusOneGens list f*diff(R_i,a)-(alpha+p-1)*a*diff(R_i,f));--the second summand
	               hodgeI = ideal mingens (RgeqPplusAlpha + ideal(flatten diffSummand) + f*hodgeMinusOne);
		      
                 );
           );

  hodgeI
  )


---------------------------------------------------------------
---------------------------------------------------------------

hodgeIdealBrieskornPham = method();

hodgeIdealBrieskornPham(List, ZZ, ZZ) :=
hodgeIdealBrieskornPham(List, QQ, ZZ) := (L,alpha,p) -> (
--calculates Hodge ideals for x_1^(b_1)+\cdots+x_n^(b_n)
--L is the list {b_1,\cdots,b_n}
--the b_i's should all be positive
--p is a non-negative integer

    n := #L;
    x := local x;
    R := QQ[x_1..x_n];
    f := sum apply(toList(0..n-1), i -> R_i^(L_i));
    w := apply(L, i -> 1/i);
    hodgeIdealWeightedHomogIsolated(f,alpha,p,w)
    )
    


---------------------------------------------------------------
---------------------------------------------------------------

higherMultiplierIdeal = method();

higherMultiplierIdeal(RingElement, ZZ, ZZ) :=
higherMultiplierIdeal(RingElement, QQ, ZZ) := (f,alpha,p) -> (
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

   pNow := 0;
   isMemberFp := hodgeCheck(f,g,alpha,pNow);
 
   while not isMemberFp do (
       pNow = pNow+1;
       isMemberFp = hodgeCheck(f,g,alpha,pNow);
       );
   
   pNow
   )


-------------------------------------------------------------- 
-------------------------------------------------------------



HRHCheck =  (f,p) -> (

    ------------------------------------------------------------------------
    -- Step 1. Basic setup and annihilator data for f^s
    ------------------------------------------------------------------------
    AnnFsf := AnnFs f;
    Ds := ring AnnFsf;
    G := gens gb AnnFsf;
    numGensDs := numgens Ds;
    ss := Ds_(numGensDs-1);
    DsF := substitute(f,Ds);

    bs := generalB({f},1_(ring f), Exponent => p+1);
    bsFac := factorBFunction bs;
    rhoFp := sort apply(toList(0..(#bsFac-1)), i -> {sub((ring bs)_0-(bsFac#i#0),QQ),bsFac#i#1});--{(root,mult)}
    --the roots are rational numbers

    ------------------------------------------------------------------------
    -- Step 3. Build J_p and the module M = D/J_p
    --(Blanco step 5 (with a different shift on Jp, consistent with b^(p+1)_f))
    ------------------------------------------------------------------------

    Gp := flatten entries G;
    Jp := ideal(Gp) + ideal(DsF^(p+1));
    M := Ds^1 / Jp;
   
    ------------------------------------------------------------------------
    -- Step 4. Compute kernels 
    ------------------------------------------------------------------------

    sLamMapsLess := {};
    sLamMapsAlpha := {};

    for i in rhoFp do (
        lam := sub(i_0, Ds);
        mult := i_1;
        a := (ss - lam)^mult;
	if i_0 < -1-p then sLamMapsLess = append(sLamMapsLess, map(M, M, matrix{{a}}))
	else if i_0 == -1-p then sLamMapsAlpha = append(sLamMapsAlpha, map(M, M, matrix{{a}}));
	);

    G1 := ideal flatten apply( sLamMapsLess, i -> flatten entries gens kernel i);--hard step
    Klams := ideal flatten apply(sLamMapsAlpha, i -> flatten entries gens kernel i);--hard step

    preVg1 := G1;--this corresponds to V^(>alpha)
    preV1 :=  (G1 + Klams);--this corresponds to V^alpha

    isSubset((ss+p+1)*preV1, preVg1+Jp)--checks if F_p(Gr^0_V(B_f))=0
    --uses that F_p(Gr^0_V(B_f))=image(N on F_p(Gr^1_V(B_f)))
    )
   


-------------------------------------------------------------- 
-------------------------------------------------------------


HRHLevel = f -> (

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


