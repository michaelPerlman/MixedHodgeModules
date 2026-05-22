--Copyright 2026 Andras Lorincz and Michael Perlman




nuAlpha = method(Options => {NuMethod => ByAnnFs})
--calculates nu_(g,alpha}(f), the number of roots beta of sufficiently high enough power-b-functions
--(up to a rescaling this is the b-function wrt powers of f), or that of
--the b-function of p_(g,alpha)gf^s, where p_(g,alpha) is the the pFunction
--strategy Malgrange only implemented for g=1

--default NuMethod: ByAnnFs

--used in: weightLength(ByNuAlpha), pFunction

nuAlpha(RingElement,RingElement,ZZ) :=
nuAlpha(RingElement,RingElement,QQ) := options -> (f,g,alpha) -> (
    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    if not (ring g) === (ring f) then error "Expected g to belong to same ring as f"
    else if options.NuMethod == PowerBFunction then (
	nuAlphaPowerBFunction(f,g,alpha))
    else if options.NuMethod == Malgrange then (
	if g != 1_(ring f) then (error "This strategy requires that g is 1 in the ring of f")
	else nuAlphaMalgrange(f,alpha))
    else nuAlphaAnnFs(f,g,alpha)
    )


nuAlpha(RingElement,RingElement,RingElement,ZZ) :=
nuAlpha(RingElement,RingElement,RingElement,QQ) := options -> (f,g,bgfs,alpha) -> (
    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    if not (ring g) === (ring f) then error "Expected g to belong to same ring as f"
    else if options.NuMethod == PowerBFunction then (
	nuAlphaPowerBFunction(f,g,bgfs,alpha))
    else if options.NuMethod == Malgrange then (
	if g != 1_(ring f) then (error "This strategy requires that g is 1 in the ring of f")
	else nuAlphaMalgrange(f,bgfs,alpha))
    else nuAlphaAnnFs(f,g,bgfs,alpha)
    )



--strategy functions:

nuAlphaPowerBFunction=method();
--formerly called nuAlpha
--uses b-function of power

nuAlphaPowerBFunction(RingElement,RingElement,RingElement,ZZ) :=
nuAlphaPowerBFunction(RingElement,RingElement,RingElement,QQ) := (f,g,bgfs,alpha) -> (
    W:=makeWeylAlgebra ring f;
    newG := substitute(g,W);
    newF := substitute(f,W);
    bg := bgfs;
    Roots := bFunctionRoots(bg);
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
    k0 := 1;
    if #AlphaInts> 0 then (
        if AlphaInts_0 > -alpha then (
            k0 =floor(AlphaInts_0 + alpha + 1)));
    PowerB := bg;
    if k0 > 1 then (
	 Anng := polynomialAnnihilator newG;
	 PowerB = (globalB(Anng, newF^(k0))).Bpolynomial);
    facPowerB := factorBFunction PowerB;
    goodFactor := ((ring(PowerB))_0)+(alpha/k0);
    preGoodFactorIndex := select(toList(0..(#facPowerB-1)), i-> facPowerB#i#0 == goodFactor);
    mult := 0;
    if #preGoodFactorIndex > 0 then (
	goodFactorIndex := (preGoodFactorIndex)_0;
	mult = facPowerB#(goodFactorIndex)#1;
	);
    mult
    )
    

nuAlphaPowerBFunction(RingElement,RingElement,ZZ) :=
nuAlphaPowerBFunction(RingElement,RingElement,QQ) := (f,g,alpha) -> (
    W:=makeWeylAlgebra ring f;
    newG := substitute(g,W);
    newF := substitute(f,W);
    Anng := polynomialAnnihilator newG;
    bg := (globalB(Anng, newF)).Bpolynomial;
    Roots := bFunctionRoots(bg);
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
    k0 := 1;
    if #AlphaInts> 0 then (
        if AlphaInts_0 > -alpha then (
            k0 =floor(AlphaInts_0 + alpha + 1)));
    PowerB := bg;
    if k0 > 1 then PowerB = (globalB(Anng, newF^(k0))).Bpolynomial;
    facPowerB := factorBFunction PowerB;
    goodFactor := ((ring(PowerB))_0)+(alpha/k0);
    preGoodFactorIndex := select(toList(0..(#facPowerB-1)), i-> facPowerB#i#0 == goodFactor);
    mult := 0;
    if #preGoodFactorIndex > 0 then (
	goodFactorIndex := (preGoodFactorIndex)_0;
	mult = facPowerB#(goodFactorIndex)#1;
	);
    mult
    )

------------------------------

nuAlphaAnnFs=method();
--default strategy



nuAlphaAnnFs(RingElement,RingElement,RingElement,ZZ) :=
nuAlphaAnnFs(RingElement,RingElement,RingElement,QQ) := (f,g,bgfs,alpha) -> (
    Annfs := local Annfs;
    B := local B;
    goodFactor := local goodFactor;
    Anng :=local Anng;
    newG := local newG;
    W := makeWeylAlgebra ring f;
    newF := sub(f,W);
    bg := bgfs;
    Roots := bFunctionRoots bg;
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
    k0:=1;
    if #AlphaInts > 0 then (
        if AlphaInts_0 > -alpha then (
            k0 = floor (AlphaInts_0+alpha+1)));
    mult := 0;
    if k0 == 1 then (
    	B = factorBFunction bg;
    	goodFactor = ((ring(bg))_0) + (alpha);
    	goodFactorIndexList := select(toList(0..(#B-1)), i-> B#i#0 == goodFactor);
    	if #goodFactorIndexList >0 then mult = B#((goodFactorIndexList)_0)#1;
	)
    else (
    	if g == 1_(ring f) then (
	    Annfs = AnnFs newF;	
	    )
    	else (
	    newG = sub(g,W);
	    Anng = polynomialAnnihilator newG;
	    Annfs = AnnIFs(Anng,newF);
	    );
   	newerF := sub(newF, ring(Annfs));
    	J := Annfs + ideal(newerF^(k0));
    	twists := toList(0..k0-1);
    	listOfTwists := apply(twists, i-> sub(bg,{(ring bg)_0 => (ring bg)_0+i}));
    	B = product listOfTwists;
    	factorB := factorBFunction B;
    	numFacsB := #factorB;
    	factorIndex := toList(0..numFacsB-1);
    	listOfReducedFacs := apply(factorIndex, i-> factorB#i#0);
    	counter := 0;
    	bool := false;
    	while counter < numFacsB and bool == false do (
	    bool = (listOfReducedFacs_(counter) == ( (ring B)_0+ alpha));
	    counter = counter + 1);
  	if bool == true then (
       	    goodFactor = factorB#(counter-1)#0;
       	    mult = factorB#(counter-1)#1;
       	    memb := true;
       	    while (memb == true) and (mult > 0) do (
	   	mult = mult - 1;
	   	B = B/ (sub(goodFactor, ring B));
	   	B = sub(B, (ring goodFactor));
	   	BB := sub(B, ring(Annfs));
	   	----BB:=sub(B, {(ring B)_0 => -(ring MalgIn)_(2*n+1)*(ring MalgIn)_n});
	   	memb = ((BB % J) == 0_(ring J)));
       	    mult = mult + 1));
     	mult
     	)


nuAlphaAnnFs(RingElement,RingElement,ZZ) :=
nuAlphaAnnFs(RingElement,RingElement,QQ) := (f,g,alpha) -> (
    Annfs := local Annfs;
    B :=local B;
    goodFactor := local goodFactor;
    bf := local bf;
    Anng := local Anng;
    newG := local newG;
    W := makeWeylAlgebra ring f;
    newF := sub(f,W);
    if g == 1_(ring f) then (
       	bf = cachedGlobalBFunction f;
	)
    else (
	newG = sub(g,W);
	Anng = polynomialAnnihilator newG;
	bf = (globalB(Anng, newF)).Bpolynomial;
	);
    Roots := bFunctionRoots bf;
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i + alpha,QQ))==1);
    k0 := 1;
    if #AlphaInts> 0 then (
        if AlphaInts_0 > -alpha then (
            k0=floor (AlphaInts_0+alpha+1)));
    mult := 0;
    if k0 ==1 then (
    	B = factorBFunction bf;
    	goodFactor = ((ring(bf))_0) + (alpha);
    	goodFactorIndexList := select(toList(0..(#B-1)), i-> B#i#0 == goodFactor);
    	if #goodFactorIndexList >0 then mult= B#((goodFactorIndexList)_0)#1;
	)
    else (
    	if g == 1_(ring f) then (
	    Annfs = AnnFs newF;	
	    )
    	else (
	    Annfs = AnnIFs(Anng,newF);
	    );
   	newerF := sub(newF, ring(Annfs));
    	J := Annfs + ideal(newerF^(k0));
    	twists := toList(0..k0-1);
    	listOfTwists := apply(twists, i-> sub(bf,{(ring bf)_0 => (ring bf)_0+i}));
    	B = product listOfTwists;
    	factorB := factorBFunction B;
    	numFacsB := #factorB;
    	factorIndex := toList(0..numFacsB-1);
    	listOfReducedFacs := apply(factorIndex, i-> factorB#i#0);
    	counter := 0;
    	bool := false;
    	while counter < numFacsB and bool == false do (
	    bool = (listOfReducedFacs_(counter) == ( (ring B)_0 + alpha));
	    counter = counter + 1);
  	if bool == true then (
       	    goodFactor = factorB#(counter-1)#0;
       	    mult = factorB#(counter-1)#1;
       	    memb :=true;
       	    while (memb == true) and (mult > 0) do (
	   	mult = mult - 1;
	   	B = B/ (sub(goodFactor, ring B));
	   	B = sub(B, (ring goodFactor));
	   	BB := sub(B, ring(Annfs));
	   	memb = ((BB % J) == 0_(ring J)));
       	    mult = mult + 1));
     	mult
     	)

        
---------------------------------------------------------
    
 
nuAlphaMalgrange=method();
--formerly called nuAlphaF
--uses initial of Malgrange ideal and product as upper bound
--only implemented for g=1


--uses: MalgrangeIdeal 


nuAlphaMalgrange(RingElement, ZZ) :=
nuAlphaMalgrange(RingElement, QQ) := (f, alpha) -> (
    bf := cachedGlobalBFunction f;
    nuAlphaMalgrange(f,bf,alpha)
    )
    

nuAlphaMalgrange(RingElement,RingElement,ZZ) :=
nuAlphaMalgrange(RingElement,RingElement,QQ) := (f,bf,alpha) -> (
    W := makeWeylAlgebra ring f;
    Roots := bFunctionRoots bf;
    AlphaInts := reverse sort select(Roots, i-> denominator(sub(i+alpha,QQ))==1);
    k0 := 1;
    if #AlphaInts> 0 then (
        if AlphaInts_0 > -alpha then (
            k0 = floor (AlphaInts_0 + alpha + 1)));
    newF := sub(f,W);
    AnnI := MalgrangeIdeal {(newF)^(k0)};
    DY := ring AnnI;
    n := numgens DY // 2 -1; -- DY = k[x_1,...,x_n,t_1,...,t_r,dx_1,...,dx_n,dt_1,...,dt_r]
    w := toList(n:0) | toList(1:1);
    MalgIn := inw(AnnI, -w|w);
    twists := toList(0..k0-1);
    listOfTwists := apply(twists, i-> sub(bf,{(ring bf)_0 => (ring bf)_0+i}));
    bigB := product listOfTwists;
    B := sub(bigB, {(ring bigB)_0 => k0 * ((ring bigB)_0)});
    factorB := factorBFunction B;
    numFacsB := #factorB;
    factorIndex := toList(0..numFacsB-1);
    listOfReducedFacs := apply(factorIndex, i-> factorB#i#0);
    counter := 0;
    bool := false;
    while counter < numFacsB and bool == false do (
	bool = (listOfReducedFacs_(counter) == ( (ring B)_0 + alpha/(k0)));
	counter = counter + 1);
    mult := 0;
    if bool == true then (
       goodFactor := factorB#(counter-1)#0;
       mult = factorB#(counter-1)#1;
       memb :=true;
       while (memb == true) and (mult > 0) do (
	   mult = mult-1;
	   B = B/ (sub(goodFactor, ring B));
	   B = sub(B, (ring goodFactor));
	   BB := sub(B, {(ring B)_0 => -(ring MalgIn)_(2*n+1)*(ring MalgIn)_n});
	   memb = ((BB % MalgIn) == 0_(ring MalgIn)));
       mult = mult+1);
     mult
     )
    


-------------------------------------------------------------- 
------------------------------------------------------------- 
--helper for pFunction below 
 
translateToAlpha = (alpha,beta) -> (
     translates := {};
     bool := (beta >= -alpha);
     counter :=beta;
     while bool do (
	 translates = append(translates, counter);   
	 counter = counter - 1;
	 bool = (counter >= -alpha));
     translates
     )
 
 
 ----------------------------------------------------------------
 
 
pFunction = method(Options => {NuMethod => ByAnnFs});
--to do: modify so don't compute high nus

--user may select NuMethod using syntax above
--default is ByAnnFs

--uses: translateToAlpha, nuAlpha


pFunction(RingElement,RingElement,ZZ) :=
pFunction(RingElement,RingElement,QQ) := options -> (f,g,alpha) -> (
    if sub(alpha,QQ) <= 0 then error "expected alpha to be a positive rational number";
    if (options.NuMethod == Malgrange) and (g != 1_(ring f)) then (
	error "This strategy requires that g is 1 in the ring of f");
    W := makeWeylAlgebra ring f;
    newG := substitute(g,W);
    newF := substitute(f,W);
    Anng := polynomialAnnihilator newG;
    bgfs := (globalB(Anng, newF)).Bpolynomial;
    Roots := bFunctionRoots bgfs;
    translatesRoots := unique flatten apply(Roots, i-> translateToAlpha(alpha,i));
    nus := apply(translatesRoots, i-> (i,nuAlpha(f,g,bgfs,-i, NuMethod => options.NuMethod)));
    pp := product apply(nus, i-> (((ring bgfs)_0)-i_0)^(i_1));
    pp
    )
