--Copyright 2026 by Andras Lorincz and Michael Perlman





--------------------------------------------------------------------
--helper functions to represent S-presentations of Gr^F_p
--------------------------------------------------------------------

Epair = (f,I,J,k) -> (
--presents Gr^F_k H^1_f(S).  Besides the image f*I_{k-1}, the numerator
--f^(k+1) represents an element of S inside S_f and is therefore zero in
--H^1_f(S)=S_f/S.  The latter relation is essential for k=0; for k<0 both
--filtration pieces are zero.

    relations := module(f*I);
    if k >= 0 then relations = relations + module ideal(f^(k+1));
    M := module(J)/relations;
    P := matrix entries presentation M;

    P
    )

--------------------------------------------------------------------

Ematrix = (f,idealList,k0) -> (

    S := ring f;

    t := #idealList;
    Epairs := apply(toList(0..t-2), i -> Epair(f,idealList_i, idealList_(i+1),k0+i));

    EpairsFixZeros := apply(Epairs, M -> if M == 0_ZZ then matrix{{0_S}} else M);
    
    Emat := directSum(EpairsFixZeros);

    Emat
    )


--------------------------------------------------------------------
--helper functions for UV matrices, which encode syzygies between
--two graded pieces via differential operators
--------------------------------------------------------------------



diSyz = (f,IJ,aa,p,i) -> (
--output {U^i,V^i} for IJ in idealPairs of UVList    

    S := ring f;

    I := IJ_1;
    J := IJ_2;

    diI := ideal(apply(flatten entries gens I, m -> diff(S_i,m)*f-(aa+p+IJ_0+1)*m*diff(S_i,f)));--double check factor

    syzMat := matrix entries syz gens (diI + J);

    UBlockRows := toList(0..IJ_3-1);


    UBlock := submatrix(syzMat, UBlockRows, );

    VBlockRows := toList(IJ_3..(IJ_3+IJ_4-1));
    

    VBlock := submatrix(syzMat, VBlockRows, );


    {UBlock, VBlock}   
    )

--------------------------------------------------------------------


UVpair = (f,IJ,aa,p,n) -> (
--output UV={{U^1,V^1}..{U^n,V^n}} for IJ in idealPairs of UVList
--each U has IJ_3 many rows and each V has IJ_4 many rows

    syzMats := apply(toList(0..n-1), i -> diSyz(f,IJ,aa,p,i));

    syzMats

    )

--------------------------------------------------------------------


UVList = (f,idealList, rList, aa, p, S) -> (
--output list of UVpair for each IJ in idealPairs


    n := numgens S;

    --idealList = {I_{p+aa-1}..I_{p+bb}}

    newIdealList := drop(idealList,1);--don't need first ideal

    --newIdealList = {I_{p+aa}..I_{p+bb}}

    t := #newIdealList;--\geq 2

    idealPairs := apply(toList(0..t-2), i -> {i,newIdealList_i, newIdealList_(i+1),rList_i, rList_(i+1)});--{i,I_i, I_{i+1}, rList_i, rList_(i+1)}
    

    UVlist := apply(idealPairs, IJ -> UVpair(f,IJ,aa,p,n));

    UVlist
    )


--------------------------------------------------------------------


scaleU = (x,A) -> (

    n := numgens A;

    xScaled := apply(toList(0..n-1), i -> {A_i*((x_i)_0), (x_i)_1});

    xScaled
    )


-------------------------------------------------------------------

zeroMat = (A,r,c) -> map(A^r, A^c, (i,j) -> 0_A)

hcat = L -> (
    if #L == 0 then error "hcat: empty list";
    M := L_0;
    for i from 1 to #L-1 do M = M | L_i;
    M
    )

vcat = L -> (
    if #L == 0 then error "vcat: empty list";
    M := L_0;
    for i from 1 to #L-1 do M = M || L_i;
    M
    )

-------------------------------------------------------------------

putTogether = (UVAscaled, rList) -> (
-- UVAscaled_k = {{xi_1 U_k^1, V_k^1}, ..., {xi_n U_k^n, V_k^n}}
-- rList = {r_{p+aa}, ..., r_{p+bb}}
-- output = the full UV block in staircase form

    numRowss := #rList;
    numSlots := #UVAscaled;

    if numSlots != numRowss-1 then error "putTogether: wrong number of UV blocks";

    XiUList := apply(UVAscaled, x -> hcat(apply(x, uv -> uv_0)));
    VList   := apply(UVAscaled, x -> hcat(apply(x, uv -> uv_1)));

    A := ring XiUList_0;
    cList := apply(XiUList, M -> numcols M);

    rows := for j from 0 to numRowss-1 list (
        blocks := for s from 0 to numSlots-1 list (
            if s == j then XiUList_s
            else if s == j-1 then VList_(j-1)
            else zeroMat(A, rList_j, cList_s)
        );
        hcat(blocks)
    );

    vcat(rows)
)


    

    
--------------------------------------------------------------------
---helper function to encode f^k torsion of local cohomology and IC
--------------------------------------------------------------------

fTorsion = (f,rList,aa,p) -> (

    diagEntries := {};

    m := #rList;

    for j from 0 to m-1 do (
	fPowerList := apply(toList(0..rList_j-1), i -> f^(max{0,j+aa+p+1}));
	diagEntries = append(diagEntries, fPowerList);
	);

    diagonalMatrix flatten diagEntries
    )

    

--------------------------------------------------------------------
---helper function to create koszul complex, tensor with N, and take
--appropriate graded piece
--------------------------------------------------------------------   
    
    
koszulSlice = (A,S,p,N) -> (
--N is a graded A module
--form the Koszul complex, tensor with N
--take the p-th graded piece
--sub back into S

    n := numgens S;

    xigens := ideal( apply(toList(0..n-1), i -> A_(i)));
    
    K :=  freeResolution xigens;

    Ktwisted := K**A^{1:{n,0}};
    
    DRN := Ktwisted**N;

    FF := freeResolution(DRN);
    
    DRp := part({p,0},FF);--a complex of QQ[e_1..e_n]-modules

    QQe := ring DRp;

    QQeToS := map(S,QQe, gens S);

    DRpS := QQeToS(DRp);

    DRpS
    )



--------------------------------------------------------------------
---Main technical function for de Rham cohomology
--------------------------------------------------------------------


deRhamInterval = method(Options => {InputType => HodgeIdeals, Weights => null})

--HodgeIdeals will give H^1_f(S)
--WeightedHodgeIdeals will give IC_f
--WeightedHomogIsolated will use hodgeIdealWeightedHomogIsolated with Weights => w

deRhamInterval(RingElement, List, ZZ) := options -> (f,B,p) -> (
--B={a,b} is a list with -n\leq a\leq b \leq 0
--p is an integer
--output is the complex (K\otimes N)_p
--which gives correct cohomology in degrees a..b
--if a=-n and b=0, then gives correct cohomology everywhere

--gives as free complex quasi-isomorphic to DR


--step 0: initialize information
     S := ring f;
     n := numgens S;
     aa := max{B_0-2,-n}+n;--shift so aa>= 0
     bb := min{B_1,0}+n;--shift so bb<= n
 

     -- step 1: calculate necessary ideals and get their numgens
     --need Gr_{p+j}(M) for all aa <= j <= bb
     --need I_{p+j} for all aa-1 <= j <= bb

     idealList := {};

     if options.InputType == HodgeIdeals then (
	 HodgeIdealTemp := (f,k) -> (if k < 0 then ideal(0_S) else hodgeIdeal(f,1,k));
	 idealList = apply(toList(aa-1..bb), i -> HodgeIdealTemp(f,p+i));
	 );

     if options.InputType == WeightedHodgeIdeals then (
	 WHodgeIdealTemp := (f,k) -> (if k <0 then ideal(0_S) else weightedHodgeIdeal(f,1,k,1));
	 idealList = apply(toList(aa-1..bb), i -> WHodgeIdealTemp(f,p+i));
	 );

     if options.InputType == WeightedHomogIsolated then (
	 w := options.Weights;
	 if w === null then error "deRhamInterval: InputType => WeightedHomogIsolated requires a Weights option";
	 if #w != n then error("deRhamInterval: Weights must have length ", toString n);
	 WHIHodgeIdealTemp := (f,k) -> (if k < 0 then ideal(0_S) else hodgeIdealWeightedHomogIsolated(f,1,k,w));
	 idealList = apply(toList(aa-1..bb), i -> WHIHodgeIdealTemp(f,p+i));
	 );

     --idealList = {I_{p+aa-1}..I_{p+bb}}
     
     rList := drop( apply(idealList, I -> numgens I), 1);
     
     -- step 2: calculate Ematrix

     E := Ematrix(f,idealList,p+aa);--a matrix with sum rList many rows
     

     -- step 3: calculate UV matrix

     UV := UVList(f,idealList,rList,aa,p,S);--{UV_(p+aa)..UV_(p+bb-1)}
     --here, UV_j={{U^1,V^1}..{U^n,V^n}}


     -- step 4: calculate f-torsion matrix

     fTors := fTorsion(f,rList,aa,p);--a matrix with sum rList many rows

     -- step 5: create A



     e := local e;
     de := local de;
     eDegs := apply(toList(0..n-1), i-> 0);
     deDegs := apply(toList(0..n-1), i -> 1);
     A := QQ[e_1..e_n, Degrees => eDegs][de_1..de_n, Degrees => deDegs];

     StoA := map(A,S, apply(toList(0..n-1), i -> e_(i+1)));

     -- step 6: sub all matrices into A
     -- scale U matrix approproately

     EA := StoA(E);
     fTorsA := StoA(fTors);


     --UV stuff

     UVA := apply(UV, x -> apply(x, uv -> {StoA(uv_0), StoA(uv_1)}));


     UVAscaled := apply(UVA, x -> scaleU(x,A));


     UVMatrix := putTogether(UVAscaled, rList);

     presN := EA|UVMatrix|fTorsA;


     --step 7: create presentation matrix N
     --make sure all matrices are stripped of degrees first

     summnds := apply(toList(0..#rList-1), i -> A^{rList_i:{ -aa-p-i,0}});
     

     freeMod := directSum(summnds);


     presNGraded := map(freeMod, ,presN);


     N := cokernel presNGraded;

     --step 8: calculate Koszul cohomology and slice
     --and sub back into S

     GRDRp := koszulSlice(A,S,p,N);

     prune GRDRp
     )
     



--------------------------------------------------------------------
-------------------------------------------------------------------- 

gradedDeRhamComplex = method(Options => {InputType => HodgeIdeals, Weights => null})
--gives the entire de Rham complex

gradedDeRhamComplex(RingElement,ZZ) := options -> (f,p) -> (

    S := ring f;
    n := numgens S;
    B := {-n,0};

    DRp := deRhamInterval(f,B,p,
	InputType => options.InputType,
	Weights => options.Weights);

    DRp
    )

--------------------------------------------------------------------
-------------------------------------------------------------------- 


gradedDeRhamComplexH1 = method();

gradedDeRhamComplexH1(RingElement, ZZ) := (f,p) -> (
--in *homological* degrees n..0
--free complex quasi-isomorphic

    DRp := gradedDeRhamComplex(f,p, InputType => HodgeIdeals);

    DRp
    )

gradedDeRhamComplexH1(RingElement, ZZ, List) := (f,p,w) -> (
--WeightedHomogIsolated variant: f is weighted homogeneous of weight 1 with
--respect to w with an isolated singularity at the origin.

    gradedDeRhamComplex(f, p,
	InputType => WeightedHomogIsolated,
	Weights => w)
    )

--------------------------------------------------------------------
-------------------------------------------------------------------- 

intersectionDuBoisComplex = method();

intersectionDuBoisComplex(RingElement, ZZ) := (f,p) -> (
--gives a free complex
--interesting homological degrees 0,-1,..,-n
--our convention makes it the same as gradedDuBoisComplex for RHM

    S := ring f;
    n := numgens S;

    DRp := gradedDeRhamComplex(f,p-n, InputType => WeightedHodgeIdeals);

    RHomDRp := Hom(DRp, S^1);
    
    RHomDRpShift := RHomDRp[-p-1];

    prune RHomDRpShift
    )

--------------------------------------------------------------------
-------------------------------------------------------------------- 


gradedDeRhamCohomologyH1 = method()

gradedDeRhamCohomologyH1(RingElement,ZZ,ZZ) := (f,p,q) -> (
--q should be an integer between -n and 0
--p should be \geq -n

    S := ring f;
    n := numgens S;

    B := {q-1,q+1};
    if q == -n then B={-n,-n+1};
    if q == 0 then B={-1,0};

    DRpq := deRhamInterval(f,B,p, InputType => HodgeIdeals);

    H := prune HH_(-q)(DRpq);

    H
    )

gradedDeRhamCohomologyH1(RingElement,ZZ,ZZ,List) := (f,p,q,w) -> (
--WeightedHomogIsolated variant.

    S := ring f;
    n := numgens S;

    B := {q-1,q+1};
    if q == -n then B={-n,-n+1};
    if q == 0 then B={-1,0};

    DRpq := deRhamInterval(f,B,p,
	InputType => WeightedHomogIsolated,
	Weights => w);

    prune HH_(-q)(DRpq)
    )


--------------------------------------------------------------------
-------------------------------------------------------------------- 


duBoisComplex = method();

duBoisComplex(RingElement, ZZ) := (f,p) -> (
--gives Omega^p_D as a free complex
--interesting homological degrees 0,-1,..,-n
--Prop 13.1 Mustata--Popa "Hodge filtration on local cohomology..."

    S := ring f;
    n := numgens S;

    DRk := gradedDeRhamComplexH1(f,p-n);

    RHomDRk := Hom(DRk, S^1);

    --do the correct cohomological shift

    RHomDRkShift := RHomDRk[-p-1];--1 codim

    prune RHomDRkShift
    )


duBoisComplex(RingElement, ZZ, List) := (f,p,w) -> (
--WeightedHomogIsolated variant.

    S := ring f;
    n := numgens S;

    DRk := gradedDeRhamComplexH1(f,p-n,w);

    RHomDRk := Hom(DRk, S^1);

    RHomDRkShift := RHomDRk[-p-1];

    prune RHomDRkShift
    )

    
--------------------------------------------------------------------
--------------------------------------------------------------------
--preDB functions
    

isPreDuBois = method();

isPreDuBois(RingElement, ZZ) := (f,m) -> (
--outputs true if f is pre m-Du Bois. False otherwise.

   S := ring f;
   n := numgens S;

   isPre := true;
   p := 0;

   while isPre and p<=m do (

       DBC := duBoisComplex(f,p);
       coho := apply(toList(1..n), i -> prune HH_(-i)(DBC));
       if any(coho, H -> H != 0) then isPre = false;
       p = p+1;
       );

   isPre
   )


isPreDuBois(RingElement, ZZ, List) := (f,m,w) -> (
--WeightedHomogIsolated variant.

   S := ring f;
   n := numgens S;

   isPre := true;
   p := 0;

   while isPre and p<=m do (

       DBC := duBoisComplex(f,p,w);
       coho := apply(toList(1..n), i -> prune HH_(-i)(DBC));
       if any(coho, H -> H != 0) then isPre = false;
       p = p+1;
       );

   isPre
   )
       

  


end

restart
load "MixedHodgeModules.m2"
installPackage "MixedHodgeModules"
load "MixedHodgeModules/deRhamDraft.m2"
viewHelp MixedHodgeModules
check "MixedHodgeModules"
uninstallPackage "MixedHodgeModules"

