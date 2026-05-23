--Copyright 2026 by Andras Lorincz and Michael Perlman

--Important: to use this package one must do:
--(1) maintain the current file structure:
-- MixedHodgeModules(folder) -> {MixedHodgeModules (folder), MixedHodgeModules.m2 (this file)} -> {aux files}
--(2) launch the package from this file, using code after "end"



--to do:
--(1) fill in descriptions and examples in documentation for nuAlpha, pFunction (and their strategies)
--(2) add examples in doc of other strategies for weightLength
--(3) create tests for every function and add to tests.m2



--change/add soon:
--(1) change weightCheck and weightLevel to input alpha\in [0,1) and k\in ZZ? (instead of beta=alpha+k)
--same for weightLength? Same for hodgeLevel, hodgeCheck?
--(2) allow alpha as a parameter in hodgeIdeal (work over QQ(alpha) where alpha is variable)?
--(3) weighted Homogeneous isolated option in de Rham and Du Bois functions?
--(4) There is an old comment on pFunction "to do: modify so don't compute high nus". 
--(5) cache generalized b-functions for speed.



--to add eventually?:
--(1) determinantal code (GLmnReps on my website can calculate I_{lambda})?
--(2) use minimal exponent formula for generation level as optional strategy
--in Hodge ideal?
--(3) allow user to input generation level as optional strategy?
--(4) get generation level on IC_f using weighted Hodge ideal?
--(5) Bernstein--Sato polynomials on singular ambient varieties X (Dirks) and microlocal b-functions on
-- singular ambient varieties?
--(8) nearby cycles and vanishing cycles? vanishing cycles are implicitly calculated in HRHCheck



newPackage(
    "MixedHodgeModules",
    Version => "1.0",
    Date => "May 18, 2026",
    Headline => "Calculations involving Hodge and weight filtrations on localizations",
    Authors => {{ Name => "András C. Lőrincz",  Email => "lorincz@ou.edu",  HomePage => "https://math.ou.edu/~alorincz/"},
	        { Name => "Michael Perlman",    Email => "mperlman@ua.edu", HomePage => "https://sites.google.com/view/michaelperlman/home"}},
    HomePage => "https://github.com/michaelPerlman/MixedHodgeModules",
    Keywords => {"D-modules"},
    AuxiliaryFiles => true,
    DebuggingMode => false,
    PackageImports => {"Polyhedra"},
    PackageExports => {"BernsteinSato", "Complexes"}
    )

importFrom(BernsteinSato, {"MalgrangeIdeal", "eliminateWA"})


export {
    --functions:
    "adjointIdeal",
    "doesGenerateNext",
    "generateNext",
    "generationLevel",
    "higherMultiplierIdeal",
    "hodgeCheck",
    "hodgeIdeal",
    "hodgeIdealBrieskornPham",
    "hodgeIdealWeightedHomogIsolated",
    "hodgeLevel",
    "hodgeOnV",
    "HRHCheck",
    "HRHLevel",
    "IHmoduleAdjoint",
    "localCohomFW",
    "monodromyWeightHodgeOnV",
    "nuAlpha",
    "pFunction",
    "weightCheck",
    "weightLength",
    "weightLevel",
    "weightedHodgeIdeal",
    "weightHodgeOnV",

    "duBoisComplex",
    "intersectionDuBoisComplex",
    "isPreDuBois",
    "gradedDeRhamComplexH1",
    "gradedDeRhamCohomologyH1",

    --symbols:
     "ByAnnFs",
     "ByNuAlpha",
     "ByWeightLevel",
     "LengthStrategy",
     "Malgrange",
     "NuMethod",
     "PowerBFunction",
     "dtBasis",
     "False",
     "sBasis",
     "True",
     "UseBasis",
     "UseGenLevel",
    }



--internal symbols for deRham
protect symbol InputType;
protect symbol WeightedHodgeIdeals;
protect symbol HodgeIdeals;

--internal symbol for caching AnnFs on a RingElement
protect symbol AnnFsCache;

--internal symbol for caching globalBFunction on a RingElement
protect symbol GlobalBCache;

--internal symbol for caching rhoFp (roots of generalized b-function with multiplicities) keyed by (f,p)
protect symbol RhoFpCache;

--internal symbol for caching hodgeOnV outputs keyed by (f, alpha, p, UseBasis) or (f, p, UseBasis)
protect symbol HodgeOnVCache;

--internal symbol for caching weightHodgeOnV outputs keyed by (f, alpha, p, m, UseBasis)
protect symbol WeightHodgeOnVCache;


--Convention: dt^p \in F_p(B_f)
-- f should be reduced


---------------------------------------------------------------
--hodgeOnV, hodge ideals, higher multiplier ideals, HRH
---------------------------------------------------------------
load "./MixedHodgeModules/hodgeFiltrations.m2"

---------------------------------------------------------------
--weightHodgeOnV, weighted hodge ideals, weightCheck, weightLevel, localCohomFW
---------------------------------------------------------------
load "./MixedHodgeModules/weightFiltrations.m2"

---------------------------------------------------------------
--nuAlpha and p-functions
---------------------------------------------------------------
load "./MixedHodgeModules/nuAlpha.m2"

---------------------------------------------------------------
--deRham and Du Bois complexes
---------------------------------------------------------------
load "./MixedHodgeModules/deRham.m2"

---------------------------------------------------------------
--files related to tests and documentation
---------------------------------------------------------------
load "./MixedHodgeModules/tests.m2"

load "./MixedHodgeModules/documentation.m2"

---------------------------------------------------------------

end



---------------------------------------------------------------
---------------------------------------------------------------
	

restart
load "MixedHodgeModules.m2"
viewHelp MixedHodgeModules
check "MixedHodgeModules"--runs all tests in tests.m2 
viewHelp hodgeIdeal
viewHelp higherMultiplierIdeal

restart
installPackage "MixedHodgeModules"
viewHelp MixedHodgeModules
check "MixedHodgeModules"
uninstallPackage "MixedHodgeModules"


--tests/ examples

S=QQ[x,y,z]
f=x^2+y^3+y*z^2
time hodgeOnV(f,2)

time hodgeOnV(f,3/4,2)


R=QQ[x,y,z,w]
f=x*w-y*z
g=1_R
alpha=1+1
weightCheck(f,g,alpha,1)

--------------------------------------------------------------
--smooth case
R=QQ[x,y]
f=y^2+x
isPreDuBois(f,1)

g=1_R
alpha=6
weightLevel(f,g,alpha)



for p from 0 to 4 do print doesGenerateNext(f,1,p)
generationLevel(f)
p=5
u={1,1/2}
testHodgeIdealforWeightedHomog(f,p,u, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,u)

p=3
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)


p=3
compareHodgeIequalHigherMultI(f,p)

testHodgeIhigherMultImodf(f,p)



---------------------------------------------------------------

--normal crossing
n=3
R=QQ[x_1..x_n]
f=product gens R
hodgeOnV(f,0)
generationLevel(f)
for p from 0 to 4 do print doesGenerateNext(f,1,p)
for p from 0 to 5 do print hodgeIdeal(f,1,p)

p=2
m=n
testWeightedHodge(f,p,m)
adjointIdeal(f)



S=QQ[x,y,z,w]
f=x*y*z
alpha=1
hodgeIdeal(f,alpha,0)
hodgeIdeal(f,alpha,1)
hodgeIdeal(f,alpha,2)
hodgeIdeal(f,alpha,3)

---------------------------------------------------------------

--2x2 determinant
S=QQ[x,y,z,w]
f=x*w-y*z
isPreDuBois(f,3)
alpha=1
hodgeIdeal(f,alpha,0)
hodgeIdeal(f,alpha,1)
hodgeIdeal(f,alpha,2)
hodgeIdeal(f,alpha,3)
hodgeIdeal(f,alpha,4)


hodgeOnV(f,2, UseBasis => sBasis)
weightHodgeOnV(f,1,2,1, UseBasis => sBasis)
ideal mingens ideal apply(apply(oo, g -> (sub(f,ring g))^2*sub(g, {(ring g)_0 => -1})), h -> sub(h,ring f))
weightHodgeOnV(f,1,2,1)
weightHodgeOnV(f,1,2,1, UseBasis => dtBasis)
hodgeOnV(f,1,1)



---------------------------------------------------------------


--cusp
R=QQ[x,y]
f=x^2+y^3
hodgeOnV(f,1) 
alpha=1
generationLevel(f,alpha)
w={1/3,1/2}
p=4
isPreDuBois(f,0)

--Davis--Yang page 4
--Zhang Conjecture E

alpha=5/6+1/10
p=2
hodgeIdeal(f,alpha,p)
higherMultiplierIdeal(f,alpha,p)



---------------------------------------------------------------

--A1 sing
R=QQ[x,y,z]
f=y^2-x*z
hodgeOnV(f,1)
w={1/2,1/2,1/2}
p=5



---------------------------------------------------------------

--cusp intersect line
R=QQ[x,y]
f=x*(y^2-x^3)
hodgeOnV(f,1)
for p from 0 to 2 do print hodgeIdeal(f,1,p)

ll=weightLength(f,1)
p=1
m=ll
adjointIdeal(f)



---------------------------------------------------------------


--MP Remark 17.13
R=QQ[x,y]
f=x*y*(x+y)
hodgeOnV(f,1)
alpha=1
for p from 0 to 2 do print hodgeIdeal(f,alpha,p)

p=0
m=1
adjointIdeal(f)


---------------------------------------------------------------

--ordinary m-fold point in PP^{n-1}
--see MP Example 20.1
n=3
m=3
m <= n/2
R=QQ[x_1..x_n]
f=sum apply(toList (1..n), i-> x_i^m)
alpha=1
time for p from 0 to 1 do print hodgeIdeal(f,alpha,p)
2*m-n

---------------------------------------------------------------

--3x3 determinant
R=QQ[x_(1,1)..x_(3,3)]
f=determinant genericMatrix(R,x_(1,1),3,3)
M=genericMatrix(R,x_(1,1),3,3)
I=minors(2,M)

HRHCheck(f,1)
HRHLevel(f)
isPreDuBois(f,2)

alpha=1
for p from 0 to 1 do print hodgeIdeal(f,alpha,p)

p=0
m=1


DB=gradedDuBoisComplex(f,1)
dd^DB
---------------------------------------------------------------


--3x3 symmetric determinant
R=QQ[x_1..x_6]
f=determinant genericSymmetricMatrix(R,x_1,3)
hodgeOnV(f,1)
alpha=1
time for p from 0 to 1 do print hodgeIdeal(f,alpha,p)
alpha=1/2
time for p from 0 to 1 do print hodgeIdeal(f,alpha,p)

p=0
m=1
adjointIdeal(f)


---------------------------------------------------------------

--Blanco example 1
R=QQ[x,y]
f=x^5+y^5+x^2*y^2
V=hodgeOnV(f,0)

for a in keys V do (
     for p from 0 to 0 do print {a,hodgeIdeal(f,a,p)})

 
p=0
m=1
adjointIdeal(f)



---------------------------------------------------------------

--Blanco example 2
R=QQ[x,y]
lambda=1/2
f=(y^2-x^3)*(y^2+lambda*x^3)
V=hodgeOnV(f,0) 
for a in keys V do (
     for p from 0 to 0 do print {a,hodgeIdeal(f,a,p)})

p=0
m=1
adjointIdeal(f)



---------------------------------------------------------------


--Blanco example 3
R=QQ[x,y,z]
f=x^3+y^3+z^3+x*y*z
hodgeOnV(f,0)--doesnt finish on my computer

w={1/3,1/3,1/3}
p=0


adjointIdeal(f)


---------------------------------------------------------------

--An singularity
--n\geq 1
R=QQ[x,y,z];
n=2
f=x^2+y^2+z^(n+1);
hodgeOnV(f,1)
w={1/2,1/2,1/(n+1)}
p=2



---------------------------------------------------------------


--Dn singularity
--n\geq 4
R=QQ[x,y,z];
n=4
f=x^2+y^(n-1)+y*z^2;
factorBFunction generalB({f},1_S,Exponent => 2)

hodgeOnV(f,1)
w={1/2,1/(n-1),(n-2)/(2*(n-1))}
p=1


---------------------------------------------------------------
    


--E6 singularity
S=QQ[x,y,z]
f=x^2+y^3+z^4
hodgeOnV(f,1) 
w={1/2,1/3,1/4}
alpha =
hodgeIdealWeightedHomogIsolated(f,alpha,p,w)

---------------------------------------------------------------
    

--E7 singularity
S=QQ[x,y,z]
f=x^2+y^3+y*z^3
globalBFunction(f)
factorBFunction oo
alpha=1
w={1/2,1/3,2/9}
p=1

---------------------------------------------------------------

--E8 singularity
S=QQ[x,y,z]
f=x^2+y^3+z^5
globalBFunction(f)
factorBFunction oo
alpha=1/2
w={1/2,1/3,1/5}
hodgeIdealWeightedHomogIsolated(f,alpha,3,w)


---------------------------------------------------------------

--Saito Example
S=QQ[x,y]
f=x^(14)+y^(14)-x^6*y^6
alpha=1
hodgeIdeal(f,alpha,0)--wont run


---------------------------------------------------------------

restart
load "MixedHodgeModules.m2"

--find counterexample to Zhang Conjecture E

S=QQ[x,y,z]
f=x^3+y^3+z^3+x*y*z
hodgeOnV(f,1)

alpha=1
p=2

for i from 0 to p do (
    I := hodgeIdeal(f,alpha,i);
    print I;
    J := higherMultiplierIdeal(f,alpha,i);
    print J;
    )



--------------------------------------------------------------
S=QQ[x_(1,1)..x_(3,2)]
M= transpose genericMatrix(S,x_(1,1),2,3)
I=minors(2,M)
H=localCohomFW(I,2,1,2+6)
for  i from 0 to 10 do print hilbertFunction(i,H)


-----------------------------------------------------

S=QQ[x,y,z]
I=ideal(x,y,z)

H=localCohomFW(I,3,1,3+3)
for  i from 0 to 10 do print hilbertFunction(i,H)


-----------------------------------------------------

--rational quartic

S = QQ[x,y,z,w]
R = QQ[s,t]
F=map(R,S,{s^4, s^3*t, s*t^3, t^4})
I=kernel F
codim I

localCohomFW(I,2,0,2+4)



--------------------------------------------------

restart
load "MixedHodgeModules.m2"
needsPackage "Complexes"


--------------------------------------------

S = QQ[x,y,z,w]
f = x*w-y*z
for p from 0 to 3 do print (prune HH_(-(4-p-1))(intersectionDuBoisComplex(f,p)))
