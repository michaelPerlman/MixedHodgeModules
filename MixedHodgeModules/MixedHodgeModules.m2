--Copyright 2026 by Andras Lorincz and Michael Perlman

--Important: to use this package one must do:
--(1) maintain the current file structure:
-- MixedHodgeModules(folder) -> {MixedHodgeModules (folder), MixedHodgeModules.m2 (this file)} -> {aux files}
--(2) launch the package from this file, using code after "end"


--recent changes:
--(1) added HRHCheck and HRHLevel
--(2) added localCohomWF
--(3) added de Rham and Du Bois functions

--to do:
--(1) monodromy weight filtration on F_p(Gr_V^{alpha}(B_f)). Should be easy to add
--(2) fill in descriptions and examples in documentation for nuAlpha, pFunction (and their strategies)
--(3) add examples in doc of other strategies for weightLength
--(4) create tests for every function and add to tests.m2
--(5) determine which strategy for weightLength should be default.
--(6) There is an old comment on pFunction "to do: modify so don't compute high nus". Action needed?
--(7) There is an old comment on weightCheck "possible bug here, use percent instead?". Action needed?
--(8) un-export weightHodgeOnV once monodromy is implemented?
--(9) document de Rham and Du Bois functions


--change/add soon:
--(1) change weightCheck and weightLevel to input alpha\in [0,1) and k\in ZZ? (instead of beta=alpha+k)
--same for weightLength? Same for hodgeLevel, hodgeCheck?
--(2) allow alpha as a parameter in hodgeIdeal (work over QQ(alpha) where alpha is variable)?
--(3) weighted Homogeneous isolated option in de Rham and Du Bois functions?


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
    Date => "April 24, 2026",
    Headline => "Calculations involving Hodge and weight filtrations on localizations",
    Authors => {{ Name => "András C. Lőrincz",  Email => "lorincz@ou.edu",  HomePage => "https://math.ou.edu/~alorincz/"},
	        { Name => "Michael Perlman",    Email => "mperlman@ua.edu", HomePage => "https://sites.google.com/view/michaelperlman/home"}},
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
    "nuAlpha",
    "pFunction",
    "weightCheck",
    "weightLength",
    "weightLevel",
    "weightedHodgeIdeal",
    "weightHodgeOnV",

    "compareHodgeIequalHigherMultI",
    "testBrieskornPham",
    "testHodgeIdealforWeightedHomog",
    "testHodgeIhigherMultImodf",
    "testWeightedHodge",

    "gradedDuBoisComplex",
    "intersectionDuBoisComplex",
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
--functions to facilitate testing
---------------------------------------------------------------
load "./MixedHodgeModules/testFunctions.m2"

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


needsPackage "SchurComplexes"
R=QQ[x,y,z]
f=y^2-x-z
DB=gradedDuBoisComplex(f,1)
M=prune HH_0 DB
F=res M
F.dd
G=schurComplex({1,1},F)
G.dd
prune HH_0(G)
prune HH_1(G)
prune HH_2(G)


R=QQ[x,y,z,w]
f=x*w-y*z
DB=intersectionDuBoisComplex(f,3)
dd^DB
prune HH_(-1) DB
for i from -1 to 10 do print hilbertFunction(i,o31)


--------------------------------------------------------------
--smooth case
R=QQ[x,y]
f=y^2+x
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


p=3
compareHodgeIequalHigherMultI(f,p)--equal by Lorincz--Yang

testHodgeIhigherMultImodf(f,p)

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
alpha=1
hodgeIdeal(f,alpha,0)
hodgeIdeal(f,alpha,1)
hodgeIdeal(f,alpha,2)
hodgeIdeal(f,alpha,3)
hodgeIdeal(f,alpha,4)

for p from 0 to 4 do print doesGenerateNext(f,1,p)
generationLevel(f)
p=5
u={1/2,1/2,1/2,1/2}
testHodgeIdealforWeightedHomog(f,p,u, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,u)

p=3
m=2
testWeightedHodge(f,p,m)
adjointIdeal(f)


hodgeOnV(f,2, UseBasis => sBasis)
weightHodgeOnV(f,1,2,1, UseBasis => sBasis)
ideal mingens ideal apply(apply(oo, g -> (sub(f,ring g))^2*sub(g, {(ring g)_0 => -1})), h -> sub(h,ring f))
weightHodgeOnV(f,1,2,1)
weightHodgeOnV(f,1,2,1, UseBasis => dtBasis)
hodgeOnV(f,1,1)


p=3
compareHodgeIequalHigherMultI(f,p)--equal by Lorincz--Yang
--see also Zhang Conjecture E

testHodgeIhigherMultImodf(f,p)

---------------------------------------------------------------


--cusp
R=QQ[x,y]
f=x^2+y^3
hodgeOnV(f,1) 
alpha=1
generationLevel(f,alpha)
w={1/3,1/2}
p=4
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)
testBrieskornPham({2,3},p)
p=3
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)


p=3
compareHodgeIequalHigherMultI(f,p)

testHodgeIhigherMultImodf(f,p)

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
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

p=3
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

p=3
compareHodgeIequalHigherMultI(f,p)--equal by Lorincz--Yang
--see also Zhang Conjecture E

testHodgeIhigherMultImodf(f,p)

---------------------------------------------------------------


--Brieskorn Pham
L={2,2,2,2,2}
p=2
testBrieskornPham(L,p, UseGenLevel => False)
testBrieskornPham(L,p)


---------------------------------------------------------------

--cusp intersect line
R=QQ[x,y]
f=x*(y^2-x^3)
hodgeOnV(f,1)
for p from 0 to 2 do print hodgeIdeal(f,1,p)

ll=weightLength(f,1)
p=1
m=ll
testWeightedHodge(f,p,m)
adjointIdeal(f)


p=2
compareHodgeIequalHigherMultI(f,p)

testHodgeIhigherMultImodf(f,p)

---------------------------------------------------------------


--MP Remark 17.13
R=QQ[x,y]
f=x*y*(x+y)
hodgeOnV(f,1)
alpha=1
for p from 0 to 2 do print hodgeIdeal(f,alpha,p)

p=0
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

p=2
compareHodgeIequalHigherMultI(f,p)

testHodgeIhigherMultImodf(f,p)

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
alpha=1
for p from 0 to 1 do print hodgeIdeal(f,alpha,p)

p=0
m=1
testWeightedHodge(f,p,m)--slow
adjointIdeal(f)

p=1
compareHodgeIequalHigherMultI(f,p)--equal by Lorincz--Yang

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
testWeightedHodge(f,p,m)--slow, but runs
adjointIdeal(f)

p=1
compareHodgeIequalHigherMultI(f,p)--equal by Lorincz--Yang

testHodgeIhigherMultImodf(f,p)
---------------------------------------------------------------

--Blanco example 1
R=QQ[x,y]
f=x^5+y^5+x^2*y^2
V=hodgeOnV(f,0)

for a in keys V do (
     for p from 0 to 0 do print {a,hodgeIdeal(f,a,p)})

 
p=0
m=1
testWeightedHodge(f,p,m)
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
testWeightedHodge(f,p,m)
adjointIdeal(f)



---------------------------------------------------------------


--Blanco example 3
R=QQ[x,y,z]
f=x^3+y^3+z^3+x*y*z
hodgeOnV(f,0)--doesnt finish on my computer

w={1/3,1/3,1/3}
p=0
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

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
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)
testBrieskornPham({2,2,n+1},p)


p=0
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)


p=2
compareHodgeIequalHigherMultI(f,p)
--see also Zhang Conjecture E

testHodgeIhigherMultImodf(f,p)

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
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

p=0
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

p=2
compareHodgeIequalHigherMultI(f,p)
--see also Zhang Conjecture E

testHodgeIhigherMultImodf(f,p)

---------------------------------------------------------------
    


--E6 singularity
S=QQ[x,y,z]
f=x^2+y^3+z^4
hodgeOnV(f,1) 
w={1/2,1/3,1/4}
p=1
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

p=0
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

testHodgeIhigherMultImodf(f,p)
---------------------------------------------------------------
    

--E7 singularity
S=QQ[x,y,z]
f=x^2+y^3+y*z^3
globalBFunction(f)
factorBFunction oo
alpha=1
w={1/2,1/3,2/9}
p=1
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

p=0
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

testHodgeIhigherMultImodf(f,p)
---------------------------------------------------------------

--E8 singularity
S=QQ[x,y,z]
f=x^2+y^3+z^5
globalBFunction(f)
factorBFunction oo
alpha=1
w={1/2,1/3,1/5}
p=1
testHodgeIdealforWeightedHomog(f,p,w, UseGenLevel => False)
testHodgeIdealforWeightedHomog(f,p,w)

p=1
m=1
testWeightedHodge(f,p,m)
adjointIdeal(f)

testHodgeIhigherMultImodf(f,p)

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


testHodgeIhigherMultImodf(f,p)


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
