end


---------------------------------------------------------------
---------------------------------------------------------------
	

restart
load "MixedHodgeModules.m2"



--tests/ examples


--------------------------------------------------------------
--smooth case
R=QQ[x,y]
f=y^2+x

hodgeOnV(f,0)
hodgeOnV(f,1)
hodgeOnV(f,2)

netList (for p from 0 to 10 list hodgeIdeal(f,1,p))



---------------------------------------------------------------

--normal crossing
n=3
R=QQ[x_1..x_n]
f=product gens R

hodgeOnV(f,0)
hodgeOnV(f,1)
hodgeOnV(f,2)
hodgeOnV(f,3)

netList (for p from 0 to 10 list hodgeIdeal(f,1,p))


---------------------------------------------------------------

--2x2 determinant
S=QQ[x,y,z,w]
f=x*w-y*z

hodgeOnV(f,0)
hodgeOnV(f,1)
hodgeOnV(f,2)
hodgeOnV(f,3)

netList (for p from 0 to 10 list hodgeIdeal(f,1,p))


---------------------------------------------------------------


--cusp
R=QQ[x,y]
f=x^2+y^3

hodgeOnV(f,0)
hodgeOnV(f,1)
hodgeOnV(f,2)
hodgeOnV(f,3)

w={1/2,1/3}

netList (for p from 0 to 10 list hodgeIdeal(f,1,p))
netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))
netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))

netList (for p from 0 to 10 list hodgeIdeal(f,1/6,p))
netList (for p from 0 to 10 list hodgeIdeal(f,1/6,p, UseGenLevel => False))
netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/6,p,w))

netList (for p from 0 to 10 list hodgeIdeal(f,5/6,p))
netList (for p from 0 to 10 list hodgeIdeal(f,5/6,p, UseGenLevel => False))
netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,5/6,p,w))




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

p=3
m=1
adjointIdeal(f)


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
alpha=1
elapsedTime hodgeOnV(f,alpha,2)
elapsedTime hodgeIdeal(f,alpha,2);--10 sec
elapsedTime hodgeIdeal(f,alpha,3);--281 sec
elapsedTime hodgeIdeal(f,alpha,4);



for p from 0 to 2 do print hodgeIdeal(f,alpha,p)
hodgeIdeal(f,alpha,3)

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
restart
loadPackage "MixedHodgeModules"
R=QQ[x,y]
f=x^5+y^5+x^2*y^2
elapsedTime hodgeOnVNew(f,1,1)
elapsedTime V=hodgeOnV(f,1,1)

for a in keys oo do (
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
V=hodgeOnV(f,1) 
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
hodgeOnV(f,1)

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
restart
loadPackage "MixedHodgeModules"
elapsedTime hodgeOnVNew(f,1,3)
elapsedTime V=hodgeOnV(f,1,2)
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
elapsedTime hodgeIdeal(f,alpha,0)--wont run


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

restart
load "MixedHodgeModules.m2"
installPackage "MixedHodgeModules"
viewHelp MixedHodgeModules
check "MixedHodgeModules"
uninstallPackage "MixedHodgeModules"


-------------------------------------


S=QQ[x,y,z]
f=y^2-x-z
p=-3
DR=gradedDeRhamComplexH1(f,p)

O0=gradedDuBoisComplex(f,0)
O1=gradedDuBoisComplex(f,1)
O2=gradedDuBoisComplex(f,2)

HH_(0) O0
prune HH_(0) O1
prune HH_(0) O2


prune HH_0 DR
prune HH_1 DR
prune HH_2 DR
prune HH_3 DR

p=-3
IC=intersectionDuBoisComplex(f,p)
prune HH_0 IC
prune HH_1 IC
prune HH_2 IC
prune HH_3 IC




S=QQ[x,y,z,w]
f=x*w-y*z
p=-3
DR=gradedDeRhamComplexH1(f,p)



prune HH_0 DR
prune HH_1 DR
prune HH_2 DR
prune HH_3 DR
prune HH_4 DR

S=QQ[x,y,z,w]
f=x*w-y*z
p=1
IC=intersectionDuBoisComplex(f,p)
prune HH_0 IC
prune HH_1 IC
prune HH_2 IC
prune HH_3 IC
prune HH_4 IC

O0=gradedDuBoisComplex(f,0)
O1=gradedDuBoisComplex(f,1)
O2=gradedDuBoisComplex(f,2)
O3=gradedDuBoisComplex(f,3)

IC0=intersectionDuBoisComplex(f,0)
IC1=intersectionDuBoisComplex(f,1)
IC2=intersectionDuBoisComplex(f,2)
IC3=intersectionDuBoisComplex(f,3)

HH_(0) O0
prune HH_(0) O1
prune HH_(0) O2
prune HH_0 O3


S=QQ[x,y]
f=x^2+y^3
p=-2
DR=gradedDeRhamComplexH1(f,p)
prune HH_0 DR
prune HH_1 DR
prune HH_2 DR

O0=gradedDuBoisComplex(f,0)
O1=gradedDuBoisComplex(f,1)
O2=gradedDuBoisComplex(f,2)

prune HH_(0) O0
prune HH_(0) O1
prune HH_0 O2

IC=intersectionDuBoisComplex(f,p)
KK= S^1/ideal(gens S)



prune HH_0 IC
prune HH_1 IC
prune HH_2 IC
prune HH_3 IC



-----------------------------------------------------------------

S=QQ[x,y]
f=x^2+y^3

for p from -2 to 0 do (
    for q from 0 to 2 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))

for p from -2 to 0 do (
    for q from 0 to 2 do (
	print {p,q,GRFdeRhamIC(f,p,q)}
	))


------------------------------------

S=QQ[x,y,z,w]
f=x*w-y*z

deRhamInterval(f,{-1},-1, InputType => HodgeIdeals)
GRFdeRhamH1f(f,0,1)

for p from -4 to 0 do (
    for q from 0 to 4 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))

for p from -4 to 0 do (
    for q from 0 to 4 do (
	print {p,q,GRFdeRhamIC(f,p,q)}
	))

------------------------------------


S=QQ[x,y,z,w]
f=x^2+y+z+w
for p from -4 to 0 do (
    for q from 0 to 4 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))

for p from -4 to 0 do (
    for q from 0 to 4 do (
	print {p,q,GRFdeRhamIC(f,p,q)}
	))


---------------------------
---

S=QQ[x,y,z]
f=x^2+y+z
for p from -3 to 4 do (
    for q from 0 to 3 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))


S=QQ[x,y,z]
f=x^2+y+z
for p from -3 to 4 do (
    for q from 0 to 3 do (
	print {p,q,GRFdeRhamSf(f,p,q)}
	))


------------------------------------

S=QQ[x,y,z]
f=y^2-x*z


for p from -3 to 0 do (
    for q from 0 to 3 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))


for p from -3 to 0 do (
    for q from 0 to 3 do (
	print {p,q,GRFdeRhamIC(f,p,q)}
	))


------------------------------------


--3x3 symmetric determinant
R=QQ[x_1..x_6]
f=determinant genericSymmetricMatrix(R,x_1,3)

for p from -6 to -3 do (
    for q from 0 to 6 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))


for p from -6 to -3 do (
    for q from 0 to 6 do (
	print {p,q,GRFdeRhamIC(f,p,q)}
	))


------------------------------------

S=QQ[x,y]
f=x^2+y
p=-1


DR = GRFdeRhamComplex (f,p)

DRp = part({p+2,0},DR)


freeResolution(DRp)



for p from -2 to 0 do (
    for q from 0 to 2 do (
	print {p,q,GRFdeRhamH1f(f,p,q)}
	))


