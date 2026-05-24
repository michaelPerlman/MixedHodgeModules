end


---------------------------------------------------------------
---------------------------------------------------------------
	

restart
load "MixedHodgeModules.m2"



--------------------------------------------------------------
--smooth case

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=y^2+x
w={1/1,1/2}

elapsedTime hodgeOnV(f,0)--.04s
elapsedTime hodgeOnV(f,1)--.04s
elapsedTime hodgeOnV(f,2)--.04s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.04s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))--.30s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--.17s


---------------------------------------------------------------

--normal crossing


restart
load "MixedHodgeModules.m2"

n=3
R=QQ[x_1..x_n]
f=product gens R

elapsedTime hodgeOnV(f,0)--.03s
elapsedTime hodgeOnV(f,1)--.05s
elapsedTime hodgeOnV(f,2)--.07s
elapsedTime hodgeOnV(f,3)--.12s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.14s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))--11.80s


---------------------------------------------------------------

--2x2 determinant

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z,w]
f=x*w-y*z
u={1/2,1/2,1/2,1/2}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.05s
elapsedTime hodgeOnV(f,2)--.08s
elapsedTime hodgeOnV(f,3)--.12s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.33s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))--31.45s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,u))
elapsedTime netList (for p from 0 to 10 list hodgeIdealDet(2,p))--.04s

---------------------------------------------------------------


--cusp

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=x^2+y^3
w={1/2,1/3}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.05s
elapsedTime hodgeOnV(f,2)--.07s
elapsedTime hodgeOnV(f,3)--.11s



elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.07s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))--22.96s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--.37s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/6,p))--.05s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/6,p, UseGenLevel => False))--2.07s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/6,p,w))--.24s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,5/6,p))--.06s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,5/6,p, UseGenLevel => False))--.81s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,5/6,p,w))--.28s




--Davis--Yang page 4
--Zhang Conjecture E

alpha=5/6+1/10
p=2
hodgeIdeal(f,alpha,p)
higherMultiplierIdeal(f,alpha,p)



---------------------------------------------------------------

--A1 sing

restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z]
f=y^2-x*z
w={1/2,1/2,1/2}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.05s
elapsedTime hodgeOnV(f,2)--.07s
elapsedTime hodgeOnV(f,3)--.10s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.13s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p, UseGenLevel => False))--6.65s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--22.78s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/2,p))--.11s
elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/2,p, UseGenLevel => False))--1.55s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/2,p,w))--14.62s



---------------------------------------------------------------


--MP Remark 17.13

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=x*y*(x+y)

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.08s
elapsedTime hodgeOnV(f,2)--.32s

elapsedTime netList (for p from 0 to 5 list hodgeIdeal(f,1,p))--.03s
elapsedTime netList (for p from 0 to 5 list hodgeIdeal(f,1,p, UseGenLevel => False))

elapsedTime netList (for p from 0 to 5 list hodgeIdeal(f,1/3,p))--.03s
elapsedTime netList (for p from 0 to 5 list hodgeIdeal(f,1/3,p, UseGenLevel => False))

elapsedTime netList (for p from 0 to 5 list hodgeIdeal(f,2/3,p))--.04s
elapsedTime netList (for p from 0 to 1 list hodgeIdeal(f,2/3,p, UseGenLevel => False))


---------------------------------------------------------------

--ordinary m-fold point in PP^{n-1}
--see MP Example 20.1

restart
load "MixedHodgeModules.m2"

--m\leq n/2 case

elapsedTime netList flatten (for n from 2 to 5 list (
	for m from 1 to floor(n/2) list (
		R:=QQ[x_1..x_n];
                f:=sum apply(toList (1..n), i-> x_i^m);
                hodgeIdeal(f,1,1)
	)
)
)--.25s


-- n/2 <= m <= n-1 case

elapsedTime netList flatten (for n from 2 to 3  list (
	for m from ceiling(n/2) to n-1 list (
		R:=QQ[x_1..x_n];
                f:=sum apply(toList (1..n), i-> x_i^m);
                hodgeIdeal(f,1,1)
	)
)
)--.08s

---------------------------------------------------------------

--Blanco example 1

restart
loadPackage "MixedHodgeModules"

R=QQ[x,y]
f=x^5+y^5+x^2*y^2

elapsedTime V0=hodgeOnV(f,0)--.14s
elapsedTime V1=hodgeOnV(f,1)--12.41s

netList apply(keys V0, k -> {k,hodgeIdeal(f,k,0)})
netList apply(keys V1, k -> {k,hodgeIdeal(f,k,1)})


---------------------------------------------------------------

--Blanco example 2

restart
loadPackage "MixedHodgeModules"

R=QQ[x,y]
lambda=1/2
f=(y^2-x^3)*(y^2+lambda*x^3)

elapsedTime V0=hodgeOnV(f,0)--.14s
elapsedTime V1=hodgeOnV(f,1)--

netList apply(keys V0, k -> {k,hodgeIdeal(f,k,0)})
netList apply(keys V1, k -> {k,hodgeIdeal(f,k,1)})



---------------------------------------------------------------


--Blanco example 3

restart
loadPackage "MixedHodgeModules"

R=QQ[x,y,z]
f=x^3+y^3+z^3+x*y*z
hodgeOnV(f,1)

elapsedTime V0=hodgeOnV(f,0)--
elapsedTime V1=hodgeOnV(f,1)--

netList apply(keys V0, k -> {k,hodgeIdeal(f,k,0)})
netList apply(keys V1, k -> {k,hodgeIdeal(f,k,1)})


---------------------------------------------------------------

--3x3 determinant [PR21]

restart
load "MixedHodgeModules.m2"

R=QQ[x_(1,1)..x_(3,3)]
f=determinant genericMatrix(R,x_(1,1),3,3)

elapsedTime hodgeOnV(f,0)--.17s
elapsedTime hodgeOnV(f,1)--.27s
elapsedTime hodgeOnV(f,2)--3.67s


elapsedTime netList (for p from 0 to 3 list hodgeIdeal(f,1,p))--110.11s
elapsedTime netList (for p from 0 to 3 list hodgeIdealDet(3,p))--.13s

---------------------------------------------------------------


--3x3 symmetric determinant
R=QQ[x_1..x_6]
f=determinant genericSymmetricMatrix(R,x_1,3)

elapsedTime hodgeOnV(f,0)--.17s
elapsedTime hodgeOnV(f,1)--.27s
elapsedTime hodgeOnV(f,2)--3.67s


elapsedTime netList (for p from 0 to 3 list hodgeIdeal(f,1,p))

elapsedTime netList (for p from 0 to 3 list hodgeIdeal(f,1/2,p))



---------------------------------------------------------------

--An singularity
--n\geq 1

restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z];
n=2
f=x^2+y^2+z^(n+1);
w={1/2,1/2,1/(n+1)}




---------------------------------------------------------------


--Dn singularity
--n\geq 4

restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z];
n=4
f=x^2+y^(n-1)+y*z^2;
w={1/2,1/(n-1),(n-2)/(2*(n-1))}



---------------------------------------------------------------
    

--E6 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+z^4
w={1/2,1/3,1/4}



---------------------------------------------------------------
    

--E7 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+y*z^3
w={1/2,1/3,2/9}



---------------------------------------------------------------

--E8 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+z^5
w={1/2,1/3,1/5}


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











