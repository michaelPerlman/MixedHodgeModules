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

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=y^2+x

elapsedTime generationLevel(f)--.00s

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=y^2+x

elapsedTime HRHLevel(f)--.04s

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=y^2+x

elapsedTime netList for p from 0 to 1 list duBoisComplex(f,p)


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

restart
load "MixedHodgeModules.m2"

n=3
R=QQ[x_1..x_n]
f=product gens R

elapsedTime generationLevel(f)--.09s

restart
load "MixedHodgeModules.m2"

n=3
R=QQ[x_1..x_n]
f=product gens R

elapsedTime HRHLevel(f)--.03s

restart
load "MixedHodgeModules.m2"

n=3
R=QQ[x_1..x_n]
f=product gens R

elapsedTime netList for p from 0 to 2 list duBoisComplex(f,p)--.33s


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

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z,w]
f=x*w-y*z

elapsedTime generationLevel(f)--.14s

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z,w]
f=x*w-y*z

elapsedTime HRHLevel(f)--.05s

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z,w]
f=x*w-y*z

elapsedTime netList for p from 0 to 3 list duBoisComplex(f,p)--.93s


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


restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=x^2+y^3

elapsedTime generationLevel(f)--.00s

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=x^2+y^3

elapsedTime HRHLevel(f)--.04s

restart
load "MixedHodgeModules.m2"

R=QQ[x,y]
f=x^2+y^3

elapsedTime netList for p from 0 to 1 list duBoisComplex(f,p)--.14s


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
elapsedTime V1=hodgeOnV(f,1)--does not run

netList apply(keys V0, k -> {k,hodgeIdeal(f,k,0)})
netList apply(keys V1, k -> {k,hodgeIdeal(f,k,1)})--does not run



---------------------------------------------------------------


--Blanco example 3

restart
loadPackage "MixedHodgeModules"

R=QQ[x,y,z]
f=x^3+y^3+z^3+x*y*z


elapsedTime V0=hodgeOnV(f,0)--.08s
elapsedTime V1=hodgeOnV(f,1)--98.71s
elapsedTime V2=hodgeOnV(f,2)--did not try yet

netList apply(keys V0, k -> {k,hodgeIdeal(f,k,0)})
netList apply(keys V1, k -> {k,hodgeIdeal(f,k,1)})
netList apply(keys V2, k -> {k,hodgeIdeal(f,k,2)})

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

restart
load "MixedHodgeModules.m2"

R=QQ[x_(1,1)..x_(3,3)]
f=determinant genericMatrix(R,x_(1,1),3,3)

elapsedTime generationLevel(f)--

restart
load "MixedHodgeModules.m2"

R=QQ[x_(1,1)..x_(3,3)]
f=determinant genericMatrix(R,x_(1,1),3,3)

elapsedTime HRHLevel(f)--.16s

restart
load "MixedHodgeModules.m2"

R=QQ[x_(1,1)..x_(3,3)]
f=determinant genericMatrix(R,x_(1,1),3,3)

elapsedTime netList for p from 0 to 1 list duBoisComplex(f,p)--

---------------------------------------------------------------


--3x3 symmetric determinant

restart
load "MixedHodgeModules.m2"

R=QQ[x_1..x_6]
f=determinant genericSymmetricMatrix(R,x_1,3)

elapsedTime hodgeOnV(f,0)--.11s
elapsedTime hodgeOnV(f,1)--.15s
elapsedTime hodgeOnV(f,2)--.94s


elapsedTime netList (for p from 0 to 3 list hodgeIdeal(f,1,p))--9.10s

elapsedTime netList (for p from 0 to 3 list hodgeIdeal(f,1/2,p))--1.94s



---------------------------------------------------------------

--An singularity
--n\geq 1
--f=x^2+y^2+z^(n+1)

restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z]

n=2
f=x^2+y^2+z^(n+1)
w={1/2,1/2,1/(n+1)}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.07s
elapsedTime hodgeOnV(f,2)--.10s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.16s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--51.98s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/3,p))--.13s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/3,p,w))--26.96s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,2/3,p))--.16s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,2/3,p,w))--36.96s



restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z]

n=3
f=x^2+y^2+z^(n+1)
w={1/2,1/2,1/(n+1)}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.07s
elapsedTime hodgeOnV(f,2)--.13s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.17s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--89.62s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/2,p))--.13s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/2,p,w))--61.38s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,3/4,p))--.18s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,3/4,p,w))--63.40s




---------------------------------------------------------------


--Dn singularity
--n\geq 4

restart
load "MixedHodgeModules.m2"

R=QQ[x,y,z];
n=4
f=x^2+y^(n-1)+y*z^2;
w={1/2,1/(n-1),(n-2)/(2*(n-1))}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.11s
elapsedTime hodgeOnV(f,2)--.15s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.27s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--125.63

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/2,p))--.07s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/2,p,w))--71.02s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/6,p))--.17s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/6,p,w))--54.80s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,5/6,p))--.24s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,5/6,p,w))--93.63s


---------------------------------------------------------------
    

--E6 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+z^4
w={1/2,1/3,1/4}

elapsedTime hodgeOnV(f,0)--.05s
elapsedTime hodgeOnV(f,1)--.23s
elapsedTime hodgeOnV(f,2)--3.67s

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--.36s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/12,p))--.25s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/12,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,5/12,p))--.35s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,5/12,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,7/12,p))--.35s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,7/12,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,11/12,p))--.37s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,11/12,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1/3,p))--.31s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1/3,p,w))--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,2/3,p))--.37s
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,2/3,p,w))--



---------------------------------------------------------------
    

--E7 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+y*z^3
w={1/2,1/3,2/9}

elapsedTime hodgeOnV(f,0)--.06s
elapsedTime hodgeOnV(f,1)--.71s
elapsedTime hodgeOnV(f,2)--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--

--add other roots

---------------------------------------------------------------

--E8 singularity

restart
load "MixedHodgeModules.m2"

S=QQ[x,y,z]
f=x^2+y^3+z^5
w={1/2,1/3,1/5}

elapsedTime hodgeOnV(f,0)--.07s
elapsedTime hodgeOnV(f,1)--.36s
elapsedTime hodgeOnV(f,2)--

elapsedTime netList (for p from 0 to 10 list hodgeIdeal(f,1,p))--
elapsedTime netList (for p from 0 to 10 list hodgeIdealWeightedHomogIsolated(f,1,p,w))--

--add other roots


---------------------------------------------------------------

--Saito Example

restart
load "MixedHodgeModules.m2"

S=QQ[x,y]
f=x^(14)+y^(14)-x^6*y^6

elapsedTime hodgeIdeal(f,1,0)--


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

