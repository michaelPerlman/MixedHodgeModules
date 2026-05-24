--Copyright 2026 by Andras Lorincz and Michael Perlman

--Important: to use this package one must do:
--(1) maintain the current file structure:
-- MixedHodgeModules(folder) -> {MixedHodgeModules (folder), MixedHodgeModules.m2 (this file)} -> {aux files}
--(2) launch the package from this file, using code after "end"



--to do:
--(1) fill in descriptions and examples in documentation for nuAlpha, pFunction (and their strategies)
--(2) add examples in doc of other strategies for weightLength
--(3) create tests for every function and add to tests.m2
--(4) Bugs to fix (discovered while writing tests):
--    (a) adjointIdeal(f) and other weight-filtration functions internally do
--        `use ring If` (the Weyl algebra), leaving the user's current ring as
--        the Weyl algebra after the call.  The output ideal lives in R (the
--        original ring), but a subsequent `ideal(b, a)` from the user is now
--        constructed in the Weyl algebra, so equality tests with the output
--        fail with "expected ideals for the same ring."  Fix: restore the
--        caller's current ring (e.g. `use R` at the end) and audit all
--        downstream code that depends on this behavior.  Affects at least:
--        adjointIdeal, weightLevel (the `use ring If` at line ~269 of
--        weightFiltrations.m2), hodgeIdealWeightedHomogIsolated, the WHI
--        helper, and probably others.
--    (b) Calling duBoisComplex(f, 0) then duBoisComplex(f, 1) in the same
--        session yields a DB1 with an extra direct summand:
--        cokernel diag(f, f) instead of cokernel | f |.  Cohomology is the
--        same up to isomorphism but the presentation is wrong.  Likely a
--        mutable-cache issue in koszulSlice's Ktwisted, freeResolution, or
--        Hom.  Reproduce: see tests #28-#33 cluster.
--    (c) Calling hodgeLevel(f, g, alpha) then hodgeCheck(f, g, alpha, p) in
--        the same session errors with "expected pair to have a method for
--        '*'".  Likely related to (a) - hodgeLevel leaves the Weyl algebra
--        as current ring, then hodgeCheck builds `g * f^k` where g is now
--        re-resolved in the wrong ring.



--change/add soon:
--(1) allow alpha as a parameter in hodgeIdeal (work over QQ(alpha) where alpha is variable)?
--(2) There is an old comment on pFunction "to do: modify so don't compute high nus". 



--to add eventually?:
--(1) use minimal exponent formula for generation level as optional strategy
--in Hodge ideal?
--(2) allow user to input generation level as optional strategy?
--(3) get generation level on IC_f using weighted Hodge ideal?
--(4) Bernstein--Sato polynomials on singular ambient varieties X (Dirks) and microlocal b-functions on
-- singular ambient varieties?
--(5) nearby cycles and vanishing cycles? vanishing cycles are implicitly calculated in HRHCheck



newPackage(
    "MixedHodgeModules",
    Version => "1.0",
    Date => "May 22, 2026",
    Headline => "Calculations involving Hodge and weight filtrations on localizations",
    Authors => {{ Name => "András C. Lőrincz",  Email => "lorincz@ou.edu",  HomePage => "https://math.ou.edu/~alorincz/"},
	        { Name => "Michael Perlman",    Email => "mperlman@ua.edu", HomePage => "https://sites.google.com/view/michaelperlman/home"}},
    HomePage => "https://github.com/michaelPerlman/MixedHodgeModules",
    Keywords => {"D-modules"},
    AuxiliaryFiles => true,
    DebuggingMode => false,
    PackageImports => {"Polyhedra", "SchurRings"},
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
    "hodgeIdealDet",
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

--Convention: dt^p \in F_p(B_f)
-- f should be reduced    



--internal symbols for deRham
protect symbol InputType;
protect symbol WeightedHodgeIdeals;
protect symbol HodgeIdeals;
protect symbol WeightedHomogIsolated;

--internal symbol for caching AnnFs on a RingElement
protect symbol AnnFsCache;

--internal symbol for caching globalBFunction on a RingElement
protect symbol GlobalBCache;

--internal symbol for caching rhoFp (roots of generalized b-function with multiplicities) keyed by (f,p)
protect symbol RhoFpCache;

--internal symbol for caching J0 = Jp \cap R[s] keyed by (f,p); populated by
--cachedJ0 and reused by hodgeOnV / weightHodgeOnV
protect symbol J0Cache;

--internal symbol for caching hodgeOnV outputs keyed by (f, alpha, p, UseBasis) or (f, p, UseBasis)
protect symbol HodgeOnVCache;

--internal symbol for caching weightHodgeOnV outputs keyed by (f, alpha, p, m, UseBasis)
protect symbol WeightHodgeOnVCache;

--internal symbol for caching hodgeIdeal outputs keyed by (f, alpha, p, UseGenLevel)
protect symbol HodgeIdealCache;

--internal symbol for caching weightedHodgeIdeal outputs keyed by (f, alpha, p, m)
protect symbol WeightedHodgeIdealCache;

--internal symbol for caching hodgeIdealWeightedHomogIsolated outputs keyed by (f, alpha, p, w)
protect symbol HodgeIdealWeightedHomogIsolatedCache;

--internal symbol for caching makeWeylAlgebra on a polynomial ring
protect symbol WeylAlgebraCache;

--internal symbol for caching polynomialAnnihilator g on a RingElement (with W from cachedWeylAlgebra)
protect symbol PolyAnnCache;





---------------------------------------------------------------
--V-filtration computation kernel and helpers:
--hodgeOnV, weightHodgeOnV, monodromyWeightHodgeOnV, HRHCheck, HRHLevel
---------------------------------------------------------------
load "./MixedHodgeModules/vFiltrations.m2"

---------------------------------------------------------------
--Hodge ideals: hodgeIdeal, generateNext, doesGenerateNext, generationLevel,
--hodgeIdealWeightedHomogIsolated, hodgeIdealBrieskornPham,
--higherMultiplierIdeal, hodgeCheck, hodgeLevel, hodgeIdealDet
---------------------------------------------------------------
load "./MixedHodgeModules/hodgeFiltrations.m2"

---------------------------------------------------------------
--Weighted Hodge ideals and weight-filtration invariants:
--weightedHodgeIdeal, adjointIdeal, IHmoduleAdjoint, weightCheck, weightLevel,
--weightLength, localCohomFW
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
installPackage ("MixedHodgeModules", RemakeAllDocumentation => true)
viewHelp MixedHodgeModules
check "MixedHodgeModules"
uninstallPackage "MixedHodgeModules"


