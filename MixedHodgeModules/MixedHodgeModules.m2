--Copyright 2026 by Andras Lorincz and Michael Perlman

--Important: to use this package one must do:
--(1) maintain the current file structure:
-- MixedHodgeModules(folder) -> {MixedHodgeModules (folder), MixedHodgeModules.m2 (this file)} -> {aux files}
--(2) launch the package from this file, using code after "end"



--to do:
--(1) Fill in examples for p-function and nuAlpha, and mention strategies




newPackage(
    "MixedHodgeModules",
    Version => "1.0",
    Date => "September 22, 2026",
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
    "localCohomFW",
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

    "clearMHMCache",

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

--internal symbols for hodgeIdealDet session-global caches (stored on detSessionRing.cache)
protect symbol DetGenericMatrixRingCache;
protect symbol ILambdaDetCache;
protect symbol HodgeIdealDetCache;





---------------------------------------------------------------
--V-filtration computation kernel and helpers:
--hodgeOnV, weightHodgeOnV, HRHCheck, HRHLevel
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



