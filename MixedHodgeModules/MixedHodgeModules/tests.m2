-- Copyright 2026 by Andras Lorincz and Michael Perlman




-* Test section *-



---------------------------------------------------------------
--1. Filtrations on graph module
---------------------------------------------------------------



TEST ///
 R = QQ[x,y];
 f = y^2+x;
 assert(keys(hodgeOnV(f,0)) == {1_QQ});
 m = ((values hodgeOnV(f,0))_0)_0;
 assert( m == 1_(ring m));
 assert(keys(hodgeOnV(f,1)) == {1_QQ});
 assert( #( (values hodgeOnV(f,1))_0 ) == 2);
 assert(HRHCheck(f,0));
 assert(HRHCheck(f,1));
 assert(HRHCheck(f,2));
 assert(HRHLevel(f) == "rational homology manifold");
 assert(sub((weightHodgeOnV(f,1,0,0))_0,R) == (R_1^2+R_0));
 assert(#monodromyWeightHodgeOnV(f,1,1,0) == 2);
///

---------------------------------------------------------------


TEST ///
 R = QQ[x,y,z];
 f = x*y*z;
 assert(keys(hodgeOnV(f,0)) == {1_QQ});
 m = ((values hodgeOnV(f,0))_0)_0;
 assert( m == 1_(ring m));
 assert(not HRHCheck(f,0));
 assert(HRHLevel(f) == -1);
 assert(#monodromyWeightHodgeOnV(f,1,1,0) == 7);
 ///

 ---------------------------------------------------------------


 TEST ///
  R = QQ[x,y,z,w];
  f = x*y-z*w;
  assert(keys(hodgeOnV(f,0)) == {1_QQ});
  m = ((values hodgeOnV(f,0))_0)_0;
  assert( m == 1_(ring m));
  assert(HRHCheck(f,0));
  assert(not HRHCheck(f,1));
  assert(HRHLevel(f) == 0);
  assert(#monodromyWeightHodgeOnV(f,1,1,0) == 5);
 ///

---------------------------------------------------------------

 TEST ///
  R = QQ[x,y,z];
  f = x*z-y^2;
  assert(keys(hodgeOnV(f,0)) == {1_QQ});
  assert(keys hodgeOnV(f,1) == {1_QQ,1/2});
  assert(#monodromyWeightHodgeOnV(f,1/2,1,0) == 2);
  assert(#monodromyWeightHodgeOnV(f,1,1,0) == 4);
  assert(HRHCheck(f,2));
  assert(HRHLevel(f) == "rational homology manifold");
 ///



---------------------------------------------------------------
--2. Hodge ideal, weighted Hodge ideal, higher multiplier ideal tests
---------------------------------------------------------------

---------------------------------------------------------------
--2.1 Smooth and normal crossing examples
---------------------------------------------------------------

TEST ///
 R = QQ[x];
 f = x;
 alpha = 1;
 w = {1};
  --hodgeIdeal tests
 assert(keys(hodgeOnV(f,0)) == {1_QQ});
 assert(generationLevel(f) == 0);
 assert(hodgeIdeal(f,alpha,0) == ideal(1_R));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,0,w) == ideal(1_R));
 assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == ideal(1_R));
 assert(hodgeIdeal(f,alpha,3) == ideal(1_R));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,w) == ideal(1_R));
 --weightedHodgeIdeal tests
 assert(adjointIdeal(f) == ideal(1_R));
 assert(weightedHodgeIdeal(f,alpha,1,0) == ideal(f^2));
 assert(weightedHodgeIdeal(f,alpha,2,0) == ideal(f^3));
 --higherMultiplierIdeal tests
 assert(higherMultiplierIdeal(f,alpha,1) == ideal(1_R));
 assert(higherMultiplierIdeal(f,alpha,2) == ideal(1_R));
///

---------------------------------------------------------------

TEST ///
 R = QQ[x,y];
 f = y^2+x;
 w={1,1/2};
 alpha = 1/2;
  --hodgeIdeal tests
 assert(keys(hodgeOnV(f,0)) == {1_QQ});
 assert(generationLevel(f,alpha) == 0);
 assert(generationLevel(f) == 0);
 assert(hodgeIdeal(f,alpha,0) == ideal(1_R));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,0,w) == ideal(1_R));
 assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == ideal(1_R));
 assert(hodgeIdeal(f,alpha,3) == ideal(1_R));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,w) == ideal(1_R));
 --weightedHodgeIdeal tests
 assert(adjointIdeal(f) == ideal(1_R));
 assert(weightedHodgeIdeal(f,alpha,1,0) == ideal(1_R));
 assert(weightedHodgeIdeal(f,alpha,2,0) == ideal(1_R));
  --higherMultiplierIdeal tests
 assert(higherMultiplierIdeal(f,alpha,1) == ideal(1_R));
 assert(higherMultiplierIdeal(f,alpha,2) == ideal(1_R));
///

---------------------------------------------------------------

TEST ///
  R=QQ[x,y,z]
  f=x*y*z;
  alpha = 1;
  --hodgeIdeal tests
  assert(generationLevel(f) == 0);
  assert(keys(hodgeOnV(f,0)) == {1_QQ});
  assert(hodgeIdeal(f,alpha,0) == ideal(1_R));
  assert(hodgeIdeal(f,alpha,1) == sub(ideal(x*y, x*z, y*z),R));
  assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == sub(ideal(y^3*z^3,x*y^2*z^3,x^2*y*z^3,x^3*z^3,x*y^3*z^2,x^2*y^2*z^2,x^3*y*z^2,x^2*y^3*z,x^3*y^2*z,x^3*y^3),R));
  assert(hodgeIdeal(f,alpha,3) == sub(ideal(y^3*z^3,x*y^2*z^3,x^2*y*z^3,x^3*z^3,x*y^3*z^2,x^2*y^2*z^2,x^3*y*z^2,x^2*y^3*z,x^3*y^2*z,x^3*y^3),R));
  --weighted Hodge ideal tests
  assert(weightedHodgeIdeal(f,alpha,0,0) == sub(ideal(f),R));
  assert(weightedHodgeIdeal(f,alpha,2,0) == sub(ideal(f^3),R));
  assert(adjointIdeal(f) == sub(ideal(x*y, x*z, y*z),R));
  assert(weightedHodgeIdeal(f,alpha,1,1) == sub(ideal(x^2*y^2, x^2*z^2, y^2*z^2),R));
  assert(weightedHodgeIdeal(f,alpha,2,1) == sub(ideal(y^3*z^3,x^3*z^3,x^3*y^3),R));
  assert(weightedHodgeIdeal(f,alpha,2,2) == sub(ideal(y^2*z^3,x*y*z^3,x^2*z^3,y^3*z^2,x^3*z^2,x*y^3*z,x^3*y*z,x^2*y^3,x^3*y^2),R));
  assert(weightedHodgeIdeal(f,alpha,2,3) == hodgeIdeal(f,alpha,2));
  --higherMultiplierIdeal tests (equal to Hodge ideals by Lorincz--Yang)
  assert(higherMultiplierIdeal(f,alpha,0) == ideal(1_R));
  assert(higherMultiplierIdeal(f,alpha,1) == sub(ideal(x*y, x*z, y*z),R));
  assert(higherMultiplierIdeal(f,alpha,3) == sub(ideal(y^3*z^3,x*y^2*z^3,x^2*y*z^3,x^3*z^3,x*y^3*z^2,x^2*y^2*z^2,x^3*y*z^2,x^2*y^3*z,x^3*y^2*z,x^3*y^3),R));  
  ///


---------------------------------------------------------------
--2.2 Weighted homogeneous isolated examples
---------------------------------------------------------------

TEST ///
 S = QQ[x,y,z,w];
 f = x*y-z*w;
 alpha = 1;
 --hodgeIdeal tests
 assert(not doesGenerateNext(f,alpha,0));
 assert(doesGenerateNext(f,alpha,1));
 assert(generationLevel(f,alpha) == 1);
 assert(hodgeIdeal(f,alpha,0) == ideal(1_S));
 assert(hodgeIdeal(f,alpha,1) == ideal(1_S));
 assert(hodgeIdeal(f,alpha,2, UseGenLevel => False) == sub(ideal(x,y,z,w),S));
 assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == (sub(ideal(x,y,z,w),S))^2);
 u={1/2,1/2,1/2,1/2};
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,0,u) == ideal(1_S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,1,u) == ideal(1_S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,2,u) == sub(ideal(x,y,z,w),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,u) == (sub(ideal(x,y,z,w),S))^2);
 --weighted Hodge ideal tests
 assert(weightedHodgeIdeal(f,alpha,0,0) == sub(ideal(f),S));
 assert(weightedHodgeIdeal(f,alpha,1,0) == sub(ideal(f^2),S));
 assert(weightedHodgeIdeal(f,alpha,0,1) == ideal(1_S));
 assert(weightedHodgeIdeal(f,alpha,1,1) == sub(ideal(x,y,z,w),S));
 assert(weightedHodgeIdeal(f,alpha,2,1) == sub(ideal(w^2,y*w,x*w,z^2,y*z,x*z,y^2,x*y+z*w,x^2),S));
 --higher multiplier ideal tests
 assert(higherMultiplierIdeal(f,alpha,0) == ideal(1_S));
 assert(higherMultiplierIdeal(f,alpha,1) == ideal(1_S));
 assert(higherMultiplierIdeal(f,alpha,2) == sub(ideal(x,y,z,w),S));
 assert(higherMultiplierIdeal(f,alpha,3) == (sub(ideal(x,y,z,w),S))^2);
///

---------------------------------------------------------------

TEST ///
 S = QQ[x,y,z];
 f = x*z-y^2;
 alpha = 1/2;
  --hodgeIdeal tests
 assert(not doesGenerateNext(f,alpha,0));
 assert(doesGenerateNext(f,alpha,1));
 assert(generationLevel(f,alpha) == 1);
 assert(hodgeIdeal(f,alpha,0) == ideal(1_S));
 assert(hodgeIdeal(f,alpha,1) == ideal(1_S));
 assert(hodgeIdeal(f,alpha,2, UseGenLevel => False) == sub(ideal(x,y,z),S));
 assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == (sub(ideal(x,y,z),S))^2);
 u={1/2,1/2,1/2};
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,0,u) == ideal(1_S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,1,u) == ideal(1_S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,2,u) == sub(ideal(x,y,z),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,u) == (sub(ideal(x,y,z),S))^2);
 --weighted Hodge ideal tests
 assert(weightedHodgeIdeal(f,alpha,0,0) == ideal(1_S));
 assert(weightedHodgeIdeal(f,alpha,1,0) == sub(ideal(x,y,z),S));
 assert(weightedHodgeIdeal(f,alpha,0,1) == ideal(1_S)); 
 assert(weightedHodgeIdeal(f,alpha,1,1) == ideal(1_S));
 assert(weightedHodgeIdeal(f,alpha,2,1) == sub(ideal(x,y,z),S));
 --higher multiplier ideal tests
 assert(higherMultiplierIdeal(f,alpha,0) == ideal(1_S));
 assert(higherMultiplierIdeal(f,alpha,1) == ideal(1_S));
 assert(higherMultiplierIdeal(f,alpha,2) == sub(ideal(x,y,z),S));
 assert(higherMultiplierIdeal(f,alpha,3) == (sub(ideal(x,y,z),S))^2);
///

---------------------------------------------------------------

TEST ///
 S = QQ[x,y];
 f = x^2+y^3;
 alpha = 1;
 --hodgeIdeal tests
 assert(generationLevel(f,alpha) == 0);
 assert(hodgeIdeal(f,alpha,0) == sub(ideal(x,y),S));
 assert(hodgeIdeal(f,alpha,1, UseGenLevel => False) == sub( ideal(x*y,x^2,y^3),S));
 assert(hodgeIdeal(f,alpha,2,UseGenLevel => False) == sub(ideal(x^3,y^4-3*x^2*y,x*y^3,x^2*y^2),S));
 assert(hodgeIdeal(f,alpha,3, UseGenLevel => False) == sub( ideal(x*y^4-x^3*y,x^2*y^3-x^4,x^3*y^2,x^4*y,y^6-5*x^4),S));
 assert(hodgeIdeal(f,alpha,1) == sub( ideal(x*y,x^2,y^3),S));
 assert(hodgeIdeal(f,alpha,2) == sub(ideal(x^3,y^4-3*x^2*y,x*y^3,x^2*y^2),S));
 assert(hodgeIdeal(f,alpha,3) == sub( ideal(x*y^4-x^3*y,x^2*y^3-x^4,x^3*y^2,x^4*y,y^6-5*x^4),S));
 w = {1/2,1/3};
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,1,w) == sub( ideal(x*y,x^2,y^3),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,2,w) == sub(ideal(x^3,y^4-3*x^2*y,x*y^3,x^2*y^2),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,w) == sub( ideal(x*y^4-x^3*y,x^2*y^3-x^4,x^3*y^2,x^4*y,y^6-5*x^4),S));
 --weighted Hodge ideal tests
 assert(weightedHodgeIdeal(f,alpha,0,0) == sub(ideal(f),S));
 assert(weightedHodgeIdeal(f,alpha,1,0) == sub(ideal(f^2),S));
 assert(weightedHodgeIdeal(f,alpha,0,1) == sub(ideal(x,y),S));
 assert(weightedHodgeIdeal(f,alpha,1,1) == sub( ideal(x*y,x^2,y^3),S));
 assert(weightedHodgeIdeal(f,alpha,2,1)  == sub(ideal(x^3,y^4-3*x^2*y,x*y^3,x^2*y^2),S));
 --higher multiplier ideal tests
 assert(higherMultiplierIdeal(f,alpha,0) == sub(ideal(y,x),S));
 assert(higherMultiplierIdeal(f,alpha,1) == sub(ideal(x*y,x^2,y^3),S));
 assert(higherMultiplierIdeal(f,alpha,2) == sub(ideal(x^2*y,x^3,x*y^3,y^5),S));
 assert(higherMultiplierIdeal(f,alpha,3) == sub(ideal(x^3*y,x^4,x^2*y^3,x*y^5,y^7),S));
 assert(higherMultiplierIdeal(f,alpha,3)+ideal(f) == hodgeIdeal(f,alpha,3)+ideal(f));
///

---------------------------------------------------------------

TEST ///
 S = QQ[x,y];
 f = x^2+y^3;
 alpha = 5/6;
 --hodgeIdeal tests
 assert(generationLevel(f,alpha) == 0);
 assert(hodgeIdeal(f,alpha,0) == ideal(1_S))
 assert(hodgeIdeal(f,alpha,1) == sub(ideal(x,y^2),S));
 assert(hodgeIdeal(f,alpha,2) == sub( ideal(3*y^3-8*x^2,x*y^2,x^2*y),S));
 assert(hodgeIdeal(f,alpha,3) == sub( ideal(9*x*y^3-8*x^3,x^3*y,x^4,3*y^5-14*x^2*y^2),S));
 w = {1/2,1/3};
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,0,w) == ideal(1_S))
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,1,w) == sub(ideal(x,y^2),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,2,w) == sub( ideal(3*y^3-8*x^2,x*y^2,x^2*y),S));
 assert(hodgeIdealWeightedHomogIsolated(f,alpha,3,w) == sub( ideal(9*x*y^3-8*x^3,x^3*y,x^4,3*y^5-14*x^2*y^2),S));
///



---------------------------------------------------------------
--2.3 Misc. examples
---------------------------------------------------------------


TEST ///
 n = 3;
 R = QQ[x_(1,1)..x_(n,n)];
 M = genericMatrix(R,x_(1,1),n,n);
 f = determinant(M);
 alpha = 1;
 assert(hodgeIdeal(f,1,0) == ideal(1_R));
/// 


---------------------------------------------------------------

---------------------------------------------------------------
--3. nuAlpha, pFunction, weightLevel tests
---------------------------------------------------------------

---------------------------------------------------------------
--3.1 Smooth and normal crossing examples
---------------------------------------------------------------

TEST ///
 R = QQ[x];
 g = 1_R;
 f = x;
 alpha = 1;
 assert(nuAlpha(f,g,alpha) == 1);
 assert(nuAlpha(f,g,alpha,NuMethod => PowerBFunction) == 1);
 assert(nuAlpha(f,g,alpha,NuMethod => Malgrange) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 1);
 assert(not weightCheck(f,g,alpha,0));
 assert(weightLevel(f,g,alpha) == 1);
 assert(first degree (pFunction(f,g,alpha)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => Malgrange)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => PowerBFunction)) == 1); 
///

---------------------------------------------------------------

TEST ///
 R = QQ[x,y];
 g = 1_R;
 f = y^2+x;
 alpha = 1/2;
 assert(nuAlpha(f,g,alpha) == 0);
 assert(nuAlpha(f,g,alpha,NuMethod => PowerBFunction) == 0);
 assert(nuAlpha(f,g,alpha,NuMethod => Malgrange) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 0);
 assert(weightLevel(f,g,alpha) == 0);
 assert(pFunction(f,g,alpha) == 1_ZZ);
 assert(pFunction(f,g,alpha, NuMethod => Malgrange) == 1_ZZ);
 assert(pFunction(f,g,alpha, NuMethod => PowerBFunction) == 1_ZZ); 
 assert(weightLevel(f,g,alpha) == 0);
///

---------------------------------------------------------------

TEST ///
  R=QQ[x,y,z];
  f=x*y*z;
  alpha = 1;
  assert(hodgeIdeal(f,1,2) == sub(ideal(y^2*z^2,x*y*z^2,x^2*z^2,x*y^2*z,x^2*y*z,x^2*y^2),R));
///




---------------------------------------------------------------
--3.2 Weighted homogeneous isolated examples
---------------------------------------------------------------

TEST ///
 S = QQ[x,y,z,w];
 f = x*y-z*w;
 g = 1_S;
 alpha = 1/2;
 assert(weightLevel(f,g,alpha) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByWeightLevel) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 0);
 assert(pFunction(f,g,alpha) == 1);
 assert(pFunction(f,g,alpha, NuMethod => Malgrange) == 1);
 assert(pFunction(f,g,alpha, NuMethod => PowerBFunction) == 1); 
 alpha = 1;
 assert(not weightCheck(f,g,alpha,0));
 assert(weightCheck(f,g,alpha,1));
 assert(weightLevel(f,g,alpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByWeightLevel) == 2);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 2);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 2);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 2);
 assert(first degree (pFunction(f,g,alpha)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => Malgrange)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => PowerBFunction)) == 1);
 alpha = 2;
 assert(not weightCheck(f,g,alpha,1));
 assert(weightCheck(f,g,alpha,2));
 assert(weightLevel(f,g,alpha) == 2);
 g = x;
 assert(not weightCheck(f,g,alpha,0));
 assert(weightCheck(f,g,alpha,1));
 assert(weightLevel(f,g,alpha) == 1);
 assert(weightLevel(f,g,alpha) == 1);
///

---------------------------------------------------------------

TEST ///
 S = QQ[x,y,z];
 f = x*z-y^2;
 g = 1_S;
 alpha = 1;
 assert(not weightCheck(f,g,alpha,0));
 assert(weightLevel(f,g,alpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByWeightLevel) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 1);
 assert(first degree (pFunction(f,g,alpha)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => Malgrange)) == 1);
 assert(first degree (pFunction(f,g,alpha, NuMethod => PowerBFunction)) == 1);
 alpha = 1/2;
 assert(weightLevel(f,g,alpha) == 0);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => Malgrange) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha, NuMethod => PowerBFunction) == 1);
 assert(weightLength(f,alpha, LengthStrategy => ByWeightLevel) == 1);
///

---------------------------------------------------------------

TEST ///
 R = QQ[x,y];
 f = x^2+y^3;
 alpha = 1;
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 1);
 assert(weightLength(f,alpha,LengthStrategy => ByWeightLevel) == 1);
 alpha = 1/6;
 assert(weightLength(f,alpha, LengthStrategy => ByNuAlpha) == 1);
 assert(weightLength(f,alpha, LengthStrategy => ByWeightLevel) == 1);
///

---------------------------------------------------------------


TEST ///
  L = {2,2,2};
  alpha = 1/2;
  I = hodgeIdealBrieskornPham(L,alpha,2);
  R = ring I;
  assert(I == ideal gens R);
///

---------------------------------------------------------------


TEST ///
   L = {1,2};
   alpha = 1/6;
   I = hodgeIdealBrieskornPham(L,1/6,2);
   R = ring I;
   assert(I == ideal(1_R));
///

---------------------------------------------------------------



---------------------------------------------------------------
--3.3 Misc. examples
---------------------------------------------------------------


TEST ///
 n = 3;
 R = QQ[x_(1,1)..x_(n,n)];
 M = genericMatrix(R,x_(1,1),n,n);
 f = determinant(M);
 alpha = 1;
 assert(weightLength(f,alpha,LengthStrategy => ByNuAlpha) == 3);
 assert(weightLength(f,alpha, LengthStrategy => ByWeightLevel) == 3);
 g = x_(1,1)*x_(2,2)-x_(1,2)*x_(2,1);
 assert(not weightCheck(f,g,alpha,0));
 assert(weightCheck(f,g,alpha,1));
 assert(weightLevel(f,g,alpha) == 1);
 alpha = 2;
 g = x_(1,1);
 assert(not weightCheck(f,g,alpha,0));
 assert(not weightCheck(f,g,alpha,1));
 assert(weightCheck(f,g,alpha,2));
 assert(weightLevel(f,g,alpha) == 2);
/// 




---------------------------------------------------------------
--3. de Rham complexes, cohomology, Du Bois complexes, pre-Du Bois tests
---------------------------------------------------------------


TEST ///
  R = QQ[x];
  f = x;
  G = gradedDeRhamComplexH1(f,0)
  assert(dim HH_0(G) == -1);
  assert(dim HH_1(G) == -1);
///


---------------------------------------------------------------
--4. Input validation tests
---------------------------------------------------------------

TEST ///
  R = QQ[x,y];
  f = y^2 + x;

  -- alpha out of (0,1]: each call should error
  assert(try (hodgeOnV(f,0,0);             false) else true);
  assert(try (hodgeOnV(f,-1/2,0);          false) else true);
  assert(try (hodgeOnV(f,3/2,0);           false) else true);
  assert(try (hodgeOnV(f,2,0);             false) else true);
  assert(try (hodgeIdeal(f,0,0);           false) else true);
  assert(try (hodgeIdeal(f,3/2,0);         false) else true);
  assert(try (weightedHodgeIdeal(f,0,0,0); false) else true);
  assert(try (weightedHodgeIdeal(f,2,0,0); false) else true);
  assert(try (weightHodgeOnV(f,0,0,0);     false) else true);
  assert(try (weightHodgeOnV(f,3/2,0,0);   false) else true);
  assert(try (generateNext(f,0,0);         false) else true);
  assert(try (generateNext(f,3/2,0);       false) else true);
  assert(try (doesGenerateNext(f,-1,0);    false) else true);
  assert(try (generationLevel(f,0);        false) else true);
  assert(try (generationLevel(f,3/2);      false) else true);
  assert(try (higherMultiplierIdeal(f,0,0);   false) else true);
  assert(try (higherMultiplierIdeal(f,3/2,0); false) else true);

  -- p negative: each call should error
  assert(try (hodgeOnV(f,1,-1);             false) else true);
  assert(try (hodgeOnV(f,-1);               false) else true);
  assert(try (hodgeIdeal(f,1,-1);           false) else true);
  assert(try (weightedHodgeIdeal(f,1,-1,0); false) else true);
  assert(try (weightHodgeOnV(f,1,-1,0);     false) else true);
  assert(try (generateNext(f,1,-1);         false) else true);
  assert(try (doesGenerateNext(f,1,-1);     false) else true);
  assert(try (higherMultiplierIdeal(f,1,-1);false) else true);

  -- m negative: should error
  assert(try (weightedHodgeIdeal(f,1,0,-1); false) else true);
  assert(try (weightHodgeOnV(f,1,0,-1);     false) else true);

  -- boundary values are accepted (alpha = 1, alpha = 1/2, p = 0, m = 0)
  assert(hodgeIdeal(f,1,0)             == ideal(1_R));
  assert(hodgeIdeal(f,1/2,0)           == ideal(1_R));
  assert(weightedHodgeIdeal(f,1,0,0)   == ideal(f));
  assert(higherMultiplierIdeal(f,1,0)  == ideal(1_R));
///


TEST ///
  R = QQ[x,y];
  f = y^2 + x;

  -- hodgeIdealBrieskornPham bad inputs
  assert(try (hodgeIdealBrieskornPham({},1,0);       false) else true);   -- empty L
  assert(try (hodgeIdealBrieskornPham({0,2},1,0);    false) else true);   -- b_i = 0
  assert(try (hodgeIdealBrieskornPham({-1,2},1,0);   false) else true);   -- b_i < 0
  assert(try (hodgeIdealBrieskornPham({1/2,2},1,0);  false) else true);   -- non-integer b_i
  assert(try (hodgeIdealBrieskornPham({2},3/2,0);    false) else true);   -- alpha out of range
  assert(try (hodgeIdealBrieskornPham({2},1,-1);     false) else true);   -- p < 0

  -- hodgeIdealWeightedHomogIsolated bad inputs
  assert(try (hodgeIdealWeightedHomogIsolated(f,1,0,{1});        false) else true); -- #w mismatch
  assert(try (hodgeIdealWeightedHomogIsolated(f,1,0,{1/2,1/2});  false) else true); -- f not w-homog
  assert(try (hodgeIdealWeightedHomogIsolated(f,3/2,0,{1,1/2});  false) else true); -- alpha
  assert(try (hodgeIdealWeightedHomogIsolated(f,1,-1,{1,1/2});   false) else true); -- p

  -- hodgeCheck bad inputs
  assert(try (hodgeCheck(f,1_R,1,-1);  false) else true);    -- p < 0
  assert(try (hodgeCheck(f,1_R,0,0);   false) else true);    -- alpha = 0
  assert(try (hodgeCheck(f,1_R,-1,0);  false) else true);    -- alpha < 0

  -- hodgeLevel bad inputs
  assert(try (hodgeLevel(f,1_R,0);   false) else true);      -- alpha = 0
  assert(try (hodgeLevel(f,1_R,-1);  false) else true);      -- alpha < 0

  -- alpha > 1 is accepted for hodgeCheck and hodgeLevel
  assert(instance(hodgeCheck(f,1_R,3/2,0), Boolean));
  assert(instance(hodgeLevel(f,1_R,3/2), ZZ));
///


TEST ///
  R = QQ[x,y];
  f = y^2 + x;
  g = 1_R;

  -- monodromyWeightHodgeOnV: alpha and p validation
  assert(try (monodromyWeightHodgeOnV(f,0,0,0);   false) else true);  -- alpha = 0
  assert(try (monodromyWeightHodgeOnV(f,3/2,0,0); false) else true);  -- alpha > 1
  assert(try (monodromyWeightHodgeOnV(f,1,-1,0);  false) else true);  -- p < 0

  -- nuAlpha: alpha > 0
  assert(try (nuAlpha(f,g,0);  false) else true);
  assert(try (nuAlpha(f,g,-1); false) else true);

  -- pFunction: alpha > 0
  assert(try (pFunction(f,g,0);  false) else true);
  assert(try (pFunction(f,g,-1); false) else true);

  -- weightCheck: alpha > 0, w >= 0
  assert(try (weightCheck(f,g,0,0);   false) else true);
  assert(try (weightCheck(f,g,-1,0);  false) else true);
  assert(try (weightCheck(f,g,1,-1);  false) else true);

  -- weightLevel: alpha > 0 (both 3-arg and 4-arg overloads)
  assert(try (weightLevel(f,g,0);  false) else true);
  assert(try (weightLevel(f,g,-1); false) else true);
  bfs = globalBFunction f;
  assert(try (weightLevel(f,g,bfs,0);  false) else true);
  assert(try (weightLevel(f,g,bfs,-1); false) else true);

  -- weightLength: alpha in (0,1]
  assert(try (weightLength(f,0);   false) else true);
  assert(try (weightLength(f,3/2); false) else true);
  assert(try (weightLength(f,-1);  false) else true);

  -- alpha > 1 is accepted for nuAlpha, weightCheck, weightLevel
  assert(instance(nuAlpha(f,g,3/2),       ZZ));
  assert(instance(weightCheck(f,g,3/2,0), Boolean));
  assert(instance(weightLevel(f,g,3/2),   ZZ));
///

