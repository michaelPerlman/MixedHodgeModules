-- Copyright 2026 by Andras Lorincz and Michael Perlman




-* Test section *-


---------------------------------------------------------------
--1. Hodge ideal tests
---------------------------------------------------------------

---------------------------------------------------------------
--1.1 Smooth and normal crossing examples
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
--1.2 Weighted homogeneous isolated examples
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

--to add: BrieskornPham test

---------------------------------------------------------------
--1.3 Semi-invariant examples
---------------------------------------------------------------


TEST ///
 n = 3;
 R = QQ[x_(1,1)..x_(n,n)];
 M = genericMatrix(R,x_(1,1),n,n);
 f = determinant(M);
 alpha = 1;

/// 



---------------------------------------------------------------
--1.4 Misc. examples
---------------------------------------------------------------







---------------------------------------------------------------

---------------------------------------------------------------
--2. nuAlpha, pFunction, weightLevel tests
---------------------------------------------------------------

---------------------------------------------------------------
--2.1 Smooth and normal crossing examples
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
  R=QQ[x,y,z]
  f=x*y*z;
  alpha = 1;
  
///




---------------------------------------------------------------
--2.2 Weighted homogeneous isolated examples
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
--2.3 Semi-invariant examples
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
--2.4 Misc. examples
---------------------------------------------------------------


