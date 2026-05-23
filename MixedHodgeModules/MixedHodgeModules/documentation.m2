--Copyright 2026 by Andras Lorincz and Michael Perlman


----------------------------------------------------------------
--------------------------------------------------------------
 -* Documentation section *-
beginDocumentation()


doc ///
  Key
    MixedHodgeModules
  Headline
    a package for computing Hodge and weight filtrations on localizations
  Description
    Text
     Let $S=\mathbb{C}[x_1,\cdots,x_n]$ and let $f\in S$ be a reduced non-constant polynomial. This package
     has functionality for calculating Hodge filtrations, weight filtrations, and related invariants on the localization $S_f$ and, more generally,
     twisted localizations $S_ff^{-\alpha}$, where $\alpha$ is a rational number in $(0,1]$.
     Many of our algorithms are based on [LY26+] and extensions of [Bla22].

     Writing $\mathcal{D}$ for the Weyl algebra of $S$, the $\mathcal{D}$-module $S_ff^{-\alpha}$ underlies
     a mixed Hodge module. This implies that it is endowed with two increasing filtrations indexed by integers:
     (1) the Hodge filtration $F_{\bullet}$, an infinite filtration by finitely-generated $S$-modules and
     (2) the weight filtration $W_{\bullet}$, a finite filtration by holonomic $\mathcal{D}$-modules.
     These two filtrations yield three important classes of ideals in $S$, namely the Hodge ideals
     [MP19, MP20], the weighted Hodge ideals [Ola23], and the higher multiplier ideals [SY25].

     The Hodge and weight filtrations on $S_ff^{-\alpha}$ are constructed using the graph module $B_f=\Gamma_+(S)$,
     obtained via $\mathcal{D}$-module push forward of $S$ along the embedding $\Gamma$ of $\mathbb{A}^n$ into
     the graph of $f$. The module $B_f$ is endowed with a decreasing filtration $V$, indexed by rational numbers,
     known as the Kashiwara–Malgrange filtration. One can obtain the Hodge-theoretic invariants of $S_ff^{-\alpha}$
     discussed above using the module $V^{\alpha}(B_f)$.

     This package has functionality for calculating the Hodge filtration $F_{\bullet}$ on $V^{\alpha}(B_f)$,
     as well as the monodromy weight filtration on a Hodge filtered piece $F_p(\operatorname{Gr}_V^{\alpha}(B_f))$,
     using algorithms based on  [Bla22]. As applications, this package computes the Hodge ideals, weighted Hodge ideals,
     and higher multiplier ideals for $\mathbb{Q}$-divisors, as well as related invariants including the generation
     level of the Hodge filtration on $S_ff^{-\alpha}$ and the HRH level.

     In another direction, this package has functionality for the $p$-functions of rational functions $(g/f^k)\cdot f^{-\alpha}$
     and the related numerical invariants $\nu_{f,g,\alpha}$, based on [LY26+]. As applications, this package can check
     the weight level of $(g/f^k)\cdot f^{-\alpha}$ in $S_ff^{-\alpha}$,
     as well as the length of the weight filtration on $S_ff^{-\alpha}$.

     As applications, the functionality above is used to calculate the Du Bois complexes $\underline{\Omega}^p_{V(f)}$ and
     the intersection Du Bois complexes $I\underline{\Omega}^p_{V(f)}$.
     
    Tree
      :Filtrations on the Graph Embedding
	@TOH "hodgeOnV"@
	@TOH "HRHCheck"@
	@TOH "HRHLevel"@
        @TOH "monodromyWeightHodgeOnV"@
	@TOH "weightHodgeOnV"@

      :Hodge Filtrations
	@TOH "doesGenerateNext"@
	@TOH "generateNext"@
	@TOH "generationLevel"@
	@TOH "higherMultiplierIdeal"@
	@TOH "hodgeCheck"@
	@TOH "hodgeIdeal"@
	@TOH "hodgeIdealBrieskornPham"@
	@TOH "hodgeIdealDet"@
	@TOH "hodgeIdealWeightedHomogIsolated"@
	@TOH "hodgeLevel"@
	@TOH "hodgeOnV"@
	@TOH "HRHCheck"@
	@TOH "HRHLevel"@

      :Weight Filtrations
	@TOH "adjointIdeal"@
	@TOH "localCohomFW"@
	@TOH "monodromyWeightHodgeOnV"@
	@TOH "nuAlpha"@
	@TOH "pFunction"@
	@TOH "weightCheck"@
	@TOH "weightedHodgeIdeal"@
	@TOH "weightHodgeOnV"@
	@TOH "weightLength"@
	@TOH "weightLevel"@

      :Graded de Rham and Du Bois Complexes
	@TOH "gradedDeRhamComplexH1"@
	@TOH "gradedDeRhamCohomologyH1"@
	@TOH "duBoisComplex"@
	@TOH "intersectionDuBoisComplex"@
	@TOH "isPreDuBois"@


  Caveat
    Throughout this package, the input polynomial $f$ is assumed to be reduced
    and to lie in a polynomial ring over $\mathbb{Q}$.
  References
    See the bibliography at @TO "Works Cited"@.
  Acknowledgement
    We thank Guillem Blanco, Bradley Dirks, Timothy Duff, Mahrud Sayrafi, and Ruijie Yang for helpful conversations. Work on this package began as part of the {\it Macaulay2 Workshop and Mini-School} held at University of Minnesota - Twin Cities,
    funded by NSF Award DMS 2302476.
  SeeAlso
    "Dmodules :: Dmodules"
    "BernsteinSato :: BernsteinSato"
  Subnodes
///



--references


doc ///
  Key
    "Works Cited"
  Headline
    bibliography for the MixedHodgeModules package
  Description
    Text
      This node collects the references used throughout the package.

  References
     [Bla22] G. Blanco, An algorithm for Hodge ideals, Math. Comp. {\bf 91} (2022), no. 338, 2955--2967; MR4473109

     [DOR] B. Dirks, S. Olano, D. Raychaudhury, A Hodge theoretic generalization of $\mathbb{Q}$-homology manifolds, arXiv preprint arXiv:2501.14065 (2025)

     [DY25] D. Davis and R. Yang, On the Hodge filtrations and V-filtrations of mixed Hodge modules, arXiv:2503.16619v4
  
     [Laz04] R. Lazarsfeld, Positivity in Algebraic Geometry II,  Ergebnisse der Mathematik und ihrer Grenzgebiete. 3. Folge. A Series of Modern Surveys in Mathematics, 49, Springer, Berlin, 2004; MR2095472

     [LY25] A. C. Lőrincz and R. Yang, Filtrations of D-modules along semi-invariant functions, arXiv:2504.19383
     
     [LY26+] A. C. Lőrincz and R. Yang, Filtrations on D-modules and multiplicities of roots of Bernstein-Sato polynomials, preprint (2026)

     [MP19] M. Mustaţă and M. Popa, Hodge ideals, Mem. Amer. Math. Soc. {\bf 262} (2019), no. 1268, v+80 pp.; MR4044463
     
     [MP20] M. Mustaţă and M. Popa, Hodge ideals for $\Bbb Q$-divisors, $V$-filtration, and minimal exponent, Forum Math. Sigma {\bf 8} (2020), Paper No. e19, 41 pp.; MR4089396

     [MP20b] M. Mustaţă and M. Popa, Hodge filtration, minimal exponent, and local vanishing, Invent. Math. {\bf 220} (2020), no. 2, 453--478; MR4081135

     [MP22] M. Mustaţă and M. Popa, Hodge filtration on local cohomology, Du Bois complex, and local cohomological dimension, Forum Math. Pi 10 (2022), Paper No. e22

     [Ola22] S. Olano, Weighted multiplier ideals of reduced divisors, Math. Ann. {\bf 384} (2022), no. 3-4, 1091--1126; MR4498468

     [Ola23] S. Olano, Weighted Hodge ideals of reduced divisors, Forum Math. Sigma {\bf 11} (2023), Paper No. e51, 28 pp.; MR4603110

     [PP25] S. G. Park and M. Popa, Hodge symmetry and Lefschetz theorems for singular varieties, arXiv preprint, arXiv:2410.15638 (2025)

     [PR21] M. Perlman and C. Raicu, Hodge ideals for the determinant hypersurface, Selecta Math. (N.S.) {\bf 27} (2021), no. 1, Paper No. 1, 22 pp.; MR4198526

     [SY25] C. Schnell and R. Yang, Higher multiplier ideals, Journal für die reine und angewandte Mathematik (Crelle's Journal), DOI: 10.1515/crelle-2025-0097
     
     [SVV] W. Shen, S. Venkatesh, and A. D. Vo, On k-Du Bois and k-rational singularities, Ann. Inst. Fourier, to appear.
     
     [Zha21] M. Zhang, Hodge filtration and Hodge ideals for $\Bbb Q$-divisors with weighted homogeneous isolated singularities, Asian J. Math. {\bf 25} (2021), no. 5, 641--664; MR4456022
///




--functions/methods:



doc ///
  Key
    adjointIdeal
    (adjointIdeal, RingElement)
  Headline
    compute the adjoint ideal of a reduced divisor
  Usage
    I = adjointIdeal(f)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
  Outputs
    I: Ideal
     an ideal of the ring of f
  Description
    Text
     Returns the adjoint ideal $\operatorname{adj}(D)\subseteq S$ of the reduced divisor $D$ defined by $f$, using the identification $\operatorname{adj}(D) = I^{W_1}_0(D)$, the weighted Hodge ideal
     of level $p=0$ and weight $m=1$. The divisor $D$ has rational singularities if and only if $\operatorname{adj}(D)$ is the unit ideal.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     adjointIdeal(f)
    Text
     We see that the $2\times 2$ determinant has rational singularities. On the other hand, the cusp does not:
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     adjointIdeal(f)
  SeeAlso
     weightedHodgeIdeal
     weightHodgeOnV
  References
      See [Theorem A, Ola22] and [Section 9.3E, Laz04] at @TO "Works Cited"@.

///



doc ///
  Key
    doesGenerateNext
    (doesGenerateNext, RingElement, QQ, ZZ)
    (doesGenerateNext, RingElement, ZZ, ZZ)
  Headline
    test whether a level of the Hodge filtration on $S_ff^{-\alpha}$ generates the next level
  Usage
    B = doesGenerateNext(f,alphaQQ,p)
    B = doesGenerateNext(f,alphaZZ,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
  Outputs
    B: Boolean
     true or false
  Description
    Text
     Returns true if the the image of $F_1(\mathcal{D})$ acting on $F_p(S_ff^{-\alpha})$ is equal to $F_{p+1}(S_ff^{-\alpha})$.
     Here $F_1(\mathcal{D})$ is the first piece of the filtration on the Weyl algebra $\mathcal{D}$ by order of differential operator.

     When $f$ is the $2\times 2$ determinant, the Hodge filtration on $S_f$ has generation level one. In particular:
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     doesGenerateNext(f,1,0)
     doesGenerateNext(f,1,1)
    Text
     When $f$ is the cusp, the Hodge filtration on $S_f$ has generation level zero. In particular:
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     doesGenerateNext(f,1,0)
     doesGenerateNext(f,1,1)
     
  SeeAlso
     generateNext
     hodgeIdeal
     generationLevel
///

doc ///
  Key
    generateNext
    (generateNext, RingElement, QQ, ZZ)
    (generateNext, RingElement, ZZ, ZZ)
  Headline
    compute the ideal generated by the action of $F_1(\mathcal{D})$ on a Hodge level of $S_ff^{-\alpha}$
  Usage
    J = generateNext(f,alphaQQ,p)
    J = generateNext(f,alphaZZ,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
  Outputs
    J: Ideal
     an ideal of the ring of f
  Description
    Text
     This function computes the ideal $J\subseteq S$ satisfying $J\cdot f^{-p-1-\alpha}=F_1(\mathcal{D})\cdot F_p(S_ff^{-\alpha})$.
     The ideal $J$ satisfies $J\subseteq I_{p+1}(\alpha D)$, where $D$ is the reduced divisor corresponding to $f$.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     alpha = 1/2;
     generateNext(f,alpha,0)
     doesGenerateNext(f,alpha,0)
    Text
     We have obtained a proper subset of $I_1(\alpha D)$.

     On the other hand, in the next example, we obtain the entire ideal $I_2(\alpha D)$.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     alpha = 5/6;
     generateNext(f,alpha,1)
     doesGenerateNext(f,alpha,1)
  SeeAlso
     hodgeIdeal
     doesGenerateNext
     generationLevel
///

doc ///
  Key
    generationLevel
    (generationLevel, RingElement, QQ)
    (generationLevel, RingElement, ZZ)
    (generationLevel, RingElement)
  Headline
    compute the generation level of the Hodge filtration on $S_ff^{-\alpha}$
  Usage
    g = generationLevel(f,alphaQQ)
    g = generationLevel(f,alphaZZ)
    g = generationLevel(f)
  Inputs
    f: RingElement
     a polynomial with rational coefficients (in at least two variables)
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer. If omitted, alpha = 1.
    alphaZZ: ZZ
     the integer 1_ZZ.
  Outputs
    g: ZZ
     a non-negative integer
  Description
    Text
     This function computes the generation level of the Hodge filtration on $S_ff^{-\alpha}$.
     It does this by searching downward from the upper bound ($n-1$, or $n-2$ when $\alpha=1$), and tests doesGenerateNext until it
     finds the smallest $p$ such that generation holds for all higher levels. The returned integer is the generation level.

     When $f$ is the $2 \times 2$ determinant, the generation level of $S_f$ is one.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     generationLevel(f)
    Text
     When $f$ is the $\mathsf{A}_1$ singularity and $\alpha = 1/2$, the generation level of $S_ff^{-1/2}$ is one.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     alpha = 1/2;
     generationLevel(f,alpha)
  SeeAlso
     doesGenerateNext
     generateNext
     hodgeIdeal
  References
     See [Theorem B, MP19] and [Theorem E, MP20b] at @TO "Works Cited"@.
///

doc ///
  Key
    higherMultiplierIdeal
    (higherMultiplierIdeal, RingElement, QQ,ZZ)
    (higherMultiplierIdeal, RingElement, ZZ,ZZ)
  Headline
    compute the higher multiplier ideal associated to $\alpha$ in $(0,1]$ and level $p$
  Usage
    J = higherMultiplierIdeal(f,alphaQQ,p)
    J = higherMultiplierIdeal(f,alphaZZ,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
  Outputs
    J: Ideal
     an ideal of the ring of f
  Description
    Text
     This function uses hodgeOnV(f,alpha,p) to calculate the higher multiplier ideal $\widetilde{I}_p(\alpha D)$,
     where $D$ is the reduced divisor associated to $f$. When $\alpha=1$, this is sometimes
     called the microlocal multiplier ideal.

     When $f$ is a semi-invariant function on a multipicity-free space, these ideals coincide with the Hodge ideals [LY25].
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     alpha = 1;
     higherMultiplierIdeal(f,alpha,2)
     hodgeIdeal(f,alpha,2)
    Text
     However, the higher multiplier ideals differ from the Hodge ideals in many cases, see [Equation (1.8), DY25]
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     alpha = 11/12;
     higherMultiplierIdeal(f,alpha,2)
     hodgeIdeal(f,alpha,2)
    Text
     The following example supports [Zha21, Conjecture E]
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     alpha = 1;
     w = {1/2,1/3};
     higherMultiplierIdeal(f,alpha,0)
     hodgeIdealWeightedHomogIsolated(f,alpha,0,w)
     higherMultiplierIdeal(f,alpha,1)
     hodgeIdealWeightedHomogIsolated(f,alpha,1,w)
     higherMultiplierIdeal(f,alpha,2)
     hodgeIdealWeightedHomogIsolated(f,alpha,2,w)
  SeeAlso
     hodgeOnV
     hodgeIdeal
  References
     See [SY25], [Theorem A, MP20], [Zha21], [DY25] at @TO "Works Cited"@.
///


doc ///
  Key
    hodgeCheck
    (hodgeCheck, RingElement, RingElement,QQ,ZZ)
    (hodgeCheck, RingElement, RingElement,ZZ,ZZ)
  Headline
    check if an element of the twisted localization $S_ff^{-\alpha}$ belongs to a particular Hodge level
  Usage
    B = hodgeCheck(f,g,alphaQQ,p)
    B = hodgeCheck(f,g,alphaZZ,w)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
    p: ZZ
     an integer, representing the level
  Outputs
    B: Boolean
     true or false
  Description
    Text
     This function checks if $g/f^{\alpha}$ is in $F_p(S_ff^{-\alpha})$ for a rational number $\alpha>0$.
     If $\beta\in (0,1]$ and $k\in \mathbb{Z}$ with $\alpha = \beta+k$ then it checks if $(g/f^k)\cdot f^{-\beta}$
     is in $F_p(S_ff^{-\alpha})$.

     The following example exhibits that, for $f$ smooth, $1/f^{p+1}$ belongs to $F_p(S_f)$ but not $F_{p-1}(S_f)$.
    Example
     S = QQ[x,y]
     f = x^2+y;
     g = 1_S;
     alpha = 1 + 6;
     hodgeCheck(f,g,alpha,4)
     hodgeCheck(f,g,alpha,5)
    Text
     So $f^{-6}$ belongs to $F_5(S_f)$ but not $F_4(S_f)$.

     We carry out another example.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     g = x*y;
     alpha = 1/2 + 3;
     hodgeCheck(f,g,alpha,1)
     hodgeCheck(f,g,alpha,2)
    Text
     So $(g/f^3)\cdot f^{-1/2}$ belongs to $F_2(S_ff^{-1/2})$ but not $F_1(S_ff^{-1/2})$.
  SeeAlso
    hodgeLevel
    hodgeIdeal
    weightLength
    weightLevel
///

doc ///
  Key
    hodgeIdeal
    (hodgeIdeal, RingElement, QQ, ZZ)
    (hodgeIdeal, RingElement, ZZ, ZZ)
  Headline
    compute the Hodge ideal $I_p(\alpha D)$ of a reduced divisor $D$
  Usage
    I = hodgeIdeal(f,alphaQQ,p)
    I = hodgeIdeal(f,alphaQQ,p, UseGenLevel => False)
    I = hodgeIdeal(f,alphaZZ,p)
    I = hodgeIdeal(f,alphaZZ,p, UseGenLevel => False)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ. 
    p: ZZ
     a non-negative integer  
  Outputs
    I: Ideal
     an ideal of the ring of f
  Description
    Text
     For $\alpha\in (0,1]$ this function computes the Hodge ideal $I_p(\alpha D)$ of a reduced divisor $D$ associated to $f\in S$.
     By default it uses generation level shortcuts (based on [Theorem A, MP20]
     and [Theorem E, MP20b]) via generateNext for large $p$. When UseGenLevel => False, it computes directly from hodgeOnV without using
     general bounds on generation level.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     alpha = 1;
     hodgeIdeal(f,alpha,2)
    Text
     This recovers [Theorem 1.1, PR21] in the case $p=2$ of the $2 \times 2$ determinant (see also [LY25]).

     The next example recovers [Remark 17.13, MP19] (see also [Zha21]).
    Example
     S = QQ[x,y]
     f = x^2+y^3;
     alpha = 1;
     hodgeIdeal(f,alpha,0)
     hodgeIdeal(f,alpha,1)
     hodgeIdeal(f,alpha,2)
    Text
     The next example recovers part of [Example 1, Bla22].
    Example
     S = QQ[x,y];
     f = x^5+y^5+x^2*y^2;
     alpha = 9/10;
     hodgeIdeal(f,alpha,0)
  SeeAlso
     hodgeOnV
     generateNext
     generationLevel
     hodgeIdealWeightedHomogIsolated
     hodgeIdealBrieskornPham
  References
     See [MP19], [MP20], [Bla22], [PR21] at @TO "Works Cited"@.
///

doc ///
  Key
    [hodgeIdeal, UseGenLevel]
  Headline
   turn off use of generation level when calculating hodgeIdeal
  Description
    Text
     Use UseGenLevel => False to turn off use of generation level and to calculate Hodge ideal
     directly using hodgeOnV.
  SeeAlso
     UseGenLevel
     True
     False
///


doc ///
  Key
    hodgeIdealBrieskornPham
    (hodgeIdealBrieskornPham, List, QQ, ZZ)
    (hodgeIdealBrieskornPham, List, ZZ, ZZ)
  Headline
    compute Hodge ideals for Brieskorn–Pham polynomials $x_1^{b_1}+...+x_n^{b_n}$
  Usage
    I = hodgeIdealBrieskornPham(L,alphaQQ,p)
    I = hodgeIdealBrieskornPham(L,alphaZZ,p)
  Inputs
    L: List
     a list {b_1,...,b_n} of positive integers
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
  Outputs
    I: Ideal
     an ideal of the Brieskorn–Pham polynomial ring QQ[x_1..x_n]
  Description
    Text
     This function constructs the Brieskorn–Pham polynomial $f = x_1^{b_1}+...+x_n^{b_n}$ and computes the Hodge ideal
     $I_p(\alpha D)$ where $D$ is the reduced divisor of $f$.
     This routine uses the weighted homogeneous isolated singularity formula [Zha21], implemented in hodgeIdealWeightedHomogIsolated.

     The following example is for the $\mathsf{A}_1$ singularity $x_1^2+x_2^2+x_3^2$.
    Example
     L = {2,2,2};
     alpha = 1/2;
     hodgeIdealBrieskornPham(L,alpha,2)
    Text
     The next example is for the $\mathsf{E}_6$ singularity $x_1^2+x_2^3+x_3^4$.
    Example
     L = {2,3,4};
     alpha = 1;
     hodgeIdealBrieskornPham(L,alpha,0)
    Text
     This recovers the fact that this polynomial has log-canonical singularities
  SeeAlso
     hodgeIdeal
     hodgeIdealWeightedHomogIsolated
  References
     See [Zha21] at @TO "Works Cited"@.
///


doc ///
  Key
    hodgeIdealDet
    (hodgeIdealDet, ZZ, ZZ)
  Headline
    compute Hodge ideals for the determinant of a generic $n\times n$ matrix
  Usage
    I = hodgeIdealDet(n,p)
  Inputs
    n: ZZ
     a positive integer (size of the generic matrix)
    p: ZZ
     a non-negative integer (Hodge filtration index)
  Outputs
    I: Ideal
     the Hodge ideal $I_p(f)$ where $f$ is the determinant of the generic
     $n\times n$ matrix, returned as an ideal in $(\mathbb{Z}/32003)[x_{1,1},\ldots,x_{n,n}]$.
  Description
    Text
     Let $f = \det(x_{i,j})$ be the determinant of the generic $n\times n$ matrix.
     This function computes the Hodge ideal $I_p(f)$ using the description of [PR21]
     as an intersection of symbolic powers of ideals of minors:
     $$I_p(f) = \bigcap_{q=1}^{n-1} I_q^{((n-q)(p-1) - \binom{n-q}{2})},$$
     where $I_q$ is the ideal of $q\times q$ minors and the symbolic powers are
     understood to be the trivial ideal when the exponent is non-positive.
     The symbolic powers themselves are built from the $\mathrm{GL}\times\mathrm{GL}$
     representation theory of the coordinate ring: for each relevant partition $\lambda$,
     ILambda assembles a $\mathrm{GL}\times\mathrm{GL}$-stable subspace of the
     symbolic power as a sum of products of minors of randomly twisted matrices.

     For efficiency, the computation is performed over $\mathbb{Z}/32003$
     throughout.

     The case $n=2$:
    Example
     hodgeIdealDet(2,0)
     hodgeIdealDet(2,1)
     hodgeIdealDet(2,2)
    Text
     The case $n=3$.  For $p=2$ the result is the ideal of $2\times 2$ minors;
     for $p=3$ it is generated by quadrics involving differences of products of
     entries:
    Example
     hodgeIdealDet(3,2)
     numgens hodgeIdealDet(3,3)
  SeeAlso
     hodgeIdeal
     hodgeIdealBrieskornPham
     hodgeIdealWeightedHomogIsolated
  References
     See [PR21] at @TO "Works Cited"@.
///


doc ///
  Key
    hodgeIdealWeightedHomogIsolated
    (hodgeIdealWeightedHomogIsolated, RingElement, QQ, ZZ, List)
    (hodgeIdealWeightedHomogIsolated, RingElement, ZZ, ZZ, List)
  Headline
    compute Hodge ideals for weighted homogeneous isolated singularities
  Usage
    I = hodgeIdealWeightedHomogIsolated(f,alphaQQ,p,w)
    I = hodgeIdealWeightedHomogIsolated(f,alphaZZ,p,w)
  Inputs
    f: RingElement
     a weighted homogeneous polynomial with isolated singularity
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
    w: List
     a list of rational weights $w_1,...,w_n$ such that $w\cdot \beta=1$ for every monomial $x^{\beta}$ of $f$
  Outputs
    I: Ideal
     an ideal of the ring of f
  Description
    Text
     Implements the recursive formula for Hodge ideals of $\mathbb{Q}$-divisors with weighted homogeneous isolated singularities, using
     [Zha21]. The polynomial $f$ should have isolated singularities, and $w$ should be normalized so that the dot product of $w$ with any
     exponent of any term of $f$ is one.
    Example
      S = QQ[x,y];
      f = x^2+y^3;
      w = {1/2,1/3};
      hodgeIdealWeightedHomogIsolated(f,1,2,w)
      hodgeIdealWeightedHomogIsolated(f,5/6,2,w)
    Text
      We carry out another example.
    Example
      S = QQ[x,y,z];
      f = y^2-x*z;
      w = {1/2,1/2,1/2};
      hodgeIdealWeightedHomogIsolated(f,1,1,w)
      hodgeIdealWeightedHomogIsolated(f,1/2,1,w)     
  SeeAlso
     hodgeIdeal
     hodgeIdealBrieskornPham
  References
      See [Zha21] at @TO "Works Cited"@.
///


doc ///
  Key
    hodgeLevel
    (hodgeLevel, RingElement, RingElement, QQ)
    (hodgeLevel, RingElement, RingElement, ZZ)
  Headline
    find the Hodge level of an element in a twisted localization $S_ff^{-\alpha}$
  Usage
    L = hodgeLevel(f,g,alphaQQ)
    L = hodgeLevel(f,g,alphaZZ)
    L = hodgeLevel(f,g,bgfs,alphaQQ)
    L = hodgeLevel(f,g,bgfs,alphaZZ)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
  Outputs
    L: ZZ
  Description
    Text
     This function finds the minimal $p$ for which $g/f^{\alpha}$ is in $F_p(S_ff^{-\alpha})$ for $\alpha>0$.
     If $\beta\in (0,1]$ and $k\in \mathbb{Z}$ with $\alpha =\beta+k$ then it finds the minimal $p$
     for which $(g/f^k)\cdot f^{-\beta}$ is in $F_p(S_ff^{-\beta})$.
     
     The following example exhibits that, for $f$ smooth, $1/f^{p+1}$ belongs to $F_p(S_f)$ but not $F_{p-1}(S_f)$.
    Example
     S = QQ[x,y]
     f = x^2+y;
     g = 1_S;
     alpha = 1 + 6;
     hodgeLevel(f,g,alpha)
    Text
     So $f^{-6}$ belongs to $F_5(S_f)$ but not $F_4(S_f)$.

     We carry out another example.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     g = x*y;
     alpha = 1/2 + 3;
     hodgeLevel(f,g,alpha)
    Text
     So $(g/f^3)\cdot f^{-1/2}$ belongs to $F_2(S_ff^{-1/2})$ but not $F_1(S_ff^{-1/2})$.
  SeeAlso
    hodgeCheck
    hodgeIdeal
    weightCheck
    weightLevel
///

doc ///
  Key
    hodgeOnV
    (hodgeOnV, RingElement, QQ, ZZ)
    (hodgeOnV, RingElement, ZZ, ZZ)
    (hodgeOnV, RingElement, ZZ)
  Headline
    compute bases for the Hodge filtration on $V^{\alpha}(B_f)$ 
  Usage
    L = hodgeOnV(f,alphaQQ,p)
    L = hodgeOnV(f,alphaZZ,p)
    L = hodgeOnV(f,alpha,p, UseBasis => dtBasis)
    L = hodgeOnV(f,alpha,p, UseBasis => sBasis)
    V = hodgeOnV(f,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer
  Outputs
    L: List
     when alpha is given, a list giving an $S$-basis for $F_p(V^{\alpha}(B_f))$
    V: HashTable
     when alpha is omitted, a hashtable with keys $\alpha$ in $(0,1]$ and values bases for $F_p(V^{\alpha}(B_f)))$
  Description
    Text
     This function uses a modification of Blanco's algorithm [Bla22] to calculate an $S$-basis of $F_p(V^{\alpha}(B_f))$. By default,
     the output is in terms of $\partial_t$, where $t$ is the $(n+1)$-st coordinate of the graph of $f$.
     Our convention is that $\partial_t^p$ is the largest power of $\partial_t$ in $F_p(B_f)$.

     If alpha is omitted, the function outputs a hash table with keys equal to rational numbers $\alpha$ such that $-\alpha-p$ is a root of the
     $(p+1)$-st generalized $b$-function of $f$, and values equal to an $S$-basis of $F_p(V^{\alpha}(B_f))$.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     hodgeOnV(f,5/6,1)
     hodgeOnV(f,1)
    Text
     The Malgrange isomorphism gives an identification of $\Gamma_+(S_f)$ with $S_f[s]f^s$. When $\alpha>0$, we get
     an isomorphism between $V^{\alpha}(B_f)$ and $V^{\alpha}(S_f[s]f^s)$, where the latter will have denominators.
     If the user selects, UseBasis => sBasis, then the result expresses $F_p(V^{\alpha}(B_f))$ as a submodule of $S_f[s]f^s$.
    Example
     S = QQ[x,y,z]
     f = y^2-x*z;
     hodgeOnV(f,1/2,1, UseBasis => sBasis)
     hodgeOnV(f,1, UseBasis => sBasis)
  SeeAlso
     hodgeIdeal
     weightHodgeOnV
     dtBasis
     sBasis
     UseBasis
  References
     See [Bla22] at @TO "Works Cited"@.
///

doc ///
  Key
    [hodgeOnV, UseBasis]
  Headline
   toggle between $\partial_t$-basis of $V^{\alpha}(B_f)$ and $s$ basis inside $S_f[s]f^s$
  Description
    Text
     The Malgrange isomorphism gives an identification between $V^{\alpha}(\Gamma_+(S_f))$ and $S_f[s]f^s$.
     If UseBasis is set to dtBasis (default), then outputs an $S$-basis of $V^{\alpha}(B_f) in terms
     of powers of $\partial_t$. If UseBasis is set to sBasis, then outputs an $S$-basis of $V^{\alpha}(B_f) as
     a submodule of $S_f[s]f^s$, which generally will have denominators.
  SeeAlso
     hodgeOnV
     monodromyWeightHodgeOnV
     weightHodgeOnV
     dtBasis
     sBasis
     UseBasis
///


doc ///
  Key
    HRHCheck
    (HRHCheck, RingElement, ZZ)
  Headline
    check if the HRH level of $f$ is at least some number
  Usage
    B = HRHCheck(f,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    p: ZZ
     an integer, representing the level
  Outputs
    B: Boolean
     true or false
  Description
    Text
     This function returns true if $\operatorname{HRH}(f)\geq p$, where $\operatorname{HRH}(f)$ is the HRH level of [DOR, PP25].

     The following example shows that the HRH level of the $2\times 2$ determinant is zero.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     HRHCheck(f,0)
     HRHCheck(f,1)
    Text
     In contrast, finite quotient singularities, such as the $\mathsf{A}_1$ singularity, are rational homology manifolds.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     HRHCheck(f,0)
     HRHCheck(f,1)
     HRHCheck(f,2)     
  SeeAlso
    hodgeCheck
    hodgeLevel
    HRHLevel
    weightCheck
    weightLevel
  References
    See [DOR] and [PP25] at @TO "Works Cited"@.
///

doc ///
  Key
    HRHLevel
    (HRHLevel, RingElement)
  Headline
    determine the HRH level of $f$
  Usage
    L = HRHLevel(f)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
  Outputs
    L: ZZ
  Description
    Text
     This function determines $\operatorname{HRH}(f)$, the HRH level of [DOR, PP25].
     If $\operatorname{HRH}(f)=\infty$ then this function outputs the string "rational homology manifold".

     The following example shows that the HRH level of the $2\times 2$ determinant is zero.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     HRHLevel(f)
    Text
     In contrast, finite quotient singularities, such as the $\mathsf{A}_1$ singularity, are rational homology manifolds.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     HRHLevel(f)   
  SeeAlso
    hodgeLevel
    HRHCheck
    weightLevel
  References
    See [DOR] and [PP25] at @TO "Works Cited"@.
///

doc ///
  Key
    IHmoduleAdjoint
    (IHmoduleAdjoint, RingElement)
  Headline
    construct a presentation of the intersection cohomology $\mathcal{D}$-module using the adjoint ideal
  Usage
    M = IHmoduleAdjoint(f)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
  Outputs
    M: Module
     a left module over the Weyl algebra of ring f, equal to D^r/K
  Description
    Text
     This function constructs a $\mathcal{D}$-module presentation of the intersection cohomology module associated to the
     hypersurface $f$ using generators of adjointIdeal(f).
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     IHmoduleAdjoint(f)
    Text
     When $f$ has rational singularities, this presents the intersection cohomology module as a cyclic module.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     IHmoduleAdjoint(f)
  SeeAlso
     adjointIdeal
     weightedHodgeIdeal
///

doc ///
  Key
    localCohomFW
    (localCohomFW, Ideal, ZZ, ZZ, ZZ)
  Headline
    compute the Hodge filtration on the weight filtration on local cohomology
  Usage
    M = localCohomFW(I,q,p,m)
  Inputs
    I: Ideal
     an ideal in the polynomial ring with rational coefficients
    q: ZZ
     an integer representing cohomological degree
    p: ZZ
     an integer representing the Hodge level
    m: ZZ
     an integer $\geq n$ representing the weight level, where $n$ is the number of generators of $S$ 
  Outputs
    M: Module
     a module over the polynomial ring
  Description
    Text
     This function computes $F_p(W_m(H^q_I(S)))$, where $I$ is an ideal in the polynomial ring $S=\mathbb{Q}[x_1,\cdots, x_n]$.
     If $I$ is prime and $c=\operatorname{codim}(I)$, then the first interesting example is $F_0(W_{c+n}(H^c_I(S)))$.

     The following example reflects that the $2\times 2$ minors of a $2\times 3$ matrix have rational singularities.
     Indeed, $F_0(W_{c+6}(H^c_I(S)))$ is the canonical module of $S/I$.
    Example
     S = QQ[x_(1,1)..x_(3,2)];
     M = transpose genericMatrix(S,x_(1,1),2,3);
     I = minors(2,M);
     M = localCohomFW(I,2,0,2+6)
    Text
     In the next example, we get the first two steps of the Hodge filtration on the IC module corresponding to the origin of $\mathbb{A}^3$.
    Example
     S = QQ[x,y,z];
     I = ideal(x,y,z);
     localCohomFW(I,3,0,3+3)
     localCohomFW(I,3,1,3+3)
  SeeAlso
     weightedHodgeIdeal
///

doc ///
  Key
    monodromyWeightHodgeOnV
    (monodromyWeightHodgeOnV, RingElement, QQ, ZZ, ZZ)
    (monodromyWeightHodgeOnV, RingElement, ZZ, ZZ, ZZ)
  Headline
    compute the monodromy weight filtration on $F_p(\operatorname{Gr}^{\alpha}_V(B_f))$
  Usage
    L = monodromyWeightHodgeOnV(f,alphaQQ,p,m)
    L = monodromyWeightHodgeOnV(f,alphaQQ,p,m,UseBasis => sBasis)
    L = monodromyWeightHodgeOnV(f,alphaZZ,p,m)
    L = monodromyWeightHodgeOnV(f,alphaZZ,p,m,UseBasis => sBasis)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in $(0,1]$. May also be given as an integer.
    alphaZZ: ZZ
     the integer $1$.
    p: ZZ
     a non-negative integer
    m: ZZ
     a non-negative integer indexing the monodromy weight filtration (centered at $0$)
  Outputs
    L: List
     a list representing an $S$-basis for a lift of $F_p(W(N)_m(\operatorname{Gr}^{\alpha}_V(B_f)))$ to $F_p(V^{\alpha}(B_f))$
  Description
    Text
     Let $s=-\partial_t t$. Multiplication by $(s+\alpha)$ is a nilpotent operator on $\operatorname{Gr}^{\alpha}_V(B_f)$.
     Writing $N = s+\alpha$, we obtain the monodromy weight filtration $W(N)_\bullet$ on $\operatorname{Gr}^{\alpha}_V(B_f)$,
     centered at $0$, given by
    Text
     $W(N)_m = \sum_{i+j=m} \operatorname{ker}(N^{i+1}) \cap \operatorname{im}(N^{-j})$.
    Text
     This is an increasing filtration by $\mathcal{D}$-modules. The function computes a lift of
     $F_p(W(N)_m(\operatorname{Gr}^{\alpha}_V(B_f)))$ to a submodule of $F_p(V^{\alpha}(B_f))$.
     The output is an $S$-basis for a submodule of $F_p(V^{\alpha}(B_f))$ whose image in
     $\operatorname{Gr}^{\alpha}_V(B_f)$ is equal to $F_p(W(N)_m(\operatorname{Gr}^{\alpha}_V(B_f)))$.
     In particular, the output always contains $F_p(V^{>\alpha}(B_f))$, and equality with
     $F_p(V^{>\alpha}(B_f))$ corresponds to the case when
     $F_p(W(N)_m(\operatorname{Gr}^{\alpha}_V(B_f))) = 0$.
    Text
     For the $\mathsf{A}_1$ singularity, the output stabilizes immediately.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     monodromyWeightHodgeOnV(f,1,1,0)
     monodromyWeightHodgeOnV(f,1,1,1)
    Text
     For the $2\times 2$ determinant, one sees a jump between $m=0$ and $m=1$.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     monodromyWeightHodgeOnV(f,1,1,0)
     monodromyWeightHodgeOnV(f,1,1,1)
  SeeAlso
     hodgeOnV
     weightHodgeOnV
     weightedHodgeIdeal
///

doc ///
  Key
    [monodromyWeightHodgeOnV, UseBasis]
  Headline
   toggle between $\partial_t$-basis of $V^{\alpha}(B_f)$ and $s$ basis inside $S_f[s]f^s$
  Description
    Text
     The Malgrange isomorphism gives an identification between $V^{\alpha}(\Gamma_+(S_f))$ and $S_f[s]f^s$.
     If UseBasis is set to dtBasis (default), then outputs an $S$-basis of $V^{\alpha}(B_f) in terms
     of powers of $\partial_t$. If UseBasis is set to sBasis, then outputs an $S$-basis of $V^{\alpha}(B_f) as
     a submodule of $S_f[s]f^s$, which generally will have denominators.
  SeeAlso
     hodgeOnV
     monodromyWeightHodgeOnV
     weightHodgeOnV
     dtBasis
     sBasis
     UseBasis
///

doc ///
  Key
    nuAlpha
    (nuAlpha, RingElement, RingElement, QQ)
    (nuAlpha, RingElement, RingElement, ZZ)
    (nuAlpha, RingElement, RingElement, RingElement, QQ)
    (nuAlpha, RingElement, RingElement, RingElement, ZZ)
  Headline
    determine $\nu_{\alpha}$ for a pair of polynomials and a rational number
  Usage
    N = nuAlpha(f,g,alphaQQ)
    N = nuAlpha(f,g,alphaZZ)
    N = nuAlpha(f,g,bgfs,alphaQQ)
    N = nuAlpha(f,g,bgfs,alphaZZ)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
    bgfs: RingElement
     the $b$-function of $g$ with respect to $f$, given as an element of $\mathbb{Q}[s]$.
  Outputs
    N: ZZ
  Description
    Text
     A description of what
  SeeAlso
     pFunction
     weightLength
     
  References
     See [LY26+] at @TO "Works Cited"@.
///

doc ///
  Key
    [nuAlpha, NuMethod]
  Headline
   toggle between the strategies for calculating nuAlpha
  Description
    Text
     The three NuMethods are Malgrange (only for g=1_R), PowerBFunction, and ByAnnFs.
  SeeAlso
     Malgrange
     PowerBFunction
     ByAnnFs
///

doc ///
  Key
    pFunction
    (pFunction, RingElement, RingElement, QQ)
    (pFunction, RingElement, RingElement, ZZ)
  Headline
    calculate the $p$-function of a pair of polynomials associated to a rational number
  Usage
    p = pFunction(f,g,alphaQQ)
    p = pFunction(f,g,alphaZZ)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
  Outputs
    p: RingElement
     a univariate polynomial
  Description
    Text
     A description of what 
  SeeAlso
     nuAlpha
     
  References
     See [LY26+] at @TO "Works Cited"@.
///

doc ///
  Key
    [pFunction, NuMethod]
  Headline
   toggle between the strategies for calculating nuAlpha
  Description
    Text
     The three NuMethods are Malgrange (only for g=1_S), PowerBFunction, and ByAnnFs.
  SeeAlso
     Malgrange
     PowerBFunction
     ByAnnFs
///


doc ///
  Key
    weightCheck
    (weightCheck, RingElement, RingElement, QQ, ZZ)
    (weightCheck, RingElement, RingElement, ZZ, ZZ)
  Headline
    check if an element of the twisted localization $S_ff^{-\alpha}$ belongs to a particular weight
  Usage
    B = weightCheck(f,g,alphaQQ,m)
    B = weightCheck(f,g,alphaZZ,m)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
    m: ZZ
     an integer, representing the level
  Outputs
    B: Boolean
     true or false
  Description
    Text
     This function checks if $g/f^{\alpha}$ is in $W_m(S_ff^{-\alpha})$ for some rational $\alpha>0$.
     If $\beta\in (0,1]$ and $k\in \mathbb{Z}$  with $\alpha = \beta +k$ then
     it checks if $(g/f^k)\cdot f^{-\beta}$ is in $W_m(S_ff^{-\beta})$.
     
     For $k>0$, when $f$ is smooth, the quotient $1/f^k$ belongs to $W_1(S_f)$ but not $W_0(S_f)$.
    Example
     S = QQ[x,y]
     f = x^2+y;
     g = 1_S;
     alpha = 1 + 6;
     weightCheck(f,g,alpha,0)
     weightCheck(f,g,alpha,1)
    Text
     so $g/f^6$ belongs to $W_1(S_f)$ but not $W_0(S_f)$.
     
     We carry out another example.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     g = x*y;
     alpha = 1/2 + 3;
     weightCheck(f,g,alpha,0)
     weightCheck(f,g,alpha,1)
    Text
     so $(g/f^3)\cdot f^{-1/2}$ belongs to $W_1(S_ff^{-1/2})$ but not $W_0(S_ff^{-1/2})$.

     We carry out another example.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     g = x;
     alpha = 1 + 3;
     weightCheck(f,g,alpha,1)
     weightCheck(f,g,alpha,2)
    Text
     so $g/f^3$ belongs to $W_2(S_f)$ but not $W_1(S_f)$.
  SeeAlso
    hodgeCheck
    hodgeLevel
    weightedHodgeIdeal
    weightLength
    weightLevel
     
  References
     See [LY26+] at @TO "Works Cited"@.
///

doc ///
  Key
    weightLength
    (weightLength, RingElement, QQ)
    (weightLength, RingElement, ZZ)
  Headline
    determine the length of the weight filtration on a twisted localization $S_ff^{-\alpha}$
  Usage
    L = weightLength(f,alphaQQ)
    L = weightLength(f,alphaZZ)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
  Outputs
    L: ZZ
  Description
    Text
     This function finds the length of the weight filtration on $S_ff^{-\alpha}$, i.e. the maximal index $m$
     for which the graded piece $\operatorname{Gr}^W_{m+n}(S_ff^{-\alpha})$ is nonzero.
     In this setup, the smallest weight $m$ such that $W_{m+n}(S_ff^{-\alpha})$ is nonzero is $m=0$.

     When $D=V(f)$ is a rational homology manifold and $\alpha=1$, the weight length is one.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     weightLength(f,1)
    Text
     On the other hand, when $f$ is the $n \times n$ determinant and $\alpha=1$, the weight length is $n$. See [PR21] or [LY25].
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     weightLength(f,1)
    Text
     We carry out another example.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     alpha = 5/6;
     weightLength(f,alpha)       
  SeeAlso
     weightCheck
     weightLevel
     weightedHodgeIdeal
     
  References
    See [LY26+] and [PR21] at @TO "Works Cited"@.
///


doc ///
  Key
    [weightLength, LengthStrategy]
  Headline
    toggle between the strategies for calculating weightLength
  Description
    Text
     The two LengthStrategies are ByNuAlpha and ByWeightLevel.
  SeeAlso
     ByNuAlpha
     ByWeightLevel
///

doc ///
  Key
    [weightLength, NuMethod]
  Headline
    when using LengthStrategy => ByNuAlpha, toggle between the strategies for calculating nuAlpha
  Description
    Text
     The three NuMethods are Malgrange (only for g=1_S), PowerBFunction, and ByAnnFs.
  SeeAlso
     Malgrange
     PowerBFunction
     ByAnnFs
///

doc ///
  Key
    weightLevel
    (weightLevel, RingElement, RingElement, QQ)
    (weightLevel, RingElement, RingElement, ZZ)
    (weightLevel, RingElement, RingElement, RingElement, QQ)
    (weightLevel, RingElement, RingElement, RingElement, ZZ)
  Headline
    find the weight level of an element in a twisted localization $S_ff^{-\alpha}$
  Usage
    L = weightLevel(f,g,alphaQQ)
    L = weightLevel(f,g,alphaZZ)
    L = weightLevel(f,g,bgfs,alphaQQ)
    L = weightLevel(f,g,bgfs,alphaZZ)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    g: RingElement
     an element of the ring of f
    alphaQQ: QQ
     a positive rational number. May also be given as an integer.
    alphaZZ: ZZ
     a positive integer.
    bgfs: RingElement
     the b-function of g with respect to f, as an element of QQ[s].
  Outputs
    L: ZZ
  Description
    Text
     This function finds the minimal weight level $m$ for which $g/f^{\alpha}$ is in $W_m(S_ff^{-\alpha})$.
     If $\beta\in (0,1]$ and $k\in \mathbb{Z}$  with $\alpha = \beta+k$ then it
     finds minimal weight level $m$ for which $(g/f^k)\cdot f^{-\beta}$ is in $W_m(S_ff^{-\beta})$.
     
     For $k>0$, if $f$ is smooth, the quotient $1/f^k$ belongs to $W_1(S_f)$ but not $W_0(S_f)$:
    Example
     S = QQ[x,y]
     f = x^2+y;
     g = 1_S;
     alpha = 1 + 6;
     weightLevel(f,g,alpha)
    Text
     We carry out another example.
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     g = x*y;
     alpha = 1/2 + 3;
     weightLevel(f,g,alpha)
    Text
     so $(g/f^3)\cdot f^{-1/2}$ belongs to $W_1(S_ff^{-1/2})$ but not $W_0(S_ff^{-1/2})$.

     We carry out another example.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     g = x;
     alpha = 1 + 3;
     weightLevel(f,g,alpha)
    Text
     so $g/f^3$ belongs to $W_2(S_f)$ but not $W_1(S_f)$.
  SeeAlso
    hodgeCheck
    hodgeLevel
    weightedHodgeIdeal
    weightCheck
    weightLength
  References
     See [LY26+] at @TO "Works Cited"@.
///

doc ///
  Key
    weightedHodgeIdeal
    (weightedHodgeIdeal, RingElement, QQ, ZZ, ZZ)
    (weightedHodgeIdeal, RingElement, ZZ, ZZ, ZZ)
  Headline
    compute the weighted Hodge ideal $I^{W_m}_p(\alpha D)$
  Usage
    I = weightedHodgeIdeal(f,alphaQQ,p,m)
    I = weightedHodgeIdeal(f,alphaZZ,p,m)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer (Hodge index)
    m: ZZ
     a non-negative integer (weight index)
  Outputs
    I: Ideal
     an ideal of the ring of f
  Description
    Text
     This function computes the weighted Hodge ideal $I^{W_m}_p(\alpha D)$, where $D$
     is the reduced divisor of $f$. It does this using a basis for $K^{\alpha}_m$ (see weightHodgeOnV).

     When $m=0$, one gets ideals generated by powers of $f$.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     alpha = 1;
     weightedHodgeIdeal(f,alpha,0,0)
     weightedHodgeIdeal(f,alpha,1,0)
     weightedHodgeIdeal(f,alpha,2,0)
    Text
     We carry out another example (see [LY25]).
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     alpha = 1/2;
     weightedHodgeIdeal(f,alpha,0,1)
     weightedHodgeIdeal(f,alpha,1,1)
     weightedHodgeIdeal(f,alpha,2,1)
  SeeAlso
     weightHodgeOnV
     adjointIdeal
     localCohomFW
  References
     See [Bla22] and [Ola23] at @TO "Works Cited"@.
///


doc ///
  Key
    weightHodgeOnV
    (weightHodgeOnV, RingElement, QQ, ZZ, ZZ)
    (weightHodgeOnV, RingElement, ZZ, ZZ, ZZ)
  Headline
    compute a basis for a weighted piece $K^{\alpha}_m$
  Usage
    L = weightHodgeOnV(f,alphaQQ,p,m)
    L = weightHodgeOnV(f,alphaQQ,p,m, UseBasis => sBasis)
    L = weightHodgeOnV(f,alphaZZ,p,m)
    L = weightHodgeOnV(f,alphaZZ,p,m, UseBasis => sBasis)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    alphaQQ: QQ
     a rational number in (0,1]. May also be given as an integer.
    alphaZZ: ZZ
     the integer 1_ZZ.
    p: ZZ
     a non-negative integer (Hodge filtration index)
    m: ZZ
     a non-negative integer (weight index)
  Outputs
    L: List
     a list of elements representing an $S$-basis for $K^{\alpha}_m$
  Description
    Text
     This function uses a modification of Blanco's algorithm [Bla22] to calculate an $S$-basis of $K^{\alpha}_m$,
     which is by definition the lift of $\operatorname{ker}((s+\alpha)^m)$ on $\operatorname{Gr}_V^{\alpha}(F_p(B_f))$
     (see the proof of [Theorem A, Ola23]). By default,
     the output is in terms of $\partial_t$, where $t$ is the $(n+1)$-st coordinate of the graph of $f$.
     Our convention is that $\partial_t^p$ is the largest power of $\partial_t$ in $F_p(B_f)$.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     weightHodgeOnV(f,5/6,1,1)
    Text
     If the user selects, UseBasis => sBasis, then the result expresses $K^{\alpha}_m$ as a submodule of $S_f[s]f^s$
     (see hodgeOnV for more information).
    Example
     S = QQ[x,y,z]
     f = y^2-x*z;
     weightHodgeOnV(f,1/2,1,1, UseBasis => sBasis)
  SeeAlso
     monodromyWeightHodgeOnV
     weightedHodgeIdeal
     weightLength
     hodgeOnV
     dtBasis
     sBasis
     UseBasis
  References
     See [Bla22] and [Ola23] at @TO "Works Cited"@.
///


doc ///
  Key
    [weightHodgeOnV, UseBasis]
  Headline
   toggle between dt-basis of $V^{\alpha}(B_f)$ and $s$ basis inside $S_f[s]f^s$
  Description
    Text
     The Malgrange isomorphism gives an identification between $V^{\alpha}(\Gamma_+(S_f))$ and $S_f[s]f^s$.
     If UseBasis is set to dtBasis (default), then outputs an $S$-basis of $V^{\alpha}(B_f)$ in terms
     of powers of $\partial_t$. If UseBasis is set to sBasis, then outputs an $S$-basis of $V^{\alpha}(B_f)$ as
     a submodule of $S_f[s]f^s$, which generally will have denominators.
  SeeAlso
     hodgeOnV
     weightHodgeOnV
     dtBasis
     sBasis
     UseBasis
///



--------------------------------------------------------------
--de Rham functions
--------------------------------------------------------------


doc ///
  Key
    duBoisComplex
    (duBoisComplex, RingElement, ZZ)
    (duBoisComplex, RingElement, ZZ, List)
  Headline
    compute the Du Bois complexes $\underline{\Omega}^p_D$ for a divisor
  Usage
    C = duBoisComplex(f,p)
    C = duBoisComplex(f,p,w)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    p: ZZ
     a non-negative integer
    w: List
     (optional) a list of positive weights $\{w_1,\ldots,w_n\}$ for which
     $f$ is weighted-homogeneous of weight $1$ with an isolated singularity
  Outputs
    C: Complex
     a free complex quasi-isomorphic to $\underline{\Omega}^p_D$
  Description
    Text
     Given a reduced polynomial $f\in S=\mathbb{Q}[x_1,\cdots,x_n]$ and $0\leq p\leq n-1$, this function computes a free complex
     quasi-isomorphic to $\underline{\Omega}^p_D$, where $D=V(f)$. As complexes in Macaulay2
     are homologically graded, the $-i$-th homology of this complex is the $i$-th
     cohomology of $\underline{\Omega}^p_D$. The homology of this complex lives in non-positive degrees.

     If $D$ is non-singular, we get $S$-free resolutions of the Kahler differentials
    Example
     S = QQ[x,y,z];
     f = y^2-x;
     DB0 = duBoisComplex(f,0)
     prune HH_0(DB0)
     DB1 = duBoisComplex(f,1)
     prune HH_0(DB1)
     DB2 = duBoisComplex(f,2)
     prune HH_0(DB2)
    Text
     The following example shows that the cusp does not have Du Bois singularities, as $H^1(\underline{\Omega}^0_D)\neq 0$.
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     DB0 = duBoisComplex(f,0)
     prune HH_(-1)(DB0)
     DB1 = duBoisComplex(f,1)
     prune HH_(-1)(DB1)
    Text
     If $f$ is weighted-homogeneous of weight $1$ with respect to a list of
     weights $w$ and has an isolated singularity at the origin, passing $w$
     uses @TO "hodgeIdealWeightedHomogIsolated"@ internally, which is
     dramatically faster than the general method.  The $\mathsf{E}_8$
     singularity:
    Example
     S = QQ[x,y,z];
     f = x^2+y^3+z^5;
     w = {1/2, 1/3, 1/5};
     DB1 = duBoisComplex(f,1,w);
     prune HH_(-1)(DB1)
  SeeAlso
     intersectionDuBoisComplex
     gradedDeRhamComplexH1
     gradedDeRhamCohomologyH1
     hodgeIdealWeightedHomogIsolated
///


doc ///
  Key
    gradedDeRhamComplexH1
    (gradedDeRhamComplexH1, RingElement, ZZ)
    (gradedDeRhamComplexH1, RingElement, ZZ, List)
  Headline
    compute the complexes $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$
  Usage
    C = gradedDeRhamComplexH1(f,p)
    C = gradedDeRhamComplexH1(f,p,w)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    p: ZZ
     a non-positive integer (Hodge filtration index)
    w: List
     (optional) a list of positive weights $\{w_1,\ldots,w_n\}$ for which
     $f$ is weighted-homogeneous of weight $1$ with an isolated singularity
  Outputs
    C: Complex
     a complex of finitely generated $S$-modules quasi-isomorphic to $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$
  Description
    Text
     Given a reduced polynomial $f\in S=\mathbb{Q}[x_1,\cdots,x_n]$ and $p\geq -n$, this function computes a complex of finitely generated $S$-modules
     quasi-isomorphic to $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$. As complexes in Macaulay2
     are homologically graded, the $-i$-th homology of this complex is the $i$-th
     cohomology of $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$.
    Example
     S = QQ[x,y,z]
     f = x^2+y^3+z^2
     C1 = gradedDeRhamComplexH1(f,-3)
     prune HH_0(C1)
     C2 = gradedDeRhamComplexH1(f,-2)
     prune HH_0(C2)
     prune HH_1(C2)
     C3 = gradedDeRhamComplexH1(f,-1)
     prune HH_2(C3)
    Text
     When $f$ is weighted-homogeneous of weight $1$ with respect to weights $w$
     and has an isolated singularity at the origin, passing $w$ uses
     @TO "hodgeIdealWeightedHomogIsolated"@ internally, which is dramatically
     faster than the general method.  The $\mathsf{E}_8$ singularity:
    Example
     S = QQ[x,y,z]
     f = x^2+y^3+z^5
     w = {1/2, 1/3, 1/5}
     C1 = gradedDeRhamComplexH1(f,-2,w)
     prune HH_0(C1)
  SeeAlso
     intersectionDuBoisComplex
     duBoisComplex
     gradedDeRhamCohomologyH1
     hodgeIdealWeightedHomogIsolated

///



doc ///
  Key
    gradedDeRhamCohomologyH1
    (gradedDeRhamCohomologyH1,RingElement,ZZ,ZZ)
    (gradedDeRhamCohomologyH1,RingElement,ZZ,ZZ,List)
  Headline
    compute cohomology of $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$
  Usage
    H = gradedDeRhamCohomologyH1(f,p,q)
    H = gradedDeRhamCohomologyH1(f,p,q,w)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    p: ZZ
     a non-positive integer (Hodge filtration index)
    q: ZZ
     a non-positive integer (Cohomological degree)
    w: List
     (optional) a list of positive weights $\{w_1,\ldots,w_n\}$ for which
     $f$ is weighted-homogeneous of weight $1$ with an isolated singularity
  Outputs
    H: Module
     the q-th cohomology of $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$
  Description
    Text
     Given a reduced polynomial $f\in S=\mathbb{Q}[x_1,\cdots,x_n]$ and $p\geq -n$ and $q\geq -n$,
     this function computes the $q$-th cohomology $\operatorname{Gr}^F_p(\operatorname{DR}(H^1_f(S)))$
     (which is a co-chain complex lying in cohomological degrees $\geq -n$). This function is more
     efficient than gradedDeRhamComplexH1, as it does not construct the entire complex.
    Example
     S = QQ[x,y,z]
     f = x^2+y^3+z^2
     gradedDeRhamCohomologyH1(f,-3,0)
     gradedDeRhamCohomologyH1(f,-2,0)
     gradedDeRhamCohomologyH1(f,-2,-1)
     gradedDeRhamCohomologyH1(f,-1,-2)
    Text
     Compare to the calculation in @TO "gradedDeRhamComplexH1"@.

     When $f$ is weighted-homogeneous of weight $1$ with respect to weights $w$
     and has an isolated singularity at the origin, passing $w$ uses
     @TO "hodgeIdealWeightedHomogIsolated"@ internally.  The $\mathsf{E}_8$
     singularity:
    Example
     S = QQ[x,y,z]
     f = x^2+y^3+z^5
     w = {1/2, 1/3, 1/5}
     gradedDeRhamCohomologyH1(f,-2,0,w)
  SeeAlso
     intersectionDuBoisComplex
     gradedDeRhamComplexH1
     duBoisComplex
     hodgeIdealWeightedHomogIsolated

///


doc ///
  Key
    isPreDuBois
    (isPreDuBois, RingElement, ZZ)
    (isPreDuBois, RingElement, ZZ, List)
  Headline
    determine if $f$ is pre $m$-Du Bois for some $m\geq 0$
  Usage
    B = isPreDuBois(f,m)
    B = isPreDuBois(f,m,w)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    m: ZZ
     a non-negative integer
    w: List
     (optional) a list of positive weights $\{w_1,\ldots,w_n\}$ for which
     $f$ is weighted-homogeneous of weight $1$ with an isolated singularity
  Outputs
    B: Boolean
  Description
    Text
     The divisor $D=V(f)$ is pre $m$-Du Bois if $\mathcal{H}^i(\underline{\Omega}^p_D)=0$ for all $p\leq m$, $i>0$.
     This notion was studied in [SVV]. It is a necessary condition for $D$ to have $m$-Du Bois singularities.

     The $2\times 2$ generic determinant is pre $m$-Du Bois for all $m\geq 0$:
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     isPreDuBois(f,0)
     isPreDuBois(f,1)
     isPreDuBois(f,2)
    Text
     On the other hand, the cusp is not pre $0$-Du Bois:
    Example
     S = QQ[x,y];
     f = x^2+y^3;
     isPreDuBois(f,0)
    Text
     When $f$ is weighted-homogeneous of weight $1$ with respect to weights $w$
     and has an isolated singularity at the origin, passing $w$ uses
     @TO "hodgeIdealWeightedHomogIsolated"@ internally.  The $\mathsf{E}_8$
     singularity is pre $1$-Du Bois:
    Example
     S = QQ[x,y,z];
     f = x^2+y^3+z^5;
     w = {1/2, 1/3, 1/5};
     isPreDuBois(f,1,w)
  SeeAlso
     duBoisComplex
     hodgeIdealWeightedHomogIsolated
  References
    See [SVV] at @TO "Works Cited"@.

///



doc ///
  Key
    intersectionDuBoisComplex
    (intersectionDuBoisComplex, RingElement, ZZ)
  Headline
    compute the intersection Du Bois complexes $I\underline{\Omega}^p_D$ for a divisor
  Usage
    C = intersectionDuBoisComplex(f,p)
  Inputs
    f: RingElement
     a polynomial with rational coefficients
    p: ZZ
     a non-negative integer 
  Outputs
    C: Complex
     a free complex quasi-isomorphic to $I\underline{\Omega}^p_D$
  Description
    Text
     Given a reduced polynomial $f\in S=\mathbb{Q}[x_1,\cdots,x_n]$ and $0\leq p\leq n-1$, this function computes a free complex
     quasi-isomorphic to $I\underline{\Omega}^p_D$, where $D=V(f)$. As complexes in Macaulay2
     are homologically graded, the $-i$-th homology of this complex is the $i$-th
     cohomology of $I\underline{\Omega}^p_D$. The homology of this complex lives in non-positive degrees.

     If $D$ is a rational homology manifold, this complex is quasi-isomorphic to $\underline{\Omega}^p_D$:
    Example
     S = QQ[x,y,z];
     f = y^2-x*z;
     DB0 = duBoisComplex(f,0)
     IDB0 = intersectionDuBoisComplex(f,0)
     prune HH_0(DB0)
     prune HH_0(IDB0)
     DB1 = duBoisComplex(f,1)
     IDB1 = intersectionDuBoisComplex(f,1)
     prune HH_0(DB1)
     prune HH_0(IDB1)
     DB2 = duBoisComplex(f,2)
     IDB2 = intersectionDuBoisComplex(f,2)
     prune HH_0(DB2)
     prune HH_0(IDB2)
    Text

     We carry out another example.
    Example
     S = QQ[x,y,z,w];
     f = x*w-y*z;
     IDB0 = intersectionDuBoisComplex(f,0)
     prune HH_0(IDB0)
     IDB1 = intersectionDuBoisComplex(f,1)
     prune HH_0(IDB1)
     prune HH_(-1)(IDB1)
     IDB2 = intersectionDuBoisComplex(f,2)
     prune HH_0(IDB2)
     prune HH_(-1)(IDB2)
  SeeAlso
     duBoisComplex
     gradedDeRhamComplexH1
     gradedDeRhamCohomologyH1
  References
     See [PR21] in @TO "Works Cited"@.
/// 

--symbols:

doc ///
  Key
    ByAnnFs
  Headline
    option value for NuMethod computing nuAlpha via annihilators
  Usage
    NuMethod => ByAnnFs
  Description
    Text
      When NuMethod is set to ByAnnFs, nuAlpha is computed using Ann(f^s)
      (and, when needed, annihilators involving g). This is the default method used by nuAlpha.
  SeeAlso
    NuMethod
    nuAlpha
    PowerBFunction
    Malgrange
///

doc ///
  Key
    ByNuAlpha
  Headline
    option value for LengthStrategy computing weightLength via nuAlpha
  Usage
    LengthStrategy => ByNuAlpha
  Description
    Text
      When LengthStrategy is set to ByNuAlpha, weightLength is computed using nuAlpha applied to the
      relevant translated roots of the b-function.
      This is the default strategy for weightLength.
  SeeAlso
    weightLength
    nuAlpha
    ByWeightLevel
///

doc ///
  Key
    ByWeightLevel
  Headline
    option value for LengthStrategy computing weightLength via weightLevel
  Usage
    Strategy => ByWeightLevel
  Description
    Text
      When LengthStrategy is set to ByWeightLevel, weightLength is computed by calling weightLevel at a
      translated exponent determined from the integer translates of roots of the b-function.
  SeeAlso
    weightLength
    weightLevel
    ByNuAlpha
///


doc ///
  Key
    dtBasis
  Headline
    option value for UseBasis selecting the dt-basis output
  Usage
    UseBasis => dtBasis
  Description
    Text
      When UseBasis is set to dtBasis, functions such as hodgeOnV, weightHodgeOnV, and monodromyWeightHodgeOnV.
      return bases expressed in a polynomial ring S[dt] (this is the default).
  SeeAlso
    UseBasis
    sBasis
    hodgeOnV
    monodromyWeightHodgeOnV
    weightHodgeOnV
///

doc ///
  Key
    False
  Headline
    option value disabling use of generation level in Hodge ideal
  Usage
    UseGenLevel => False
  Description
    Text
      The option UseGenLevel => False turns off use of general bounds for generation level
      when calculating Hodge ideals. This is useful for testing and debugging.
  SeeAlso
    True
    UseGenLevel
///



doc ///
  Key
    NuMethod
  Headline
    choose the method used to compute nuAlpha-related quantities
  Usage
    NuMethod => ByAnnFs
    NuMethod => PowerBFunction
    NuMethod => Malgrange
  Description
    Text
      This option is supported by nuAlpha, pFunction, and weightLength.
  SeeAlso
    nuAlpha
    pFunction
    weightLength
    ByAnnFs
    PowerBFunction
    Malgrange
///


doc ///
  Key
    LengthStrategy
  Headline
    choose an algorithmic strategy (used by weightLength)
  Usage
    LengthStrategy => ByNuAlpha
    LengthStrategy => ByWeightLevel
  Description
    Text
      This option is currently used by weightLength.
  SeeAlso
    weightLength
    ByNuAlpha
    ByWeightLevel
///

doc ///
  Key
    Malgrange
  Headline
    option value for NuMethod computing nuAlpha using the Malgrange ideal (g = 1_S)
  Usage
    NuMethod => Malgrange
  Description
    Text
      When NuMethod is set to Malgrange, nuAlpha is computed using the Malgrange ideal and an
      initial ideal membership test. This strategy is implemented only when g = 1_S.
  SeeAlso
    NuMethod
    nuAlpha
    ByAnnFs
    PowerBFunction
///

doc ///
  Key
    PowerBFunction
  Headline
    option value for NuMethod computing nuAlpha using b-functions of powers
  Usage
    NuMethod => PowerBFunction
  Description
    Text
      When NuMethod is set to PowerBFunction, nuAlpha is computed by passing to a sufficiently high
      power f^{k0} and reading the multiplicity of the relevant linear factor from the b-function
      of that power.
  SeeAlso
    NuMethod
    nuAlpha
    ByAnnFs
    Malgrange
///

doc ///
  Key
    sBasis
  Headline
    option value for UseBasis selecting the s-basis output
  Usage
    UseBasis => sBasis
  Description
    Text
      When UseBasis is set to sBasis, functions such as hodgeOnV, weightHodgeOnV, and monodromyWeightHodgeOnV.
      return bases in the s-variable model S_f[s]f^s, which may introduce denominators.
  SeeAlso
    UseBasis
    dtBasis
    hodgeOnV
    monodromyWeightHodgeOnV
    weightHodgeOnV
///


doc ///
  Key
    True
  Headline
    option value enabling UseGenLevel in hodgeIdeal
  Usage
    UseGenLevel => True
  Description
    Text
      The option UseGenLevel => True enables calculation of Hodge ideal using
      general generation level bounds. This is the default.
  SeeAlso
    False
    UseGenLevel
///


doc ///
  Key
    UseBasis
  Headline
    choose a basis in Brieskorn lattice computations
  Usage
    UseBasis => dtBasis
    UseBasis => sBasis
  Description
    Text
      This option is supported by hodgeOnV, weightHodgeOnV, and monodromyWeightHodgeOnV.
  SeeAlso
    hodgeOnV
    monodromyWeightHodgeOnV
    weightHodgeOnV
    dtBasis
    sBasis
///


doc ///
  Key
    UseGenLevel
  Headline
    toggle generation-level shortcuts in Hodge ideal computations
  Usage
    UseGenLevel => True
    UseGenLevel => False
  Description
    Text
      This option is supported by hodgeIdeal.
  SeeAlso
    hodgeIdeal
///
