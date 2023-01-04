<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Bayesian linear regression with unknown noise
  precision>|<\doc-date>
    Janurary 2, 2023
  </doc-date>>

  For Bayesian linear regression with fixed noise precision, we relied on a
  theorem for linear Gaussian systems to go from the prior and likelihood to
  the posterior. Here, we first state and prove a new theorem that's very
  similar to that theorem, which will be crucially for helping us go from the
  prior and likelihood to the posterior for Bayesian linear regression with
  unknown noise precision.

  <\theorem>
    Given the following prior and likelihood

    <\eqnarray*>
      <tformat|<table|<row|<cell|p<around*|(|\<sigma\><rsup|2>|)>>|<cell|=>|<cell|IG<around*|(|a<rsub|\<sigma\>>,b<rsub|\<sigma\>>|)>>>|<row|<cell|p<around*|(|<with|font-series|bold|x>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|\<mu\>><rsub|x>,\<sigma\><rsup|2><with|font-series|bold|\<Sigma\>><rsub|x>|)>,>>|<row|<cell|p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|x>,\<sigma\><rsup|2>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|A><with|font-series|bold|x>+<with|font-series|bold|b>,<with|font-series|bold|\<Sigma\>><rsub|y>|)>>>>>
    </eqnarray*>

    the posterior is given by

    <\eqnarray*>
      <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|x>,\<sigma\><rsup|2>\<mid\><with|font-series|bold|y>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|\<mu\>><rsub|x\<mid\>y>,\<sigma\><rsup|2><with|font-series|bold|\<Sigma\>><rsub|x\<mid\>y>|)>
      IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|\<sigma\>\<mid\>y>,b<rsub|\<sigma\>\<mid\>y>|)>>>|<row|<cell|<with|font-series|bold|\<mu\>><rsub|x\<mid\>y>>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|x\<mid\>y><around*|(|<with|font-series|bold|\<Sigma\>><rsub|x><rsup|-1><with|font-series|bold|\<mu\>><rsub|x>+<with|font-series|bold|A><rsup|T><with|font-series|bold|y>|)>>>|<row|<cell|<with|font-series|bold|\<Sigma\>><rsub|x\<mid\>y>>|<cell|=>|<cell|<around*|(|<with|font-series|bold|\<Sigma\>><rsub|x><rsup|-1>+<with|font-series|bold|A><rsup|T><with|font-series|bold|A>|)><rsup|-1>>>|<row|<cell|a<rsub|\<sigma\>\<mid\>y>>|<cell|=>|<cell|a<rsub|\<sigma\>>+D<rsub|y>/2>>|<row|<cell|b<rsub|\<sigma\>\<mid\>y>>|<cell|=>|<cell|b<rsub|\<sigma\>>+<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|x><rsup|T><with|font-series|bold|\<Sigma\>><rsub|x><rsup|-1><with|font-series|bold|\<mu\>><rsub|x>+<with|font-series|bold|y><rsup|T><with|font-series|bold|y>-<with|font-series|bold|\<mu\>><rsub|x\<mid\>y><rsup|T>
      <with|font-series|bold|\<Sigma\>><rsub|x\<mid\>y><rsup|-1>
      <with|font-series|bold|\<mu\>><rsub|x\<mid\>y>|)>>>>>
    </eqnarray*>
  </theorem>

  <\proof>
    The proof should be analagous to the proof of the Bayes rule for linear
    Gaussian systems presented in Section 4.4.3 of Murphy.
  </proof>

  <with|font-series|bold|Likelihood. >

  <\equation*>
    p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X>,<with|font-series|bold|w>,\<sigma\><rsup|2>|)>=<big|prod><rsub|i=1><rsup|N><with|font|cal|N><around*|(|y<rsub|i>\<mid\>\<mu\><rsub|i>,\<sigma\><rsup|2>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X><with|font-series|bold|w>,\<sigma\><rsup|2><with|font-series|bold|I>|)>
  </equation*>

  with <math|<with|font-series|bold|\<mu\>>=<with|font-series|bold|X><with|font-series|bold|w>>.\ 

  <with|font-series|bold|Prior.>

  <\equation*>
    p<around*|(|<with|font-series|bold|w>,\<sigma\><rsup|2>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|0>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|0>|)>
    IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|0>,b<rsub|0>|)>
  </equation*>

  <with|font-series|bold|Posterior.> Applying Theorem 1, we obtain:

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|w>,\<sigma\><rsup|2>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|N>|)>
    IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>>>|<row|<cell|<with|font-series|bold|w><rsub|N>>|<cell|=>|<cell|<with|font-series|bold|V><rsub|N><around*|(|<with|font-series|bold|V><rsub|0><rsup|-1><with|font-series|bold|w><rsub|0>+<with|font-series|bold|X><rsup|T><with|font-series|bold|y>|)>>>|<row|<cell|<with|font-series|bold|V><rsub|N>>|<cell|=>|<cell|<around*|(|<with|font-series|bold|V><rsub|0><rsup|-1>+<with|font-series|bold|X><rsup|T><with|font-series|bold|X>|)><rsup|-1>>>|<row|<cell|a<rsub|N>>|<cell|=>|<cell|a<rsub|0>+N/2>>|<row|<cell|b<rsub|N>>|<cell|=>|<cell|b<rsub|0>+<frac|1|2><around*|(|<with|font-series|bold|w><rsub|0><rsup|T><with|font-series|bold|V><rsub|0><rsup|-1><with|font-series|bold|w><rsub|0>+<with|font-series|bold|y><rsup|T><with|font-series|bold|y>-<with|font-series|bold|w><rsub|N><rsup|T><with|font-series|bold|V><rsub|N><rsup|-1><with|font-series|bold|w><rsub|N>|)>,>>>>
  </eqnarray*>

  which is exactly the same as Equation 7.69 to 7.73 from Murphy.

  <with|font-series|bold|Marginal posterior.>\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>|)>>|<cell|=>|<cell|<big|int>p<around*|(|<with|font-series|bold|w>,\<sigma\><rsup|2>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>|)>
    d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|N>|)>
    IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,<around*|(|1/\<sigma\><rsup|2>|)><with|font-series|bold|V><rsub|N>|)>
    Ga<around*|(|1/\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<with|font|cal|T><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,<around*|(|b<rsub|N>/a<rsub|N>|)><with|font-series|bold|V><rsub|N>,2a<rsub|N>|)><space|1em><around*|(|<text|by
    Eq 11.61>|)>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|\<sigma\><rsup|2>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>|)>>|<cell|=>|<cell|<big|int>p<around*|(|<with|font-series|bold|w>,\<sigma\><rsup|2>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>|)>
    d<with|font-series|bold|w>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|N>|)>
    IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    d<with|font-series|bold|w>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|N>|)>
    \ d<with|font-series|bold|w> IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>>>|<row|<cell|>|<cell|=>|<cell|IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>>>>>
  </eqnarray*>

  <with|font-series|bold|Posterior predictive.>\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~>,<with|font-series|bold|X>,<with|font-series|bold|y>|)>>|<cell|=>|<cell|<big|int><big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w>,\<sigma\><rsup|2><with|font-series|bold|I><rsub|m>|)>
    NIW<around*|(|<with|font-series|bold|w>,\<sigma\><rsup|2>\<mid\><with|font-series|bold|><with|font-series|bold|w><rsub|N>,<with|font-series|bold|V><rsub|N>,a<rsub|N>,b<rsub|N>|)>
    d<with|font-series|bold|w> d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|int><around*|[|<big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w>,\<sigma\><rsup|2><with|font-series|bold|I><rsub|m>|)>
    <with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|V><rsub|N>|)>
    d<with|font-series|bold|w>|]> IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    \ d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><around*|(|<with|font-series|bold|I><rsub|m>+<wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|V><rsub|N><wide|<with|font-series|bold|X>|~>|)>|)>
    IG<around*|(|\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    \ d\<sigma\><rsup|2><space|1em><around*|(|<text|by Eq
    4.126>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w><rsub|N>,<around*|(|1/\<sigma\><rsup|2>|)><around*|(|<with|font-series|bold|I><rsub|m>+<wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|V><rsub|N><wide|<with|font-series|bold|X>|~>|)>|)>
    Ga<around*|(|1/\<sigma\><rsup|2>\<mid\>a<rsub|N>,b<rsub|N>|)>
    \ d\<sigma\><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<with|font|cal|T><around*|(|<wide|<with|font-series|bold|y>|~><mid|\|><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w><rsub|N>,<frac|b<rsub|N>|a<rsub|N>><around*|(|<with|font-series|bold|I><rsub|m>+<wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|V><rsub|N><wide|<with|font-series|bold|X>|~>|)>,2a<rsub|N>|)><space|1em><around*|(|<text|by
    Eq 11.61>|)>>>>>
  </eqnarray*>

  \;

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|page-screen-margin|false>
  </collection>
</initial>