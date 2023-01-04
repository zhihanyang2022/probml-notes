<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Bayesian linear regression with fixed noise
  precision>|<doc-date|Jan 2, 2023>>

  <section|Preliminaries>

  <\definition>
    (Linear Gaussian system; originally 4.124 on p119 of Murphy)

    Let <math|<with|font-series|bold|x>\<in\>\<bbb-R\><rsup|D<rsub|x>>> and
    <math|<with|font-series|bold|y>\<in\>\<bbb-R\><rsup|D<rsub|y>>> be two
    random variables. The following generative model

    <\eqnarray*>
      <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|x>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|\<mu\>><rsub|x>,<with|font-series|bold|\<Sigma\>><rsub|x>|)>>>|<row|<cell|p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|x>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|A><with|font-series|bold|x>+<with|font-series|bold|b>,<with|font-series|bold|\<Sigma\>><rsub|y>|)>>>>>
    </eqnarray*>

    is the linear Gaussian system, where <math|<with|font-series|bold|A>\<in\>\<bbb-R\><rsup|D<rsub|y>\<times\>D<rsub|x>>>
    and <math|<with|font-series|bold|b>\<in\>\<bbb-R\><rsup|D<rsub|y>>>.
  </definition>

  <\theorem>
    <label|th-bayes>(Bayes rule for linear Gaussian sytems, originally 4.125
    on pp119 of Murphy)

    Given a linear Gaussian system, the posterior is given by

    <\eqnarray*>
      <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|y>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|\<mu\>><rsub|post>,<with|font-series|bold|\<Sigma\>><rsub|post>|)>>>|<row|<cell|<with|font-series|bold|\<Sigma\>><rsub|post><rsup|-1>>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|x><rsup|-1>+<with|font-series|bold|A><rsup|T><with|font-series|bold|\<Sigma\>><rsub|y><rsup|-1><with|font-series|bold|A>>>|<row|<cell|<with|font-series|bold|\<mu\>><rsub|post>>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|post><around*|[|<with|font-series|bold|A><rsup|T><with|font-series|bold|\<Sigma\>><rsub|y><rsup|-1><around*|(|<with|font-series|bold|y>-<with|font-series|bold|b>|)>+<with|font-series|bold|\<Sigma\>><rsub|x><rsup|-1><with|font-series|bold|\<mu\>><rsub|x>|]>>>>>
    </eqnarray*>

    Note that the posterior uncertainty (i.e., the covariance is independent
    on <math|<with|font-series|bold|y>>) and the posterior mean is a linear
    function of <math|<with|font-series|bold|y>>.
  </theorem>

  <section|Core>

  <with|font-series|bold|Likelihood. >

  <\equation*>
    p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X>,<with|font-series|bold|w>,\<beta\>|)>=<big|prod><rsub|i=1><rsup|N><with|font|cal|N><around*|(|y<rsub|i>\<mid\>\<mu\><rsub|i>,\<sigma\><rsup|2>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X><with|font-series|bold|w>,\<sigma\><rsup|2><with|font-series|bold|I>|)>
  </equation*>

  with <math|<with|font-series|bold|\<mu\>>=<with|font-series|bold|X><with|font-series|bold|w>>.\ 

  <with|font-series|bold|Prior.>\ 

  <\equation*>
    p<around*|(|<with|font-series|bold|w>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|0>,<with|font-series|bold|V><rsub|0>|)>
  </equation*>

  <with|font-series|bold|Posterior.>

  The key insight is that <math|<with|font-series|bold|X>> transforms
  <math|<with|font-series|bold|w>> into the mean of
  <math|p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X>,<with|font-series|bold|w>,\<beta\>|)>>.
  So we can apply Theorem <reference|th-bayes> to obtain the posterior:

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>,\<sigma\><rsup|2>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|w><rsub|N>,<with|font-series|bold|V><rsub|N>|)>>>|<row|<cell|<with|font-series|bold|w><rsub|N>>|<cell|=>|<cell|<with|font-series|bold|V><rsub|N><around*|(|<with|font-series|bold|X><rsup|T><around*|(|\<sigma\><rsup|2><with|font-series|bold|I>|)><rsup|-1><with|font-series|bold|y>+<with|font-series|bold|V><rsub|0><rsup|-1><with|font-series|bold|w><rsub|0>|)>>>|<row|<cell|>|<cell|=>|<cell|<with|font-series|bold|V><rsub|N><around*|(|<with|font-series|bold|V><rsub|0><rsup|-1><with|font-series|bold|w><rsub|0>+<around*|(|1/\<sigma\><rsup|2>|)><with|font-series|bold|X><rsup|T><with|font-series|bold|y>|)>>>|<row|<cell|<with|font-series|bold|V><rsub|N>>|<cell|=>|<cell|<around*|(|<with|font-series|bold|V><rsub|0><rsup|-1>+<with|font-series|bold|X><rsup|T><around*|(|\<sigma\><rsup|2><with|font-series|bold|I>|)><rsup|-1><with|font-series|bold|X>|)><rsup|-1>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<with|font-series|bold|V><rsub|0><rsup|-1>+<around*|(|1/\<sigma\><rsup|2>|)><with|font-series|bold|X><rsup|T><with|font-series|bold|X>|)><rsup|-1>>>>>
  </eqnarray*>

  <with|font-series|bold|Posterior predictive.>

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~>,<with|font-series|bold|X>,<with|font-series|bold|y>,\<sigma\><rsup|2>|)>>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w>,\<sigma\><rsup|2><with|font-series|bold|I><rsub|m>|)>
    <with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|><with|font-series|bold|w><rsub|N>,<with|font-series|bold|V><rsub|N>|)>
    d<with|font-series|bold|w>>>|<row|<cell|>|<cell|=>|<cell|<big|int><with|font|cal|N><around*|(|<wide|<with|font-series|bold|y>|~>\<mid\><wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|w><rsub|N>,\<sigma\><rsup|2><with|font-series|bold|I><rsub|m>+<wide|<with|font-series|bold|X>|~><rsup|T><with|font-series|bold|V><rsub|N><wide|<with|font-series|bold|X>|~>|)><space|1em><around*|(|<text|By
    Eq 4.126>|)>>>>>
  </eqnarray*>
</body>

<\initial>
  <\collection>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|1>>
    <associate|th-bayes|<tuple|2|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Preliminaries>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Core>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>