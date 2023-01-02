<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Bayesian linear regression with fixed noise
  precision>|<doc-date|Jan 2, 2023>>

  <section|Preliminaries>

  <\definition>
    (Linear Gaussian system; originally 4.124 on p119)

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
    on pp119)

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
    p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X>,<with|font-series|bold|w>,\<beta\>|)>=<big|prod><rsub|i=1><rsup|N><with|font|cal|N><around*|(|y<rsub|i>\<mid\>\<mu\><rsub|i>,\<beta\><rsup|-1>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X><with|font-series|bold|w>,\<beta\><rsup|-1><with|font-series|bold|I>|)>
  </equation*>

  with <math|<with|font-series|bold|\<mu\>>=<with|font-series|bold|X><with|font-series|bold|w>>.
  Let's call <math|\<beta\>> the noise precision.

  <with|font-series|bold|Prior.>\ 

  <\equation*>
    p<around*|(|<with|font-series|bold|w>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|0>,\<alpha\><rsup|-1><with|font-series|bold|I>|)>
  </equation*>

  Let's call <math|\<alpha\>> the prior precision. Of course, we could also
  have used a Gaussian with an arbitrary mean and covariance. But let's keep
  it simple for now.

  <with|font-series|bold|Posterior.>

  The key insight is that <math|<with|font-series|bold|X>> transforms
  <math|<with|font-series|bold|w>> into the mean of
  <math|p<around*|(|<with|font-series|bold|y>\<mid\><with|font-series|bold|X>,<with|font-series|bold|w>,\<beta\>|)>>.
  So we can apply Theorem <reference|th-bayes> to obtain the posterior:

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|w>\<mid\><with|font-series|bold|y>,<with|font-series|bold|X>,\<beta\><rsup|-1>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|x>\<mid\><with|font-series|bold|\<mu\>><rsub|post>,<with|font-series|bold|\<Sigma\>><rsub|post>|)>>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|<with|font-series|bold|\<Sigma\>><rsub|post><rsup|-1>>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|w><rsup|-1>+<with|font-series|bold|X><rsup|T><with|font-series|bold|\<Sigma\>><rsub|y><rsup|-1><with|font-series|bold|X>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|\<alpha\><rsup|-1><with|font-series|bold|I>|)><rsup|-1>+<with|font-series|bold|X><rsup|T><around*|(|\<beta\><rsup|-1><with|font-series|bold|I>|)><rsup|-1><with|font-series|bold|X>>>|<row|<cell|>|<cell|=>|<cell|\<alpha\><with|font-series|bold|I>+\<beta\><with|font-series|bold|X><rsup|T><with|font-series|bold|X>>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|<with|font-series|bold|\<mu\>><rsub|post>>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|post><around*|[|<with|font-series|bold|X><rsup|T><with|font-series|bold|\<Sigma\>><rsub|y><rsup|-1><around*|(|<with|font-series|bold|y>-<with|font-series|bold|0>|)>+<with|font-series|bold|\<Sigma\>><rsub|w><rsup|-1><with|font-series|bold|\<mu\>><rsub|w>|]>>>|<row|<cell|>|<cell|=>|<cell|<with|font-series|bold|\<Sigma\>><rsub|post><around*|[|<with|font-series|bold|X><rsup|T><around*|(|\<beta\><rsup|-1><with|font-series|bold|I>|)><rsup|-1><with|font-series|bold|y>+<around*|(|\<alpha\><rsup|-1><with|font-series|bold|I>|)><rsup|-1><with|font-series|bold|0>|]>>>|<row|<cell|>|<cell|=>|<cell|\<beta\><with|font-series|bold|\<Sigma\>><rsub|post><with|font-series|bold|X><rsup|T><with|font-series|bold|y>>>>>
  </eqnarray*>

  <section|Experiment>

  \;
</body>

<\initial>
  <\collection>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
    <associate|th-bayes|<tuple|2|?>>
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

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Experiment>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>