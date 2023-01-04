<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Mean-field Variational Inference for Bayesian Mixture
  of Gaussian>|<doc-date|December 27, 2022>>

  <abstract-data|<abstract|Notes for the Bayesian Mixture of Gaussians
  section in Bishop's PRML. Compared to the book's treatment, this derivation
  more detailed. When choosing priors, we rely on Theorem 2.2 in Beal's PhD
  thesis, i.e., the priors should be conjugate to the complete-data log
  likelihood.>>

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|1<space|2spc>Complete-data
    log likelihood> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|2<space|2spc>Prior
    for <with|mode|math|\<pi\>>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|3<space|2spc>Joint
    prior for <with|mode|math|<with|font-series|bold|\<mu\>><rsub|k>> and
    <with|mode|math|<with|font-series|bold|\<Sigma\>><rsub|k>>:>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|4<space|2spc>Computing
    <with|mode|math|q<around*|(|Z|)>>: variational E-step>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|5<space|2spc>Computing
    <with|mode|math|q<around*|(|\<pi\>,\<mu\><rsub|1:C>,\<Lambda\><rsub|1:C>|)>>:
    variational M-step> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5><vspace|0.5fn>

    <with|par-left|1tab|5.1<space|2spc><with|mode|math|q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>>

    <with|par-left|1tab|5.2<space|2spc><with|mode|math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|6<space|2spc>Quantities
    required for computing <with|mode|math|q<around*|(|Z|)>>>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-8><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|7<space|2spc>Quantities
    required for computing <with|mode|math|><with|mode|math|q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>
    and <with|mode|math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-9><vspace|0.5fn>
  </table-of-contents>

  <section|Complete-data log likelihood>

  Just write out the mixture of Gaussians assuming that the latent variables
  are observed:

  <\equation*>
    <big|prod><rsub|i=1><rsup|n>p<around*|(|<with|font-series|bold|x><rsub|i>,<with|font-series|bold|z><rsub|i>\<mid\><with|font-series|bold|\<theta\>>|)>=
    <big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>> <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)><rsup|z<rsub|i
    k>>
  </equation*>

  where <math|<with|font-series|bold|\<theta\>>=<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1:k>,<with|font-series|bold|\<Lambda\>><rsub|1:k>|)>>.

  <section|Prior for <math|\<pi\>>>

  Dirichlet:

  <\equation*>
    Dir<around*|(|<with|font-series|bold|\<pi\>>\<mid\><with|font-series|bold|\<alpha\>><rsub|0>|)>\<propto\><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>
  </equation*>

  Relevant part of the complete-data log likelihood:

  <\equation*>
    <big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>>
  </equation*>

  Proof that they are conjugate:

  <\equation*>
    <around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>|)><around*|(|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>>|)>=<around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>|)><around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<big|sum><rsub|i=1><rsup|n>z<rsub|i
    k>>|)>=<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>+<big|sum><rsub|i=1><rsup|n>z<rsub|i
    k>-1>
  </equation*>

  <section|Joint prior for <math|<with|font-series|bold|\<mu\>><rsub|k>> and
  <math|<with|font-series|bold|\<Sigma\>><rsub|k>>:>

  Gaussian-Wishart:

  <\equation*>
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>\<mid\><with|font-series|bold|m><rsub|0>,\<beta\><rsub|0>,v<rsub|0>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|m><rsub|0>,<around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><rsup|-1>|)>
    Wi<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>\<mid\><with|font-series|bold|L><rsub|0>,v<rsub|0>|)>
  </equation*>

  Relevant part of the complete-data log likelihood:

  <\equation*>
    <big|prod><rsub|i=1><rsup|n><with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)><rsup|z<rsub|i
    k>>=<big|prod><rsub|i<text| s.t. >z<rsub|i,k=1>><with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)>
  </equation*>

  Proof that they are conjugate to each other: omitted, see Section 4.6.3.3
  of Murphy.

  <section|Computing <math|q<around*|(|Z|)>>: variational E-step>

  For notationaly clarity, we let <math|<with|font-series|bold|\<mu\>>=<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>>
  and <math|<with|font-series|bold|\<Lambda\>>=<around*|{|<with|font-series|bold|\<Lambda\>><rsub|k>|}>>.

  Starting from the CAVI update rule:

  <\eqnarray*>
    <tformat|<cwith|8|10|3|3|font-base-size|9>|<table|<row|<cell|>|<cell|>|<cell|log
    q<around*|(|<with|font-series|bold|Z>|)>>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>,<with|font-series|bold|Z>,<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>+log
    p<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log \<pi\><rsub|k>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>><around*|[|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|<with|font-series|bold|\<Lambda\>>><rsub|k><rsup|-1>|)><rsup|>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<frac|1|2>
    ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> <around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> <around*|[|<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k><around*|(|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>|)>+const>>|<row|<cell|>|<cell|\<triangleq\>>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log \<rho\><rsub|i k>>>>>
  </eqnarray*>

  Therefore

  <\equation*>
    <tabular|<tformat|<table|<row|<cell|q<around*|(|<with|font-series|bold|Z>|)>>|<cell|\<propto\>>|<cell|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<rho\><rsub|i
    k><rsup|z<rsub|i k>>>>>>>
  </equation*>

  Requring that the distributuon is normalized on each row:

  <\equation*>
    <tabular|<tformat|<table|<row|<cell|q<around*|(|<with|font-series|bold|Z>|)>=>|<cell|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>r<rsub|i
    k><rsup|z<rsub|i k>>>>>>>
  </equation*>

  where <math|r<rsub|i k>> is the <math|\<rho\><rsub|i k>> divided by row
  sum.\ 

  <section|Computing <math|q<around*|(|\<pi\>,\<mu\><rsub|1:C>,\<Lambda\><rsub|1:C>|)>>:
  variational M-step>

  Starting from the CAVI update rule:

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|ln
    q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>,<with|font-series|bold|Z>,<with|font-series|bold|X>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|>,<with|font-series|bold|\<Lambda\>>|)>+log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>,<with|font-series|bold|\<mu\>>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|>,<with|font-series|bold|\<Lambda\>><rsub|>|)>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>,<with|font-series|bold|\<mu\>><rsub|>,<with|font-series|bold|\<Lambda\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<around*|(|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>|)>+<big|sum><rsub|k=1><rsup|C><around*|(|p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> log <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)>|)>+const,<eq-number><label|logqparam-decompose>>>>>
  </eqnarray*>

  which contains separate terms for <math|<with|font-series|bold|\<pi\>>> and
  each <math|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>.
  This implies that the approximate posterior further factors into

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1>,\<ldots\>,.<with|font-series|bold|\<mu\>><rsub|C>,<with|font-series|bold|\<Lambda\>><rsub|1>,\<ldots\>,<with|font-series|bold|\<Lambda\><rsub|>><rsub|C>|)>>|<cell|=>|<cell|log
    q<around*|(|<with|font-series|bold|\<pi\>>|)>+<big|sum><rsub|k=1><rsup|C>ln
    q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>>|<row|<cell|q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1>,\<ldots\>,.<with|font-series|bold|\<mu\>><rsub|C>,<with|font-series|bold|\<Lambda\>><rsub|1>,\<ldots\>,<with|font-series|bold|\<Lambda\><rsub|>><rsub|C>|)>>|<cell|=>|<cell|q<around*|(|<with|font-series|bold|\<pi\>>|)><big|prod><rsub|k=1><rsup|C>q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>>>>
  </eqnarray*>

  so that we can optimize over each component distribution separately.

  <subsection|<math|q<around*|(|<with|font-series|bold|\<pi\>>|)>>>

  Isolating the terms containing <math|<with|font-series|bold|\<pi\>>> from
  Equation <reference|logqparam-decompose>, we see that

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>|<cell|=>|<cell|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|C><around*|(|\<alpha\><rsub|0>-1|)>
    ln \<pi\><rsub|k>+<big|sum><rsub|k=1><rsup|C><around*|(|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k>|)> ln \<pi\><rsub|k>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|k=1><rsup|C><around*|(|\<alpha\><rsub|0>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k>-1|)> ln \<pi\><rsub|k>+<big|sum><rsub|k=1><rsup|C> ln
    \<pi\><rsub|k>+const>>|<row|<cell|>|<cell|\<Rightarrow\>>|<cell|>>|<row|<cell|q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>|<cell|=>|<cell|Dir<around*|(|<with|font-series|bold|\<pi\>>\<mid\><with|font-series|bold|\<alpha\>>|)><space|1em><text|with><space|1em>\<alpha\><rsub|k>=\<alpha\><rsub|0>+N<rsub|k>>>>>
  </eqnarray*>

  <subsection|<math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>>

  Isolating the terms containing <math|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>
  from Equation <reference|logqparam-decompose>, we see that

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>|<cell|=>|<cell|log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>+log
    p<around*|(|<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> log <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)>+const>>>>
  </eqnarray*>

  By the Gaussian-Wishart prior we specified earlier, we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>>|<cell|=>|<cell|<with|font|cal|N><around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|m><rsub|0>,<around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><rsup|-1>|)>>>|<row|<cell|>|<cell|\<propto\>>|<cell|<around*|\||\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|1/2>
    exp<around*|(|-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>|)>>>|<row|<cell|>|<cell|\<propto\>>|<cell|<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|1/2>
    exp<around*|(|-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>|)>>>|<row|<cell|>|<cell|\<Rightarrow\>>|<cell|>>|<row|<cell|log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>>|<cell|=>|<cell|<frac|1|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>+const>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|p<around*|(|<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>>|<cell|=>|<cell|<with|font|cal|W><around*|(|<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>\<mid\><with|font-series|bold|W><rsub|0>,\<upsilon\><rsub|0>|)>>>|<row|<cell|>|<cell|\<propto\>>|<cell|<around*|\||<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|\|><rsup|<around*|(|\<upsilon\><rsub|0>-C-1|)>/2>
    exp<around*|(|-<frac|1|2>Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>|)>>>|<row|<cell|>|<cell|\<Rightarrow\>>|<cell|>>|<row|<cell|log
    p<around*|(|<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>>|<cell|=>|<cell|<frac|<around*|(|\<upsilon\><rsub|0>-C-1|)>|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>-<frac|1|2>
    Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)>+const.>>>>
  </eqnarray*>

  The normal likelihoods for all the observations can be written as

  <\equation*>
    <big|sum><rsub|i=1><rsup|N>r<rsub|i k> log
    <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k><rsup|-1>|)>=<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<frac|1|2> log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|)>+const.
  </equation*>

  Now, we can assemble all terms together:

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>|<cell|=>|<cell|<frac|1|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|<around*|(|\<upsilon\><rsub|0>-C-1|)>|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>-<frac|1|2>
    Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)>>>|<row|<cell|>|<cell|>|<cell|+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<frac|1|2> log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|)>+const.>>>>
  </eqnarray*>

  It may not be immediate clear what we should do here, so it's good to
  remind ourselves that by Theorem 2 of Beal tells us that the approximate
  posterior over <math|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>
  would also be Gaussian-Wishart.\ 

  Firstly, we see that

  <\equation*>
    <frac|<around*|(|\<upsilon\><rsub|0>-C-1|)>|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>+<frac|<around*|(|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k>|)>|2> log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>=<frac|<around*|(|\<upsilon\><rsub|0>+N<rsub|k>-C-1|)>|2>
    log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>,
  </equation*>

  which means that the Wishart posterior distribution on
  <math|<with|font-series|bold|\<Lambda\>><rsub|k>> have DOF
  <math|\<upsilon\><rsub|k>=\<upsilon\><rsub|0>+N-C>.

  Secondly, we need manipulate these terms

  <\equation*>
    -<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><space|1em>-<frac|1|2>
    Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)><space|2em><big|sum><rsub|i=1><rsup|N>r<rsub|i
    k><around*|(|-<frac|1|2><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|)>
  </equation*>

  into two terms of the form

  <\equation*>
    -<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-?<rsub|1>|)><rsup|T><around*|(|?<with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-?<rsub|1>|)><space|1em>-<frac|1|2>
    Tr<around*|(|?<rsub|2><with|font-series|bold|\<Lambda\>><rsub|k>|)>
  </equation*>

  where the <math|?<rsub|2>> in the second term shouldn't contain
  <math|<with|font-series|bold|\<mu\>><rsub|k>> because it would be a
  parameter of the Wishart posterior. For simplicity, we will drop the shared
  <math|-1/2> term shortly.

  To begin, we do an expansion with a trick:

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>-<around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>-<around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>-2
    \ <around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><wide*|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>|\<wide-underbrace\>><rsub|0>+N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>>>|<row|<cell|>|<cell|=>|<cell|N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>>>|<row|<cell|>|<cell|=>|<cell|N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>+
    Tr<around*|[|<around*|(|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>,>>>>
  </eqnarray*>

  which is great because the first term here contains
  <math|<with|font-series|bold|\<mu\>><rsub|k>> while the second term here
  doesn't and can be directly merged with
  <math|Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)>>.\ 

  We can merge the first term here with

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|\<beta\><rsub|0><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>+N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<beta\><rsub|0><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>+N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<beta\><rsub|0><around*|(|<with|font-series|bold|\<mu\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k><rsup|T>-<with|font-series|bold|\<mu\>><rsub|k><with|font-series|bold|m><rsub|0><rsup|T>-<with|font-series|bold|m><rsub|0><with|font-series|bold|\<mu\>><rsub|k><rsup|T>+<with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>|)>+N<rsub|k><around*|(|<with|font-series|bold|\<mu\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k><rsup|T>-<with|font-series|bold|\<mu\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k><rsup|T>+<wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><with|font-series|bold|\<mu\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k><rsup|T>-<with|font-series|bold|\<mu\>><rsub|k><around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>-<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><with|font-series|bold|\<mu\>><rsub|k><rsup|T>+\<beta\><rsub|0><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k><rsup|T>-<with|font-series|bold|\<mu\>><rsub|k><frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|\<beta\><rsub|0>+N<rsub|k>><with|font-series|bold|\<mu\>><rsub|k><rsup|T>|)>+\<beta\><rsub|0><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><rsup|T>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|><around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>+\<beta\><rsub|0><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>
    <around*|(|complete the squre|)>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|><around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>+\<beta\><rsub|0><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>-<frac|\<beta\><rsub|0><rsup|2><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|\<beta\><rsub|0>N<rsub|k><with|font-series|bold|m><rsub|0><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|\<beta\><rsub|0>N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><with|font-series|bold|m><rsub|0><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|N<rsub|k><rsup|2><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>+\<beta\><rsub|0><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>-<frac|\<beta\><rsub|0><rsup|2><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|\<beta\><rsub|0>N<rsub|k><with|font-series|bold|m><rsub|0><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|\<beta\><rsub|0>N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><with|font-series|bold|m><rsub|0><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|N<rsub|k><rsup|2><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>+<frac|\<beta\><rsub|0><rsup|2>+\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+<frac|\<beta\><rsub|0>N<rsub|k>+N<rsub|k><rsup|2>|\<beta\><rsub|0>+N<rsub|k>><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>-<frac|\<beta\><rsub|0>N<rsub|k><with|font-series|bold|m><rsub|0><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>-<frac|\<beta\><rsub|0>N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><with|font-series|bold|m><rsub|0><rsup|T>|\<beta\><rsub|0>+N<rsub|k>>+<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>><with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>+<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>+<around*|(|<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<with|font-series|bold|m><rsub|0><with|font-series|bold|m><rsub|0><rsup|T>-<with|font-series|bold|m><rsub|0><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k><with|font-series|bold|m><rsub|0><rsup|T>+<wide|<with|font-series|bold|x>|\<bar\>><rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k><rsup|T>|)>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>+<around*|(|<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|\<cdots\>+<around*|(|<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|Tr<around*|[|<around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><rsup|T>+<around*|(|<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><rsup|T><around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|<around*|(|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|>|\<beta\><rsub|0>+N<rsub|k>>|)><rsup|T>+Tr<around*|[|<around*|(|<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k>|]>,>>>>
  </eqnarray*>

  where the first term here has the required form of
  <math|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>-?<rsub|1>|)><rsup|T><around*|(|?<with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-?<rsub|1>|)>>
  and the second term can be directly merged with
  <math|Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)>>.

  Finally, we have

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|log
    q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>>|<row|<cell|>|<cell|=>|<cell|-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-
    <frac|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><with|font-series|bold|<wide|x|\<bar\>>><rsub|k>|\<beta\><rsub|0>+N>|)><rsup|T><around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<frac|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><with|font-series|bold|<wide|x|\<bar\>>><rsub|k>|\<beta\><rsub|0>+N>|)>+>>|<row|<cell|>|<cell|>|<cell|<frac|<around*|(|v<rsub|0>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k>-1|)>|2> log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>-<frac|1|2>
    Tr<around*|(|<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1>+N<rsub|k>
    <around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><around*|(|<with|font-series|bold|x><rsub|i>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>+<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>
    <around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>|)><with|font-series|bold|\<Lambda\>><rsub|k>|)>+const>>|<row|<cell|>|<cell|=>|<cell|log
    <with|font|cal|N><around*|(|<with|font-series|bold|\<mu\>><rsub|k><mid|\|><frac|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><with|font-series|bold|<wide|x|\<bar\>>><rsub|k>|\<beta\><rsub|0>+N<rsub|k>>,<around*|(|<around*|(|\<beta\><rsub|0>+N<rsub|k>|)><with|font-series|bold|\<Lambda\>><rsub|k>|)><rsup|-1>|)>+log
    <with|font|cal|W><around*|(|<with|font-series|bold|\<Lambda\>><rsub|k><mid|\|><around*|(|<with|font-series|bold|W><rsub|0><rsup|-1>+N<rsub|k>
    <around*|(|<with|font-series|bold|x><rsub|n>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><around*|(|<with|font-series|bold|x><rsub|n>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>+<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>
    <around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>|)><rsup|-1>,v<rsub|0>+N<rsub|k>|)>.>>>>
  </eqnarray*>

  Or we can define

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|font-series|bold|m><rsub|k>>|<cell|=>|<cell|<frac|\<beta\><rsub|0><with|font-series|bold|m><rsub|0>+N<rsub|k><with|font-series|bold|<wide|x|\<bar\>>><rsub|k>|\<beta\><rsub|0>+N<rsub|k>>>>|<row|<cell|\<beta\><rsub|k>>|<cell|=>|<cell|\<beta\><rsub|0>+N>>|<row|<cell|<with|font-series|bold|W><rsub|k><rsup|-1>>|<cell|=>|<cell|<with|font-series|bold|W><rsub|0><rsup|-1>+N<rsub|k>
    <around*|(|<with|font-series|bold|x><rsub|n>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><around*|(|<with|font-series|bold|x><rsub|n>-<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>|)><rsup|T>+<frac|\<beta\><rsub|0>N<rsub|k>|\<beta\><rsub|0>+N<rsub|k>>
    <around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><around*|(|<wide|<with|font-series|bold|x>|\<bar\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T>>>|<row|<cell|\<upsilon\><rsub|k>>|<cell|=>|<cell|v<rsub|0>+N<rsub|k>>>>>
  </eqnarray*>

  This results turn out to be very similar to the standard Gaussian-Wishart
  posterior without soft assignment. The approach here should be directly
  applicable to that, too.

  <section|Quantities required for computing <math|q<around*|(|Z|)>>>

  Standard results:\ 

  <\equation*>
    \<rho\><rsub|i k>=\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>
  </equation*>

  Please refer to Bishop's PRML from computing these expectations.

  <section|Quantities required for computing
  <math|><math|q<rsup|><around*|(|<with|font-series|bold|\<pi\>>|)>> and
  <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>>

  Please refer to Bishop's PRML from computing the associated expectations.
</body>

<\initial>
  <\collection>
    <associate|page-medium|papyrus>
    <associate|page-orientation|portrait>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|3|1>>
    <associate|auto-4|<tuple|4|1>>
    <associate|auto-5|<tuple|5|1>>
    <associate|auto-6|<tuple|5.1|1>>
    <associate|auto-7|<tuple|5.2|1>>
    <associate|auto-8|<tuple|6|1>>
    <associate|auto-9|<tuple|7|?>>
    <associate|logqparam-decompose|<tuple|1|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Complete-data
      log likelihood> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Prior
      for <with|mode|<quote|math>|\<pi\>>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Joint
      prior for <with|mode|<quote|math>|<with|font-series|<quote|bold>|\<mu\>><rsub|k>>
      and <with|mode|<quote|math>|<with|font-series|<quote|bold>|\<Sigma\>><rsub|k>>:>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Computing
      <with|mode|<quote|math>|q<around*|(|Z|)>>: variational E-step>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Computing
      <with|mode|<quote|math>|q<around*|(|\<pi\>,\<mu\><rsub|1:C>,\<Lambda\><rsub|1:C>|)>>:
      variational M-step> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <with|par-left|<quote|1tab>|5.1<space|2spc><with|mode|<quote|math>|q<around*|(|<with|font-series|<quote|bold>|\<pi\>>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|5.2<space|2spc><with|mode|<quote|math>|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Quantities
      required for computing <with|mode|<quote|math>|q<around*|(|Z|)>>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>Quantities
      required for computing <with|mode|<quote|math>|><with|mode|<quote|math>|q<rsup|><around*|(|<with|font-series|<quote|bold>|\<pi\>>|)>>
      and <with|mode|<quote|math>|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>