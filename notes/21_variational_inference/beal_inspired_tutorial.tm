<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|Mean-field Variational Inference for
  Conjugate-Exponential Models: A Tutorial>|<doc-subtitle|Revisiting Chapter
  21: Variational Inference of <with|font-shape|italic|Machine Learning: A
  Probabilistic Perspective>>|<doc-author|<author-data|<author-name|Zhihan
  Yang>|<\author-affiliation>
    Department of Mathematics and Statistics
  </author-affiliation>|<\author-affiliation>
    Carleton College
  </author-affiliation>>>>

  <abstract-data|<\abstract>
    All models encounterd in Murphy begs a question: why?\ 

    First show that the model actually belongs to the family

    Using this fact, we pick the approximate posterior that allows for
    tractable VI

    Then we derive the coordinate descent VI algorithm
  </abstract>>

  <section|Conjugate-exponential models>

  <subsection|Exponential family>

  <section|Variational inference in general>

  <section|Examples>

  <subsection|Univariate Gaussian>

  \;

  <subsection|Linear regression>

  \;

  <subsection|Linear regression with ARD>

  \;

  <subsection|Mixture of Gaussians>

  Complete-data likelihood:

  <\equation*>
    <big|prod><rsub|i=1><rsup|n>p<around*|(|<with|font-series|bold|x><rsub|i>,<with|font-series|bold|z><rsub|i>\<mid\><with|font-series|bold|\<theta\>>|)>=
    <big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>> <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Sigma\>><rsub|k>|)><rsup|z<rsub|i
    k>>
  </equation*>

  where <math|<with|font-series|bold|\<theta\>>=<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1:k>,<with|font-series|bold|\<Sigma\>><rsub|1:k>|)>>.

  Prior for <math|<with|font-series|bold|\<pi\>>>:

  Dirichlet

  <\equation*>
    Dir<around*|(|<with|font-series|bold|\<pi\>>\<mid\><with|font-series|bold|\<alpha\>><rsub|0>|)>\<propto\><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>
  </equation*>

  Relevant part of the complete data log likelihood

  <\equation*>
    <big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>>
  </equation*>

  Proof that they are conjugate (this is not actually useful)

  <\equation*>
    <around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>|)><around*|(|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|z<rsub|i
    k>>|)>=<around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>-1>|)><around*|(|<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<big|sum><rsub|i=1><rsup|n>z<rsub|i
    k>>|)>=<big|prod><rsub|k=1><rsup|c>\<pi\><rsub|k><rsup|<around*|(|\<alpha\><rsub|0>|)><rsub|k>+<big|sum><rsub|i=1><rsup|n>z<rsub|i
    k>-1>
  </equation*>

  The joint prior for <math|<with|font-series|bold|\<mu\>><rsub|k>> and
  <math|<with|font-series|bold|\<Sigma\>><rsub|k>>:

  Gaussian wishart

  <\equation*>
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>\<mid\><with|font-series|bold|m><rsub|0>,\<beta\><rsub|0>,v<rsub|0>|)>=<with|font|cal|N><around*|(|<with|font-series|bold|\<mu\>><rsub|k>\<mid\><with|font-series|bold|m><rsub|0>,<around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><rsup|-1>|)>
    Wi<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>\<mid\><with|font-series|bold|L><rsub|0>,v<rsub|0>|)>
  </equation*>

  Relevant part of the complete data log likelihood

  <\equation*>
    <big|prod><rsub|i=1><rsup|n><with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Sigma\>><rsub|k>|)><rsup|z<rsub|i
    k>>=<big|prod><rsub|i<text| s.t. >z<rsub|i,k=1>><with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Sigma\>><rsub|k>|)>
  </equation*>

  Proof that they are conjugate to each other

  Omitted, see section 4.6.3.3

  Instantiate the CAVI algorithm

  <subsection|<math|q<around*|(|Z|)>>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<around*|(|<with|font-series|bold|Z>|)>>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>,<with|font-series|bold|Z>,<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>+log
    p<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>,<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>\<comma\><around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log \<pi\><rsub|k>|]>+\<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log <with|font|cal|N><around*|(|<with|font-series|bold|x><rsub|i>\<mid\><with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Sigma\>><rsub|k>|)><rsup|>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|<frac|1|2>
    ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> <around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> <around*|[|<frac|1|2> \<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k><around*|(|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<frac|1|2> \<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<around*|{|<with|font-series|bold|\<mu\>><rsub|k>|}>,<around*|{|<with|font-series|bold|\<Sigma\>><rsub|k>|}>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>|)>>>|<row|<cell|>|<cell|\<triangleq\>>|<cell|<big|sum><rsub|i=1><rsup|n><big|sum><rsub|k=1><rsup|c>z<rsub|i
    k> log \<rho\><rsub|i k>>>|<row|<cell|q<around*|(|<with|font-series|bold|Z>|)>>|<cell|\<propto\>>|<cell|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>\<rho\><rsub|i
    k><rsup|z<rsub|i k>>>>>>
  </eqnarray*>

  Requring that the distributuon is normalized on each row:

  <\equation*>
    <tabular|<tformat|<table|<row|<cell|q<around*|(|<with|font-series|bold|Z>|)>=>|<cell|<big|prod><rsub|i=1><rsup|n><big|prod><rsub|k=1><rsup|c>r<rsub|i
    k><rsup|z<rsub|i k>>>>>>>
  </equation*>

  where <math|r<rsub|i,k>> is the rho divided by row sum. But we need values
  of rho to compute this distribution.

  <subsection|<math|q<around*|(|\<pi\>,\<mu\><rsub|1:C>,\<Lambda\><rsub|1:C>|)>>>

  Starting from the same formula:

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|ln
    q<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>|)>>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|\<pi\>>,<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>,<with|font-series|bold|Z>,<with|font-series|bold|X>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>|)>+log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>,<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+log
    p<around*|(|<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>|)>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|X>\<mid\><with|font-series|bold|Z>,<with|font-series|bold|\<mu\>><rsub|1:C>,<with|font-series|bold|\<Lambda\>><rsub|1:C>|)>|]>+const>>|<row|<cell|>|<cell|=>|<cell|<around*|<left|(|3>|log
    p<around*|(|<with|font-series|bold|\<pi\>>|)>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|Z>|)>><around*|[|log
    p<around*|(|<with|font-series|bold|Z>\<mid\><with|font-series|bold|\<pi\>>|)>|]>|<right|)|3>>+<big|sum><rsub|k=1><rsup|C><around*|(|p<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\><rsub|>><rsub|k>|)>+<big|sum><rsub|i=1><rsup|N>r<rsub|i
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

  <subsection|<math|q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>>

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

  Now, we can put everything together:

  <\eqnarray*>
    <tformat|<table|<row|<cell|log q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>|<cell|=>|<cell|<frac|1|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)><rsup|T><around*|(|\<beta\><rsub|0><with|font-series|bold|\<Lambda\>><rsub|k>|)><around*|(|<with|font-series|bold|\<mu\>><rsub|k>-<with|font-series|bold|m><rsub|0>|)>>>|<row|<cell|>|<cell|>|<cell|+<frac|<around*|(|\<upsilon\><rsub|0>-C-1|)>|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>-<frac|1|2>
    Tr<around*|(|<with|font-series|bold|W><rsub|0><rsup|-1><with|font-series|bold|\<Lambda\>><rsub|k>|)>>>|<row|<cell|>|<cell|>|<cell|+<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k> <around*|(|<frac|1|2> log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>-<frac|1|2><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|)>+const.>>>>
  </eqnarray*>

  It may not be immediate clear what we should do here, so it's good to
  remind ourselves that by Theorem TODO the approximate posterior over
  <math|<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>>
  would also be Gaussian-Wishart.\ 

  Firstly, we see that

  <\equation*>
    <frac|<around*|(|\<upsilon\><rsub|0>-C-1|)>|2>
    log<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|><rsup|>+<frac|<around*|(|<big|sum><rsub|i=1><rsup|N>r<rsub|i
    k>|)>|2> log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>=<frac|<around*|(|\<upsilon\><rsub|0>+N<rsub|k>-C-1|)>|2>
    log <around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>,
  </equation*>

  which means that the Wishart posterior distribution on
  <math|<with|font-series|bold|\<Lambda\>><rsub|k>> have
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
  <math|<with|font-series|bold|\<mu\>><rsub|k>>. For simplicity, we will drop
  the shared <math|-1/2> term shortly.

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

  <subsection|Quantities required for computing <math|q<around*|(|Z|)>>>

  Standard results:\ 

  <\equation*>
    \<rho\><rsub|i r>=\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>+<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>-<frac|D|2>
    ln 2\<pi\>-<frac|1|2> \<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><around*|(|<with|font-series|bold|x><rsub|i>-<with|font-series|bold|\<mu\>><rsub|k>|)>|]>>|<cell|=>|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|x><rsub|i><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|x><rsub|i>-2<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|x><rsub|i>+<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|<with|font-series|bold|x><rsub|i><rsup|T>\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><with|font-series|bold|x><rsub|i>-2\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k>|]><with|font-series|bold|x><rsub|i>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|<with|font-series|bold|x><rsub|i><rsup|T><with|font-series|bold|W><rsub|k><with|font-series|bold|x><rsub|i>-2\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k>|]><with|font-series|bold|x><rsub|i>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|<with|font-series|bold|x><rsub|i><rsup|T><with|font-series|bold|W><rsub|k><with|font-series|bold|x><rsub|i>-2<with|font-series|bold|m><rsub|k>v<rsub|k><with|font-series|bold|W><rsub|k><with|font-series|bold|x><rsub|i>+\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<mu\>><rsub|k>,<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|<with|font-series|bold|\<mu\>><rsub|k><rsup|T><with|font-series|bold|\<Lambda\>><rsub|k><with|font-series|bold|\<mu\>><rsub|k>|]>>>|<row|<cell|>|<cell|=>|<cell|?>>|<row|<cell|>|<cell|>|<cell|>>|<row|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<Lambda\>><rsub|k>|)>><around*|[|ln<around*|\||<with|font-series|bold|\<Lambda\>><rsub|k>|\|>|]>>|<cell|=>|<cell|>>|<row|<cell|\<bbb-E\><rsub|q<around*|(|<with|font-series|bold|\<pi\>>|)>><around*|[|log
    \<pi\><rsub|k>|]>>|<cell|=>|<cell|>>>>
  </eqnarray*>

  <subsection|Quantities required for computing
  <math|><math|q<rsup|\<ast\>><around*|(|<with|font-series|bold|\<pi\>>|)>>
  and <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>>\ 

  \;

  <subsection|Algorithm>

  Initialize approximate posteriors <math|q<around*|(|<with|font-series|bold|\<pi\>>|)>>,
  <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>> directly from the
  priors.

  q(pi) to be dirichlet with a0=10e-3; <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>

  for normalized data, maybe we can start off with wishart of mean I = Wv; if
  v=2, then w=I/2

  m0 is set to 0, and we choose beta to be 1

  \;

  Compute the statistics required for computing
  <math|q<around*|(|<with|font-series|bold|Z>|)>>

  \;

  E_mu_k_Lambda_k_of_quadratic_form = D * (1/beta_k) + ups_k * (x_i - m_k).T
  @ W_k @ (x_i - m_k)

  E_Lambda_k_of_log_det_Lambda_k = D * log(2) + log det W_k

  for i in range(D):

  <space|2em>E_Lambda_k_of_log_det_Lambda_k += gamma((ups_k + 1 - i) / 2)

  for k in range(C):

  <space|2em>E_pi_k_of_log_pi_k = gamma(alpha_t[k]) - gamma(sum(alpha_t))

  \;

  Compute <math|q<around*|(|<with|font-series|bold|Z>|)>>

  try to think of a matrix formula for this?

  r_nk_unnorm =\ 

  \;

  Compute the statistics required for computing
  <math|q<around*|(|<with|font-series|bold|\<pi\>>|)>>,
  <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>

  easy, done already

  \;

  Compute <math|q<around*|(|<with|font-series|bold|\<pi\>>|)>>,
  <math|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>

  10.58

  10.60 to 10.63

  \;

  Until convergence

  <subsection|Experiment>

  \;

  <subsection|Logistic regression*>

  \;

  <subsection|Logistic regression with ARD*>

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|page-medium|papyrus>
    <associate|page-orientation|landscape>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|3.6|3>>
    <associate|auto-11|<tuple|3.7|3>>
    <associate|auto-12|<tuple|3.8|3>>
    <associate|auto-13|<tuple|3.9|3>>
    <associate|auto-14|<tuple|3.10|3>>
    <associate|auto-15|<tuple|3.11|?>>
    <associate|auto-16|<tuple|3.12|?>>
    <associate|auto-17|<tuple|3.13|?>>
    <associate|auto-18|<tuple|3.14|?>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|2|1>>
    <associate|auto-4|<tuple|3|1>>
    <associate|auto-5|<tuple|3.1|1>>
    <associate|auto-6|<tuple|3.2|1>>
    <associate|auto-7|<tuple|3.3|1>>
    <associate|auto-8|<tuple|3.4|1>>
    <associate|auto-9|<tuple|3.5|2>>
    <associate|logqparam-decompose|<tuple|1|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Conjugate-exponential
      models> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Exponential family
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Variational
      inference in general> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Examples>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Univariate Gaussian
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Linear regression
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Linear regression with ARD
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1tab>|3.4<space|2spc>Mixture of Gaussians
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1tab>|3.5<space|2spc><with|mode|<quote|math>|q<around*|(|Z|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.6<space|2spc><with|mode|<quote|math>|q<around*|(|\<pi\>,\<mu\><rsub|1:C>,\<Lambda\><rsub|1:C>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|3.7<space|2spc><with|mode|<quote|math>|q<rsup|\<ast\>><around*|(|<with|font-series|<quote|bold>|\<pi\>>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|1tab>|3.8<space|2spc><with|mode|<quote|math>|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|1tab>|3.9<space|2spc>Quantities required for
      computing <with|mode|<quote|math>|q<around*|(|Z|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|1tab>|3.10<space|2spc>Quantities required for
      computing <with|mode|<quote|math>|><with|mode|<quote|math>|q<rsup|\<ast\>><around*|(|<with|font-series|<quote|bold>|\<pi\>>|)>>
      and <with|mode|<quote|math>|q<around*|(|\<mu\><rsub|k>,\<Lambda\><rsub|k>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|3.11<space|2spc>Algorithm
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|3.12<space|2spc>Experiment
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1tab>|3.13<space|2spc>Logistic regression*
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|1tab>|3.14<space|2spc>Logistic regression with
      ARD* <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>
    </associate>
  </collection>
</auxiliary>