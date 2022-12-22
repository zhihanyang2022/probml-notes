<TeXmacs|2.1.1>

<style|generic>

<\body>
  <doc-data|<doc-title|PCA and SVD>>

  <section|What does PCA try to achieve?>

  Maximum variance perspective

  Project remainig data to obtain maximum variance

  \;

  Reconstruction perspective

  <\theorem>
    Suppose we want to find an orthonormal set of <math|L> basis vectors
    <math|<with|font-series|bold|w><rsub|j>\<in\>\<bbb-R\><rsup|D>> (they
    form columns of <math|<with|font-series|bold|W>\<in\>\<bbb-R\><rsup|D\<times\>L>>)
    such that the following reconstruction objective is minimized:

    <\equation*>
      J<around*|(|<with|font-series|bold|W>,<with|font-series|bold|Z>|)>=<around*|\<\|\|\>|<with|font-series|bold|X>-<with|font-series|bold|W><with|font-series|bold|Z>|\<\|\|\>><rsub|F><rsup|2>.
    </equation*>

    The solution is to set <math|<with|font-series|bold|<wide|W|^>>> to be
    the <math|L> eigenvectors of the empirical covariance matrix
    <math|<wide|<with|font-series|bold|\<Sigma\>>|^>=<with|font-series|bold|X><rsup|T><with|font-series|bold|X>>.
    </theorem>

  One can show that the errors can be decomposed due to triangles, then it's
  again the maximum variance perspective kicks in.

  Again, in either cases, the principal components are by definition
  orthonormal.y

  <section|How to compute the principal components?>

  Eigendecomposition of covariance matrix

  <section|Computing principal compoents using SVD>

  \;

  <section|SVD computes more than just the principal components>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?|../.TeXmacs/texts/scratch/no_name_6.tm>>
    <associate|auto-2|<tuple|2|?|../.TeXmacs/texts/scratch/no_name_6.tm>>
    <associate|auto-3|<tuple|3|?|../.TeXmacs/texts/scratch/no_name_6.tm>>
    <associate|auto-4|<tuple|4|?|../.TeXmacs/texts/scratch/no_name_6.tm>>
  </collection>
</references>