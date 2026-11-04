

For Brown's method, the Fisher statistic was defined as

\[
X_F=-2\sum_{i=1}^{3}\ln(p_i),
\]

where \(p_i\) denotes the P value from the \(i\)-th component test. Under independence, \(X_F\) follows a \(\chi^2\) distribution with six degrees of freedom. To account for dependence, Brown's method approximates its null distribution using a scaled chi-square distribution, \(X_F\sim c\chi^2_{\nu}\), where

\[
c=\frac{\operatorname{Var}(X_F)}
{2\operatorname{E}(X_F)},\qquad
\nu=\frac{2[\operatorname{E}(X_F)]^2}
{\operatorname{Var}(X_F)}.
\]

The variance was estimated empirically from the covariance between the Fisher-transformed component P values across ambient barcodes. Under uniform marginal null P values, \(\operatorname{E}(X_F)=6\).

For ACAT, the three component P values were combined using the Cauchy transformation:

\[
T_{\mathrm{ACAT}}=
\sum_{i=1}^{3}w_i
\tan\left[\left(0.5-p_i\right)\pi\right],
\]

where \(w_i=1/3\) for equal weighting. The combined P value was calculated as

\[
p_{\mathrm{ACAT}}=
\frac{1}{2}-\frac{1}{\pi}
\arctan(T_{\mathrm{ACAT}}).
\]

ACAT provides an asymptotically valid combination procedure under arbitrary dependence among component tests.
