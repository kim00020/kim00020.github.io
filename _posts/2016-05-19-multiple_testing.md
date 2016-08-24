---
layout: post
title: "Basics of Model Selection: Likelihood Testing and Multple Testing"
output: html_document
comments: true
---
 
 I am going to talk about model selection from a statistical hypothesis perspective. A model can be defined by setting a likelihood. When two models are nested (one is a special case of the other one), the asymptotic distribution of -2 times log likelihood ratio converges to a chi-square distribution, so the statistical hypothesis can be constructed. Even though the likelihood ratio test is a simple and strong tool for model selection, it can be only applied to nested models, and it is not appropriate to use for the comparision of multiple models.
  
  Another popular testing-based approach for model selection is mutiple testing procedures. In linear models, by marginally considering each covariate, we can get $$p$$ number of t-statistics and its p-values, and we can select the variables whose p-value is small. Since Bonferroni's correction is too conservative, False Discovery Rate (FDR) control [(Benjamini and Hochberg, 1995)](http://www.math.tau.ac.il/~ybenja/MyPapers/benjamini_hochberg1995.pdf) has been popular, which controls the expected proportion of Type I errors among the rejected hypotheses.
  
   However, these mutiple testing approaches could have a critical problem when test statistics are correlated. Let me explain it with a simple example. Consider the model as 
   
   $$Y = X\beta_0 + \epsilon,$$   
   
   where $$\epsilon \sim N(0,0.1 I_n)$$ and the true regression coefficients $$\beta_0 = (1,1,-1,-1,-1,0,\dots,0)$$ with $$n=200$$ and $$p=20$$, So only first five variables are involved in the data-generating model. Also, suppose that the covariate is generated from the Gaussian distribution with a confound symmetry covariance with correlation 0.5; i.e., $$\Sigma_{i,j}=0.5$$, if $$i\neq j$$, and $$1$$, otherwise for $$1\leq i,j\leq p$$. Similar settings with this are commonly used in simulation studies of many variable selection papers such as  [nonlocal prior paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3867525/) , [the reciprocal lasso](http://www.tandfonline.com/doi/abs/10.1080/01621459.2014.984812), and [the BASAD](https://arxiv.org/pdf/1405.6545.pdf). The below is the R code to generate the data. 

{% highlight r %}
n=200; p =20
Sig = matrix(0.5,p,p)
diag(Sig)=1
ch = chol(Sig)
X = matrix(rnorm(n*p),n,p)%*%ch
beta0 = rep(0,p)
beta0[1:5] = c(1,1,-1,-1,-1)
y = X%*%beta0 + rnorm(n)*sqrt(0.1)
{% endhighlight %}

{% highlight r %}
as.vector(round(cor(X,y),3))
{% endhighlight %}



{% highlight text %}
##  [1]  0.018  0.069 -0.477 -0.532 -0.544 -0.231 -0.328 -0.386 -0.306 -0.117
## [11] -0.168 -0.226 -0.245 -0.313 -0.320 -0.304 -0.236 -0.238 -0.271 -0.210
{% endhighlight %}
  I calculated the marginal sample correlation between the response $$y$$ and each individual covariate $$X_j$$ for $$j=1,\dots,p$$. 

   Intuitively, the first five sample correlation coefficients shoud be more significant than the others, since the true model includes the first five, but the first two correlation coefficients are almost zero: 0.018 and 0.069. Since the t-statistic is a monotone transformation of the abolute value of the correlation coefficient, multiple testing procedures for the data set, including the FDR control, will never select the first two variables that are in the true model. As suggested by [Fan and Lv (2008)](http://orfe.princeton.edu/~jqfan/papers/06/SIS.pdf), *Sure Independence Screening* (SIS) has been highlightened in many areas, but the above result shows that SIS may fail under some  settings where the covariates are correlated. 

###### *What's going on here?* 

To see the details of the problem, let me write a covariance structure of the model. Roughly speaking, the correlation coefficients betweeen the response and the covariates can be identified by the following quantity:

$$\begin{eqnarray*} 
Cov({\bf X},y) &\approx&  Cov\left( {\bf X} ,  X_1 + X_2 - X_3 - X_4 - X_5 \right)\\
&=&  \begin{bmatrix}
1 & 0.5 & \dots & 0.5 \\
0.5 & 1 & \dots & 0.5 \\
0.5 & \dots& \dots &\dots \\
0.5& \dots  &  \dots & \dots\\
0.5 &\dots &\dots & 1\\
\end{bmatrix}
\begin{bmatrix}
1\\
1\\
-1\\
-1\\
-1\\
0\\
\dots\\
0
\end{bmatrix}
\\
&=&
\begin{bmatrix}
0\\
0\\
-1\\
-1\\
-1\\
-0.5\\
\dots\\
-0.5
\end{bmatrix}
,
\end{eqnarray*}
$$

where $${\bf X} = (X_1,\dots,X_p)^T$$.

   The first two covariance values are zero, and  even when the signal-to-ratio and the sample size are large enough, the resulting sample correlation between important variables and the response could be close to zero, and the sample correlation between the response and the noise variable that do not involve in the true model could be significantly far from zero. The above extreme example says that when the covariates are correlated the marginal correlations could force us  to select a wrong model.

   To overcome this issue, [Fan et al. (2012)](https://orfe.princeton.edu/~jqfan/papers/12/FDP-JASA.pdf) introduced a multiple testing procedure that takes account for the correlation structure of the test statistics, but its estimation of  the  covariance structure itself is tricky under high-dimensional settings. In biostatistics fields, [*Surrogate Variable Analysis*](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030161) (SVA) has been quite popular, which is a multiple-staged singular value decomposition for multiple testing problems. It really reduces the correlation between p-values, and the resulting distribution of the p-values becomes alomost like a uniform distribution. However, somtimes its performance is unstable, since the methodology is somewhat heuristic. Also, the inventor of the FDR control, Hochberg, published [a paper](https://projecteuclid.org/euclid.aos/1013699998) that asserts that  the FDR control is robust to some arbitrary correlation between test statistics exists. But their claim cannot be applied to the above example, and any multiple testing procedure based on marginal test statistics fail to select the true model in the example.
   
##### *Conclusion*
  
  Even though multiple testing procedures such as the FDR control are useful in many applications, when the test statistics are correlated, we really need to be careful to use the procedures. Sadly, a lot of people just use FDR or multiple testing for the high-dimensional model selection without considering the correlation structure. It's no good ! 
  
 
 
 

