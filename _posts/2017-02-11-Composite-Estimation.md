---
layout: post 
title: "회귀 분석의 응용"
author: JKKim
modified: 2017-02-15
tags: [통계학 개론]
comments: true
share: true
---

# 예제 


2014년에 콩고 공화국의 어느 지역에는 에볼라 바이러스 증상을 보이는 환자가 나타났습니다. 그래서 콩고 공화국에서는 급히 국제 보건 기구의 도움을 받아 그 주변 지역에서 백신 접종 캠페인을 실시 했습니다. 이 캠페인은 열흘 동안 열리는데 그 캠페인을 통해 바이러스 백신 접종률을  높이는게 아주 중요하다고 합니다. (만약 접종률이 낮으면 전염병이 다른 지역으로 퍼질 확률이 높아지게 된다고 합니다.) 따라서 백신 접종률을 높히기 위해서 다양한 노력을 하는 것도 중요하지만 실제 접종률을 정확히 파악하는 것 역시 중요합니다. 

실제 접종률을 파악하기 위해 그 지역의 마을을 다니며 서베이를 했는데 현지 사정상 제 4일차, 5일차, 6일차, 8일차에 서베이를 해서 마을 주민들의 백신 접종률을 아래와 같이 추정했다고 합니다. 



| 날(day) | 추정 접종율(%) | 표준 오차(%) |  접종자수 | 
|---------|----------------|--------------|-----------|
|  4      |       68.8     |     15.3     |   240,624 |
|  5      |       75.6     |     14.1     |   274,701 |
|  6      |       87.5     |     10.6     |   305,807 |
|  8      |       93.8     |      7.6     |   341,780 |
|---------|----------------|--------------|-----------|
| 10      |                |              |   366,783 |
|  |  | | | 




위의 테이블에서 표준오차(Standard error)는 서베이를 통해서 얻어지는 추정치의 신뢰도를 나타내는 지표입니다. 위 테이블의 접종자수는 서베이를 통해서 얻어지는 것이 아니라 캠페인 기간동안 접종을 받은 사람들 숫자를 집계한 것입니다. 이 접종자수에는 해당 마을 주민 뿐만 아니라 인근 마을 주민들도 접종을 했기에 이를 바탕으로 접종률 통계를 산출할수는 없다고 합니다. 위의 정보를 이용해서 결국 마지막날인 제 10일째에 실제 접종률이 얼마였는지를 추정하는 것이 자료분석의 목표입니다. 어떻게 해야 보다 정확한 추정을 할수 있을까요? 



# 풀이 

이러한 자료 구조는 실제 문제에서 많이 나타납니다. 우리가 원하는 변수는 실제 접종률인데 (이를 $Y$로 표현합시다) 이를 직접 관측하지는 못하고 $Y$의 추정량인 $\hat{Y}$만을 관측합니다. 다행히 확률 표본 추출로 얻어진 서베이에서는 추정량의 표준오차를 얻을수 있습니다. 위의 자료에서 접종자수는 $Y$와 관련있는 보조정보로써 $X$라고 표현할수 있습니다. 우리는 상식적으로 $Y$와 $X$에 선형관계가 있을 것이라고 모형을 세울수 있습니다. 이러한 자료에서 흔히 가정되는 모형은 ratio model 입니다. 
\begin{equation}
Y_i = X_i \gamma + e_i, \   \   \  e_i \sim N (0, X_i \sigma^2) .        \ \ \  \  \ (1) 
\end{equation}
이러한 모형에서 $\gamma$의 최적 추정량은 
$ \hat{\gamma} = (\sum_i Y_i )/( \sum_i  X_i )$ 
으로 계산됩니다. 이를 바탕으로 $Y_{10}$의 예측을 
$$\hat{Y}_{10} = X_{10} \hat{\gamma}$$ 
으로 계산할수 있습니다. 이러한 추정을 비추정(ratio estimation)이라고 부릅니다. 

그러나 문제는 우리가 $Y_i$를 관측하는 것이 아니라 $$\hat{Y}_{i}$$을 관측하므로 위의 비추정을 바로 적용할수 없다는 것입니다. 만약 $$u_i = \hat{Y}_i - Y_i$$이라고 표현하면 이 분포는 평균이 $0$ 이고 분산이  $v_i$인 정규분포로 가정할수 있으므로 $$\hat{Y}_i=Y_i + u_i$$의 주변분포는 
$$ \hat{Y}_i  \sim N( X_i \gamma , X_i \sigma^2 + v_i )
 $$으로 표현됩니다. 
우리 자료에서는 4개의 관측치가 있고 2개의 모수가 있으므로 ($v_i$는 알려져 있음)  모수 추정이 가능해 집니다.  이렇게 해서 모수가 추정되면 $Y_{10}$의 예측치는 $$\hat{Y}_{10}=X_{10}\hat{\gamma}$$으로 구해집니다. 위의 자료를 가지고 계산하면 $$\hat{Y}_{10}=1.024$$로 나오게 됩니다.


# 토론 

1. 
모수 추정을 위해서 최대우도 추정법을 사용하거나  Iterative reweighted least squares method 를 사용할수 있는데 후자는 $\gamma$와 $\sigma^2$을 번갈아가면서 추정하는 것입니다. 이를 위해서 $$W_i (\sigma^2)=1/( X_i \sigma^2 + v_i )$$으로 정의하면 $\gamma$의 추정은 
$$ \sum_i W_i ( \hat{\sigma}^2 )  ( \hat{Y}_i - X_i \gamma)^2 $$
을 최소화 하는 값인 
\begin{equation}
 \hat{\gamma}= \frac{\sum_i W_i ( \hat{\sigma}^2) X_i \hat{Y}_i  }{ \sum_i W_i ( \hat{\sigma}^2) X_i^2 } 
 \label{1} 
\end{equation}
으로 구하면 됩니다. $\sigma^2$의 추정을 위해서는 $Z_i (\gamma)= (\hat{Y}_i - X_i \gamma)^2$으로 정의하면  $Z_i(\gamma) \sim [  X_i \sigma^2 + v_i , 2(  X_i \sigma^2 + v_i )^2 ]  $임을 이용하여 
$$ \sum_i  \left\{    \frac{ Z_i (\hat{\gamma}) }{ X_i \sigma^2 + v_i  } - 1 
  \right\}^2  
$$ 
을 최소화하는 $\sigma^2$을 찾으면 됩니다. 이는 
\begin{equation}
 \hat{\sigma}^{2(t+1)}= \frac{\sum_i \{ W_i ( \hat{\sigma}^{2(t)}  ) \}^2  X_i ( \hat{Z}_i - v_i)  }{ \sum_i  \{ W_i ( \hat{\sigma}^{2(t)}  ) \}^2 X_i^2   } 
\label{2}
\end{equation}
으로 계산됩니다. 따라서 위의 두 식을 번갈아가면서 계산하면  그 수렴값이 최적 모수 추정값으로 사용될수 있습니다. 



2. 위의 방법대로 추정을 하면 ${Y}_{10}$의 예측값이 1보다 큰 값을 가지게 될수 있습니다. 이는 모형 (1)이 $$Y_i \in (0,1)$$의 range restriction 을 부여하지 않았기 때문에 발생하는 것으로 이해할수 있습니다. 이를 해결하기 위해서는 다음의 모형을 세울수 있을 것입니다. 
\\
\begin{equation}
 Y_i = p(\gamma_0+ X_i \gamma_1) + e_i, \   \   \  e_i \sim  (0, p_i (1-p_i) \sigma^2) . 
\end{equation}
이때 
$$ p(x)=  \exp(x)/\{1+ \exp(x)\} 
$$
으로 표현되고 $$p_i = p(\gamma_0 + X_i \gamma_1)$$입니다.  따라서 
$$ \hat{Y}_i \sim \left[  p(\gamma_0 + X_i \gamma_1) , v_i + p_i (1-p_i) \sigma^2 \right] 
$$
이 얻어지므로 토론 1에서 소개한 Iterative reweighted least squares method  방법을 사용하여 모수를 추정할수 있습니다. 이렇게 해서 모수 추정을 하게 되면 $$Y_{10}$$의  예측은 $$\hat{Y}_{10} = p(\hat{\gamma}_0 + X_{10} \hat{\gamma}_1)$$으로 계산됩니다. 
위의 자료를 적용하면 $$\hat{Y}_{10}$$=95.2%이 얻어집니다. 



3. 또한 이러한 모수 추정을 바탕으로 $$\hat{Y}_i$$가 관측된 날에 대하여 $Y_i$에 대한 보다 개선된 예측이 가능해질 것입니다. $Y_i$에 대한 추정량으로써 $$\hat{Y}_i$$와 $$\hat{p}_i = p(\hat{\gamma}_0 + X_{i} \hat{\gamma}_1)$$의 두가지를 생각할수 있는데 전자는 직접 조사에 의한 추정량으로 직접추정량(direct estimator)이라고 부르고 후자를 모형에 의하여 얻어진 추정량으로 모형기반 추정량이라고 부를수 있습니다. 이 두개를 최적 결합한 예측값은 
\begin{equation}
 \hat{Y}_i^*  = \hat{\alpha}_i \hat{Y}_i + (1- \hat{\alpha}_i)  \hat{p}_i
 \label{5}
 \end{equation}
으로 계산되고  이때 최적계수 $\hat{\alpha}_i$는 
\begin{equation}
 \hat{\alpha}_i = \frac{ \hat{p}_i (1-\hat{p}_i) \hat{\sigma}^2 }{v_i +   \hat{p}_i (1-\hat{p}_i) \hat{\sigma}^2} 
 \end{equation}
으로 계산됩니다. 


4. 위의 예제는 제가 했던 컨설팅에서 얻어진 실제 데이터를 사용한 것입니다. 위의 방법론을  바탕으로 각 날짜별로 백신 접종률의 추정값을 구하면 다음과 같습니다. 괄호 안의 값은 해당 추정량의 표준오차에 해당하는 값이 됩니다. 모든 경우에 대해서 최적추정의 표준 오차가 제일 작은 것을 확인할수 있을 것입니다. 이는 모형을 반영하여 일종의 shrinkage estimation 을 구현했기 때문입니다. 

| 날(day) |  직접 추정(%) |  모형기반추정(%) |  최적추정(%) | 
|---------|----------------|--------------|-----------|
|  4      |       68.8 (15.3) |     66.6 (15.0) |  68.0 (10.7)   |
|  5      |       75.6 (14.1) |     78.9 (12.9) | 77.2 (9.6)    |
|  6      |       87.5 (10.6) |     87.0 (10.7) |    87.2 (7.6)   |
|  8      |       93.8  (7.6) |      93.0 (8.1) |    93.1 (5.7) |


.
.





.
.







