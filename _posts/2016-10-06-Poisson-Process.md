---
layout: post 
title: "포아송 프로세스"
author: "JKKim"
date: "October 6, 2016"


comments: true
share: true


---




## 예제 

어느 시골 마을  입구 큰길 앞에서는 택시가 1시간당 10대가 지나간다고 합니다. 그 10대중 2대만 빈차라고 하고 합승을 해주지 않는다고 합니다. 당신이 택시 타는 줄의 두번째에 있다면 당신이 택시를 잡기 위해 1시간 이상 기다리게 될 확률은 얼마일까요? 

## 풀이  

택시가 오는 시각이 서로 독립적인 경우 한 시간동안 지나가는 택시의 갯수는 포아송 분포를 따릅니다. 위의 예에서는 한 시간 동안 오는 빈 택시의 갯수는 평균이 2인 포아송 분포, 즉 Poisson(2),가 되는 것입니다. 그런데 이를 t 시간으로 확장하면 t 시간동안 오는 빈택시의 빈도수 $N_t$ 는 t 에 따라 증가하는 확률 변수로써 그 분포가 Poisson (2t) 를 따르고 이를 Poisson process 라고 합니다. 


포아송 프로세스는 사건 발생 횟수에 대한 확률 변수라고 한다면 거꾸로 사건이 발생하는 시간에 대한 확률 변수를 생각할수 있습니다. 즉, $T_1$을 첫번째 택시가 오는데 걸리는  시간, $T_2$를 두번째 택시가 오는데 걸리는 시간, 이런 식으로 정의하면 이 문제는 
$P( T_2 >1)$을 계산하고자 하는 문제입니다. 그런데 포아송 프로세스의 정의상 $\{ T_2 >t \}$ 이라는 사건은 $\{N_t \le 1\}$이라는 사건과 동등하므로 (왜냐하면 두번째 택시가 t 시간 이후에 온다는 이야기는 t 시간 까지 빈택시 횟수가 1개 이하라는 이야기이므로) 
우리가 원하는 확률은 
$$P( T_2 >1) = P ( N_1 \le 1) = P(N_1 =0) + P(N_1=1) = e^{-2} (1+2)=3e^{-2}=0.41$$
으로 계산됩니다. 참고로 $X$ 가 Poisson(m) 를 따른다는 것은 
$P(X=k)= \frac{ m^k}{k!} e^{-m}$
임을 의미합니다. 







## 토론 

1. 포아송 프로세스 $N_t$는 처음에는 0부터 시작해서 $t$가 늘어남에 따라 하나씩 증가합니다. (단위 시간당 발생 비율은 $\lambda$ 값을 가집니다.) 이 포아송 프로세스가 처음으로 1 의 값을 갖게 되는 사건이 발생하는 시간 (제 1회 사건 발생시간)의 확률분포가 바로 지수분포가 됩니다. 왜냐하면 $\{N_t=0\}$ 이라는 사건은 시간 $t$까지는 관심 사건이 아직 일어나지 않았다는 이야기이고 이는 처음 발생은 t 이후에 일어난 다는 의미이므로 
$\{T_1>t\}$과 동등합니다. 따라서 
$$ P( T_1>t)=P(N_t=0) = \frac{(\lambda t)^0}{0 !}\exp (-\lambda t) $$
으로 되므로 $T_1$ 의 누적 분포 함수 (cumulative distribution function)은 $1-\exp (-\lambda t)$이 되어서 평균이 $1/\lambda$ 인 지수 분포를 따르는 것을 확인할수 있습니다. 


2. 만약 첫번째 사건이 발생하는 시간에 관심이 있는 것이 아니라 첫 $r(>1)$ 번째 사건이 발생하는 시간에 관심이 있는 경우에는 서로 독립인 지수분포가 $r$개 합해진 형태인 감마 분포를 따르게 됩니다. 앞의 정의를 이용하면 $T_r$은 첫 $r$-번째 사건이 발생하기까지 걸리는 시간이 되는데  $\{T_r > t\}$라는 것은 그 소요 시간이 $t$보다 크다는 것이고 따라서 시간 $t$에는 아직 $r$-번째 사건이 발생하지 않았다는 의미이므로 $\{N_t \le r-1\}$과 동일해 집니다.  따라서 
$$ P (Y_r \le t ) = 1- P ( Y_r > t) = 1 - P ( N_t  \le r-1 ) = 1 - \sum_{k=0}^{r-1} \frac{ (\lambda t)^k}{ k!} e^{-\lambda t } = \sum_{k=r}^{\infty} \frac{ (\lambda t)^k}{ k!} e^{-\lambda t } $$  
이 얻어지고 이를 $t$에 대해 미분하면 
$$ f_Y( t ) = \frac{\lambda^k }{(k-1)! }  t^{k-1}e^{- \lambda t }  $$    
이 되어서 감마 분포의 확률 밀도 함수와 동일해 집니다. 위의 택시 예에서 두번째 줄에 있으므로 택시를 잡는 사건이 두번째 일어나는 시간에 관심이 있는 것입니다. 그 시간이 감마 분포를 따르게 되는 것입니다. 

3. 지수 분포의 특징은 memoryless property 라는 것인데 이를 수식으로 표현하면 다음과 같습니다. 
$$ P(T_1 > t+s \mid T_1 >t) = P(T_1>s). $$
즉, 그 전에 얼마 동안 기다렸다는 것이 아무런 상관이 없다는 것입니다. 이는 결국 포아송 프로세스의 정의에서 연유하는 것입니다. 어떤 사건이 발생하는 사건(arrival time)이 서로 독립적으로 분포하기 때문에 이런 성질을 갖는 것입니다. 


.
.










