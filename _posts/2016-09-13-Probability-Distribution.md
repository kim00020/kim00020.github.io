---
layout: post 
title: "확률 분포"
author: "JKKim"
date: "September 13, 2016"



comments: true
share: true

---




# 예제 

컴퓨터용 팬(fan)을 만드는 회사에서 제품이 고장나면 무상으로 교체해주는 워런티 기간을 3년으로 잡았다고 합시다. 컴퓨터용 팬과 같은 기계적 제품의 생존 시간을 X 라고 할때 그 확률 값은 
$$ P( X < x ) = 1- \exp \left( - \frac{x}{\theta}\right),  \ \  x>0 $$
을 따르는 것이 알려져 있고 여기서 모수 $\theta=E(X)$값은 회사 제품별로 다른데 그 회사의 과거 자료를 분석한 결과 그 제품의 평균 수명은 5만 시간이라는 것이 밝혀졌다고 합니다. 그렇다면이 3년 워런티에 따라 발생하는 추가 비용은 제품값의 몇 \% 으로 계산해야 할까요? 



# 풀이  

3년은 최대 365 day $\times$ 24 hours $\times$ 3 years = 26,280 hours 로 계산될수 있으므로 3년 내에 고장이 날 확률은 $P(X<26,280)$을 계산하는 문제로 귀결됩니다. 따라서 위의 식에 $\theta=50,000$을 대입하면 
$$ P(X<26,280)= 1- \exp (-\frac{26280}{50000}) =0.4088$$
으로 계산되어 40.88\% 의 확률로 3년내에 고장이 발생할 것입니다. 따라서 원 제품값에 추가로 40.88\% 에 해당되는 비용을 워런티의 값으로 계산해 주어야 손해가 발생하지 않습니다. 즉, 제품 한대값이 원래 1만원이라면 워런티를 포함한 판매가는 1만 4천 88원보다 크게 잡아야 합니다. 만약 1년 워런티로 한다면 추가비용은 16\%로 줄어들 것입니다. 

# 토론

1. 우리가 어떤 확률 변수 $X$에 대해 특정 범위에서의  확률값을 계산하기 위해서는 확률 분포(probability distribution)를 알아야 합니다. 위의 예제에서는 평균이 $\theta$ 인 지수 분포(exponential distribution)를 사용했는데 생존시간에 대한 확률 분포로써 지수 분포가 아닌 다른 분포 (예: Weibull 분포)도 사용될수 있습니다. 


2. 확률 분포란 확률 변수의 실현치들의 분포가 어떤 규칙을 따른다고 보고 그 법칙을 수리적으로 기술한 것입니다. 현실 세계에서 이를 정확하게 파악하는 것은 매우 어려운 일이고 대신 이를 좀더 단순한 방식으로 근사적으로 표현하고자 하는 것이 통계적 모형의 핵심입니다. 

3. 즉 통계 모형이란 어떤 확률 변수에 대한 확률 분포를 좀더 이해하기 쉬운 형태로 근사적으로 표현하여 문제를 단순화하여 이해하고자 하는 것입니다. 통계모형은 확률분포를 범주화하여 몇개의 모수(parameter)로 표현해 줌으로써 확률 분포를 결정하는 문제를 모형 선택과 모수 추정의 문제로 단순화 시키게 됩니다. 모형 선택은 해당 분야의 특정 도메인 지식으로 알려져 있는 경우가 종종 있고 모수 추정은 수집된 데이터를 통해서 얻는 경우가 많이 있습니다. 


4. 위의 예에서 흥미로운 것은 평균 수명이 5만시간 이지만 그의 반 정도에 해당되는 2만 6천 시간이내에 고장날 확률이 40%나 된다는 것입니다. 평균 수명의 반도 안되서 고장나는 사건이 발생하는 빈도를 결코 무시할수 없다는 것입니다. 그런 면에서는 평균이 주는 착시 효과라고도 부를수 있을  것입니다. 이를 제대로 계산하지 못해서 워런티 값을 낮게 책정하면 많이 팔고도 손해를 볼수 있게 됩니다. 


