---
title: NOIP2021 方差
date: 2021-11-21 13:06:35
tags: [差分,dp]
categories: [NOIP2021]
mathjax: true
---


$$
\begin{align*}
n\overline{a}&=\sum_{i=1}^na_i\\
\\
n^2D
&=n\sum_{i=1}^n(a_i-\overline{a})^2\\
&=n\sum_{i=1}^n a_i^2-2n\overline{a}\sum_{i=1}^n a_i+n^2\overline{a}^2\\
&=n\sum_{i=1}^n a_i^2-2(\sum_{i=1}^n a_i)^2+(\sum_{i=1}^n a_i)^2\\
&=n\sum_{i=1}^n a_i^2-(\sum_{i=1}^n a_i)^2\\
\end{align*}
$$


观察题目中的式子，$a'_i\leftarrow a_{i-1}+a_{i+1}-a_i$，根据 CF1110E 的套路，可以差分，令 $d_i=a_{i+1}-a_i(1\leq i\lt n)$，一次操作 $(2\leq i\lt n)$ 即：

$$
d'_{i-1}=a'_{i}-a'_{i-1}=(a_{i-1}+a_{i+1}-a_i)-a_{i-1}=a_{i+1}-a_i=d_i\\
d'_i=a'_{i+1}-a'_{i}=a_{i+1}-(a_{i-1}+a_{i+1}-a_i)=a_i-a_{i-1}=d_{i-1}
$$

相当于交换 $d_{i-1},d_i(2\leq i\lt n)$，故 $d$ 可以通过若干次操作，变为任意 $d$ 的排列。



由于 $a_1$ 不变，那么 $a$ 和 $d$ 是一一对应的。现在要求一个 $d$ 的排列使得 $n^2D$ 最小，继续推式子：
$$
\begin{align*}
n^2D
&=n\sum_{i=1}^n a_i^2-(\sum_{i=1}^n a_i)^2\\
&=n\sum_{i=1}^n a_i^2-\sum_{i=1}^n\sum_{j=1}^na_ia_j\\
&=\frac{1}{2}(n\sum_{i=1}^n a_i^2-2\sum_{i=1}^n\sum_{j=1}^na_ia_j+n\sum_{j=1}^n a_j^2)\\
&=\frac{1}{2}(\sum_{i=1}^n\sum_{j=1}^n(a_i-a_j)^2)\\
&=\sum_{i=1}^{n-1}\sum_{j=i}^{n-1}(a_{j+1}-a_i)^2\\
&=\sum_{i=1}^{n-1}\sum_{j=i}^{n-1}(d_i+d_{i+1}+\cdots+d_j)^2\\
\end{align*}
$$
$n^2D$ 取最小值时，$d$ 一定是先递减后递增的。



考虑在分界位置（即递减到递增）从小到大往两边加数，由于
$$
\begin{align*}
n^2D
&=n\sum_{i=1}^n a_i^2-(\sum_{i=1}^n a_i)^2\\
\end{align*}
$$
维护 $dp[k][s]$ 表示当前已经加入了 $k$ 个数，现在的 $a$ 之和为 $s$ 的最小 $\displaystyle \sum_{i=1}^n a_i^2$。

初始 $dp[1][s]=0$。

转移时考虑加到左边还是右边：

- 左边：原来是 $a_1,a_2,\cdots,a_k$，现在变为 $d,a_1+d,a_2+d,\cdots,a_k+d$，新增的贡献为

$$
\begin{align*}
\Delta
&=d^2+\sum_{i=1}^k (d+a_i)^2-\sum_{i=1}^k a_i^2\\
&=(k+1)d^2+2d\sum_{i=1}^ka_i\\
&=(k+1)d^2+2ds\\
\end{align*}
$$

- 右边：原来是 $a_1,a_2,\cdots,a_k$，现在变为 $a_1,a_2,\cdots,a_k,a_k+d$，新增的贡献为
  $$
  \begin{align*}
  \Delta
  &=(a_k+d)^2\\
  \end{align*}
  $$
  其实可以发现 $a_k$ 是固定的，为之前所有的 $d$ 之和，不需要再记录。



答案为 $\displaystyle \min_s\{n\times dp[n][s]-s^2\}$。



分析一下时间复杂度，第一维是 $O(n)$ 的，第二维是 $O(nV)$ 的，转移 $O(1)$。

但是可以发现 $d$ 为 $0$ 的转移可以忽略，第一维是 $O(\min(n,V))$ 的，总复杂度为 $O(nV^2)$，可以通过。



```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll INF=1ll<<60;
int n,a[11000],d[11000]; ll dp[510000];
ll sqr(ll x){ return x*x;}
int main(){
	scanf("%d",&n);
	for (int i=1;i<=n;i++) scanf("%d",&a[i]);
	for (int i=1;i<n;i++) d[i]=a[i+1]-a[i];
	sort(d+1,d+n);
	
	for (int i=0;i<=500000;i++) dp[i]=INF;
	dp[0]=0; int lim=0,sum=0;
	for (int i=1;i<n;i++){
		if (!d[i]) continue;
		for (int s=lim;s>=0;s--){
			if (dp[s]==INF) continue;
			dp[s+sum+d[i]]=min(dp[s+sum+d[i]],dp[s]+sqr(sum+d[i]));
			dp[s+i*d[i]]=min(dp[s+i*d[i]],dp[s]+2*s*d[i]+i*sqr(d[i]));
			dp[s]=INF;
		}
		lim+=i*d[i]; sum+=d[i];
	}
	ll ans=INF;
	for (int i=0;i<=lim;i++)
		if (dp[i]!=INF) ans=min(ans,n*dp[i]-1ll*i*i);
	printf("%lld\n",ans);
	return 0;
}

```

