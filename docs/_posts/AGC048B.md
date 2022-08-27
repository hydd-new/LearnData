---
title: AGC048B
date: 2021-10-01 16:15:23
tags: [贪心]
categories: [AtCoder]
mathjax: true
---

- 设一个由小括号和中括号组成的串 $S$ 的权值为 $\sum_{S_i\in\{\texttt{'(',')'}\}} A_i+\sum_{S_i\in\{\texttt{'[',']'}\}} B_i$。
- 求所有的合法的由小括号和中括号组成的长度为 $N$ 的串中最大的权值是多少。
- $2\leq N\leq 10^5,2\mid N,1\leq A_i,B_i\leq 10^9$。
<!-- more -->

### sol

![sol](https://cdn.jsdelivr.net/gh/hydd-new/image-hosting@master/agc048b.png)

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
int n,a[110000],b[2][51000];
int main(){
	scanf("%d",&n);
	ll ans=0;
	for (int i=0;i<n;i++){
		scanf("%d",&a[i]);
		ans+=a[i];
	}
	int x;
	for (int i=0;i<n;i++){
		scanf("%d",&x); x-=a[i];
		b[i&1][i>>1]=x;
	}
	sort(b[0],b[0]+(n>>1),greater<int>());
	sort(b[1],b[1]+(n>>1),greater<int>());
	for (int i=0;i<(n>>1);i++)//not n
		if (b[0][i]+b[1][i]>0) ans+=b[0][i]+b[1][i];
		else break;
	printf("%lld\n",ans);
	return 0;
}
```
#### Note
  - 将代码里注释那一行的 `n>>1` 写成 `n` 可能导致RE
  - `ans` 需要开 `long long`