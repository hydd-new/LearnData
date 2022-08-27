---
title: AGC048C
date: 2021-10-01 16:46:04
tags: [差分]
categories: [AtCoder]
mathjax: true
---

- 给你 $N,L$ 和序列 $0=A_0<A_1<A_2<\cdots<A_n\leq A_{n+1}=L+1,0=B_0<B_1<B_2<\cdots<B_n\leq B_{n+1}=L+1$。
- 每次可以选择一个数 $x(1\leq x\leq N)$，使得 $A_x=A_{x-1}+1$ 或 $A_x=A_{x+1}-1$。
- 问把 $A$ 变成 $B$ 的最少次数，无解输出 $-1$。
- $1\leq N\leq 10^5,N\leq L\leq 10^9$。
<!-- more -->

### sol

![sol](https://cdn.jsdelivr.net/gh/hydd-new/image-hosting@master/agc048c.png)

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
int n,m,a[110000],b[110000],x[110000],y[110000];
void read(int *a,int *b,int &cnt){
	a[0]=0; a[n+1]=m+1;
	for (int i=1;i<=n;i++) scanf("%d",&a[i]);
	for (int i=n+1;i>=1;i--) a[i]-=a[i-1]+1;
	for (int i=1;i<=n+1;i++)
		if (a[i]) b[++cnt]=i;
}
int main(){
	scanf("%d%d",&n,&m);
	int cnt1=0,cnt2=0;
	read(a,x,cnt1); read(b,y,cnt2);
	ll ans=0; int j=1;
	for (int i=1;i<=cnt2;i++){
		int sum=0,l=x[j],r;
		while (j<=cnt1&&sum<b[y[i]]) sum+=a[x[j++]];
		r=x[j-1]; ans+=r-l;
		if (sum!=b[y[i]]){ puts("-1"); return 0;}
		if (y[i]<l) ans+=l-y[i];
		else if (r<y[i]) ans+=y[i]-r;
	}
	if (j<=cnt1) puts("-1");
	else printf("%lld\n",ans);
	return 0;
}
```
#### Note
  - `ans` 需要开 `long long`