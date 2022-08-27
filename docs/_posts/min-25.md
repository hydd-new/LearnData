---
title: min_25筛
date: 2021-10-11 20:00:53
tags: min_25筛
categories: 算法
---
### 是什么

Min_25 筛是一种求积性函数前缀和的算法。

具体来说，对于一个积性函数 $f$，求：

$$
ans=\sum_{i=1}^n f(i)
$$
需要快速计算素数的幂的值。计算素数个数函数 $\pi(n)$ 的时间复杂度为 $O(\dfrac{n^{\frac{3}{4}}}{\log n})$（它只用到第一步），计算其它函数时间复杂度为 $O(n^{1-\epsilon})$（一般能跑过的范围内，也可以按照 $O(\dfrac{n^{\frac{3}{4}}}{\log n})$ 计算）。



### 怎么做

#### 第一步

思考下面这个问题：
$$
G(n)=\sum_{i=1}^n [i\in \mathbb P]f(i)
$$
其中 $\mathbb P$ 是素数集。



考虑我们是怎么筛素数的，有一种方法是埃氏筛，就是枚举 $p$，然后将最小素因子是 $p$ 的合数筛去。

当然，我们不能直接像埃氏筛一样枚举最小素因子是 $p$ 的合数，这样还不如直接暴力。



先找到完全积性函数 $f'$ 满足对于任意 $x\in \mathbb P$，$f'(x)=f(x)$，显然满足条件的 $f'$ 只有一个（完全积性函数由素数位置的点值唯一确定）。

设 $g(n,j)$ 表示 $2,3,\cdots,n$ 这些数经过 $j$ 轮埃氏筛（即筛了前 $j$ 个素数）后，剩下的 $f'$ 的值和，$1$ 没有贡献，则 $G(n)=g(n,|P|)$。

即 <u>$i$ 是质数</u>或者 $i$ 的最小素因子大于 $P_j$ 的 $f'(i)$ 之和，具体地：
$$
g(n,j)=\sum_{i=2}^n [i\in \mathbb P\space \operatorname{or}\space \operatorname{minp}(i)\gt P_j]f'(i)
$$
其中 $\operatorname{minp}(i)$ 表示 $i$ 的最小素因子，$P_j$ 表示第 $j$ 小的素数。



考虑使用 $dp$ 来计算 $g(n,j)$，那么要在 $g(n,j-1)$ 中剩下的数中把最小素因子是 $P_j$ 的合数筛去。

$P_j^2$ 为满足最小素因子是 $P_j$ 的最小合数，分情况讨论：

- $P_j^2\gt n$，不存在最小素因子是 $P_j$ 的合数，$g(n,j)=g(n,j-1)$；
- $P_j^2\leq n$，$f'$ 满足完全积性，且会被筛掉的数为 $2,3,\cdots,\lfloor\dfrac{n}{P_j}\rfloor$ 中任何满足 $\operatorname{minp}(x)>P_{j-1}$ 的数 $x\times P_j$，$g(n,j)=g(n,j-1)-f'(P_j)[g(\lfloor\dfrac{n}{P_j}\rfloor,j-1)-g(P_j-1,j-1)]$（因为 $g(\lfloor\dfrac{n}{P_j}\rfloor,j-1)$ 中还计算了 $2,3,\cdots,p_j-1$ 种的素数的值所以扣去）。

 $g(P_j-1,j-1)$ 可以直接求出，不需要单独计算（$P_j^2\leq n$，就是一个素数位置的值的前缀和）。

由于 $g$ 的第一维是每次是除一个数下取整，而 $\lfloor\dfrac{\lfloor\frac{n}{a}\rfloor}{b}\rfloor=\lfloor\dfrac{n}{ab}\rfloor$，所以第一维可以只记录形如 $\lfloor\dfrac{n}{x}\rfloor$ 的值。



#### 第二步

求出了 $G$ 之后如何计算答案。

类似地，设 $S(n,j)$ 表示 $2,3,\cdots,n$ 这些数中，$i$ 的最小素因子大于 $P_j$ 的 $f(i)$ 之和，将 $1$ 的贡献单独算，则 $ans=S(n,1)$。

具体地：
$$
S(n,j)=\sum_{i=2}^n[\operatorname{minp}(i)\gt P_j]f(i)
$$


考虑仍然使用 $dp$ 计算，
$$
S(n,j)=G(n)-\sum_{i=1}^{j-1}f(P_i)+\sum_{k\ge j}\sum_{e\geq 1,P_k^e\leq n}f(P_k^e)(S(\lfloor\frac{n}{P_k^e}\rfloor,k+1)+[e\gt 1])
$$


### 应用

#### 计算素数个数函数 $\pi(n)$ 

只需要用到第一步，$f$ 函数满足素数位置的值为 $1$。

直接按照第一步的步骤做即可。

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
ll n,l,g[2100000];
int tot,p[2100000]; bool vis[2100000];
int cnt,pos1[2100000],pos2[2100000]; ll num[2100000];
void Sieve(int n){
    vis[1]=1;
    for (int i=2;i<=n;i++){
        if (!vis[i]) p[++tot]=i;
        for (int j=1;j<=tot&&i*p[j]<=n;j++){
            vis[i*p[j]]=true;
            if (!(i%p[j])) break;
        }
    }
}
int calc(ll x){ return x<=l?pos1[x]:pos2[n/x];}
void upd(ll x){
    num[++cnt]=x; g[cnt]=x-1;
    if (x<=l) pos1[x]=cnt;
    else pos2[n/x]=cnt;
}
int main(){
    scanf("%lld",&n); l=1;
    while (l*l<=n) l++;
    l--; Sieve(l);
    for (ll i=1,j;i<=n;i=j+1){
        j=n/(n/i); upd(n/i);
    }
    for (int j=1;j<=tot;j++)
        for (int i=1;1ll*p[j]*p[j]<=num[i];i++)
            g[i]-=(g[calc(num[i]/p[j])]-(j-1));
    printf("%lld\n",g[1]);
    return 0;
}
```

