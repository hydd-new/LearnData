---
title: '模板'
date: 2021-10-01 14:00:00
tags: []
published: true
hideInList: false
feature: 
isTop: true
---

- 普通幂转第二类斯特林数：$i^k=\sum\limits_{j=1}^{i}S_{k,j}\times C_{i,j}\times j! $。
- 第二类斯特林数求自然数幂和：$\sum\limits^n_{i = 0}i ^ k =  \sum \limits _ {j = 1} ^ n S_{k, j}\times j! \times C_{n + 1, j + 1}$。
- $\rm min-max$容斥：$\max(S)=\sum_{T \subseteq S}(-1)^{|T|+1}\min(T)$，$\min(S)=\sum_{T \subseteq S}(-1)^{|T|+1}\max(T)$
- $\gcd-\mathrm{lcm}$容斥：$\mathrm{lcm}(T)=\prod_{S\in T}\gcd(S)^{(-1)^{|S|+1}}$
<!-- more -->
- 多项式开根：若$B'^2(x)\equiv A(x) (mod\ x^{\frac{n}{2}})$，$B^2(x)\equiv A(x) (mod\ x^{n})$，则有以下等式成立：$B(x)\equiv \dfrac{A(x)+B'^2(x)}{2B'(x)} (mod\ x^n)$。

# 0. IO/数学

## 0.1 读入优化

```cpp
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
char readop(){
    char ch=Getchar();
    while (ch<'A'||ch>'Z') ch=Getchar();
    return ch;
}
```



## 0.2 输出优化

```cpp
char pbuf[100000],*pp=pbuf;
void pc(const char c) {
    if (pp-pbuf==100000) fwrite(pbuf,1,100000,stdout),pp=pbuf;
    *pp++=c;
}
void write(int x) {
    static int sta[35];
    if (x<0){
        pc('-');
        x=-x;
    }
    int top=0;
    do{
        sta[top++]=x%10;
    x/=10;
    } while(x);
    while (top) pc(sta[--top]+'0');
}
void myfflush(){
    fwrite(pbuf,1,pp-pbuf,stdout);
}
```


## 0.3 取模优化
```cpp
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
```


# 1. 数论

## 1.1 预处理阶乘和阶乘逆元

```cpp
fac[0]=1; for (int i=1;i<=10000000;i++) fac[i]=fac[i-1]*i%Mod;
inv[1]=1; for (int i=2;i<=10000000;i++) inv[i]=(Mod-Mod/i)*inv[Mod%i]%Mod;
inv[0]=1; for (int i=1;i<=10000000;i++) inv[i]=inv[i-1]*inv[i]%Mod;
ll C(int n,int m){
    if (n<m) return 0;
    return 1ll*fac[n]*inv[m]%Mod*inv[n-m]%Mod;
}
```

## 1.2 快速幂

```cpp
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
```



## 1.3 gcd&lcm

```cpp
int gcd(int a,int b){
    if (!b) return a;
    return gcd(b,a%b);
}
int lcm(int a,int b){ return 1ll*a/gcd(a,b)*b;}
```


## 1.4 CRT

```cpp
#include<cstdio>
#include<algorithm>
using namespace std;
typedef long long ll;
int n;
ll a[11],b[11];
ll mult(ll a,ll b,ll mod){
    if (a<b) swap(a,b);
    ll ans=0;
    while (b){
        if (b&1) ans=(ans+a)%mod;
        a=(a<<1)%mod; b>>=1;
    }
    return ans;
}
void exgcd(ll a,ll b,ll &x,ll &y){
    if (!b){
        x=1; y=0;
        return;
    }
    exgcd(b,a%b,x,y);
    ll tmp=x; x=y; y=tmp-a/b*y;
}
ll CRT(int n){
    ll x,y;
    ll ans=0,lcm=1;
    for (int i=1;i<=n;i++) lcm*=b[i];
    for (int i=1;i<=n;i++){
        ll tmp=lcm/b[i];
        exgcd(tmp,b[i],x,y);
        x=(x%b[i]+b[i])%b[i];
        ans=(ans+mult(mult(tmp,x,lcm),a[i],lcm))%lcm;
    }
    return (ans+lcm)%lcm;
}

int main(){
    scanf("%d",&n);
    for (int i=1;i<=n;i++) scanf("%lld",&a[i]);
    for (int i=1;i<=n;i++){
        scanf("%lld",&b[i]);
        a[i]=(a[i]%b[i]+b[i])%b[i];
    }
    printf("%lld\n",CRT(n));
    return 0;
}
```

## 1.5 exgcd&excrt

- 考虑合并两个方程 $x\equiv a_1(\operatorname{mod} m_1),x\equiv a_2(\operatorname{mod} m_2)$。
- 那么，有 $x=a_1+k_1\times m_1,x=a_2+k_2\times m_2$。
- $k_2\times m_2-k_1\times m_1=a_1-a_2$。
- 令 $c=a_1-a_2$，$g=gcd(m_2,m_1)$。
- 那么，若 $g\nmid c$，则方程无解。否则 $g\mid c,\frac{c}{g}\in \mathcal{Z}$。
- 我们可以用 $exgcd$ 求出 $k_2'\times m_2+k_1'\times m_1=g$ 中满足条件的 $k_1'$。
- 将两边同时乘 $\frac{c}{g}$，即 $k_1',k_2'$ 乘上 $\frac cg$，右边变成 $c$ 。
- 我们得到了 $k_2'\times m_2+k_1'\times m_1=c$ 的解 $k_1'$。可以将 $k_1' \mod m_2$。
- 但是注意，我们的原方程是 $k_2\times m_2-k_1\times m_1=a_1-a_2$，故 $k_1=-k_1'$。
- 也就是 $x=a_1-k_1'\times m_1$。我们就求出了一个解 $x_0$。那么通解是 $x\equiv x_0(\operatorname{mod} \operatorname{lcm}(m_1,m_2))$。

```cpp
typedef long long ll;
ll exgcd(ll a,ll b,ll &x,ll &y){
    if (!b){
        x=1; y=0;
        return a;
    }
    ll g=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return g;
}
ll a[N],m[N];
ll excrt(){
    ll M=m[1],A=a[1],x,y,d;
    for (int i=2;i<=n;i++){
        d=exgcd(M,m[i],x,y);
        if ((A-a[i])%d) return -1;
        x=(A-a[i])/d*x%m[i];
        A-=M*x;
        M=M/d*m[i];
        A%=M;
    }
    return (A%M+M)%M;
}
```



## 1.6 Pollard Rho

- 输出 $n$ 的最大质因子（$1\leqslant n\leqslant 10^{18}$）

![tmp1.png](/wp-content/uploads/2021/04/tmp1.png)

```cpp
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<algorithm>
#define ctz __builtin_ctzll
using namespace std;
typedef long long ll;
ll ans;
/*inline ll mul(ll x,ll y,ll Mod){
    ll ret=x*y-((ll)((long double)x/Mod*y+0.5))*Mod;
    return (ret%Mod+Mod)%Mod;
}*/
inline ll mul(ll x,ll y,ll p){
    return (__int128)x*y%p;
}
inline ll qpow(ll x,ll a,ll Mod){
    ll res=1;
    while (a){
        if (a&1) res=mul(res,x,Mod);
        x=mul(x,x,Mod); a>>=1;
    }
    return res;
}
inline ll rdm(){ return 1ll*rand()<<31|rand();}
ll gcd(ll x,ll y){
    if (!y) return x;
    return gcd(y,x%y);
}
bool Miller_Rabin(ll n){
    if (n==2||n==3) return true;
    if (!(n&1)||(n==1)||!(n%3)) return false;
    ll p=n-1,m=0;
    while (!(p&1)) p>>=1,++m;
    int Case=8;
    while (Case--){
        ll lst=qpow(rdm()%(n-1)+1,p,n),now=lst;
        for (int i=1;i<=m;i++){
            now=mul(now,now,n);
            if (now==1&&lst!=1&&lst!=n-1) return false;
            lst=now;
        }
        if (now!=1) return false;
    }
    return true;
}
ll Pollard_Rho(ll n,int c){
    ll i=1,k=2,x=rand()%(n-1)+1,y=x,sum=1;
    while (true){
        i++; x=(mul(x,x,n)+c)%n;
        sum=mul(sum,(y-x+n)%n,n);
        if (x==y||!sum) return n;
        if (i==k||i%127==0){
            int d=gcd(sum,n);
            if (d!=1) return d; 
            if (i==k){ y=x; k<<=1;}
        }
    }
}
void work(ll n){
    if (n<=ans) return;
    if (Miller_Rabin(n)){
        ans=max(ans,n);
        return;
    }
    ll tmp=n;
    while (tmp==n) tmp=Pollard_Rho(n,rand()%n);
    while (n%tmp==0) n/=tmp;
    work(tmp); work(n);
}
int main(){
    srand(time(NULL));
    int T; ll x;
    scanf("%d",&T);
    while (T--){
        scanf("%lld",&x);
        if (Miller_Rabin(x)){ puts("Prime"); continue;}
        ans=0; work(x);
        printf("%lld\n",ans);
    }
    return 0;
}

//贴个超级短的：
LL pollardRho(LL n, int a){
    LL x=2,y=2,d=1;
    while(d==1){
        x=(x*x+a)%n;
        y=(y*y+a)%n;y=(y*y+a)%n;
        d=gcd(abs(x-y),n);
    }
    if(d==n) return pollardRho(n,a+1);
    return d;
}
```



## 1.7 欧拉筛（线性筛）最小质因子

```cpp
//注意1不为质数 
for (int i=2;i<=n;i++){
    if (!lst[i]){
        lst[i]=i;
        prime[++cnt]=i;
    }
    for (int j=1;j<=cnt&&i*prime[j]<=n;j++){
        lst[i*prime[j]]=prime[j];
        if (i%prime[j]==0) break;
    }
}
```

## 1.8 BSGS

```cpp
int BSGS(ll a,ll b){
    int blo=(int)(sqrt(p)+1);
    ll base=b;
    for(int i=0;i<blo;i++){
        myhash[base]=i;
        base=base*a%Mod;
    }
    ll tmp=1;
    base=qpow(a,blo);
    for(int i=1;i<=blo+1;i++){
        tmp=tmp*base%Mod;
        int t=tmp;
        if (myhash.count(t)) return i*blo-myhash[t];
    }
    return -1;
}
```

## 1.9 欧拉筛（线性筛）mu/phi

```cpp
int cnt,mu[110000],p[110000],phi[110000];
bool vis[110000];
void init(int n){
    vis[1]=1;//注意，1不是质数
    phi[1]=0; mu[1]=1;
    for (int i=2;i<=n;i++){
        if (!vis[i]){
            p[++cnt]=i;
            phi[i]=i-1; mu[i]=-1;
        }
        for (int j=1;j<=cnt&&i*p[j]<=n;j++){
            vis[i*p[j]]=true;
            if (i%p[j]==0){
                phi[i*p[j]]=phi[i]*p[j];
                mu[i*p[j]]=0;
                break;
            }
            phi[i*p[j]]=phi[i]*(p[j]-1);
            mu[i*p[j]]=-mu[i];
        }
    }
}
```

## 1.10 求原根

![tmp2.png](/wp-content/uploads/2021/04/tmp2.png)

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
int p,phi;
int gcd(int a,int b){
	if (!b) return a;
	return gcd(b,a%b);
}
ll qpow(ll x,ll a,ll Mod){
	ll res=1;
	while (a){
		if (a&1) res=res*x%Mod;
		x=x*x%Mod; a>>=1;
	}
	return res;
}
int cnt,fac[110000];
void getfac(int x){
	for (int i=2;i*i<=x;i++)
		if (x%i==0){
			fac[++cnt]=i;
			while (x%i==0)x/=i;
		}
	if (x!=1) fac[++cnt]=x;
}
int primitive_root(){
	getfac(phi);
	for (int i=1;i<p;i++){
		if (gcd(i,p)!=1) continue; 
		int flag=true;
		for (int j=1;j<=cnt;j++)
			if (qpow(i,phi/fac[j],p)==1){
				flag=false;
				break;
			}
		if (flag) return i;
	}
	return -1;
}
int getphi(int x){
	int ans=x;
	for (int i=2;i*i<=x;i++)
		if (x%i==0){
			ans=1ll*ans*(i-1)/i;
			while (x%i==0) x/=i;
		}
	if (x>1) ans=1ll*ans*(x-1)/x;
	return ans;
}

int main(){
	scanf("%d",&p);
	phi=getphi(p);
	int g=primitive_root();
	printf("%d\n",g);
	return 0;
}
```



# 2.数据结构

## 2.1 树状数组

### 单点加，区间求和

```cpp
void add(int *tree,int x,int y){
    for (;x<=N;x+=x&-x) tree[x]+=y;
}
int getsum(int *tree,int x){
    int sum=0;
    for (;x;x-=x&(-x)) sum+=tree[x];
    return sum;
}
```

### 区间加，区间求和

```cpp
/*
    b[i]=a[i]-a[i-1]
    a[1] + a[2] + ... + a[n]
    = (b[1]) + (b[1]+b[2]) + ... + (b[1]+b[2]+...+b[n])
    = n*b[1] + (n-1)*b[2] +... +b[n]
    = n * (b[1]+b[2]+...+b[n]) - (0*b[1]+1*b[2]+...+(n-1)*b[n])
    sum1[i] = \sum b[i],sum2[i] = \sum b[i]*(i-1)
*/
int n,a[N];
int sum1[N],sum2[N];
void change(int x,int y){
    for (int i=x;i<=n;i+=(i&-i)){
        sum1[i]+=y;
        sum2[i]+=y*(x-1);
    }
}
void range_change(int l,int r,int y){
    change(l,y);
    change(r+1,-y);
}
int query(int x){
    int res=0;
    for (int i=x;i>=1;i-=(i&-i)) res+=x*sum1[i]-sum2[i];
    return res;
}
int range_query(int l,int r){ return query(r)-query(l-1);}
 
int main(){
    for (int i=1;i<=n;i++) change(i,a[i]-a[i-1]);
    //[l,r] +x        range_change(l,r,x); 
    //sum(l,r)        range_query(l,r)
    return 0;
}
```



## 2.2 线段树（ 区间加，区间求和）

```cpp
#include<cstdio>
#define ls ((now)<<1)//左儿子 
#define rs ((now)<<1|1)//右儿子 
using namespace std;
const int MAXN=110000;
typedef long long ll;
int n,m;
ll a[MAXN];//n->区间长度，m->询问个数，a->区间初始值 
struct node{
    ll val,add;//线段树的值，加法标记 
} tree[MAXN<<2];
void pushup(int now){//上传标记 
    tree[now].val=tree[ls].val+tree[rs].val; 
}
void pushdown(int now,int l,int r){//下传标记 
    int mid=(l+r)>>1;
    //更新实际值 
    tree[ls].val+=tree[now].add*(mid-l+1);
    tree[rs].val+=tree[now].add*(r-mid);
    //下传加法标记 
    tree[ls].add+=tree[now].add;
    tree[rs].add+=tree[now].add;
    //清空当前标记 
    tree[now].add=0;
}
void buildtree(int now,int l,int r){//建树 
    //零标记 
    tree[now].add=0;
    if (l==r){
        tree[now].val=a[l];//初始化 
        return;
    }
    int mid=(l+r)>>1;
    buildtree(ls,l,mid); buildtree(rs,mid+1,r);//建子树 
    pushup(now);//上传更新 
}
void change(int now,int l,int r,int x,int y,ll k){
    if (l>r||x>y) return;
    if (l==x&&r==y){
        tree[now].add+=k;//更新加法标记 
        tree[now].val+=k*(r-l+1);//更新值
        return;
    }
    int mid=(l+r)>>1; pushdown(now,l,r);//先下传 
    if (y<=mid) change(ls,l,mid,x,y,k);//更新左子树 
    else if (x>mid) change(rs,mid+1,r,x,y,k);//更新右子树 
    else change(ls,l,mid,x,mid,k),change(rs,mid+1,r,mid+1,y,k);//更新左右子树 
    pushup(now);//上传标记 
}
ll query(int now,int l,int r,int x,int y){
    if (l>r||x>y) return 0;
    if (l==x&&r==y) return tree[now].val;//返回实际值 
    int mid=(l+r)>>1; pushdown(now,l,r);//先下传 
    if (y<=mid) return query(ls,l,mid,x,y);//询问左子树 
    else if (x>mid) return query(rs,mid+1,r,x,y);//询问右子树 
    else return query(ls,l,mid,x,mid)+query(rs,mid+1,r,mid+1,y);//将左右子树的询问答案相加 
    pushup(now);//上传标记 
}
int main(){
    scanf("%d%d",&n,&m);
    for (int i=1;i<=n;i++) scanf("%lld",&a[i]);//读入初始值 
    int op,x,y; ll k;
    buildtree(1,1,n);//建树 
    while (m--){
        scanf("%d%d%d",&op,&x,&y);
        if (op==1){//1代表加 
            scanf("%lld",&k);
            change(1,1,n,x,y,k);
        } else printf("%lld\n",query(1,1,n,x,y));//其它数代表询问 
    }
    return 0;
}
```

## 2.3 splay

```cpp
#include<cstdio>
#include<algorithm>
using namespace std;
const int INF=1e9;
int n,m,op,x,a[110000];
int root,sz,tag[110000],ch[110000][2],f[110000],cnt[110000],key[110000],size[110000];
inline void clear(int x){
    ch[x][0]=ch[x][1]=f[x]=cnt[x]=key[x]=size[x]=0;
}
inline int get(int x){
    return ch[f[x]][1]==x;
}
inline void update(int x){
    if (x){
        size[x]=cnt[x];
        if (ch[x][0]) size[x]+=size[ch[x][0]];
        if (ch[x][1]) size[x]+=size[ch[x][1]];
    }
}
inline void pushdown(int x){
    if (x&&tag[x]){
        tag[ch[x][0]]^=1;
        tag[ch[x][1]]^=1;
        swap(ch[x][0],ch[x][1]);
        tag[x]=0;
    }
}

inline void rotate(int x){
    pushdown(f[x]); pushdown(x);
    int fa=f[x],gfa=f[fa],dir=get(x);
    ch[fa][dir]=ch[x][dir^1]; f[ch[fa][dir]]=fa;
    f[fa]=x; ch[x][dir^1]=fa;
    f[x]=gfa;
    if (gfa) ch[gfa][ch[gfa][1]==fa]=x;
    update(fa); update(x);
}
void splay(int x,int goal=0){
    for (int fa;(fa=f[x])!=goal;rotate(x))
        if (f[fa]!=goal)
            rotate((get(x)==get(fa))?fa:x);
    if (!goal) root=x;
}
inline int build(int l,int r,int fa){
    if (l>r) return 0;
    int mid=(l+r)>>1;
    int now=++sz;
    key[now]=a[mid]; f[now]=fa;
    cnt[now]=1; tag[now]=0;
    ch[now][0]=build(l,mid-1,now);
    ch[now][1]=build(mid+1,r,now);
    update(now);
    return now;
}
inline void insert(int v){
    if (!root){
        root=++sz;
        ch[sz][0]=ch[sz][1]=f[sz]=0;
        key[sz]=v; cnt[sz]=size[sz]=1;
        return;
    }
    int now=root,fa=0;
    while (true){
        if (key[now]==v){
            cnt[now]++;
            update(now); update(fa);
            splay(now);
            break;
        }
        fa=now;
        now=ch[now][key[now]<v];
        if (!now){
            sz++;
            ch[sz][0]=ch[sz][1]=0; key[sz]=v; 
            size[sz]=1; cnt[sz]=1;
            f[sz]=fa; ch[fa][key[fa]<v]=sz;
            update(fa);
            splay(sz);
            break;
        }
    }
}
inline int Rank(int v){
    int ans=0,now=root;
    while (true){
        if (v<key[now]) now=ch[now][0];
        else{
            ans+=(ch[now][0]?size[ch[now][0]]:0);
            if (v==key[now]){
                splay(now);
                return ans+1;
            }
            ans+=cnt[now];
            now=ch[now][1];
        }
    }
}
inline int find(int x){
    int now=root;
    while (true){
        pushdown(now);
        if (ch[now][0]&&x<=size[ch[now][0]]) now=ch[now][0];
        else{
            int temp=(ch[now][0]?size[ch[now][0]]:0)+cnt[now];
            if (x<=temp) return now;//key[now];
            x-=temp; now=ch[now][1];
        }
    }
}
inline int pre(){
    int now=ch[root][0];
    while (ch[now][1]) now=ch[now][1];
    return now;
}
 
inline int next(){
    int now=ch[root][1];
    while (ch[now][0]) now=ch[now][0];
    return now;
}
inline void del(int x){
    Rank(x);
    if (cnt[root]>1){
        cnt[root]--;
        update(root);
        return;
    }
    if (!ch[root][0]&&!ch[root][1]){
        clear(root);
        root=0;
        return;
    }
    if (!ch[root][0]){
        int oldroot=root;
        root=ch[root][1];
        f[root]=0;
        clear(oldroot);
        return;
    }
    else if (!ch[root][1]){
        int oldroot=root;
        root=ch[root][0];
        f[root]=0;
        clear(oldroot);
        return;
    }
    int leftbig=pre(),oldroot=root;
    splay(leftbig);
    f[ch[oldroot][1]]=root; ch[root][1]=ch[oldroot][1];
    clear(oldroot);
    update(root);
    return;
}
inline void print(int now){
    pushdown(now);
    if (ch[now][0]) print(ch[now][0]);
    if (key[now]!=-INF&&key[now]!=INF)
      printf("%d ",key[now]);
    if (ch[now][1]) print(ch[now][1]);
}
void reverse(int x,int y){//翻转[x,y] 
    int a=find(x);
    int b=find(y+2);
    splay(a); splay(b,a);
    tag[ch[ch[root][1]][0]]^=1;
}
int main(){
    scanf("%d",&n);
    insert(-INF); insert(INF);
    for (int i=1;i<=n;i++){
        scanf("%d%d",&op,&x);
        switch (op){
            case 1:
                insert(x);
                break;
            
            case 2:
                del(x);
                break;
            
            case 3:
                printf("%d\n",Rank(x)-1);
                break;
            
            case 4:
                printf("%d\n",key[find(x+1)]);
                break;
            
            case 5:
                insert(x);
                printf("%d\n",key[pre()]);
                del(x);
                break;
            
            case 6:
                insert(x);
                printf("%d\n",key[next()]);
                del(x);
                break;
            
        }
    }
    return 0;
}
```



## 2.4 左偏树
```cpp
struct pq{
    int v[110000];
    int l[110000],r[110000],d[110000];
    int merge(int x,int y){
        if (!x||!y) return x|y;
        if (v[x]<v[y]) swap(x,y);
        r[x]=merge(r[x],y);
        if (d[r[x]]>d[l[x]]) swap(l[x],r[x]);
        d[x]=d[r[x]]+1;
        return x;
    }
    void pop(int &x){ x=merge(l[x],r[x]);}
    int top(int x){ return v[x];}
} que;
```



## 2.5 LCT（Link-Cut Tree）

- 这部分代码敲了好几遍了，易错点
  1. 把fa[x]和f[x]混起来写。。。
  2. access中加注释的地方写成nroot(x)。
  3. link的时候写成f[y]=x。

```cpp
#include<cstdio>
#include<algorithm> 
#define ls c[x][0]
#define rs c[x][1]
using namespace std;
int n,m,type,x,y,f[310000],c[310000][2],v[310000],s[310000],st[310000];
bool r[310000];
inline bool nroot(int x){
    return x==c[f[x]][0]||x==c[f[x]][1];
}
inline void pushup(int x){
    s[x]=s[ls]^s[rs]^v[x];
}
inline void pushr(int x){
    swap(ls,rs);
    r[x]^=1;
}
inline void pushdown(int x){
    if (r[x]){
        if (ls) pushr(ls);
        if (rs) pushr(rs);
        r[x]=0;
    }
}
inline int dir(int x){ return c[f[x]][1]==x;}
void rotate(int x){
    int y=f[x],z=f[y],k=dir(x),w=c[x][!k];
    if (nroot(y)) c[z][dir(y)]=x;
    c[x][!k]=y; c[y][k]=w;
    if (w) f[w]=y; 
    f[y]=x; f[x]=z;
    pushup(y);
}
void splay(int x){
    int y=x,z=0;
    st[++z]=y;
    while (nroot(y)) st[++z]=y=f[y];
    while (z) pushdown(st[z--]);
    while (nroot(x)){
        y=f[x];
        if (nroot(y)) rotate(dir(x)!=dir(y)?x:y);
        rotate(x);
    }
    pushup(x);
}
void access(int x){
    for (int y=0;x;y=x,x=f[y]){//注意此处，x是打通到整棵树的根的路径，所以是;x;而不是;nroot(x);。
        splay(x);
        rs=y; pushup(x);
    } 
}
inline void makeroot(int x){ access(x); splay(x); pushr(x);}
int findroot(int x){
    access(x); splay(x);
    while (ls){
        pushdown(x); x=ls;
    }
    splay(x);
    return x;
}
inline void split(int x,int y){ makeroot(x); access(y); splay(y);}
void link(int x,int y){
    makeroot(x);
    if (findroot(y)!=x) f[x]=y;
}
void cut(int x,int y){
    makeroot(x);
    if (findroot(y)==x&&f[y]==x&&!c[y][0]){
        f[y]=0; c[x][1]=0;
        pushup(x);
    }
}
int main(){
    scanf("%d%d",&n,&m);
    for (int i=1;i<=n;i++) scanf("%d",&v[i]);
    while (m--){
        scanf("%d%d%d",&type,&x,&y);
        switch (type){
            case 0: split(x,y); printf("%d\n",s[y]); break;
            case 1: link(x,y); break;
            case 2: cut(x,y); break;
            case 3: splay(x); v[x]=y; break;
        }
    }
    return 0;
}
```

- 新写的板子

```cpp
#include<bits/stdc++.h>
int n,q,u[110000],v[110000];
template <typename T>
class LCT{
	struct node{
		node *ch[2],*fa,*pathfa;
		T val,sum,max;
		bool rev;
		node(T v=0){ rev=false; val=sum=max=v; ch[0]=ch[1]=fa=pathfa=NULL;}
		bool dir(){ return this==fa->ch[1];}
		void pushup(){
			sum=val;
			if (ch[0]!=NULL) sum+=ch[0]->sum;
			if (ch[1]!=NULL) sum+=ch[1]->sum;
			
			max=val;
			if (ch[0]!=NULL) max=std::max(max,ch[0]->max);
			if (ch[1]!=NULL) max=std::max(max,ch[1]->max);
		}
		void pushdown(){
			if (rev){
				std::swap(ch[0],ch[1]);
				if (ch[0]!=NULL) ch[0]->rev^=1;
				if (ch[1]!=NULL) ch[1]->rev^=1;
				rev=false;
			}
		}
		void rotate(){
			node *y=fa,*z=fa->fa;
			if (z!=NULL) z->pushdown();
			y->pushdown(); pushdown();
			std::swap(pathfa,y->pathfa);//!!!!!!
			
			int k=dir();
			if (z!=NULL) z->ch[y->dir()]=this;
			fa=z;
			
			y->ch[k]=ch[k^1];
			if (ch[k^1]) ch[k^1]->fa=y;
			
			ch[k^1]=y; y->fa=this;
			y->pushup(); pushup();
		}
		void splay(){
			while (fa!=NULL){
				if (fa->fa!=NULL){
					fa->fa->pushdown(); fa->pushdown();
					if (dir()==fa->dir()) fa->rotate();
					else rotate();
				}
				rotate();
			}
		}
		void expose(){//使当前点是重链上最下面一个点
			splay();
			pushdown();
			if (ch[1]){
				ch[1]->fa=NULL;
				ch[1]->pathfa=this;
				ch[1]=NULL;
				pushup();
			}
		}
		bool splice(){//将当前这条路径和上面一条路径拼起来
			splay();
			if (pathfa==NULL) return false;
			pathfa->expose();
			pathfa->ch[1]=this; fa=pathfa;
			pathfa=NULL; fa->pushup();
			return true;
		}
		void access(){//将这个点到根的路径存到一个splay中
			expose();
			while (splice());
		}
		void evert(){//把一个点作为整棵树的根
			access();
			splay(); rev^=1;
		}
		node* findroot(){
			access(); splay();
			node *x=this;
			x->pushdown();
			while (x->ch[0]) x=x->ch[0],x->pushdown();
			x->splay(); return x;
		}
		T querymax(){
			access(); splay();
			return max;
		}
		T querysum(){
			access(); splay();
			return sum;
		}
	} *p[110000];
	public:
	void evert(int u,T k){
		p[u]=new node(k);
	}
	bool link(int u,int v){
		p[v]->evert();
		if (p[u]->findroot()==p[v]) return false;
		p[v]->pathfa=p[u];
		return true;
	}
	bool cut(int u,int v){
		p[u]->evert();
		p[v]->access(); p[v]->splay();
		// p[v].pushdown();????
		if (p[v]->ch[0]!=p[u]) return false;
		p[v]->ch[0]->fa=NULL; p[v]->ch[0]=NULL;
		p[v]->pushup(); return true;
	}
	T querymax(int u,int v){
		p[u]->evert();
		return p[v]->querymax();
	}
	T querysum(int u,int v){
		p[u]->evert();
		return p[v]->querysum();
	}
	void change(int u,T k){
		p[u]->splay();
		p[u]->val=k;
		p[u]->pushup();
	}
};
LCT<int> tree;
int main(){
	scanf("%d",&n); int x;
	for (int i=1;i<n;i++) scanf("%d%d",&u[i],&v[i]);
	for (int i=1;i<=n;i++){
		scanf("%d",&x);
		tree.evert(i,x);
	}
	for (int i=1;i<n;i++) tree.link(u[i],v[i]);
	scanf("%d",&q);
	char op[11]; int u,v;
	while (q--){
		scanf("%s%d%d",op,&u,&v);
		if (op[0]=='C') tree.change(u,v);
		else if (op[1]=='S') printf("%d\n",tree.querysum(u,v));
		else if (op[1]=='M') printf("%d\n",tree.querymax(u,v));
	}
	return 0;
}
```

## 2.6 虚树

- 现在改为luogu2495的某篇题解

```cpp
#include<bits/stdc++.h>
using namespace std;
#define il inline
#define re register
#define debug printf("Now is Line : %d\n",__LINE__)
#define file(a) freopen(#a".in","r",stdin);freopen(#a".out","w",stdout)
#define int long long
#define inf 123456789000000000
#define mod 1000000007
il int read()
{
    re int x = 0, f = 1; re char c = getchar();
    while(c < '0' || c > '9') { if(c == '-') f = -1; c = getchar();}
    while(c >= '0' && c <= '9') x = x * 10 + c - 48, c = getchar();
    return x * f;
}
#define maxn 250005
struct edge
{
	int v, w, next;
}e[maxn << 1];
int n, m, head[maxn], cnt, is[maxn], mi[maxn], dfn[maxn], col, t;
int size[maxn], fa[maxn], top[maxn], son[maxn], dep[maxn], s[maxn];
vector<int>v[maxn];
il void add(int u, int v, int w)
{
	e[++ cnt] = (edge){v, w, head[u]};
	head[u] = cnt;
}
il bool cmp(int a, int b){return dfn[a] < dfn[b];}
il void dfs1(int u, int fr)
{
	size[u] = 1, fa[u] = fr, dep[u] = dep[fr] + 1;
	for(re int i = head[u]; i; i = e[i].next)
	{
		int v = e[i].v;
		if(v == fr) continue;
		mi[v] = min(mi[u], e[i].w);
		dfs1(v, u), size[u] += size[v];
		if(size[son[u]] < size[v]) son[u] = v;
	}
}
il void dfs2(int u, int fr)
{
	top[u] = fr, dfn[u] = ++ col;
	if(!son[u]) return;
	dfs2(son[u], fr);
	for(re int i = head[u]; i; i = e[i].next)
	{
		int v = e[i].v;
		if(v != fa[u] && v != son[u]) dfs2(v, v);
	}
}
il int lca(int a, int b)
{
	while(top[a] != top[b]) dep[top[a]] > dep[top[b]] ? a = fa[top[a]] : b = fa[top[b]];
	return dep[a] < dep[b] ? a : b;
}
il void push(int x)
{
	if(t == 1) {s[++ t] = x;return;}
	int l = lca(x, s[t]);
	if(l == s[t]) return;
	while(t > 1 && dfn[s[t - 1]] >= dfn[l]) v[s[t - 1]].push_back(s[t]), --t;
	if(s[t] != l) v[l].push_back(s[t]), s[t] = l;
	s[++ t] = x;
}
il int dp(int u)
{
	if(v[u].size() == 0) return mi[u];
	int temp = 0;
	for(re int i = 0; i < v[u].size(); ++ i) temp += dp(v[u][i]);
	v[u].clear();
	return min(mi[u], temp);
}
signed main()
{
	file(a);
	n = read();
	for(re int i = 1; i < n; ++ i)
	{
		int u = read(), v = read(), w = read();
		add(u, v, w), add(v, u, w);
	}
	mi[1] = inf, dfs1(1, 0), dfs2(1, 1);
	int T = read();
	while(T --)
	{
		m = read();
		for(re int i = 1; i <= m; ++ i) is[i] = read();
		sort(is + 1, is + m + 1, cmp);
		s[t = 1] = 1;
		for(re int i = 1; i <= m; ++ i) push(is[i]);
		while(t > 0) v[s[t - 1]].push_back(s[t]), --t;
		printf("%lld\n", dp(1));
	}
	return 0;
}
```



# 3.网络流&图论

## 3.1 最大流
### dinic（建议）

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int INF=0x3f3f3f3f;
int n,m,s,t,vis[11000],dep[11000]; ll maxflow;
int edgenum=1,F[210000],V[210000],Next[210000],Head[11000],cur[11000];
void addedge(int u,int v,int f){
    V[++edgenum]=v; F[edgenum]=f;
    Next[edgenum]=Head[u];
    Head[u]=edgenum;
}
void Add(int u,int v,int f){
    addedge(u,v,f);
    addedge(v,u,0);
}
bool bfs(){
    for (int i=1;i<=n;i++) vis[i]=false,dep[i]=INF;
    dep[s]=0;
    vis[s]=true;
    queue<int> que;
    que.push(s);
    while (!que.empty()){
        int u=que.front(); que.pop();
        vis[u]=false;
        for (int e=Head[u];e;e=Next[e]){
            int v=V[e];
            if (dep[v]>dep[u]+1&&F[e]){
                dep[v]=dep[u]+1;
                if (!vis[v]){
                    que.push(v);
                    vis[v]=true;
                }
            }
        }
    }
    return dep[t]!=INF;
}
int dfs(int u,int flow){
    if (u==t){ 
        maxflow+=flow;
        return flow;
    }
    int used=0;
    for (int &e=cur[u];e;e=Next[e]){
        int v=V[e];
        if (F[e]&&dep[v]==dep[u]+1){
            int minflow=dfs(v,min(flow-used,F[e]));
            F[e]-=minflow; F[e^1]+=minflow; used+=minflow;
            if (used==flow) break;
        }
    }
    return used;
}
int dinic(){
    while (bfs()){
        for (int i=1;i<=n;i++) cur[i]=Head[i];
        while (dfs(s,INF));
    }
    return maxflow;
}
int main(){
    scanf("%d%d%d%d",&n,&m,&s,&t);
    int u,v,f;
    for (int i=1;i<=m;i++){
        scanf("%d%d%d",&u,&v,&f);
        Add(u,v,f);
    }
    dinic();
    printf("%lld\n",maxflow);
    return 0;
}

```

### ISAP
```cpp
#include<cstdio>
#include<algorithm>
#include<cstring>
using namespace std;
const int INF=0x3f3f3f3f;
int n,m,s,t,gap[11000],cur[11000],d[11000],que[11000];
int edgenum,vet[210000],val[210000],Next[210000],Head[11000];
void addedge(int u,int v,int cost){
    vet[++edgenum]=v; val[edgenum]=cost;
    Next[edgenum]=Head[u]; Head[u]=edgenum;
}
void bfs(){
    memset(gap,0,sizeof(gap));
    memset(d,0,sizeof(d));
    d[t]=1; gap[1]++;
    int head=1,tail=1; que[1]=t;
    while (head<=tail){
        int u=que[head++];
        for (int e=Head[u];e;e=Next[e])
            if (!d[vet[e]]){
                d[vet[e]]=d[u]+1;
                gap[d[vet[e]]]++;
                que[++tail]=vet[e];
            }
    }
}
int dfs(int u,int flow){
    if (u==t) return flow;
    int used=0;
    for (int &e=cur[u];e;e=Next[e])
        if (d[u]==d[vet[e]]+1){
            int tmp=dfs(vet[e],min(flow,val[e]));
            used+=tmp; flow-=tmp;
            val[e]-=tmp; val[e^1]+=tmp;
            if (!flow) return used;
        }
    if (!(--gap[d[u]])) d[s]=n+1;
    d[u]++; gap[d[u]]++; cur[u]=Head[u];
    return used;
}
int maxflow(){
    bfs();
    for (int i=1;i<=n;i++) cur[i]=Head[i];
    int ans=dfs(s,INF);
    while (d[s]<=n) ans+=dfs(s,INF);
    return ans;
}
int main(){
    scanf("%d%d%d%d",&n,&m,&s,&t);
    int u,v,cost; edgenum=1;
    for (int i=1;i<=m;i++){
        scanf("%d%d%d",&u,&v,&cost);
        addedge(u,v,cost); addedge(v,u,0);
    }
    printf("%d\n",maxflow());
    return 0;
}
```




## 3.2 最小费用最大流
### EK+SPFA（建议）
```cpp
#include<bits/stdc++.h>
using namespace std;
const int INF=0x3f3f3f3f;
int n,m,s,t,ans,maxflow;
int vis[5100],dis[5100],lst[5100]; 
int edgenum=1,F[110000],V[110000],Next[110000],W[110000],Head[5100],cur[5100];
void addedge(int u,int v,int f,int w){
	V[++edgenum]=v; W[edgenum]=w; F[edgenum]=f;
	Next[edgenum]=Head[u];
	Head[u]=edgenum;
}
void Add(int u,int v,int f,int w){
	addedge(u,v,f,w);
	addedge(v,u,0,-w);
}
bool spfa(){
	for (int i=1;i<=n;i++) vis[i]=false,dis[i]=INF;
	dis[s]=0; vis[s]=true;
	queue<int> que;
	que.push(s);
	while (!que.empty()){
		int u=que.front(); que.pop();
		vis[u]=false;
		for (int e=Head[u];e;e=Next[e]){
			int v=V[e];
			if (dis[v]>dis[u]+W[e]&&F[e]){
				dis[v]=dis[u]+W[e];
				lst[v]=e;
				if (!vis[v]){
					que.push(v);
					vis[v]=true;
				}
			}
		}
	}
	return dis[t]!=INF;
}
int mincostmaxflow(){
	while (spfa()){
		int tmp=INF;
		for (int i=t;i!=s;i=V[lst[i]^1]) tmp=min(tmp,F[lst[i]]);
		for (int i=t;i!=s;i=V[lst[i]^1]) F[lst[i]]-=tmp,F[lst[i]^1]+=tmp;
		ans+=dis[t]*tmp; maxflow+=tmp;		
	}
	return maxflow;
}
int main(){
	scanf("%d%d%d%d",&n,&m,&s,&t);
	int u,v,f,w;
	for (int i=1;i<=m;i++){
		scanf("%d%d%d%d",&u,&v,&f,&w);
		Add(u,v,f,w);
	}
	mincostmaxflow();
	printf("%d %d\n",maxflow,ans);
	return 0;
}
```

### EK+1次SPFA+dijkstra（建议）
```cpp
#include<bits/stdc++.h>
using namespace std;
typedef pair<int,int> pii;
const int INF=0x3f3f3f3f;
int n,m,s,t,ans,maxflow;
int vis[5100],dis[5100],lst[5100],h[5100]; 
int edgenum=1,F[110000],V[110000],Next[110000],W[110000],Head[5100],cur[5100];
void addedge(int u,int v,int f,int w){
    V[++edgenum]=v; W[edgenum]=w; F[edgenum]=f;
    Next[edgenum]=Head[u];
    Head[u]=edgenum;
}
void Add(int u,int v,int f,int w){
    addedge(u,v,f,w);
    addedge(v,u,0,-w);
}
bool dijkstra(){
    for (int i=1;i<=n;i++) vis[i]=false,dis[i]=INF;
    dis[s]=0;
	priority_queue<pii> que; que.push(pii(0,s));
    while (!que.empty()){
        int u=que.top().second; que.pop();
        if (vis[u]) continue;
        vis[u]=true;
        for (int e=Head[u];e;e=Next[e]){
            int v=V[e];
            if (dis[v]>dis[u]+W[e]+h[u]-h[v]&&F[e]){
                dis[v]=dis[u]+W[e]+h[u]-h[v]; lst[v]=e;
                que.push(pii(-dis[v],v));
            }
        }
    }
    return dis[t]!=INF;
}
bool spfa(){
    for (int i=1;i<=n;i++) vis[i]=false,dis[i]=INF;
    dis[s]=0; vis[s]=true;
    queue<int> que; que.push(s);
    while (!que.empty()){
        int u=que.front(); que.pop();
        vis[u]=false;
        for (int e=Head[u];e;e=Next[e]){
            int v=V[e];
            if (dis[v]>dis[u]+W[e]&&F[e]){
                dis[v]=dis[u]+W[e]; lst[v]=e;
                if (!vis[v]){
                    que.push(v);
                    vis[v]=true;
                }
            }
        }
    }
    return dis[t]!=INF;
}

int mincostmaxflow(){
	spfa();
    while (dijkstra()){
        int tmp=INF;
        for (int i=t;i!=s;i=V[lst[i]^1]) tmp=min(tmp,F[lst[i]]);
        for (int i=t;i!=s;i=V[lst[i]^1]) F[lst[i]]-=tmp,F[lst[i]^1]+=tmp;
        ans+=(h[t]-h[s]+dis[t])*tmp; maxflow+=tmp;
    	for (int i=1;i<=n;i++) h[i]+=dis[i]; 
    }
    return maxflow;
}
int main(){
    scanf("%d%d%d%d",&n,&m,&s,&t);
    int u,v,f,w;
    for (int i=1;i<=m;i++){
        scanf("%d%d%d%d",&u,&v,&f,&w);
        Add(u,v,f,w);
    }
    mincostmaxflow();
    printf("%d %d\n",maxflow,ans);
    return 0;
}
```

### dinic
```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<queue>
using namespace std;
const int INF=0x3f3f3f3f;
int n,m,s,t;
int ans,maxflow;
int vis[5100];//是否到达过该点 
int dist[5100]; 
int edgenum=1,F[110000],V[110000],Next[110000],W[110000],Head[5100];
//f:最大流量（flow），w:单位费用 
void addedge(int u,int v,int f,int w){
    V[++edgenum]=v; W[edgenum]=w; F[edgenum]=f;
    Next[edgenum]=Head[u];
    Head[u]=edgenum;
}
void Add(int u,int v,int f,int w){
    addedge(u,v,f,w);
    addedge(v,u,0,-w);
}
bool spfa(){
    memset(vis,false,sizeof(vis));
    memset(dist,0x3f,sizeof(dist));
    dist[s]=0;
    vis[s]=true;
    queue<int> que;
    que.push(s);
    while (!que.empty()){
        int u=que.front(); que.pop();
        vis[u]=false;
        for (int e=Head[u];e;e=Next[e]){
            int d=V[e];
            if (dist[d]>dist[u]+W[e]&&F[e]){
                dist[d]=dist[u]+W[e];
                if (!vis[d]){
                    que.push(d);
                    vis[d]=true;
                }
            }
        }
    }
    return dist[t]<INF;
}
int dfs(int u,int flow){
    if (u==t){ 
        vis[t]=true; maxflow+=flow;
        return flow;
    }
    int used=0;
    vis[u]=true;
    for (int e=Head[u];e;e=Next[e]){
        int d=V[e];
        if ((!vis[d]||d==t)&&F[e]&&dist[d]==dist[u]+W[e]){
            int minflow=dfs(d,min(flow-used,F[e]));
            if (minflow!=0) ans+=W[e]*minflow,F[e]-=minflow,F[e^1]+=minflow,used+=minflow;
            if (used==flow) break;
        }
    }
    return used;
}
int mincostmaxflow(){
    while (spfa()){
        vis[t]=1;
        while (vis[t]){
            memset(vis,false,sizeof(vis));
            dfs(s,INF); 
        }
    }
    return maxflow;
}
int main(){
    n=read(),m=read(),s=read(),t=read();
    int u,v,f,w;
    for (register int i=1;i<=m;i++){
        u=read(),v=read(),f=read(),w=read();
        Add(u,v,f,w);
    }
    mincostmaxflow();
    printf("%d %d\n",maxflow,ans);
    return 0;
}
```


### EK+dijkstra（最小费用流！！！！不是最小费用最大流！！！）
```cpp
/*
CF164C
有n个任务，每个任务给定一个开始时间si、持续时间ti、所得收益ci；
有k台机器，每台机器一个时刻只能处理一个任务；
求如何安排，使得收益和最大。
*/
#include<bits/stdc++.h>
using namespace std;
typedef pair<int,int> pii;
const int INF=0x3f3f3f3f;
int n,k,s,t,ans,maxflow,u[5100],v[5100],c[5100];
int cnt,num[5100],h[5100],vis[5100],dis[5100],lst[5100]; 
int edgenum=1,F[110000],V[110000],Next[110000],W[110000],Head[5100],cur[5100];
void addedge(int u,int v,int f,int w){
	V[++edgenum]=v; W[edgenum]=w; F[edgenum]=f;
	Next[edgenum]=Head[u];
	Head[u]=edgenum;
}
void Add(int u,int v,int f,int w){
	if (w<0){
		num[u]-=f; num[v]+=f;
		addedge(v,u,f,-w);
//		cerr<<"#"<<v<<' '<<u<<' '<<f<<' '<<-w<<endl;
		addedge(u,v,0,w); ans+=w*f;
	} else{
		addedge(u,v,f,w);
//		cerr<<"#"<<u<<' '<<v<<' '<<f<<' '<<w<<endl;
		addedge(v,u,0,-w);
	}
}
bool dijkstra(){
    for (int i=s;i<=t;i++) vis[i]=false,dis[i]=INF;
    dis[s]=0;
	priority_queue<pii> que; que.push(pii(0,s));
    while (!que.empty()){
        int u=que.top().second; que.pop();
        if (vis[u]) continue;
        vis[u]=true;
        for (int e=Head[u];e;e=Next[e]){
            int v=V[e];
            if (dis[v]>dis[u]+W[e]+h[u]-h[v]&&F[e]){
                dis[v]=dis[u]+W[e]+h[u]-h[v]; lst[v]=e;
                que.push(pii(-dis[v],v));
            }
        }
    }
    return dis[t]!=INF;
}
int mincostmaxflow(){
    while (dijkstra()){
        int tmp=INF;
        for (int i=t;i!=s;i=V[lst[i]^1]) tmp=min(tmp,F[lst[i]]);
        for (int i=t;i!=s;i=V[lst[i]^1]) F[lst[i]]-=tmp,F[lst[i]^1]+=tmp;
        ans+=(h[t]-h[s]+dis[t])*tmp; maxflow+=tmp;
    	for (int i=s;i<=t;i++) h[i]+=dis[i]; 
    }
    return maxflow;
}
int main(){
	scanf("%d%d",&n,&k);
	for (int i=1;i<=n;i++){
		scanf("%d%d%d",&u[i],&v[i],&c[i]); v[i]+=u[i];
		num[++cnt]=u[i]; num[++cnt]=v[i];
	}
	sort(num+1,num+cnt+1); cnt=unique(num+1,num+cnt+1)-num-1;
	s=0; t=cnt+1;
	for (int i=1;i<=n;i++){
		u[i]=lower_bound(num+1,num+cnt+1,u[i])-num;
		v[i]=lower_bound(num+1,num+cnt+1,v[i])-num;
	}
	for (int i=1;i<=cnt;i++) num[i]=0;
	num[1]+=k; num[cnt]-=k;
	for (int i=1;i<=n;i++) Add(u[i],v[i],1,-c[i]);
	for (int i=1;i<=cnt;i++){
		if (num[i]>0) Add(s,i,num[i],0);
		if (num[i]<0) Add(i,t,-num[i],0);
	}
	for (int i=1;i<cnt;i++) Add(i,i+1,k,0);
	mincostmaxflow();
//	printf("%d\n",-ans);
	for (int i=1;i<=n;i++)
		printf("%d ",F[i<<1]);
	return 0;
}

```
## 3.3 无向图欧拉回路
```cpp
// 求出欧拉图中的一条欧拉回路
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
int head[100010], ver[1000010], Next[1000010], tot; // 邻接表
int stack[1000010], ans[1000010]; // 模拟系统栈，答案栈
bool vis[1000010];
int n, m, top, t;

void add(int x, int y) {
    ver[++tot] = y, Next[tot] = head[x], head[x] = tot;
}

void euler() {
    stack[++top] = 1;
    while (top > 0) {
        int x = stack[top], i = head[x];
        // 找到一条尚未访问的边
        while (i && vis[i]) i = Next[i];
        // 沿着这条边模拟递归过程，标记该边，并更新表头
        if (i) {
            stack[++top] = ver[i];
            head[x] = Next[i];
            vis[i] = vis[i ^ 1] = true;
        }        
        // 与x相连的所有边均已访问，模拟回溯过程，并记录于答案栈中
        else {
            top--;
            ans[++t] = x;
        }
    }
}

int main() {
    cin >> n >> m;
    tot = 1;
    for (int i = 1; i <= m; i++) {
        int x, y;
        scanf("%d%d", &x, &y);
        add(x, y), add(y, x);
    }
    euler();
    for (int i = t; i; i--) printf("%d\n", ans[i]);
}
```

## 3.4 无向图Tarjan算法

- 注意!的部分

```cpp
// tarjan算法求无向图的桥、边双连通分量并缩点
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
using namespace std;
const int SIZE = 100010;
int head[SIZE], ver[SIZE * 2], Next[SIZE * 2];
int dfn[SIZE], low[SIZE], c[SIZE];
int n, m, tot, num, dcc, tc;
bool bridge[SIZE * 2];
int hc[SIZE], vc[SIZE * 2], nc[SIZE * 2];

void add(int x, int y) {
    ver[++tot] = y, Next[tot] = head[x], head[x] = tot;
}

void add_c(int x, int y) {
    vc[++tc] = y, nc[tc] = hc[x], hc[x] = tc;
}

void tarjan(int x, int in_edge) {
    dfn[x] = low[x] = ++num;
    for (int i = head[x]; i; i = Next[i]) {
        int y = ver[i];
        if (!dfn[y]) {
            tarjan(y, i);
            low[x] = min(low[x], low[y]);
            if (low[y] > dfn[x])
                bridge[i] = bridge[i ^ 1] = true;
        }
        else if (i != (in_edge ^ 1))
            low[x] = min(low[x], dfn[y]);
    }
}

void dfs(int x) {
    c[x] = dcc;
    for (int i = head[x]; i; i = Next[i]) {
        int y = ver[i];
        if (c[y] || bridge[i]) continue;
        dfs(y);
    }
}

int main() {
    cin >> n >> m;
    tot = 1;
    for (int i = 1; i <= m; i++) {
        int x, y;
        scanf("%d%d", &x, &y);
        add(x, y), add(y, x);
    }
    for (int i = 1; i <= n; i++)
        if (!dfn[i]) tarjan(i, 0);
    for (int i = 2; i < tot; i += 2)
        if (bridge[i])
            printf("%d %d\n", ver[i ^ 1], ver[i]);

    for (int i = 1; i <= n; i++)
        if (!c[i]) {
            ++dcc;
            dfs(i);
        }
    printf("There are %d e-DCCs.\n", dcc);
    for (int i = 1; i <= n; i++)
        printf("%d belongs to DCC %d.\n", i, c[i]);

    tc = 1;
    for (int i = 2; i <= tot; i++) {
        int x = ver[i ^ 1], y = ver[i];
        if (c[x] == c[y]) continue;
        add_c(c[x], c[y]);
    }
    printf("缩点之后的森林，点数 %d，边数 %d\n", dcc, tc / 2);
    for (int i = 2; i < tc; i += 2)
        printf("%d %d\n", vc[i ^ 1], vc[i]);
}


// tarjan算法求无向图的割点、点双连通分量并缩点
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<vector>
using namespace std;
const int SIZE = 100010;
int head[SIZE], ver[SIZE * 2], Next[SIZE * 2];
int dfn[SIZE], low[SIZE], stack[SIZE], new_id[SIZE], c[SIZE];
int n, m, tot, num, root, top, cnt, tc;
bool cut[SIZE];
vector<int> dcc[SIZE];
int hc[SIZE], vc[SIZE * 2], nc[SIZE * 2];

void add(int x, int y) {
    ver[++tot] = y, Next[tot] = head[x], head[x] = tot;
}

void add_c(int x, int y) {
    vc[++tc] = y, nc[tc] = hc[x], hc[x] = tc;
}

void tarjan(int x) {
    dfn[x] = low[x] = ++num;
    stack[++top] = x;
    if (x == root && head[x] == 0) { // 孤立点 !!!!!!!!!!!!!!!!!!!!!!
        dcc[++cnt].push_back(x);
        return;
    }
    int flag = 0;
    for (int i = head[x]; i; i = Next[i]) {
        int y = ver[i];
        if (!dfn[y]) {
            tarjan(y);
            low[x] = min(low[x], low[y]);
            if (low[y] >= dfn[x]) {
                flag++;
                if (x != root || flag > 1) cut[x] = true;
                cnt++;
                int z;
                do {
                    z = stack[top--];
                    dcc[cnt].push_back(z);
                } while (z != y);
                dcc[cnt].push_back(x);
            }
        }
        else low[x] = min(low[x], dfn[y]);
    }
}

int main() {
    cin >> n >> m;
    tot = 1;
    for (int i = 1; i <= m; i++) {
        int x, y;
        scanf("%d%d", &x, &y);
        if (x == y) continue;
        add(x, y), add(y, x);
    }
    for (int i = 1; i <= n; i++)
        if (!dfn[i]) root = i, tarjan(i);
    for (int i = 1; i <= n; i++)
        if (cut[i]) printf("%d ", i);
    puts("are cut-vertexes");
    for (int i = 1; i <= cnt; i++) {
        printf("v-DCC #%d:", i);
        for (int j = 0; j < dcc[i].size(); j++)
            printf(" %d", dcc[i][j]);
        puts("");
    }
    // 给每个割点一个新的编号(编号从cnt+1开始)
    num = cnt;
    for (int i = 1; i <= n; i++)
        if (cut[i]) new_id[i] = ++num;
    // 建新图，从每个v-DCC到它包含的所有割点连边
    tc = 1;
    for (int i = 1; i <= cnt; i++)
        for (int j = 0; j < dcc[i].size(); j++) {
            int x = dcc[i][j];
            if (cut[x]) {
                add_c(i, new_id[x]);
                add_c(new_id[x], i);
            }
            else c[x] = i; // 除割点外，其它点仅属于1个v-DCC
        }
    printf("缩点之后的森林，点数 %d，边数 %d\n", num, tc / 2);
    printf("编号 1~%d 的为原图的v-DCC，编号 >%d 的为原图割点\n", cnt, cnt);
    for (int i = 2; i < tc; i += 2)
        printf("%d %d\n", vc[i ^ 1], vc[i]);
}

```



## 3.5 有向图Tarjan算法（强连通分量）

- 注意!的部分

```cpp
void tarjan(int u)
{
     int v;
     dfn[u]=low[u]=++idx;//每次dfs，u的次序号增加1
     s.push(u);//将u入栈
     ins[u]=1;//标记u在栈内
     for(int i=head[u];i!=-1;i=e[i].next)//访问从u出发的边
     {
         v=e[i].v;
         if(!dfn[v])//如果v没被处理过
         {
             tarjan(v);//dfs(v)
             low[u]=min(low[u],low[v]);//u点能到达的最小次序号是它自己能到达点的最小次序号和连接点v能到达点的最小次序号中较小的
         }
         else if(ins[v])low[u]=min(low[u],dfn[v]);//如果v在栈内，u点能到达的最小次序号是它自己能到达点的最小次序号和v的次序号中较小的
     }     
     if(dfn[u]==low[u])//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     {
         Bcnt++;
         do
         {
             v=s.top();
             s.pop();
             ins[v]=0;
             Belong[v]=Bcnt;
         }while(u != v);
     }
}
```

## 3.6 二分图匹配（匈牙利算法）

```cpp
bool dfs(int x) {
    for (int i = head[x], y; i; i = next[i])
        if (!visit[y = ver[i]]) {
            visit[y] = 1;
            if (!match[y] || dfs(match[y])) {
                match[y]=x;
                return true;
            }
        }
    return false;
}

for (int i = 1; i <= n; i++) {
    memset(visit, 0, sizeof(visit));
    if (dfs(i)) ans++;
}
```

## 3.7 二分图最大权匹配（KM算法）
```cpp
#include<bits/stdc++.h>
using namespace std;
const int INF=0x3f3f3f3f;
int n,m,w[110][110],label_x[110],label_y[110],match[110],slack[110],lst[110]; bool vis[110];
int main(){
//	freopen("matrix.in","r",stdin);
//	freopen("matrix.out","w",stdout);
	scanf("%d%d",&n,&m);
	for (int i=1;i<=n;i++)
		for (int j=1;j<=n;j++){
			scanf("%d",&w[i][j]);
			label_x[i]=max(label_x[i],w[i][j]);
		}
	for (int u=1;u<=n;u++){
		memset(slack,0x3f,sizeof(slack));
		memset(vis,false,sizeof(vis));
		match[0]=u; int i,j;
		for (i=0;match[i];i=j){
			vis[i]=true; int x=match[i],delta=INF;
			for (int y=1;y<=n;y++)
				if (!vis[y]){
					int tmp=label_x[x]+label_y[y]-w[x][y];
					if (tmp<slack[y]){
						slack[y]=tmp;
						lst[y]=i;
					}
					if (delta>slack[y]){
						delta=slack[y];
						j=y;
					}
				}
			for (int y=0;y<=n;y++)
				if (vis[y]){
					label_x[match[y]]-=delta;
					label_y[y]+=delta;
				} else slack[y]-=delta;
		}
		while (i) match[i]=match[lst[i]],i=lst[i];
	}
	int ans=0;
	for (int i=1;i<=n;i++) ans+=label_x[i]+label_y[i];
	printf("%d\n",ans);
//	printf("%d\n",ans*n);
//	if (m)
//		for (int i=1;i<=n;i++){
//			for (int j=1;j<=n;j++)
//				printf("%d ",label_x[i]+label_y[j]);
//			putchar('\n');
//		}
	return 0;
}

```

# 4. 矩阵&多项式

## 4.1 多项式基础模板

```cpp
struct poly{
    vector<int> v;
    poly(){ v.clear();}
    inline int size(){ return v.size();}
    friend poly operator+(poly A,poly B){
        int lenA=A.v.size()-1,lenB=B.v.size()-1;
        int len=max(lenA,lenB);
        A.v.resize(len+1); B.v.resize(len+1);
        poly C; C.v.resize(len+1);
        for (int i=0;i<=len;i++) C.v[i]=add(A.v[i],B.v[i]);
        return C;
    }
    friend poly operator*(poly A,int v){
        for (int i=0;i<A.size();i++) A.v[i]=mul(A.v[i],v);
        return A;
    }
    friend poly operator*(const poly &A,const poly &B){
        int len=A.v.size()+B.v.size()-2;
        poly C; C.v.resize(len+1);
        for (int i=0;i<(int)A.v.size();i++)
            for (int j=0;j<(int)B.v.size();j++)
                C.v[i+j]=add(C.v[i+j],mul(A.v[i],B.v[j]));
        return C;
    }
    friend poly operator%(poly A,const poly &B){
        int lenA=(int)A.v.size()-1,lenB=(int)B.v.size()-1;
        int inv=getinv(B.v[lenB]);
        if (lenA<lenB) return A;
        for (int i=lenA;i>=lenB;i--)
            if (A.v[i]){
                int t=mul(A.v[i],inv);
                for (int j=i;j>=i-lenB;j--) A.v[j]=dec(A.v[j],mul(t,B.v[lenB-i+j]));
            }
        A.v.resize(lenB+1); return A;
    }
};
poly qpow(poly x,char *a,const poly &M){//a为二进制数
    poly res; res.v.clear(); res.v.push_back(1);
    for (int i=strlen(a)-1;i>=0;i--){
        if (a[i]=='1') res=(res*x)%M;
        x=(x*x)%M;
    }
    return res;
}
```



## 4.2 FWT

```cpp
// 均类似于FFT/NTT的三重循环。
// Or卷积，正变换是大的+小的，逆变换是大的-小的。（两个数Or起来会变大）。
// And卷积，正变换是小的+大的，逆变换是小的-大的。（两个数And起来会变小）。
// Xor卷积，正变换和FFT/NTT一样，小的=(X+Y)，大的=(X-Y)。逆变换把结果都/2。
void FWT_or(int *a,int op){
    for (int i=1;i<N;i<<=1)
        for (int p=i<<1,j=0;j<N;j+=p)
            for (int k=0;k<i;++k)
                if (op==1) a[i+j+k]=(a[j+k]+a[i+j+k])%Mod;
                else a[i+j+k]=(a[i+j+k]+Mod-a[j+k])%Mod;
}
void FWT_and(int *a,int op){
    for (int i=1;i<N;i<<=1)
        for (int p=i<<1,j=0;j<N;j+=p)
            for (int k=0;k<i;++k)
                if (op==1) a[j+k]=(a[j+k]+a[i+j+k])%Mod;
                else a[j+k]=(a[j+k]+Mod-a[i+j+k])%Mod;
}
void FWT_xor(int *a,int op){
    for (int i=1;i<N;i<<=1)
        for (int p=i<<1,j=0;j<N;j+=p)
            for (int k=0;k<i;k++){
                int X=a[j+k],Y=a[i+j+k];
                a[j+k]=(X+Y)%Mod; a[i+j+k]=(X+Mod-Y)%Mod;
                if (op==-1){
                    a[j+k]=1ll*a[j+k]*inv2%Mod;
                    a[i+j+k]=1ll*a[i+j+k]*inv2%Mod;//inv2为2的逆元 
                }
            }
}
```



## 4.3 FFT

- 首先，重要的一点是```rev[x]=(rev[x>>1]>>1)|((x&1)<<(bit-1))```。相当于把它的最后一个位置删去，然后反过来，此时所有位置都往左移了一位，那么右移一位，最后把最后一个位置的数放在最高位。
- 然后，将原序列`x<rev[x]`的位置，交换`x,rev[x]`。（否则会交换两次，相当于没交换）。
- 考虑正变换
  - 枚举mid，表示当前长度的一半。此时，$\omega=\cos \frac{2\pi}{2mid}+sin\frac{2\pi}{2mid}=\cos\frac{\pi}{mid}+sin \frac{\pi}{mid}$。
  - 然后枚举起始位置$i$，必须是 $2mid$ 的倍数。
  - 初始时，$w=1$，然后，枚举段中的位置 $j$ ，则 $i+j$ 为左边段的真实位置，$i+j+mid$ 是右边段的真实位置。
  - $a[i + j] = (X + Y), a[i + j + mid] = w * (X - Y)$，最后$w*=\omega$。
- 接着考虑逆变换，其实就是改变了$\omega$，此时$\omega=cos \frac{\pi}{mid}-i\sin\frac{\pi}{mid}$。最后要把每个数都/n。

```cpp
#include<cstdio>
#include<cmath>
#include<algorithm>
using namespace std;
const int MAXN=1e7+10;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
const double Pi=acos(-1.0);
struct complex{
    double x,y;
    complex (double xx=0,double yy=0){ x=xx; y=yy;}
} a[MAXN],b[MAXN];
complex operator+(complex a,complex b){ return complex(a.x+b.x,a.y+b.y);}
complex operator-(complex a,complex b){ return complex(a.x-b.x,a.y-b.y);}
complex operator*(complex a,complex b){ return complex(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);} 
int N,M;
int l,r[MAXN];
int len=1;
void FFT(complex *A,int type){
    for (int i=0;i<len;i++) 
        if (i<r[i]) swap(A[i],A[r[i]]);
    for (int mid=1;mid<len;mid<<=1){
        complex Wn(cos(Pi/mid),type*sin(Pi/mid));
        for (int R=mid<<1,j=0;j<len;j+=R){
            complex w(1,0);
            for (int k=0;k<mid;k++,w=w*Wn){
                complex x=A[j+k],y=w*A[j+mid+k];
                A[j+k]=x+y; A[j+mid+k]=x-y;
            }
        }
    }
}
int main(){
    int N=read(),M=read();
    for (int i=0;i<=N;i++) a[i].x=read();
    for (int i=0;i<=M;i++) b[i].x=read();
    while (len<=N+M) len<<=1,l++;
    for (int i=0;i<len;i++) r[i]=(r[i>>1]>>1)|((i&1)<<(l-1));
    FFT(a,1); FFT(b,1);
    for (int i=0;i<=len;i++) a[i]=a[i]*b[i];
    FFT(a,-1);
    for (int i=0;i<=N+M;i++) printf("%d ",(int)(a[i].x/len+0.5));
    return 0;
}
```



## 4.4 NTT

- 正变换：改变了FFT中的$\omega$，此时$\omega = g ^ \dfrac{\varphi(Mod)}{mid * 2} = g ^ \dfrac{Mod - 1}{mid * 2}$。
- 逆变换：改变了FFT中的$\omega$，此时$\omega = \dfrac 1{g ^ \dfrac{\varphi(Mod)}{mid * 2}}=\dfrac 1{g ^ \dfrac{Mod - 1}{mid * 2}}$。

```cpp
#include<cstdio>
#include<algorithm>
using namespace std;
typedef long long ll;
const int g=3;
const int Mod=998244353;
const int MAXN=2100000;//>n+m上取到2的幂次 
int n,m,len,rev[MAXN];
ll a[MAXN],b[MAXN],c[MAXN];
inline ll qpow(ll x,int a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline int getinv(int x){ return qpow(x,Mod-2);}
void NTT(ll *a,int inv){
    for (int i=0;i<n;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<n;mid<<=1){
        int tmp=qpow(g,(Mod-1)/(mid<<1));
        if (inv==-1) tmp=getinv(tmp);
        for (int i=0;i<n;i+=mid*2){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=omega*tmp%Mod){
                int x=a[i+j],y=omega*a[i+j+mid]%Mod;
                a[i+j]=(x+y)%Mod,a[i+j+mid]=(x-y+Mod)%Mod;
            }
        }
    }
}
int main(){
    scanf("%d%d",&n,&m); n++; m++;
    for (int i=0;i<n;i++) scanf("%lld",&a[i]);
    for (int i=0;i<m;i++) scanf("%lld",&b[i]);
    len=n+m-1;
    int bit=0;
    while ((1<<bit)<len) bit++;
    n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(a,1); NTT(b,1);
    for (int i=0;i<n;i++) a[i]=a[i]*b[i]%Mod;
    NTT(a,-1);
    ll inv=getinv(n);
    for (int i=0;i<len;i++){
        c[i]=a[i]*inv%Mod;
        printf("%lld ",c[i]);
    }
    return 0;
}
```



## 4.5 多项式求逆

- 就一个式子，将变换后的按位乘改为：$b[i] = (2 - c[i] * b[i]) * b[i]$。

```cpp
void polymul(vec &a,vec &b,vec &c,int bit){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(n,a,1); NTT(n,b,1);
    for (int i=0;i<n;i++) a[i]=a[i]*b[i]%Mod;
    NTT(n,a,-1);
    ll inv=getinv(n); c.resize(n);
    for (int i=0;i<n;i++) c[i]=a[i]*inv%Mod;
    a.clear(); b.clear();
}
void polyinv(int len,int bit,vec &a,vec &b){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    vec c=a; c.resize(len);
    NTT(n,c,1); NTT(n,b,1);
    for (int i=0;i<n;i++) b[i]=1ll*(Mod+2-1ll*c[i]*b[i]%Mod)%Mod*b[i]%Mod;
    NTT(n,b,-1); b.resize(len);
    ll inv=getinv(n);
    for (int i=0;i<len;i++) b[i]=b[i]*inv%Mod;
}
int getbit(int x){
    int bit=1;
    while ((1<<bit)<x) bit++;
    return bit;
}
void getinv(int n,vec &a,vec &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    getinv((n+1)>>1,a,b);
    int bit=getbit(n<<1);
    polyinv(n,bit,a,b);
}
```

### 非递归版
```cpp
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int Mod=998244353;
const int inv2=(Mod+1)/2;
int rev[810000];
int f[810000],g[810000];
ll fac[110000],inv[110000];
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
void NTT(int *a,int n,int op){
    for (int i=0;i<n;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<n;mid<<=1){
        int tmp=qpow(3,(Mod-1)/(mid<<1));
        if (op==-1) tmp=qpow(tmp,Mod-2);
        for (int i=0;i<n;i+=(mid<<1)){
            ll omega=1;
            for (int j=0;j<mid;j++,omega=omega*tmp%Mod){
                int x=a[i+j],y=omega*a[i+j+mid]%Mod;
                a[i+j]=(x+y)%Mod; a[i+j+mid]=(x-y+Mod)%Mod;
            }
        }
    }
    if (op==1) return;
    ll inv=qpow(n,Mod-2);
    for (int i=0;i<n;i++) a[i]=a[i]*inv%Mod;
}
void polysqr(int *a,int n){
    int len; for (len=1;len<(n<<1);len<<=1);
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)?(len>>1):0);
    NTT(a,len,1);
    for (int i=0;i<len;i++) a[i]=1ll*a[i]*a[i]%Mod;
    NTT(a,len,-1);
    for (int i=n;i<len;i++) a[i]=0;
}
void polyinv(int *a,int *b,int n){
    static int A[810000],B[810000];
    b[0]=qpow(a[0],Mod-2);
    for (int mid=2;mid<(n<<1);mid<<=1){
        int len=mid<<1;
        for (int i=0;i<mid;i++) A[i]=a[i],B[i]=b[i];
        for (int i=mid;i<len;i++) A[i]=0,B[i]=0;
        for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)?mid:0);
        
        NTT(A,len,1); NTT(B,len,1);
        for (int i=0;i<len;i++) b[i]=1ll*(Mod+2-1ll*A[i]*B[i]%Mod)*B[i]%Mod;
        NTT(b,len,-1);
        for (int i=mid;i<len;i++) b[i]=0;
    }
}
void polysqrt(int *a,int *b,int n){
    static int A[810000],B[810000];
    b[0]=qpow(a[0],Mod-2);
    for (int mid=2;mid<(n<<1);mid<<=1){
        int len=mid<<1;
        for (int i=0;i<mid;i++) A[i]=a[i];
        for (int i=mid;i<len;i++) A[i]=0;
        polyinv(b,B,mid);
        for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)?mid:0);
        
        NTT(A,len,1); NTT(B,len,1);
        for (int i=0;i<len;i++) A[i]=1ll*A[i]*B[i]%Mod;
        NTT(A,len,-1);
        
        for (int i=0;i<mid;i++) b[i]=1ll*(b[i]+A[i])%Mod*inv2%Mod;
        for (int i=mid;i<len;i++) b[i]=0;
    }
}
int T,n;
int main(){
    fac[0]=1; for (int i=1;i<=100000;i++) fac[i]=fac[i-1]*i%Mod;
    inv[1]=1; for (int i=2;i<=100000;i++) inv[i]=(Mod-Mod/i)*inv[Mod%i]%Mod;
    inv[0]=1; for (int i=1;i<=100000;i++) inv[i]=inv[i-1]*inv[i]%Mod;
    for (int i=0;i<=100000;i++) f[i]=qpow(qpow(3,1ll*i*(i-1)/2),Mod-2)*inv[i]%Mod;
    polysqr(f,100001);
    for (int i=0;i<=100000;i++) f[i]=qpow(3,1ll*i*(i-1)/2)*f[i]%Mod;
    polysqrt(f,g,100001);
    
    scanf("%d",&T);
    while (T--){
        scanf("%d",&n);
        printf("%lld\n",1ll*g[n]*fac[n]%Mod);
    }
    return 0;
}

```

## 4.6 拉格朗日插值

- $\displaystyle f(x)=\sum_{i=1}^n y_i\prod_{i\not= j} \frac{x-x_j}{x_i-x_j} $

```cpp
    scanf("%lld%lld",&n,&k);//求f(k) 
    for (int i=1;i<=n;i++) scanf("%lld%lld",&x[i],&y[i]);//f(x[i])=y[i]
    ll s1,s2,ans=0;
    for (int i=1;i<=n;i++){
        s1=y[i]%Mod; s2=1ll;
        for (int j=1;j<=n;j++)
            if (i!=j){
                s1=s1*(k-x[j])%Mod;
                s2=s2*(x[i]-x[j])%Mod;
            }
        s1=(s1+Mod)%Mod; s2=(s2+Mod)%Mod;
        ans=(ans+s1*inv(s2))%Mod;
    }
    printf("%lld\n",ans);
    
    
//1..n
int Lagrange(int k){
    if (k<=n) return f[k];
    int res=0;
    for (int i=1;i<=n;i++){
        int s=1;
        for (int j=0;j<=n;j++){
            if (i==j) continue;
            s=1ll*s*((k-j)%Mod+Mod)%Mod;
            if (i>j) s=1ll*s*inv[i-j]%Mod;
            else s=(Mod-1ll*s*inv[j-i]%Mod)%Mod;
        }
        res=(res+1ll*f[i]*s%Mod)%Mod;
    }
    return res;
}
```

## 4.7 下降幂多项式乘法

```cpp
#include<cstdio>
#include<algorithm>
using namespace std;
typedef long long ll;
const int Mod=998244353;
const int G=3;
int n,m,len;
ll f[610000],g[610000],s[610000];
ll fac[610000],invfac[610000];
int rev[610000];
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
char readop(){
    char ch=Getchar();
    while (ch<'A'||ch>'Z') ch=Getchar();
    return ch;
}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
inline ll getinv(int x){ return qpow(x,Mod-2);}
void NTT(ll *a,int inv){
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<len;mid<<=1){
        int tmp=qpow(G,(Mod-1)/(mid<<1));
        if (inv==-1) tmp=getinv(tmp);
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }
    if (inv==-1){
        ll inv=getinv(len);
        for (int i=0;i<n;i++) a[i]=mul(a[i],inv);
    }
}
void FDT(ll *f,int inv){
    if (inv==1)
        for (int i=0;i<n;i++) s[i]=invfac[i];
    else
        for (int i=0;i<n;i++)
            if (i&1) s[i]=Mod-invfac[i];
            else s[i]=invfac[i]; 
    
    for (int i=n;i<len;i++) s[i]=0;
    
    NTT(f,1); NTT(s,1);
    for (int i=0;i<len;i++) f[i]=mul(f[i],s[i]);
    NTT(f,-1);
    
    for (int i=n;i<len;i++) f[i]=0;
}
int main(){
    n=read()+1; m=read()+1;
    for (int i=0;i<n;i++) f[i]=read();
    for (int i=0;i<m;i++) g[i]=read();
    
    n=n+m-1;
    int bit=0; while ((1<<bit)<n+n) bit++;
    len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    fac[0]=1; for (int i=1;i<=len;i++) fac[i]=fac[i-1]*i%Mod;
    invfac[1]=1; for (int i=2;i<=len;i++) invfac[i]=(Mod-Mod/i)*invfac[Mod%i]%Mod;
    invfac[0]=1; for (int i=1;i<=len;i++) invfac[i]=invfac[i-1]*invfac[i]%Mod;
    
    FDT(f,1); FDT(g,1);
    for (int i=0;i<n;i++) f[i]=mul(mul(f[i],g[i]),fac[i]);
    FDT(f,-1);
    for (int i=0;i<n;i++) printf("%lld ",f[i]);
    return 0;
}
```

## 4.8 多项式除法

- $A(x)/B(x)$
- $A$的最高次项为 $x^n$，$B$ 的最高次项为 $x^m$。
- 把 $A,B$ 翻转，$B$ 的最高此项改为 $x^{n-m}$。
- 求出 $B$ 的逆，和 $A$ 乘起来，再取 $x^0...x^{n-m}$。再翻转 ，输出。

```cpp
#include<cstdio>
#include<vector>
#include<algorithm>
//#define int long long
#define ls (now<<1)
#define rs (now<<1|1)
using namespace std;
typedef long long ll;
typedef vector<int> poly;
const int Mod=998244353;
const int G=3; const int invG=332748118;
int n,m,a[310000],rev[310000],ans[310000];
poly f,g;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(int x){ return qpow(x,Mod-2);}
void NTT(poly &a,int len,int inv){
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<len;mid<<=1){
        int tmp=qpow(G,(Mod-1)/(mid<<1));
        if (inv==-1) tmp=getinv(tmp);
//        printf("%d  ",(Mod-1)/(mid<<1));
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }
    if (inv==-1){
        ll inv=getinv(len);
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
void mult(poly a,poly b,poly &c,int n){
    int bit=0; while ((1<<bit)<n+n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(a,len,1); NTT(b,len,1);
    c.resize(len);
    for (int i=0;i<len;i++) c[i]=mul(a[i],b[i]);
    NTT(c,len,-1);
}
void polyinv(int len,int bit,poly &a,poly &b){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    poly c=a; c.resize(len); b.resize(len);
    NTT(c,n,1); NTT(b,n,1);
    for (int i=0;i<n;i++) b[i]=1ll*(Mod+2-1ll*c[i]*b[i]%Mod)%Mod*b[i]%Mod;
    NTT(b,n,-1); b.resize(len);
}
int getbit(int x){
    int bit=1;
    while ((1<<bit)<x) bit++;
    return bit;
}
void getinv(int n,poly &a,poly &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    getinv((n+1)>>1,a,b);
    int bit=getbit(n<<1);
    polyinv(n,bit,a,b);
}
poly c,d;
poly getdiv(poly a,poly b){
    int n=a.size()-1,m=b.size()-1;
    if (n<m){
        d.resize(1); d[0]=0;
        return d;
    }
    reverse(a.begin(),a.end());
    reverse(b.begin(),b.end());
    b.resize(n-m+1); c.clear();
    getinv(n-m+1,b,c); d.clear();
    mult(c,a,d,n+(n-m+1)); d.resize(n-m+1);
    reverse(d.begin(),d.end());
    return d;
}
poly ans1,tmp,ans2;
void getmod(poly &a,poly b,poly &c){
    int n=a.size()-1,m=b.size()-1;
    if (n<m) return;
    ans1=getdiv(a,b);
    mult(ans1,b,tmp,n+m);
    c.resize(m);
    for (int i=0;i<m;i++) c[i]=dec(a[i],tmp[i]);
}
int main(){
    n=read(); m=read();
    f.resize(n+1); g.resize(m+1);
    for (int i=0;i<=n;i++) f[i]=read();
    for (int i=0;i<=m;i++) g[i]=read();
    getmod(f,g,ans2); ans1.resize(n-m);
    for (int i=0;i<=n-m;i++) printf("%d ",ans1[i]); putchar('\n');
    for (int i=0;i<m;i++) printf("%d ",ans2[i]);
    return 0;
}
/*
5 1
15466465 9465465 2154654 64546 4546466 8664456
1545 75468
*/
```

## 4.9 多点求值
```cpp
#include<cstdio>
#include<vector>
#include<algorithm>
#define ls (now<<1)
#define rs (now<<1|1)
using namespace std;
typedef long long ll;
typedef vector<int> poly;
const int Mod=998244353;
const int G=3; const int invG=332748118;
int n,m,a[310000],rev[310000],ans[310000];
ll Inv[2100000],powG[2100000],powInvG[2100000];
poly f;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(const int &x){ return qpow(x,Mod-2);}
void NTT(poly &a,const int &len,const int &inv){
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<len;mid<<=1){
        int tmp;
        if (inv==1) tmp=powG[mid];
        else tmp=powInvG[mid];
//        int tmp=qpow(G,(Mod-1)/(mid<<1));
//        if (inv==-1) tmp=getinv(tmp);
//        printf("%d  ",(Mod-1)/(mid<<1));
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }
    if (inv==-1){
        ll inv=Inv[len];
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
/*void NTT(poly &a,int len,int inv){
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<len;mid<<=1){
        int tmp=qpow(G,(Mod-1)/(mid<<1));
        if (inv==-1) tmp=getinv(tmp);
//        printf("%d  ",(Mod-1)/(mid<<1));
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }
    if (inv==-1){
        ll inv=getinv(len);
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}*/
void mult(poly a,poly b,poly &c,const int &n){
    int bit=0; while ((1<<bit)<n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(a,len,1); NTT(b,len,1);
    c.resize(len);
    for (int i=0;i<len;i++) c[i]=mul(a[i],b[i]);
    NTT(c,len,-1);
}
void polyinv(const int &len,const int &bit,poly &a,poly &b){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    poly c=a; c.resize(len); b.resize(len);
    NTT(c,n,1); NTT(b,n,1);
    for (int i=0;i<n;i++) b[i]=1ll*(Mod+2-1ll*c[i]*b[i]%Mod)%Mod*b[i]%Mod;
    NTT(b,n,-1); b.resize(len);
}
inline int getbit(const int &x){
    int bit=1;
    while ((1<<bit)<x) bit++;
    return bit;
}
void getinv(const int &n,poly &a,poly &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    getinv((n+1)>>1,a,b);
    int bit=getbit(n<<1);
    polyinv(n,bit,a,b);
}
poly c,d;
poly getdiv(poly a,poly b){
    int n=a.size()-1,m=b.size()-1;
    if (n<m){
        d.resize(1); d[0]=0;
        return d;
    }
    reverse(a.begin(),a.end());
    reverse(b.begin(),b.end());
    b.resize(n-m+1); c.clear();
    getinv(n-m+1,b,c); d.clear();
    mult(c,a,d,n+(n-m+1)); d.resize(n-m+1);
    reverse(d.begin(),d.end());
    return d;
}
poly ans1,tmp,ans2;
void getmod(poly a,poly b,poly &c){
    while (!b.empty()&&!b.back()) b.pop_back();
    int n=a.size()-1,m=b.size()-1;
    if (n<m){ c=a; return;}
    ans1=getdiv(a,b);
    mult(ans1,b,tmp,n+m);
    c.resize(m);
    for (int i=0;i<m;i++) c[i]=dec(a[i],tmp[i]);
}
poly p[310000];
void getpoly(const int &now,const int &l,const int &r){
    if (l==r){
        p[now].push_back(dec(0,a[l]));
        p[now].push_back(1); return;
    }
    int mid=(l+r)>>1;
    getpoly(ls,l,mid); getpoly(rs,mid+1,r);
    mult(p[ls],p[rs],p[now],r-l+2);
}
void solve(const poly &A,const int &now,const int &l,const int &r){
    if (r-l<=1000){
        for (int i=l;i<=r;i++){
            int s=0;
            for (int j=A.size()-1;j>=0;j--) s=add(mul(s,a[i]),A[j]);
            ans[i]=s;
        }
        return;
    }
    poly B;
    int mid=(l+r)>>1;
    getmod(A,p[ls],B);
    solve(B,ls,l,mid);
    getmod(A,p[rs],B);
    solve(B,rs,mid+1,r);
}
int main(){
    for (int mid=1;mid<=2000000;mid<<=1){
        powG[mid]=qpow(G,(Mod-1)/(mid<<1));
        powInvG[mid]=qpow(powG[mid],Mod-2);
    }
    Inv[1]=1; for (int i=2;i<=2000000;i++) Inv[i]=(Mod-Mod/i)*Inv[Mod%i]%Mod;
    n=read(); m=read(); f.resize(n+1);
    for (int i=0;i<=n;i++) f[i]=read();
    for (int i=1;i<=m;i++) a[i]=read();
//    n=64000; m=64000; f.resize(n+1);
//    for (int i=0;i<=n;i++) f[i]=1;
//    for (int i=1;i<=m;i++) a[i]=i;
    getpoly(1,1,m);
    getmod(f,p[1],f);
    solve(f,1,1,m);
    for (int i=1;i<=m;i++) printf("%d\n",ans[i]);
    return 0;
}
```
## 4.10 多项式ln/exp

```cpp
#include<cstdio>
#include<vector>
using namespace std;
typedef vector<int> poly;
typedef long long ll;
const int Mod=998244353;
const int G=3;
int n,m,a[310000],rev[310000],ans[310000];
poly f,g;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(const int &x){ return qpow(x,Mod-2);}

void NTT(poly &a,int len,int inv){
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int mid=1;mid<len;mid<<=1){
        int tmp=qpow(G,(Mod-1)/(mid<<1));
        if (inv==-1) tmp=getinv(tmp);
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }
    if (inv==-1){
        ll inv=getinv(len);
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
void mult(poly a,poly b,poly &c,const int &n){
    int bit=0; while ((1<<bit)<n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(a,len,1); NTT(b,len,1);
    c.resize(len);
    for (int i=0;i<len;i++) c[i]=mul(a[i],b[i]);
    NTT(c,len,-1);
}
void polyinv(const int &len,const int &bit,poly &a,poly &b){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    poly c=a; c.resize(len); b.resize(len);
    NTT(c,n,1); NTT(b,n,1);
    for (int i=0;i<n;i++) b[i]=1ll*(Mod+2-1ll*c[i]*b[i]%Mod)%Mod*b[i]%Mod;
    NTT(b,n,-1); b.resize(len);
}
inline int getbit(const int &x){
    int bit=1;
    while ((1<<bit)<x) bit++;
    return bit;
}
void getinv(const int &n,poly &a,poly &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    getinv((n+1)>>1,a,b);
    int bit=getbit(n<<1);
    polyinv(n,bit,a,b);
}
void Dao(poly &A,poly &B,int len){
    for (int i=1;i<len;i++) B[i-1]=1ll*i*A[i]%Mod;
    B[len-1]=0;
}
void Jifen(poly &A,poly &B,int len){
    for (int i=1;i<len;i++) B[i]=1ll*A[i-1]*getinv(i)%Mod;
    B[0]=0;
}
poly A,B,C;
void getln(poly &f,poly &g,int n){
    f.resize(n); g.resize(n);
    A.resize(n); B.resize(n);
    Dao(f,A,n); getinv(n,f,B);
    mult(A,B,C,n+n);
    Jifen(C,g,n);
}
poly t;
void getexp(const int &n,poly &a,poly &b){
    if (n==1){
        b.resize(1); b[0]=1;
        return;
    }
    getexp(n>>1,a,b);
    getln(b,t,n);
    for (int i=0;i<n;i++) t[i]=dec(a[i],t[i]);
    t[0]=add(1,t[0]);
    mult(t,b,b,n+n); t.clear(); b.resize(n);
}
int main(){
    n=read(); f.resize(n);
    for (int i=0;i<n;i++) f[i]=read();
    for (m=1;m<=n;m<<=1);
    f.resize(m);
    getexp(m,f,g);
    for (int i=0;i<n;i++) printf("%d ",g[i]);
    return 0;
}
```

## 4.11 普通多项式转下降幂多项式
```cpp
#include<cstdio>
#include<algorithm>
#include<vector>
#define ls (now<<1)
#define rs (now<<1|1)
using namespace std;
typedef long long ll;
typedef vector<int> poly;
const int Mod=998244353;
const int G=3; const int invG=332748118;
int n,m,a[510000],rev[510000],ans[510000];
ll Inv[1100000];
int GPow[2][19][510000];
ll fac[1100000],invfac[1100000];
poly f;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(const int &x){ return qpow(x,Mod-2);}
inline void InitG(){
    for (int p=1;p<=18;p++) {
        int buf1=qpow(G,(Mod-1)/(1<<p));
        int buf0=qpow(invG,(Mod-1)/(1<<p));
        GPow[1][p][0]=GPow[0][p][0]=1;
        for (int i=1;i<(1<<p);i++){
            GPow[1][p][i]=1LL*GPow[1][p][i-1]*buf1%Mod;
            GPow[0][p][i]=1LL*GPow[0][p][i-1]*buf0%Mod;
        }
    }
}
void NTT(poly &a,const int &len,int inv){
    if (inv==-1) inv=0;
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    /*for (int mid=1;mid<len;mid<<=1){
        int tmp;
        if (inv==1) tmp=powG[mid];
        else tmp=powInvG[mid];
        for (int i=0;i<len;i+=(mid<<1)){
            ll omega=1;
            for (ll j=0;j<mid;j++,omega=mul(omega,tmp)){
                int x=a[i+j],y=mul(omega,a[i+j+mid]);
                a[i+j]=add(x,y),a[i+j+mid]=dec(x,y);
            }
        }
    }*/
    for (int l=2,cnt=1;l<=len;l<<=1,cnt++){
        int m=l>>1;
        for (int i=0;i<len;i+=l){
            int *buf=GPow[inv][cnt];
            for (int j=0;j<m;j++,buf++) {
                int x=a[i+j],y=1LL*(*buf)*a[i+j+m]%Mod;
                a[i+j]=add(x,y),a[i+j+m]=dec(x,y);
            }
        }
    }
    if (inv!=1){
        ll inv=Inv[len];
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
void mult(poly a,poly b,poly &c,const int &n){
    int bit=0; while ((1<<bit)<n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(a,len,1); NTT(b,len,1);
    c.resize(len);
    for (int i=0;i<len;i++) c[i]=mul(a[i],b[i]);
    NTT(c,len,-1);
}
void polyinv(const int &len,const int &bit,poly &a,poly &b){
    int n=1<<bit;
    for (int i=0;i<n;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    poly c=a; c.resize(len); b.resize(len);
    NTT(c,n,1); NTT(b,n,1);
    for (int i=0;i<n;i++) b[i]=mul(dec(2,mul(c[i],b[i])),b[i]);
    NTT(b,n,-1); b.resize(len);
}
inline int getbit(const int &x){
    int bit=1;
    while ((1<<bit)<x) bit++;
    return bit;
}
void getinv(const int &n,poly &a,poly &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    getinv((n+1)>>1,a,b);
    int bit=getbit(n<<1);
    polyinv(n,bit,a,b);
}
poly c,d;
poly getdiv(poly a,poly b){
    int n=a.size()-1,m=b.size()-1;
    if (n<m){
        d.resize(1); d[0]=0;
        return d;
    }
    reverse(a.begin(),a.end());
    reverse(b.begin(),b.end());
    b.resize(n-m+1); c.clear();
    getinv(n-m+1,b,c); d.clear();
    mult(c,a,d,n+(n-m+1)); d.resize(n-m+1);
    reverse(d.begin(),d.end());
    return d;
}
poly ans1,tmp,ans2;
void getmod(poly a,poly b,poly &c){
    while (!b.empty()&&!b.back()) b.pop_back();
    int n=a.size()-1,m=b.size()-1;
    if (n<m){ c=a; return;}
    ans1=getdiv(a,b);
    mult(ans1,b,tmp,n+m);
    c.resize(m);
    for (int i=0;i<m;i++) c[i]=dec(a[i],tmp[i]);
}
poly p[510000];
void getpoly(const int &now,const int &l,const int &r){
    if (l==r){
        p[now].push_back(dec(0,a[l]));
        p[now].push_back(1); return;
    }
    int mid=(l+r)>>1;
    getpoly(ls,l,mid); getpoly(rs,mid+1,r);
    mult(p[ls],p[rs],p[now],r-l+2);
}
void solve(const poly &A,const int &now,const int &l,const int &r){
    if (r-l<=700){
        for (int i=l;i<=r;i++){
            int s=0;
            for (int j=A.size()-1;j>=0;j--) s=add(mul(s,a[i]),A[j]);
            ans[i]=s;
        }
        return;
    }
    poly B;
    int mid=(l+r)>>1;
    getmod(A,p[ls],B);
    solve(B,ls,l,mid);
    getmod(A,p[rs],B);
    solve(B,rs,mid+1,r);
}
poly s;
void FDT(poly &f,int len,int inv){
    f.resize(len); s.resize(len);
    if (inv==1)
        for (int i=0;i<n;i++) s[i]=invfac[i];
    else
        for (int i=0;i<n;i++)
            if (i&1) s[i]=Mod-invfac[i];
            else s[i]=invfac[i]; 
    
    for (int i=n;i<len;i++) s[i]=0;
    
    NTT(f,len,1); NTT(s,len,1);
    for (int i=0;i<len;i++) f[i]=mul(f[i],s[i]);
    NTT(f,len,-1);
    
    for (int i=n;i<len;i++) f[i]=0;
}
int main(){
    Inv[1]=1; for (int i=2;i<=1000000;i++) Inv[i]=(Mod-Mod/i)*Inv[Mod%i]%Mod;
    InitG();
    n=read(); f.resize(n);
//    n=100000; f.resize(n); 
//    for (int i=0;i<n;i++) f[i]=i;
    for (int i=0;i<n;i++) f[i]=read();
    m=n; for (int i=1;i<=m;i++) a[i]=i-1;
    getpoly(1,1,m);
    getmod(f,p[1],f);
    solve(f,1,1,m);
    
    int bit=0; while ((1<<bit)<n+n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    fac[0]=1; for (int i=1;i<=len;i++) fac[i]=fac[i-1]*i%Mod;
    invfac[1]=1; for (int i=2;i<=len;i++) invfac[i]=(Mod-Mod/i)*invfac[Mod%i]%Mod;
    invfac[0]=1; for (int i=1;i<=len;i++) invfac[i]=invfac[i-1]*invfac[i]%Mod;
    f.resize(n);
    for (int i=0;i<n;i++) f[i]=mul(ans[i+1],invfac[i]);
    FDT(f,len,-1);
    for (int i=0;i<n;i++) printf("%d ",f[i]);    
    return 0;
}
```

## 4.12 下降幂多项式转普通多项式
```cpp
#include<cstdio>
#include<vector>
using namespace std;
const int G=3,invG=332748118;
const int Mod=998244353;
typedef long long ll;
typedef vector<int> poly;
struct node{
    poly g,f;
};
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(int x){ return qpow(x,Mod-2);}
int n,GPow[2][20][1100000];
int rev[1100000]; ll fac[1100000],inv[1100000],invfac[1100000];
poly f,tmp,ans;
inline void InitG(){
    fac[0]=1; for (int i=1;i<=1000000;i++) fac[i]=fac[i-1]*i%Mod;
    inv[1]=1; for (int i=2;i<=1000000;i++) inv[i]=(Mod-Mod/i)*inv[Mod%i]%Mod;
    invfac[0]=1; for (int i=1;i<=1000000;i++) invfac[i]=invfac[i-1]*inv[i]%Mod;
    for (int p=1;p<=19;p++) {
        int buf1=qpow(G,(Mod-1)/(1<<p));
        int buf0=qpow(invG,(Mod-1)/(1<<p));
        GPow[1][p][0]=GPow[0][p][0]=1;
        for (int i=1;i<(1<<p);i++){
            GPow[1][p][i]=1LL*GPow[1][p][i-1]*buf1%Mod;
            GPow[0][p][i]=1LL*GPow[0][p][i-1]*buf0%Mod;
        }
    }
}
void NTT(poly &a,const int &len,int inv){
    if (inv==-1) inv=0;
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int l=2,cnt=1;l<=len;l<<=1,cnt++){
        int m=l>>1;
        for (int i=0;i<len;i+=l){
            int *buf=GPow[inv][cnt];
            for (int j=0;j<m;j++,buf++){
                int x=a[i+j],y=1LL*(*buf)*a[i+j+m]%Mod;
                a[i+j]=add(x,y),a[i+j+m]=dec(x,y);
            }
        }
    }
    if (inv!=1){
        ll inv=getinv(len);
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
poly A,B;
void mult(int n,int m,const poly &a,const poly &b,poly &c){
    n++; m++;
    int len=1,bit=0;
    while (len<(n+m)) len<<=1,bit++;
    A.resize(n); B.resize(m);
    for (int i=0;i<n;i++) A[i]=a[i];
    for (int i=0;i<m;i++) B[i]=b[i];
    c.resize(len);
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(A,len,1); NTT(B,len,1);
    for (int i=0;i<len;i++) c[i]=1ll*A[i]*B[i]%Mod;
    NTT(c,len,-1);
}
poly s;
void ADD(poly &a,const poly &b){
    if (b.size()>a.size()) a.resize(b.size());
    for (int i=0;i<(int)b.size();i++) a[i]=add(a[i],b[i]);    
}
void FDT(poly &f,int len,int inv){
    f.resize(len); s.resize(len);
    if (inv==1)
        for (int i=0;i<n;i++) s[i]=invfac[i];
    else
        for (int i=0;i<n;i++)
            if (i&1) s[i]=Mod-invfac[i];
            else s[i]=invfac[i]; 
    
    for (int i=n;i<len;i++) s[i]=0;
    
    NTT(f,len,1); NTT(s,len,1);
    for (int i=0;i<len;i++) f[i]=mul(f[i],s[i]);
    NTT(f,len,-1);
    f.resize(n);
}
node interpolation(int l,int r,const poly &A,const poly &B){
    if (l==r){
        node tmp;
        tmp.g.clear(); tmp.f.clear();
        tmp.g.push_back(dec(0,l)); tmp.g.push_back(1);
        tmp.f.push_back(mul(A[l],B[l]));
        return tmp;
    }
    int mid=(l+r)>>1;
    node a=interpolation(l,mid,A,B);
    node b=interpolation(mid+1,r,A,B);
    poly c; mult(mid-l+1,r-mid,a.g,b.g,c);
    poly d; mult(mid-l+1,r-mid-1,a.g,b.f,d);
    poly e; mult(mid-l,r-mid,a.f,b.g,e);
    ADD(d,e);
    return (node){c,d};
}
int main(){
    InitG();
    n=read(); f.resize(n); tmp.resize(n);
    for (int i=0;i<n;i++) f[i]=read();
    int bit=0; while ((1<<bit)<n+n) bit++;
    int len=1<<bit;
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));f.resize(n); tmp.resize(n);
    FDT(f,len,1);    
    for (int i=0;i<n;i++) f[i]=mul(f[i],fac[i]);
    for (int i=0;i<n;i++) tmp[i]=mul(invfac[i],((n-i)&1)?invfac[n-i-1]:Mod-invfac[n-i-1]);
    ans=interpolation(0,n-1,f,tmp).f;
    ans.resize(n);
    for (int i=0;i<n;i++) printf("%d ",ans[i]);
    return 0;
}
```

## 4.13 常系数齐次线性递推
```cpp
#include<cstdio>
#include<vector>
#include<algorithm>
using namespace std;
const int G=3,invG=332748118;
const int Mod=998244353;
typedef long long ll;
typedef vector<int> poly;
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
inline int add(int x,int y){return x+y>=Mod?x+y-Mod:x+y;}
inline int dec(int x,int y){return x-y<0?x-y+Mod:x-y;}
inline int mul(int x,int y){return 1ll*x*y%Mod;}
ll qpow(ll x,ll a){
    ll res=1;
    while (a){
        if (a&1) res=res*x%Mod;
        x=x*x%Mod; a>>=1;
    }
    return res;
}
inline ll getinv(int x){ return qpow(x,Mod-2);}
int GPow[2][19][410000];
int rev[410000],f[410000],g[410000];
poly A,B,c,a,b;
inline void InitG(){
    for (int p=1;p<=18;p++) {
        int buf1=qpow(G,(Mod-1)/(1<<p));
        int buf0=qpow(invG,(Mod-1)/(1<<p));
        GPow[1][p][0]=GPow[0][p][0]=1;
        for (int i=1;i<(1<<p);i++){
            GPow[1][p][i]=1LL*GPow[1][p][i-1]*buf1%Mod;
            GPow[0][p][i]=1LL*GPow[0][p][i-1]*buf0%Mod;
        }
    }
}
void NTT(poly &a,const int &len,int inv){
    if (inv==-1) inv=0;
    a.resize(len);
    for (int i=0;i<len;i++)
        if (i<rev[i]) swap(a[i],a[rev[i]]);
    for (int l=2,cnt=1;l<=len;l<<=1,cnt++){
        int m=l>>1;
        for (int i=0;i<len;i+=l){
            int *buf=GPow[inv][cnt];
            for (int j=0;j<m;j++,buf++) {
                int x=a[i+j],y=1LL*(*buf)*a[i+j+m]%Mod;
                a[i+j]=add(x,y),a[i+j+m]=dec(x,y);
            }
        }
    }
    if (inv!=1){
        ll inv=getinv(len);
        for (int i=0;i<len;i++) a[i]=mul(a[i],inv);
    }
}
void mult(int n,int m,const poly &a,const poly &b,poly &c){
    int len=1,bit=0;
    while (len<(n+m)) len<<=1,bit++;
    A.resize(n); B.resize(m);
    for (int i=0;i<n;i++) A[i]=a[i];
    for (int i=0;i<m;i++) B[i]=b[i];
    c.resize(len);
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    NTT(A,len,1); NTT(B,len,1);
    for (int i=0;i<len;i++) c[i]=1ll*A[i]*B[i]%Mod;
    NTT(c,len,-1);
}
void PolyInv(int n,poly &a,poly &b){
    if (n==1){
        b.resize(1);
        b[0]=getinv(a[0]);
        return;
    }
    PolyInv((n+1)>>1,a,b);
    int len=1,bit=0;
    while (len<(n<<1)) len<<=1,bit++;
    c.resize(n);
    for (int i=0;i<n;i++) c[i]=a[i];
    for (int i=0;i<len;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    
    NTT(c,len,1); NTT(b,len,1);
    for (int i=0;i<len;i++) b[i]=1ll*(2-1ll*c[i]*b[i]%Mod+Mod)%Mod*b[i]%Mod;
    NTT(b,len,-1); b.resize(n);
}
poly fr,gr,tmp,q,R,Q;
void PolyMod(int n,int m,poly &f,poly &g,poly &r){
    fr.clear(); tmp.clear(); q.clear(); R.clear();
    f.resize(n+1); fr.resize(n+1);
    for (int i=0;i<=n;i++) fr[n-i]=f[i];
    mult(n-m+1,n,gr,fr,q);
    q.resize(n-m+1);
    reverse(q.begin(),q.end());
    mult(m+1,n-m+1,g,q,R);
    for (int i=0;i<m;i++) r[i]=(f[i]-R[i]+Mod)%Mod;
    r.resize(m);
}
int n,k; poly res,x;
int val[410000];
int main(){
    InitG();
    n=read(); k=read();
    a.resize(k+1); res.resize(k+1); x.resize(k+1);
    for (int i=1;i<=k;i++) a[k-i]=(Mod-read()%Mod)%Mod; 
    for (int i=0;i<=k-1;i++) val[i]=(read()%Mod+Mod)%Mod;
    a[k]=res[0]=1;
    x[1]=1;
    gr.resize(k+1);
    for (int i=0;i<=k;i++) gr[k-i]=a[i];
    PolyInv(k+1,gr,tmp);
    for (int i=0;i<=k;i++) gr[i]=tmp[i];
    while (n){
        if (n&1){
            mult(k,k,res,x,res);
            PolyMod(2*k+1,k,res,a,res);
        }
        mult(k,k,x,x,x);
        PolyMod(2*k+1,k,x,a,x);
        n>>=1;
    }
    int ans=0; res.resize(k);
    for (int i=0;i<k;i++) ans=(ans+1ll*res[i]*val[i])%Mod;
    printf("%d\n",ans);
    return 0;
}
```

# 5.字符串

## 5.1 后缀自动机（SAM）

```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int S=1;
int z,n,cnt,sz[2100000],a[2100000],c[2100000],len[2100000];
int link[2100000],trans[2100000][27];
long long ans;
char s[1100000];
inline void add(int loc){
    int c=s[loc]-'a',v=z;
    z=++cnt; len[z]=loc;
    for (;v&&!trans[v][c];v=link[v]) trans[v][c]=z;
    if (!v) link[z]=S;
    else{
        int x=trans[v][c];
        if (len[v]+1==len[x]) link[z]=x;
        else
        {
            int y=++cnt; len[y]=len[v]+1;
            for (int i=0;i<=26;i++) trans[y][i]=trans[x][i];
            link[y]=link[x]; link[x]=y; link[z]=y;
            for (;trans[v][c]==x;v=link[v]) trans[v][c]=y;
        }
    }
    sz[z]=1;
}
void calc(){
    for (int i=1;i<=cnt;i++) c[len[i]]++;
    for (int i=1;i<=cnt;i++) c[i]+=c[i-1];
    for (int i=1;i<=cnt;i++) a[c[len[i]]--]=i;
    for (int i=cnt;i;i--){
        int u=a[i];
        sz[link[u]]+=sz[u];
        ans+=1ll*(len[u]-len[link[u]])*sz[u]*(n-sz[u]);
    }
}
int main(){
    scanf("%s",s+1); n=strlen(s+1);
    z=1; cnt=1;
    for (int i=1;i<=n;i++) add(i);
    calc();
    printf("%lld\n",ans);
    return 0;
}

```

## 5.2 回文自动机（PAM）

```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
typedef long long ll;
const int MAXN=310000;
const int Mod=1e9+7;
const int inv2=5e8+4;
char s[MAXN];
class Palindromic_tree{
    /*
        两棵树，记录长度为奇数和偶数的回文串
        0为长度为偶数的根,1为长度为奇数的根

        *1.    偶数根节点字符长度为0，不存在，将它看做一个不在字符集中的字符-1
        *2.    奇数根节点字符长度为-1，作用相当于backspace，方便处理单字符形成的回文串
        *3.    只更改了最长的满足条件回文串的个数，
            最后要从1->tot，把fail,fail[fail],...都也加上，这样显然会T
            所以要倒过来，从tot->1，cnt[fail]+=cnt[now];
    */
    private:
        int now,tot,s[MAXN];
        int fail[MAXN],cnt[MAXN],len[MAXN],ch[MAXN][26];
        inline int newnode(int x){                               //新建一个节点，长度为x
            len[++tot]=x;
            return tot;
        }
        inline int getfail(int x,int r){                         //从l节点开始找s[r]的fail
            while (s[r-len[x]-1]!=s[r]) x=fail[x];               //一直跳直到找到后缀回文
            return x;
        }
    public:
        void build(const char *t,int n){
            for (int i=1;i<=n;i++) s[i]=t[i]-'a';
            s[0]=-1; fail[0]=1; len[0]=0;                        //*1.
            len[1]=-1;                                           //*2.
            now=0; tot=1;                                        //当前在0号节点，用了<=1的位置
            int x,y;
            for (int i=1;i<=n;i++){
                x=getfail(now,i);                                //找到回文的位置
                if (!ch[x][s[i]]){                               //没转移过
                    y=newnode(len[x]+2);                         //前后都加上这个字符，长度加2
                    fail[y]=ch[getfail(fail[x],i)][s[i]];        //记录新的fail
                    ch[x][s[i]]=y;
                }
                now=ch[x][s[i]]; cnt[now]++;                     //记录出现次数（只记录了最长的满足条件的）*3.
            }
        }
        ll getans(){
            ll res=0;
            for (int i=tot;i>=1;i--){
                cnt[fail[i]]+=cnt[i];
                res=(res+1ll*cnt[i]*(cnt[i]-1)%Mod*inv2)%Mod;
            }
            return res;
        }
} tree;
int main(){
    scanf("%s",s+1);
    tree.build(s,strlen(s+1));
    printf("%lld\n",tree.getans());
    return 0;
}

```

## 5.3 后缀数组（SA）

```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
char s[1100000];
int y[1100000],x[1100000],c[1100000],sa[1100000],rk[1100000],height[1100000],wt[30];
int n,m;
void get_SA(){
    for (int i=1;i<=n;++i) ++c[x[i]=s[i]];
    //c数组是桶 
    //x[i]是第i个元素的第一关键字 
    for (int i=2;i<=m;++i) c[i]+=c[i-1]; 
    //做c的前缀和，我们就可以得出每个关键字最多是在第几名 
    for (int i=n;i>=1;--i) sa[c[x[i]]--]=i; 
    for (int k=1;k<=n;k<<=1){
        int num=0;
        for (int i=n-k+1;i<=n;++i) y[++num]=i;
        //y[i]表示第二关键字排名为i的数，第一关键字的位置 
        //第n-k+1到第n位是没有第二关键字的 所以排名在最前面 
        for (int i=1;i<=n;++i) if (sa[i]>k) y[++num]=sa[i]-k;
        //排名为i的数 在数组中是否在第k位以后
        //如果满足(sa[i]>k) 那么它可以作为别人的第二关键字，就把它的第一关键字的位置添加进y就行了
        //所以i枚举的是第二关键字的排名，第二关键字靠前的先入队 
        for (int i=1;i<=m;++i) c[i]=0;
        //初始化c桶 
        for (int i=1;i<=n;++i) ++c[x[i]];
        //因为上一次循环已经算出了这次的第一关键字 所以直接加就行了 
        for (int i=2;i<=m;++i) c[i]+=c[i-1];//第一关键字排名为1~i的数有多少个 
        for (int i=n;i>=1;--i) sa[c[x[y[i]]]--]=y[i],y[i]=0;
        //因为y的顺序是按照第二关键字的顺序来排的 
        //第二关键字靠后的，在同一个第一关键字桶中排名越靠后 
        //基数排序 
        swap(x,y);
        //这里不用想太多，因为要生成新的x时要用到旧的，就把旧的复制下来，没别的意思 
        x[sa[1]]=1;num=1;
        for (int i=2;i<=n;++i)
            x[sa[i]]=(y[sa[i]]==y[sa[i-1]] && y[sa[i]+k]==y[sa[i-1]+k]) ? num : ++num;
        //因为sa[i]已经排好序了，所以可以按排名枚举，生成下一次的第一关键字 
        if (num==n) break;
        m=num;
        //这里就不用那个122了，因为都有新的编号了 
    }
    for (int i=1;i<=n;++i) printf("%d ",sa[i]);
}
void get_height(){
    int k=0;
    for (int i=1;i<=n;++i) rk[sa[i]]=i;
    for (int i=1;i<=n;++i){
        if (rk[i]==1) continue;//第一名height为0 
        if (k) --k;//h[i]>=h[i-1]-1;
        int j=sa[rk[i]-1];
        while (j+k<=n && i+k<=n && s[i+k]==s[j+k]) ++k;
        height[rk[i]]=k;//h[i]=height[rk[i]];
    }
    putchar(10);
    for (int i=1;i<=n;++i) printf("%d ",height[i]);
}
int main(){
    gets(s+1);
    n=strlen(s+1); m=122;
    //因为这个题不读入n和m所以要自己设
    //n表示原字符串长度，m表示字符个数，ascll('z')=122 
    //我们第一次读入字符直接不用转化，按原来的ascll码来就可以了 
    //因为转化数字和大小写字母还得分类讨论，怪麻烦的 
    get_SA(); //get_height();
    return 0;
}

```

### 新的板子
```cpp
struct SA{
	int f[N],g[N],a[N],h[N],val[20][N];
	int sa[N],rk[N],oldrk[N];
	void build(){
		memset(f,0,sizeof(f)); int t;
		for (int i=1;i<=n;i++) f[s[i]-'a'+1]=1;
		for (int i=1;i<=26;i++) f[i]+=f[i-1]; t=f[26];
		for (int i=1;i<=n;i++) rk[i]=f[s[i]-'a'+1];
		for (int i=1;i<n;i<<=1){
			for (int j=1;j<=n;j++) oldrk[j]=rk[j];
			for (int j=1;j<=n;j++) a[j]=j+i<=n?rk[j+i]:0;//后i位的rank
			
			memset(f,0,sizeof(f));
			for (int j=1;j<=n;j++) f[a[j]]++;
			for (int j=1;j<=t;j++) f[j]+=f[j-1];
			for (int j=1;j<=n;j++) g[f[a[j]]--]=j;//以后i位的rank为关键字的排序结果
			
			memset(f,0,sizeof(f));
			for (int j=1;j<=n;j++) f[rk[j]]++;
			for (int j=1;j<=t;j++) f[j]+=f[j-1];
			for (int j=n;j>=1;j--) sa[f[rk[g[j]]]--]=g[j];//以前i位rank为第一关键字，后i位rank为第二关键字的排序结果
			
			t=0;
			for (int j=1;j<=n;j++){
				t+=oldrk[sa[j]]>oldrk[sa[j-1]]||(oldrk[sa[j]]==oldrk[sa[j-1]]&&a[sa[j]]>a[sa[j-1]]);
				rk[sa[j]]=t;
			}
		}
		for (int j=1;j<=n;j++) sa[rk[j]]=j;
		int j=0;
		for (int i=1;i<=n;i++){
			if (rk[i]==n) j=0;
			else{
				if (j) j--;
				while (s[i+j]==s[sa[rk[i]+1]+j]) j++;
			}
			h[rk[i]]=j;
		}
		
		for (int i=1;i<=n;i++) val[0][i]=h[i];
		for (int j=1;(1<<j)<=n;j++)
			for (int i=1;i+(1<<j)-1<=n;i++)
				val[j][i]=min(val[j-1][i],val[j-1][i+(1<<(j-1))]);
	}
	int lcp(int x,int y){
		x=rk[x]; y=rk[y];
		if (x>y) swap(x,y);
		y--; int k=Log2[y-x+1];
		return min(val[k][x],val[k][y-(1<<k)+1]);
	}
};
```

## 5.4 Manacher

```cpp
#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
int p[21000000];
char s[11000000];
char t[21000000];//存添加字符后的字符串 

int init(){//形成新的字符串 
    int len=strlen(s);//len是输入字符串的长度
    t[0]='$';//处理边界，防止越界 
    t[1]='#';
    int j=2; 
    for (int i=0;i<len;i++){
        t[j++]=s[i];
        t[j++]='#';
    } 
    t[j]='\0';//处理边界，防止越界（容易忘记） 
    return j;// 返回t的长度 
}
int Manacher(){//返回最长回文串 
    int len=init();//取得新字符串的长度， 完成向t的转换
    int max_len=-1;//最长回文长度
    int id;
    int mx=0;
    for (int i=1;i<=len;i++){
        if (i<mx)
            p[i]=min(p[2*id-i],mx-i);//上面图片就是这里的讲解 
        else p[i]=1;
        while (t[i-p[i]]==t[i+p[i]])//不需边界判断，因为左有'$'，右有'\0'标记；
            p[i]++;//mx对此回文中点的贡献已经结束，现在是正常寻找扩大半径
        if (mx<i+p[i]){//每走移动一个回文中点，都要和mx比较，使mx是最大，提高p[i]=min(p[2*id-i],mx-i)效率 
            id=i;//更新id 
            mx=i+p[i];//更新mx 
        }
        max_len=max(max_len,p[i]-1); 
    } 
    return max_len; 
}
int main(){
    scanf("%s",s);
    printf("%d",Manacher());
    return 0;
}
```

## 5.5 0/1可持久化trie

```cpp
void change(int &x,int dep,int v){
    tree[++cnt]=tree[x]; x=cnt;
    tree[x].sz++;
    if (dep==-1) return;
    bool k=v&(1<<dep);
    if (k) change(tree[x].ch[1],dep-1,v);
    else change(tree[x].ch[0],dep-1,v);
}
void build(){
    cnt=0;
    for (int i=1;i<=n;i++){
        rt[i]=rt[i-1];
        change(rt[i],30,a[i]);
    }
}
```

## 5.6 离线可删线性基

```cpp
#include<cstdio>
#include<algorithm>
#include<unordered_map>
using namespace std;
int n,m,Q,ans,sz,op[2100000],x[2100000],lst[2100000],nxt[2100000];
bool vis[35];
struct node{
    int x,y;
    bool operator>(const node &a) const{
        return nxt[y]>nxt[a.y];
    }
} c[35];
char Getchar(){
    static char now[1<<20],*S,*T;
    if (T==S){
        T=(S=now)+fread(now,1,1<<20,stdin);
        if (T==S) return EOF;
    }
    return *S++;
}
int read(){
    int x=0,f=1;
    char ch=Getchar();
    while (ch<'0'||ch>'9'){
        if (ch=='-') f=-1;
        ch=Getchar();
    }
    while (ch<='9'&&ch>='0') x=x*10+ch-'0',ch=Getchar();
    return x*f;
}
unordered_map<int,int> num;
void ins(node x){
    for (int i=n;i>=0;i--)
        if (x.x&(1<<i)){
            if (!vis[i]){
                vis[i]=true; sz++;
                c[i]=x; break;
            } else{
                if (x>c[i]) swap(c[i],x);
                x.x^=c[i].x;
            }
        }
}
void del(int x){
    for (int i=n;i>=0;i--)
        if (vis[i]&&c[i].y==x){
            vis[i]=false; sz--;
            break;
        }
    
}
int main(){
    n=read(); Q=read(); m=(1<<n);
    for (int i=1;i<=Q;i++){
        op[i]=read(); x[i]=read();
        nxt[num[x[i]]]=i; lst[i]=num[x[i]];
        num[x[i]]=i;
    }
    for (int i=1;i<=Q;i++)
        if (!nxt[i]) nxt[i]=Q+1;
    for (int i=1;i<=Q;i++){
        if (op[i]==1) ins((node){x[i],i});
        else del(lst[i]);
        ans^=(m/(1<<sz));
    }
    printf("%d\n",ans);
    return 0;
}

```

## 5.7 AC自动机

```cpp
#include<cstdio>
#include<cstring>
using namespace std;
int T,n,q,a[510000],pos[510000],len[510000],l[510000];
int tot,sum[510000],fail[510000],ch[510000][27];
int head,tail,que[510000];
char st[510000],s[510000],pre[510000],suf[510000];
int ins(){
    int now=0,len;
    len=strlen(suf);
    for (int i=0;i<len;i++){
        if (!ch[now][suf[i]-'a']) ch[now][suf[i]-'a']=++tot;
        now=ch[now][suf[i]-'a'];
    }
    
    if (!ch[now][26]) ch[now][26]=++tot;
    now=ch[now][26];

    len=strlen(pre);
    for (int i=0;i<len;i++){
        if (!ch[now][pre[i]-'a']) ch[now][pre[i]-'a']=++tot;
        now=ch[now][pre[i]-'a'];
    }
    l[now]=strlen(suf)+strlen(pre);
    return now;
}
void build(){
    int head=1,tail=0;
    for (int i=0;i<=26;i++)
        if (ch[0][i]) que[++tail]=ch[0][i];
    while (head<=tail){
        int now=que[head++];
        for(int i=0;i<=26;i++){
            if (ch[now][i]){
                fail[ch[now][i]]=ch[fail[now]][i];
                que[++tail]=ch[now][i];
            }
            else ch[now][i]=ch[fail[now]][i];
        }
    }
}
void query(int x,int len){
    int now=0;
    for (int i=x+1;i<=x+len;i++) now=ch[now][st[i]-'a'];
    now=ch[now][26];
    for (int i=x+1;i<=x+len;i++){
        now=ch[now][st[i]-'a'];
        int tmp=now;
        while (l[tmp]>len) tmp=fail[tmp];
        sum[tmp]++;
    }
}

int main(){
    scanf("%d",&T);
    while (T--){
        memset(fail,0,sizeof(fail));
        memset(sum,0,sizeof(sum));
        memset(ch,0,sizeof(ch));
        scanf("%d%d",&n,&q);
        int cnt=0; tot=0;
        for (int i=1;i<=n;i++){
            scanf("%s",s+1); len[i]=strlen(s+1);
            a[i]=cnt;
            for (int j=1;j<=len[i];j++) st[++cnt]=s[j];
        }
        for (int i=1;i<=q;i++){
            scanf("%s%s",pre,suf);
            pos[i]=ins();
        }

        build();
        for (int i=1;i<=n;i++) query(a[i],len[i]);
        for (int i=tot;i>=1;i--) sum[fail[que[i]]]+=sum[que[i]];
        for (int i=1;i<=q;i++) printf("%d\n",sum[pos[i]]);
    }
    return 0;
}
```



# 6.分治

## 6.1 整体二分

```cpp
#include<cstdio>
#include<vector>
using namespace std;
const int MAXN=1100000;
typedef long long ll;
int n,m,x,k,t[MAXN],ans[MAXN];
ll tree[MAXN];
vector<int> vec[MAXN];
struct node{
    int l,r,num;
    bool isq;
    ll k;
} a[MAXN],Left[MAXN],Right[MAXN];
inline void add(int x,ll y){
    for (;x<=m;x+=x&-x) tree[x]+=y;
}
inline ll getsum(int x){
    ll sum=0;
    for (;x;x-=x&-x) sum+=tree[x];
    return sum;
}
void solve(int L,int R,int l,int r){
    if (l>r) return;
    if (L==R){
        for (int i=l;i<=r;i++)
            if (a[i].isq) ans[a[i].num]=L;
        return;
    }
    int mid=(L+R)>>1,n=0,m=0;
    for (int i=l;i<=r;i++){
        if (a[i].isq){
            ll tmp=0;
            for (int j=0;j<vec[a[i].num].size();j++){
                tmp+=getsum(vec[a[i].num][j]);
                if (tmp>=a[i].k) break;
            }
            if (tmp>=a[i].k) Left[++n]=a[i];
            else a[i].k-=tmp,Right[++m]=a[i];
        }
        else{
            if (a[i].num<=mid){
                if (a[i].l<=a[i].r){
                    add(a[i].l,a[i].k);
                    add(a[i].r+1,-a[i].k);
                }
                else{
                    add(1,a[i].k);
                    add(a[i].r+1,-a[i].k);
                    add(a[i].l,a[i].k);
                }
                Left[++n]=a[i];
            }
            else Right[++m]=a[i];
        }
    }
    for (int i=l;i<=r;i++)
        if (a[i].num<=mid&&!a[i].isq){
            if (a[i].l<=a[i].r){
                add(a[i].l,-a[i].k);
                add(a[i].r+1,a[i].k);
            }
            else{
                add(1,-a[i].k);
                add(a[i].r+1,a[i].k);
                add(a[i].l,-a[i].k);
            }
        }
    for (int i=1;i<=n;i++) a[l+i-1]=Left[i];
    for (int i=1;i<=m;i++) a[l+n+i-1]=Right[i];
    solve(L,mid,l,l+n-1); solve(mid+1,R,l+n,r);
}
int main(){
    scanf("%d%d",&n,&m);
    for (int i=1;i<=m;i++){
        scanf("%d",&x);
        vec[x].push_back(i);
    }
    for (int i=1;i<=n;i++) scanf("%d",&t[i]);
    scanf("%d",&k);
    for (int i=1;i<=k;i++){
        scanf("%d%d%lld",&a[i].l,&a[i].r,&a[i].k);
        a[i].isq=false; a[i].num=i;
    }
    for (int i=1;i<=n;i++){
        a[i+k].k=t[i]; a[i+k].isq=true;
        a[i+k].num=i;
    }
    solve(1,k+1,1,k+n);
    for (int i=1;i<=n;i++)
        if (ans[i]!=k+1) printf("%d\n",ans[i]);
        else puts("NIE");
    return 0;
}
```

## 6.2 CDQ分治（三维偏序）

```cpp
#include<cstdio>
#include<algorithm>
using namespace std;
const int MAXN=110000;
const int MAXM=210000;
struct node{
    int s,c,m,num,tot;
} t[MAXN],a[MAXN];
int n,k,cnt,tree[MAXM],tot[MAXN],num[MAXN],ans[MAXN];
inline void add(int x,int y){
    for (;x<=k;x+=x&(-x)) tree[x]+=y;
}
inline int getsum(int x){
    int sum=0;
    for (;x;x-=x&(-x)) sum+=tree[x];
    return sum;
}
inline int cmp1(node a,node b){
    return (a.s<b.s)||(a.s==b.s&&a.c<b.c)||(a.s==b.s&&a.c==b.c&&a.m<b.m);
}
inline int cmp2(node a,node b){
    return (a.c<b.c)||(a.c==b.c&&a.m<b.m);
}
void CDQ(int l,int r){
    if (l==r){
        tot[a[l].num]+=a[l].tot-1;
        return;
    }
    int mid=(l+r)>>1;
    CDQ(l,mid); CDQ(mid+1,r);
    sort(a+l,a+mid+1,cmp2); sort(a+mid+1,a+r+1,cmp2);
    int i,j;
    for (i=l,j=mid+1;j<=r;j++){
        for (;i<=mid&&a[i].c<=a[j].c;i++) add(a[i].m,a[i].tot);
        tot[a[j].num]+=getsum(a[j].m);
    }
    while ((--i)>=l) add(a[i].m,-a[i].tot);
}
int main(){
    scanf("%d%d",&n,&k);
    for (int i=1;i<=n;i++) scanf("%d%d%d",&t[i].s,&t[i].c,&t[i].m);
    sort(t+1,t+n+1,cmp1);
    for (int i=1;i<=n;i++)
        if (i==1||t[i-1].s!=t[i].s||t[i-1].m!=t[i].m||t[i-1].c!=t[i].c){
            a[++cnt]=t[i];
            a[cnt].tot=1;
            a[cnt].num=cnt;
        } else a[cnt].tot++;
    for (int i=1;i<=cnt;i++) num[i]=a[i].tot;
    CDQ(1,cnt);
    for (int i=1;i<=cnt;i++) ans[tot[i]]+=num[i];
    for (int i=0;i<n;i++) printf("%d\n",ans[i]);
    return 0;
}

```

# 7. 计算几何
## 7.1 基础模板
```cpp
struct Point{
    double x,y;
    Point(double x=0,double y=0):x(x),y(y){}
};
typedef Point Vector;
Vector operator+(Vector A,Vector B){ return Vector(A.x+B.x,A.y+B.y);}
Vector operator-(Point A,Point B){ return Vector(A.x-B.x,A.y-B.y);}
Vector operator*(Vector A,double p){ return Vector(A.x*p,A.y*p);}
Vector operator/(Vector A,double p){ return Vector(A.x/p,A.y/p);}
const double eps=1e-6;
int sgn(double x){
    if (fabs(x)<eps) return 0;
    if (x<0) return -1;
    return 1;
}
bool operator==(const Point& a,const Point& b){
    if (sgn(a.x-b.x)==0&&sgn(a.y-b.y)==0) return true;
    return false;
}
double Dot(Vector A,Vector B){ return A.x*B.x+A.y*B.y;}
double Length(Vector A){ return sqrt(Dot(A, A));}
double Angle(Vector A,Vector B){ return acos(Dot(A,B)/Length(A)/Length(B));}
double Cross(Vector A,Vector B){ return A.x*B.y-A.y*B.x;}
double Area(Point A,Point B,Point C){//计算两向量构成的三角形形有向面积
    return Cross(B-A,C-A)/2;
}
double Area2(Point A,Point B,Point C){//计算两向量构成的平行四边形有向面积
    return Cross(B-A,C-A);
}
Vector Rotate(Vector A,double rad){//rad为弧度 且为逆时针旋转的角
    return Vector(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
}
Vector Normal(Vector A){//向量A左转90°的单位法向量
    double L=Length(A);
    return Vector(-A.y/L,A.x/L);
}
bool ToLeftTest(Point a,Point b,Point c){ return Cross(b-a,c-b)>0;}

double DistanceToLine(Point P,Point A,Point B){//点P到直线AB距离公式
    Vector v1=B-A,v2=P-A;
    return fabs(Cross(v1,v2)/Length(v1));//不去绝对值，得到的是有向距离
}
double DistanceToSegment(Point P,Point A,Point B){//点P到线段AB距离公式
    if (A==B) return Length(P-A);
    Vector v1=B-A,v2=P-A,v3=P-B;
    if (dcmp(Dot(v1,v2))<0) return Length(v2);
    if (dcmp(Dot(v1,v3))>0) return Length(v3);
    return DistanceToLine(P,A,B);
}
```


# 8.矩阵
```cpp
## 4.6 矩阵基础模板

​```cpp
struct matrix{
    int v[55][55];
    matrix(){ memset(v,0,sizeof(v));}
    int det(){
        int res=1;
        for (int i=1;i<=K;i++){
            if (!v[i][i]){
                for (int j=i+1;j<=K;j++)
                    if (v[j][i]){ swap(v[i],v[j]); break;}
                res=dec(0,res); 
                if (!v[i][i]) return 0;
            }
            ll inv=getinv(v[i][i]);
            for (int j=i+1;j<=K;j++){
                int tmp=mul(v[j][i],inv);
                for (int k=i;k<=K;k++) v[j][k]=dec(v[j][k],mul(tmp,v[i][k]));
            }
        }
        for (int i=1;i<=K;i++) res=mul(res,v[i][i]);
        return res;
    }
    friend matrix operator+(const matrix &a,const matrix &b){
        matrix res;
        for (int i=1;i<=K;i++)
            for (int j=1;j<=K;j++)
                res.v[i][j]=add(a.v[i][j],b.v[i][j]);
        return res;
    }
    friend matrix operator-(const matrix &a,const matrix &b){
        matrix res;
        for (int i=1;i<=K;i++)
            for (int j=1;j<=K;j++)
                res.v[i][j]=dec(a.v[i][j],b.v[i][j]);
        return res;
    }
    friend matrix operator*(const matrix &a,const int b){
        matrix res;
        for (int i=1;i<=K;i++)
            for (int j=1;j<=K;j++)
                res.v[i][j]=mul(a.v[i][j],b);
        return res;
    }
    friend matrix operator*(const matrix &a,const matrix &b){
        matrix res;
        for (int i=1;i<=K;i++)
            for (int k=1;k<=K;k++)
                if (a.v[i][k])
                    for (int j=1;j<=K;j++)
                        res.v[i][j]=add(res.v[i][j],mul(a.v[i][k],b.v[k][j]));
        return res;
    }
    friend matrix operator^(matrix x,int a){
        matrix res;
        for (int i=1;i<=K;i++) res.v[i][i]=1;
        while (a){
            if (a&1) res=res*x;
            x=x*x; a>>=1;
        }
        return res;
    }
};
```

ToDo List：多项式多点插值

