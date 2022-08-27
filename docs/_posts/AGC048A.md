---
title: AGC048A
date: 2021-10-01 15:33:10
tags: []
categories: [AtCoder]
mathjax: true
---

- 给定一个小写字母组成的字符串 $S$，每次操作可以交换相邻两个字符。
- 问最少需要多少次操作，使得 S 的字典序严格大于 `atcoder`，无解输出 `-1`。
- 有 $T$ 组数据，$1\leq T\leq 100,1\leq |S|\leq 1000$。
<!-- more -->

### sol

![sol](https://cdn.jsdelivr.net/gh/hydd-new/image-hosting@master/agc048a.png)

```cpp
#include<bits/stdc++.h>
using namespace std;
const char t[8]={'a','t','c','o','d','e','r',0};
int T,n; char s[1100];
int main(){
	scanf("%d",&T);
	while (T--){
		for (int i=0;i<8;i++) s[i]=0;
		scanf("%s",s); n=strlen(s);
		bool flag=false;
		for (int i=0;i<8;i++)
			if (t[i]!=s[i]){
				flag=(t[i]<s[i]);
				break;
			}
		if (flag){ puts("0"); continue;}
		int pos=-1;
		for (int i=0;i<n;i++)
			if (s[i]!='a'){
				pos=i;
				break;
			}
		if (pos==-1){ puts("-1"); continue;}
		printf("%d\n",pos-1+(s[pos]<='t'));
	}
	return 0;
}
```