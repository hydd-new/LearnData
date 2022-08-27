---
title: 'trees'
date: 2022-01-17 14:12:53
tags: [trees]
categories: [CodeForces2100-2300]
mathjax: true
---


#### 600E
- 线段树合并模板，直接维护max和sum。
- 不过应该也可以用 dsu on tree 做，把小的子树先dfs，大子树dfs后不清空，再加入小子树。


#### 208E
- 可能不连通，不妨加入 $0$ 号点变成一棵树。
- 相当于求 $p$ 级祖先有多少个 $p$ 的儿子。
- 先离线求出 $p$ 级祖先是谁，把询问挂到 $p$ 级祖先上。
- 再把某个点上有的询问，记录进入子树前的值，出来再减去。


#### 620E
- 变成 dfs 序，区间修改，区间询问数量。
- 线段树每个区间维护某种颜色的数量即可。


#### 1555E
- **求 $\max-\min$ 可以枚举其中一个**。
- 假设我枚举 $\min$，就要求最小的 $\max$ 使得只用 $[\min,\max]$ 中的区间可以覆盖每一个位置。
- 线段树维护区间最小值，区间+1/-1，当最小值非 $0$ 的时候说明所有都被覆盖了。枚举 $\min$，双指针枚举 $\max$ 即可。


#### 1009F
- 暴力的想法是 dp 来转移子树距离为 $x$ 的数量。
- 长链剖分，维护 $d$ 最大时最小的深度。先递归轻儿子，清空，再递归重儿子，不清空，再把轻儿子加入。


#### 587C
- 倍增维护每个点往上的边前 $A=\max\{a\}$ 小的。
- 询问的时候也暴力合并这些值，复杂度为 $O(n\log n\times A)$


#### 282E
- 就算前后缀有交问题也不大，可以调整成没交的。
- 所以把所有后缀异或和加入 trie，每个前缀在 trie 上找异或它最大的值即可。


#### 540E
- 求出交换过的关键位置现在的值。
- 非关键位置之间一定没有，关键位置之间可以用树状数组。
- 非关键位置和关键位置之间，考虑 $a_x=y$，那么 $y$ 到 $a_x$ 之间的数都有贡献，由于两个一定都是关键点，所以可以轻易求出中间有多少关键点，就能计算了。


#### 1181D
- 如果初始举办次数都一样，那么会按照 $1,2,\cdots ,m$ 的顺序举办。
- 现在相当于是，之前某个每举办一次比赛，就删去它第一次出现的位置。
- 那么你可以求出每次删去的是哪个位置，二分求出现在应该到哪个位置，mod 一下就知道在哪里举办了。


#### 1380E
- 有一种简单的构造就是 $\{1\}\rightarrow 2,\{1,2\}\rightarrow 3,\{1,2,3\}\rightarrow 4,\cdots$，这样的方案是 $i$ 和 $i+1$ 不在同一堆的个数之和。
- 而且，$i$ 和 $i+1$ 不在同一堆，一定需要有一次移动，且移动不会同时满足其它的 $i$。
- 用并查集+启发式来优化合并，把小的集合的贡献和影响重新计算（前一个/后一个在那个塔里，减少答案，vector 维护集合里的点）。


#### 960D
- 发现无论怎么移动，深度为 $i$ 的（从 $0$ 开始）从左到右一定是 $2^i+x,2^i+x+1,\cdots,2^{i+1}-1,2^i,2^i+1,\cdots,2^i+x-1$。
- 1 操作就是把一行每个数往右移 $k$，也就是 $x$ 减少 $k$。
- 2 操作首先把这行每个数往右移 $k$，对于下一行它应该右移 $2k$，再下一行应该右移 $4k$，以此类推（只有 $60$ 层，可以暴力）。
- 3 操作首先要找到它在那一行的位置，然后不停 $/2$，找到上一行对应位置，输出它的数。


#### 817F
- 如果只有 1,2 操作可以直接set维护线段，因为复杂度是可以均摊的。
- 不过有 3 操作不行，但是发现还是可以离散化，离散化后就变成区间赋 0/1，区间取反，求全局第一个 0。
- 线段树维护区间翻和不翻的第一个 0 和标记即可。


#### 1615D
- 这什么翻译，还以为是点带权，原来是边带权。
- 那就好做了，先变成两个点到根路径上的异或和，然后就要求两个点值相同或不同，不过不需要 2-SAT，这个是双向的，直接并查集合并（相同合并 (x,y),(x',y')，不同合并 (x,y'),(x',y)）。


#### 165D
- 相当于询问简单路径上有没有白边，有就是 -1，否则就是两点简单路径长度。
- 树状数组维护每个点到根路径上白边数量，也就是改成白边的时候给子树 +1，改成黑边给子树 -1。
- 询问就找到 LCA，query(x)+query(y)-2query(lca)看白边数量，如果非 0 就输出 -1，否则输出 dep[x]+dep[y]-2dep[lca]。


#### 893F
- 一眼可持久化线段树合并。
- 每个点维护不超过某个的绝对深度的点权最小值（维护相对深度到上一层合并的时候要移位，不好做）。
- 询问的时候就在那个点的线段树里询问区间最小值即可。


#### 406D
- 所谓最右边的一个能到达的山坡就是凸包上下一个点（在它之后都会被它拦住）。
- 每个点连向右边的这个点，然后能共同到达的就是求 LCA。


#### 894D
- 每个子树求出到它距离不超过 $x$ 的数量。
- 然后枚举LCA，查询子树里贡献即可。


#### 690C3
- 根据经典套路，新直径端点一定是原来直径的两个点和新的点中选两个点。
- 把新点加入的时候更新倍增数组，就可以求 LCA 和距离了。


#### 191E
- 先二分答案，然后判断有多少个区间和 $\geq mid$。
- 具体的方法是固定一个端点，要求前缀和不超过某个值的个数，可以离散化后树状数组求得。复杂度是 $\log^2$ 的。


#### 1575I
- 观察式子，发现两个数符号没有关系，$\max(|a_i-a_j|,|a_i+a_j|)=|a_i|+|a_j|$。
- 那么 $i,j$ 间除了 $i,j$ 其它点都贡献两次，问题就变成路径上的权值之和，维护每个点到根路径和，即修改时给子树加/减（树状数组），询问求 LCA 即可。


#### 85C
- 只能错恰好一次，你错了之后，如果是 0 走成了 1，那么之后都会走 0；反之是 1 走成了 0，那么之后都会走 1。
- 怎么找一步不错的叶子节点？在原序列upper_bound一下，看这个点的位置，如果是叶子节点说明就是它，如果不是，这一步会左走，而之后就没有比它大的点了，一直往右，记录一下就知道是哪个了。
- 而确定了目标点后，往上的答案都可以处理出来（还是要记录一直往右和往左走达到的点），从父亲转移过来。


#### 1252B
- 容易想到设 $dp_{u,0/1}$ 表示 $u$ 的子树覆盖完，有没有路径以 $u$ 为端点的方案数。
- 但是发现这不行，原因是你有可能子树里原先不合法，现在合法，因为根节点的路径可以接到父亲。
- 改为设 $dp_{u,0/1/2}$ 表示 $u$ 的子树覆盖完，没有以 $u$ 为端点/有以 $u$ 为端点且不必须和父亲连（子树本身合法）/有以 $u$ 为端点且必须和父亲连（子树在根节点处不合法）的方案数。
- $dp_{u,0}$ 要从子树选择两条连上来，那么选择两个 $x,y$，$dp_{x,1/2},dp_{y,1/2}$（能连上来），剩下的就只能是 $dp_{v,0/1}$（能不连上来），乘起来累加。
- $dp_{u,1}$ 要从子树选择一条连上来且子树没有其它能连上来的，那么选择一个 $x$，$dp_{x,1/2}$（能连上来），剩下的只能是 $dp_{v,0}$（不能连上来），乘起来累加。
- $dp_{u,2}$ 要从子树选择一条连上来且子树没有其它必须连上来的，那么选择一个 $x$，$dp_{x,1/2}$（能连上来），剩下的只能是 $dp_{v,0/1}$（能不连上来），乘起来累加（计算的并不是必须和父亲连的，而是总方案，即有端点在根上的），最后减去 $dp_{u,1}$。


- 但是 $dp_{u,0}$ 需要枚举两个变量，而其它只需要枚举一个变量，逆元算一下就行。其实也可以做，枚举两个中靠后的一个，维护之前选了另一个的答案，把当前的能不连上来的除掉，乘上能连上来的即可。


#### 916D
- 可以撤销撤销操作看上去很吓人，不过仔细想一想，撤销了之前的 $k$ 次操作，相当于回到 $k$ 步之前操作后，和前一步的状态并没有关系，到达当前状态的操作已经被撤销了。
- 所以应该用可持久化数据结构，可以考虑两个线段树，一个权值线段树维护 $a_i$ 的优先级，另一个权值线段树维护优先级小于 $x$ 的数的个数。一个是单点修改单点查询，一个是区间修改单点查询，都可以用可持久化。


#### 1403B
- 容易想到，对于一个子树，你留下来给父亲匹配的一定不超过两个点，否则你选两个点自行匹配，剩下的留给父亲也可行。留下两个时，到父亲的边会计算两次（你往上保留的奇偶性一定不变，在之后和别的合并一定不如现在就合并）。
- 发现留下的只和子树里叶子节点的个数的奇偶性有关，也就是子树叶子节点个数是偶数的需要多增加 $1$（边多走了一次）。
- 树链剖分可以维护，然后再暴力撤销。

