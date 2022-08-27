import{_ as n}from"./plugin-vue_export-helper.21dcd24c.js";import{o as s,c as a,d as p}from"./app.e99e8fad.js";const t={},o=p(`<ul><li>\u7ED9\u5B9A\u4E00\u4E2A 01 \u4E32 $s$\uFF0C\u8981\u6C42\u4F60\u627E\u5230\u4E00\u4E2A\u4E0E $s$ \u957F\u5EA6\u76F8\u7B49\u7684 01 \u4E32 $t$ \u5E76\u6EE1\u8DB3\u4EE5\u4E0B\u6761\u4EF6\uFF1A <ol><li>\u5BF9\u4E8E\u4EFB\u610F\u7684 $l,r(1\\leq l\\leq r\\leq n)$\uFF0C$s[l:r]$ \u548C $t[l:r]$ \u7684\u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u76F8\u540C\uFF1B</li><li>$t$ \u4E2D $0$ \u7684\u6570\u91CF\u5C3D\u53EF\u80FD\u591A\u3002</li></ol></li><li>\u6709\u591A\u89E3\u8F93\u51FA\u4EFB\u610F\u4E00\u79CD\u3002$|S|\\leq 10^5$\u3002</li></ul><h3 id="sol" tabindex="-1"><a class="header-anchor" href="#sol" aria-hidden="true">#</a> sol</h3><ul><li><p>\u8FD9\u4E2A\u9898\u975E\u5E38\u6709\u610F\u601D\u3002\u5BB9\u6613\u60F3\u5230\u6309\u7167\u4E00\u4E2A\u987A\u5E8F\u4F9D\u6B21\u786E\u5B9A\u6BCF\u4E2A\u4F4D\u7F6E\u7684\u503C\u3002</p></li><li><p>\u4F46\u662F\u524D\u9762\u7684\u4F4D\u7F6E\u786E\u5B9A\u4E86\uFF0C\u5F53\u524D\u4F4D\u7F6E\u7684\u503C\u4E0D\u4EC5\u548C\u4E4B\u524D\u7684\u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u6709\u5173\uFF0C\u8FD8\u548C $1$ \u7684\u4E2A\u6570\u6709\u5173\u3002</p></li><li><p>\u53D1\u73B0 <code>10</code> \u662F\u4E0D\u80FD\u88AB\u5176\u5B83\u7B49\u957F\u7684\u4E32\u6240\u66FF\u6362\u7684\uFF0C\u5B83\u7684\u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u4E3A $1$\uFF0C\u800C\u5176\u5B83\u7B49\u957F\u7684\u4E32\u90FD\u4E3A $2$\u3002</p></li><li><p>\u6709\u4E00\u4E2A\u601D\u8DEF\uFF08\u731C\u6D4B\uFF09\uFF0C\u5C31\u662F\u628A\u4E32\u4E2D\u7684 <code>10</code> \u4E0D\u65AD\u5220\u9664\uFF0C\u76F4\u5230\u4E0D\u5B58\u5728\u3002\u4F46\u8FD9\u6837\u771F\u7684\u884C\u5417\uFF1F</p></li><li><p>\u5220\u9664\u67D0\u4E2A <code>10</code> \u540E\u5BF9\u6240\u6709 $l,r$ \u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u7684\u5F71\u54CD\uFF1A</p><ol><li>$l,r$ \u5B8C\u5168\u5305\u542B\u88AB\u5220\u9664\u7684 <code>10</code>\uFF0C\u90A3\u4E48\u539F\u5148\u65E0\u8BBA\u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u662F\u600E\u4E48\u9009\u7684\uFF0C\u90FD\u4E00\u5B9A\u4F1A\u9009\u62E9 <code>10</code> \u4E2D\u7684\u6070\u597D\u4E00\u4E2A\u6570\uFF0C\u73B0\u5728\u5B83\u7684\u957F\u5EA6\u5FC5\u5B9A\u51CF\u5C11 1\uFF1B</li><li>$l,r$ \u548C \u88AB\u5220\u9664\u7684 <code>10</code> \u65E0\u4EA4\uFF0C\u90A3\u4E48\u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u4E0D\u5F71\u54CD\uFF1B</li><li>$l,r$ \u548C \u88AB\u5220\u9664\u7684 <code>10</code> \u6709\u4EA4\u4F46\u4E0D\u5305\u542B\uFF0C\u5373 $l$ \u6216 $r$ \u5728 <code>10</code> \u4E4B\u95F4\uFF0C\u5982\u679C\u662F $l$ \u5219\u539F\u5B50\u5E8F\u5217\u4E00\u5B9A\u4F1A\u9009\u62E9\u6700\u5DE6\u8FB9\u7684 $0$\uFF0C\u5982\u679C\u662F $r$ \u5219\u539F\u5B50\u5E8F\u5217\u4E00\u5B9A\u4F1A\u9009\u62E9\u6700\u53F3\u8FB9\u7684 $1$\uFF0C\u5220\u53BB\u540E\u957F\u5EA6\u5FC5\u5B9A\u51CF\u5C11 1\u3002</li></ol></li><li><p>\u5220\u53BB\u540E\u7684\u95EE\u9898\uFF0C\u5982\u679C $s&#39;$ \u548C $t&#39;$ \u6EE1\u8DB3\u6761\u4EF6\uFF0C\u90A3\u4E48\u52A0\u56DE\u8FD9\u4E2A <code>10</code> \u540E\u4ECD\u7136\u6EE1\u8DB3\u6761\u4EF6\u3002</p></li><li><p>\u6240\u4EE5\u53EF\u4EE5\u4E00\u76F4\u5220\u53BB <code>10</code>\uFF0C\u76F4\u5230\u6CA1\u6709 <code>10</code>\u3002\u8003\u8651\u73B0\u5728\u7684\u4E32\u4E00\u5B9A\u4E3A <code>00...011...1</code>\uFF0C\u53EF\u4EE5\u628A\u6240\u6709\u6570\u5B57\u90FD\u53D8\u6210 $0$\uFF0C\u6240\u6709 $l,r$ \u6700\u957F\u4E0D\u4E0B\u964D\u5B50\u5E8F\u5217\u957F\u5EA6\u4E0D\u53D8\u3002</p></li><li><p>\u5177\u4F53\u5B9E\u73B0\u53EF\u4EE5\u4ECE\u5DE6\u5230\u53F3\u626B\u63CF\uFF0C\u7EF4\u62A4\u4E00\u4E2A\u6808\u8BB0\u5F55\u4E4B\u524D\u8FD8\u6CA1\u88AB\u5220\u53BB\u7684 $1$ \u7684\u4F4D\u7F6E\uFF0C\u5982\u679C\u5F53\u524D\u4E3A $0$ \u4E14\u6808\u4E0D\u4E3A\u7A7A\u5219\u53EF\u4EE5\u548C\u524D\u9762\u7684\u4E00\u8D77\u5220\u53BB\uFF0C\u5426\u5219\u52A0\u5165\u6808\u3002\u6700\u540E\u6808\u4E2D\u7684 $1$ \u5373\u6CA1\u6709\u88AB\u5220\u53BB\u7684 $1$\uFF0C\u90FD\u53D8\u4E3A $0$ \u5373\u53EF\u3002</p></li></ul><div class="language-cpp ext-cpp line-numbers-mode"><pre class="language-cpp"><code><span class="token macro property"><span class="token directive-hash">#</span><span class="token directive keyword">include</span><span class="token string">&lt;bits/stdc++.h&gt;</span></span>
<span class="token keyword">using</span> <span class="token keyword">namespace</span> std<span class="token punctuation">;</span>
string s<span class="token punctuation">;</span> <span class="token keyword">int</span> top<span class="token punctuation">,</span>t<span class="token punctuation">[</span><span class="token number">110000</span><span class="token punctuation">]</span><span class="token punctuation">;</span>
<span class="token keyword">int</span> <span class="token function">main</span><span class="token punctuation">(</span><span class="token punctuation">)</span><span class="token punctuation">{</span>
	cin<span class="token operator">&gt;&gt;</span>s<span class="token punctuation">;</span> <span class="token keyword">int</span> n<span class="token operator">=</span>s<span class="token punctuation">.</span><span class="token function">length</span><span class="token punctuation">(</span><span class="token punctuation">)</span><span class="token punctuation">;</span>
	<span class="token keyword">for</span> <span class="token punctuation">(</span><span class="token keyword">int</span> i<span class="token operator">=</span><span class="token number">0</span><span class="token punctuation">;</span>i<span class="token operator">&lt;</span>n<span class="token punctuation">;</span>i<span class="token operator">++</span><span class="token punctuation">)</span>
		<span class="token keyword">if</span> <span class="token punctuation">(</span>s<span class="token punctuation">[</span>i<span class="token punctuation">]</span><span class="token operator">==</span><span class="token char">&#39;1&#39;</span><span class="token punctuation">)</span> t<span class="token punctuation">[</span><span class="token operator">++</span>top<span class="token punctuation">]</span><span class="token operator">=</span>i<span class="token punctuation">;</span>
		<span class="token keyword">else</span> <span class="token keyword">if</span> <span class="token punctuation">(</span>top<span class="token punctuation">)</span> top<span class="token operator">--</span><span class="token punctuation">;</span>
	<span class="token keyword">while</span> <span class="token punctuation">(</span>top<span class="token punctuation">)</span> s<span class="token punctuation">[</span>t<span class="token punctuation">[</span>top<span class="token operator">--</span><span class="token punctuation">]</span><span class="token punctuation">]</span><span class="token operator">=</span><span class="token char">&#39;0&#39;</span><span class="token punctuation">;</span>
	cout<span class="token operator">&lt;&lt;</span>s<span class="token operator">&lt;&lt;</span>endl<span class="token punctuation">;</span>
	<span class="token keyword">return</span> <span class="token number">0</span><span class="token punctuation">;</span>
<span class="token punctuation">}</span>
</code></pre><div class="line-numbers" aria-hidden="true"><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div><div class="line-number"></div></div></div>`,4),e=[o];function c(l,i){return s(),a("div",null,e)}var d=n(t,[["render",c],["__file","CF1204D.html.vue"]]);export{d as default};
