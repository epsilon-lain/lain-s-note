"""Draw the hand-sketch model. Only drawing uses floats; bounds use rational arithmetic."""
from fractions import Fraction
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.gridspec import GridSpec

folder = Path(__file__).resolve().parent
cjk = "Noto Sans CJK JP"
if any(f.name == "Noto Sans CJK SC" for f in font_manager.fontManager.ttflist):
    cjk = "Noto Sans CJK SC"
plt.rcParams.update({"font.family": [cjk, "DejaVu Sans"], "font.size": 10.8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor":"#bfc8cd", "axes.labelcolor":"#324755",
                     "xtick.color":"#435764","ytick.color":"#435764",
                     "axes.unicode_minus":False})
navy, teal, amber, soft = "#1c3545", "#008279", "#ad5827", "#687c84"
fig = plt.figure(figsize=(10.7, 5.25), layout="constrained", facecolor="white")
gs = GridSpec(1, 2, figure=fig, width_ratios=[0.92, 1.16])
left, right = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
alpha = 1/(2+1/(3+(math.sqrt(5)-1)/2))
q, p = [1,2,7,9,16], [0,1,3,4,7]
h = [a-alpha*b for a,b in zip(p,q)]
left.fill([h[0],h[1],h[2]],[q[0],q[1],q[2]],color="#e8efef")
marks = [1/(abs(h[i])+abs(h[i+1])) for i in range(3)]
left.fill([0,h[1],0],[marks[0],q[1],marks[1]],color="#f3c69c")
left.axvline(0,color=soft,ls="--",lw=1,zorder=1)
for parity in (0,1):
    left.plot(h[parity::2],q[parity::2],color=teal,lw=2.1,zorder=2)
left.plot([0,h[1]],[0,q[1]],color=amber,lw=1.6,ls=":",zorder=3)
left.scatter([h[0]+j*h[1] for j in (1,2)], [q[0]+j*q[1] for j in (1,2)],
             s=28,marker="D",color=teal,edgecolor="white",lw=.6,zorder=6)
left.plot(h,q,color=navy,lw=1.95,zorder=3)
left.scatter(h,q,s=23,color=navy,zorder=5)
left.scatter([0]*3,marks,s=38,edgecolor=amber,facecolor="white",lw=1.7,zorder=6)
left.scatter([0],[0],color=amber,s=23,zorder=6)
for i,j in enumerate((0,1,2)):
    left.annotate(f"$P_{j}$",(h[i],q[i]),xytext=((-2,8) if i==0 else ((8,1) if i==1 else (-23,3))),textcoords="offset points",color=navy)
left.annotate("$O$",(0,0),xytext=(-17,-2),textcoords="offset points",color=amber)
left.text(-.433,8.8,"绿线连接同侧顶点\n空心点是中线穿越点",color=teal,fontsize=10)
left.set(xlim=(-.49,.19),ylim=(-.35,9.8),xlabel="横向偏移  h",ylabel="沿中线的位置  q",
         xticks=[-.4,-.2,0],yticks=[0,2,4,6,8])
left.set_title("同一次折返，量两种面积",loc="left",color=navy,fontweight="bold",pad=12)
left.text(.03,-.215,"整个折角三角形面积为 3/2，橙色为波与中线围出的面积。",transform=left.transAxes,fontsize=9,color=soft)

def digit(j):
    if j >= 4 and j & (j-1) == 0:
        return j
    return 1

def ratio(n):
    lo, hi = Fraction(0), Fraction(1)
    for j in range(115, n+1, -1):
        lo, hi = 1/(digit(j)+hi), 1/(digit(j)+lo)
    a = digit(n+1)
    lower, upper = 1+a/(1+hi), 1+a/(1+lo)
    assert upper-lower < Fraction(1,10**18)
    return float((lower+upper)/2)

ns = list(range(1,34))
rs = [ratio(n) for n in ns]
right.axhline((1+math.sqrt(5))/2,ls="--",color=teal,lw=2,label="全部部分商为 1")
right.plot(ns,rs,".-",color=amber,label="在第 4、8、16、32、… 项出现大部分商",lw=1.45,ms=5.1)
right.set_title("同样无限次返回，间距比可以无界",loc="left",color=navy,fontweight="bold",pad=12)
right.set(xlim=(0,34),ylim=(1.25,24),xlabel="折返序号  n",ylabel="相邻穿越位置之比",xticks=[1,7,15,23,31])
right.set_yscale("log",base=2)
right.set_yticks([1.5,2,4,8,16])
right.set_yticklabels(["1.5","2","4","8","16"])
right.grid(axis="y",color="#e4e9ea",lw=.8)
right.legend(loc="upper left",frameon=False,fontsize=9.1)
right.text(.02,-.215,"位置从固定原点测量。图中只画前 33 段，无界性由生成规则证明。",transform=right.transAxes,fontsize=9,color=soft)
fig.savefig(folder/"crossing-wave.png",dpi=170,bbox_inches="tight",pad_inches=.18)
