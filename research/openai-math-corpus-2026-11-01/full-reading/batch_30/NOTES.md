# Queue 53 — 完整正文首读

- 稿件：Abelian covers, Gale correspondences, and the Hodge conjecture for powers（October 6, 2026；family 032）。
- 固定来源：openai/math@4eaee87548250b2682affda7d05d947c05554bae；directory tree 10b943586f6ea65e2eab389f1587fa844c0cbd5c（未截断）。
- README 将实际 PDF 定位为本目录 paper.pdf；单一实际正文为 build/paper.tex。完整读取 README 1–21、paper.tex 1–3779，合计3800行、178723字节。TeX 按300行以下连续块显示，全部输出无截断。
- 两个缓存按原字节（含末尾换行）计算 Git blob SHA，分别等于 d9b256307bca27bc8c1391a5d06b8d070369c0c5、0088171f5929a6438c6c0f6de82de781fb929b5d。
- 全目录仅这两个文本 blob；无递归 input/include、附录、外部 bibliography、figure/caption 或 includegraphics。交换图与表格全在正文。3627–3778行为全部22条 inline bibliography。无未读局部数学源。成品 PDF、README 所链前版和被引文献正文未读取，不计入局部源闭包。

原稿声称在每个连通 topological abelian curve cover 的全标记 tensor Hodge-generic locus，所有 cover-Jacobian 与任意 CM companion 的混合幂满足有理 Hodge 猜想；并通过 perfect Gale tensor 及 whole-degree rational correspondence supply，将结论转移到至多两方程的 full labelled smooth diagonal complete intersections 的全部幂。CM Hodge theorem、saturated tensor criterion、correspondence transfer 等 companion 输入只按本篇准确接口读取，没有独立核验被引证明。

本次保留三条条件性构造候选（完整假设、机制、桥梁和缺步见 records.json）：

1. 互补外幂 coproduct 的基项一一配对，作为保留符号的满秩证书；锚点 thm:gale（2817）、Lagrange系数（2859）、exact exterior coproduct（2903）。须先找到实际 Q(Z) 行/列和共同基的同一 exterior model，并承担插值分母及高度。
2. 精确 support quotient 加字符级 perfect contraction，借有理 whole-degree maps 推出零 cokernel；锚点 prop:characters（2748）、prop:supply（2964）、rational source degree（3014）。须有实际 finite symmetry/quotient，不能对任意 polynomial matrix 用连通性推常秩。
3. fixed-hyperplane 零交与 rank-one 方向非正交图给不可约性；共享 ideal 另用整个作用识别，避免误算块独立性；锚点 p4:twist-model（1448）、propagation（1518）、full Γ links（2245）、supporting ideals（2456）。须在实际算术行空间上构造这些精确保留结构的 operators。

阅读不等于独立有效性。三条均不提供 π 的 BA/nonBA 结论；同一整数行列式的非零与解析小量、实际 Q(Z) row rank、共同基、列高及二次占用 cΣn_a² 的联合成本仍待证明。本批只写 source cache、records.json、NOTES.md，未修改 ledger 或发布。

