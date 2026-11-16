# 原生数学正文与算法来源

本目录保留研究产物本身；不发布聊天trace、运行控制脚本、凭证、主机路径或私有节点/任务标识。原始 `.tex` 包保持未修改，公开文件与原来源的哈希映射在上一级 [PROVENANCE.json](../PROVENANCE.json)。

## 数学正文

| 公开文件 | 原生条目 | 提取范围 |
|---|---|---|
| [interval-criterion.md](interval-criterion.md) | `M-05ac7d1d25ae1353@1` | 完整数学正文 |
| [near-double-family.md](near-double-family.md) | `M-eea68ed0edda0837@1` | 完整数学正文 |
| [ac-resource-and-macro.md](ac-resource-and-macro.md) | `M-0d87e4c929e155d1@1` | 证据/provenance小节以前的原样数学段落 |
| [resource-signature.md](resource-signature.md) | `M-5606eba96406db45@1` | 有限运行回执与可执行代码以前的原样数学段落 |

这些文件是原有数学条目正文或段落的直接摘取；未将本次复核写回原数学叙述。英文摘录保留原有公式文本格式，中文主文另用GitHub数学格式排版。历史正文中的“尚无独立审阅”是发表时的状态，不覆盖后来实际发生的同伴审阅。

## 真正的原生研究算法

以下五份Python文件从权威 `.tex` 的raw正文提取，**源码未经修改**，包括原有结尾换行。提取后的Python SHA256与原始生成程序哈希对应；准确数值见provenance。它们保留原算法身份，没有以本次外部复核代码替换原研究代码。

| 公开文件 | 身份与作用 |
|---|---|
| [native-theory-v2.py](native-theory-v2.py) | 理论作者原生v2：分拆谱、直接图检查、两类族、分量界；`M-76595f02ff6f5258@2` |
| [native-application.py](native-application.py) | 应用作者原生：完整B4球、全部8路径、资源与单宏；`M-7984e3d81670a37f@1` |
| [native-independent-theory-verifier.py](native-independent-theory-verifier.py) | 原生不同作者的逐切分及直接异型图核验；`M-8551c6027291f241@1` |
| [native-independent-ac-verifier.py](native-independent-ac-verifier.py) | 原生不同作者的字符串/半路径树/宏DP核验；`M-809143f7c563abe4@1` |
| [native-resource-signature-audit.py](native-resource-signature-audit.py) | 原生一般资源签名有限审计；`M-5606eba96406db45@1` |

这些程序由原生隔离研究任务执行并产生既有产物；本次发布整理只做文本提取、隐私检查及语法解析，**没有执行它们**。源码中的 `/output` 是原隔离任务的虚拟输出接口，不是私人主机路径。历史相对源文件名只是种子引用；算法不据此读取这些文件。独立运行原程序需要提供其原输出接口，或自行明确修改输出路径；修改后的程序当然会有不同的源码哈希。

如只需在普通目录重现本次轻量数学复核，使用上一级 [verify_native_results.py](../verify_native_results.py)。它是后来另写、可直接运行的程序，身份与上述五份原生算法明确分开。
