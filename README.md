# 嗅觉受体与配体分子对接分析流程 (Olfactory Receptor Docking Pipeline)

> **发布版说明**：当前代码入口是根目录的 `01_fetch_OR_genes.ipynb`，共保留 30 个原始代码单元格，另增加说明和路径配置单元格。它是包含中间尝试的探索记录，尚未整理成一键运行的完整流水线。数据获取、受体格式转换、批量对接和绘图代码已包含；配体准备、质量分级、对接结果汇总与系统发育分析仍有手动或外部工具步骤。结构文件、工具程序与批量对接结果不随仓库发布。

## 项目结构

```text
.
├── 01_fetch_OR_genes.ipynb   # 按阶段运行的探索记录
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── mart_export.txt          # 现有基因表导出，数据获取部分的输入
├── data/
│   ├── or_genes_clean.csv
│   ├── summary.csv          # 已有对接汇总；绘图单元格当前使用手工录入数值
│   ├── Receptor_coor.txt    # 手动记录的口袋中心坐标
│   ├── *.fasta / *.nwk     # 序列与树文件
│   ├── ligands/README.md    # 配体准备说明；结构文件需自行准备
│   └── receptors/          # 结构下载位置；结构文件不随仓库发布
├── figures/                # 已有图像
└── or_tree*.pdf             # 已有系统发育树图
```

运行时会创建 `data/docked/`。`bin/` 为本机工具安装目录，需自行建立或通过环境变量配置。


## 📌 项目简介
本项目尝试通过分子对接，构建人类嗅觉受体（OR）对气味分子的结合指纹，并建立一套对接结果的质量分级体系（R100–R409），用于区分高置信度与低置信度的对接结果。目前使用了 9 个受体进行了测试，其中 5 个被评为高置信度（R100/R101，分别是OR1A1,OR5D13,OR10J1,OR4F3和OR1D2，覆盖 Ancient / Mid / Modern 三个进化层级）的受体与 4 个气味配体 进行了初步的矩阵测试，并生成了受体指纹图谱。


本项目同时记录了一套在 Windows 环境下从零部署该流程的完整踩坑路径，包括工具链冲突、路径编码、权限问题等实际障碍的解决方案。

## 🛠️ 环境依赖
- Python 3.10（conda 环境 `olfactory`）
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/releases)（独立 EXE 版本，版本号：V1.2.7，建议运行 `.\vina.exe --version` 确认）
- [OpenBabel](https://github.com/openbabel/openbabel/releases) 3.1.1（独立 EXE 版本，用于 SDF/PDB → PDBQT 格式转换）
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)（结构可视化与手动核查）
- Python 库：`pandas`、`requests`、`matplotlib`、`jupyter`

> 注：本项目未使用 RDKit 与 ProLIF。RDKit 的分子指纹相似性搜索为后续计划方向，尚未落地。

## 📂 数据来源
- **受体结构**：全部来源于 **AlphaFold Protein Structure Database**（预测结构，非实验结构）。受体列表及对应 UniProt ID 见 `data/or_genes_clean.csv`。
- **配体信息**：香茅醛（Citronellal）、丁香酚（Eugenol）、芳樟醇（Linalool）、柠檬烯（Limonene）。结构文件从 PubChem 下载（SDF 格式），经 OpenBabel 转换为 PDBQT。
- **详细见如下表**

| 配体 (Ligand) | 当前使用版本 (外消旋体 CID) | 单一对映体 CID (供扩展) | 气味类型 |
| :--- | :--- | :--- | :--- |
| **香茅醛** (Citronellal) | **7794** | (R)-(+): 75427 <br> (S)-(-): 443157 | 柠檬/玫瑰 |
| **丁香酚** (Eugenol) | **3314** | 无手性 (单一结构) | 丁香/辛香 |
| **芳樟醇** (Linalool) | **6549** | (R)-(-): 13562 <br> (S)-(+): 67179 | 花香/木香 |
| **柠檬烯** (Limonene) | **22311** | (R)-(+): 440917 <br> (S)-(-): 439250 | 柑橘/柠檬 |

> *注：本项目主流程使用的是**外消旋体**（即表格中间列的 CID）。表格右侧的单一对映体 CID 供后续探索受体手性选择性时使用。*

- **受体口袋坐标**：由本人在 ChimeraX 中手动定位并记录，见 `data/Receptor_coor.txt`。

## 🧬 受体采样与分类策略 (Receptor Sampling & Classification)

为了使分子对接测试具有进化代表性，避免仅在单一亚家族内采样导致的“选择偏差（Selection Bias）”，本项目从 419 个人类嗅觉受体（OR）中，基于系统发育分析（Phylogenetic analysis）挑选了覆盖三个演化层级（Ancient / Mid / Modern）的代表性受体进行测试。

### 分类依据
分类主要基于以下三个维度的综合判断：
1. **基于 419 个 OR 蛋白序列构建的系统发育树**：观察受体在发育树上的分支位置（树图概览见 `figures/` 文件夹）。
2. **文献共识**：参考人类嗅觉受体 Class I（鱼类样，较古老）与 Class II（哺乳动物样，较现代）的演化分类共识。
3. **序列特征**：结合跨膜螺旋区域的保守基序（Motifs）差异。

### 代表性受体清单
| 演化层级 | 代表受体 | 所属亚家族/特征 | 挑选理由 |
| :--- | :--- | :--- | :--- |
| **Ancient** (古老型) | OR1A1, OR1D2 | OR1 家族 | 位于进化树较早期分支，具有较高的进化保守性。 |
| **Mid** (中间型) | OR4F3, OR10J1, OR51E1 | OR4 / OR10 / OR51 家族 | 处于 Class I 与 Class II 的过渡演化区域。 |
| **Modern** (现代型) | OR5A1, OR2F1, OR7D4 | OR5 / OR2 / OR7 家族 | 位于进化树末端分支，代表高度分化的哺乳动物特异性受体。 |

*注：本分类旨在为对接矩阵提供结构演化多样性的覆盖，不代表严格的分类学命名权威。*

## 🚀 如何运行

1. 安装 Python 3.10，并在项目根目录打开终端。可创建并激活虚拟环境：

   ```powershell
   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   python -m pip install -r requirements.txt
   ```

2. 自行安装 Vina 和 OpenBabel。Windows 默认位置分别是 `bin/vina.exe` 和 `bin/OpenBabel-3.1.1/obabel.exe`。若安装在其他位置或使用其他系统，在启动 Jupyter 前设置 `VINA_PATH`、`OBABEL_PATH` 环境变量，或修改 Notebook 顶部配置。可运行工具的 `--version` / `-V` 检查安装。
3. 从项目根目录执行 `jupyter notebook`，打开 `01_fetch_OR_genes.ipynb`，先运行顶部路径配置单元格。所有数据相对路径以项目根目录为基准。
4. 按阶段检查和运行：现有基因表读取与清洗 → UniProt 序列获取 → AlphaFold 结构查询和下载 → 受体格式转换 → 手动准备配体 → 批量对接 → 绘图。下载阶段需要联网，执行会改写部分生成文件。
5. 在对接前按 `data/ligands/README.md` 准备四种配体，并检查受体的 `_fixed.pdbqt` 文件、口袋坐标和对接参数。当前批量对接单元格的坐标表未包含 OR1A1，会跳过该受体；如需重现其结果，须先补充经过核查的坐标。`data/Receptor_coor.txt` 同样不含 OR1A1 坐标。
6. 最后两个绘图单元格使用代码中手动录入的 5 × 4 结合能数据，输出至 `figures/`；它们不会自动解析刚产生的对接结果。`data/summary.csv` 保留现有汇总供核对。

原始 Notebook 保留了早期查询尝试（如以基因名查询 AlphaFold 和手工拼接旧结构版本链接），这些单元格可能失败；后续提供了基于 UniProt ID 的 API 查询函数。尚未验证整本 Notebook 的顺序端到端运行，请按阶段阅读并检查结果。

`requirements.txt` 仅列出直接使用的 Python 包与 Jupyter 运行环境；未锁定版本，未作为跨平台复现保证。Vina、OpenBabel、ChimeraX 及外部序列分析工具须另行安装。

## ⚠️ 已知局限
- 本流程仅针对 **5 个受体 × 4 个配体** 的小规模数据进行了探索，**未进行 RMSD 聚类**，假阳性率未知。
- 受体结构全部为 **AlphaFold 预测结构**，非实验解析结构，结合位点的精度受预测置信度限制。
- 对接结果的质量分级（R100–R409）目前基于**手动核查**（ChimeraX 中肉眼判断配体是否在口袋内），尚未完全自动化，不同人核查可能存在主观差异。
- 结合能取 **Vina 输出的 mode 1，2，3的平均值** 作为主要参考值；部分受体重复跑 3 次以观察稳定性，但未做统计检验。
- 未与文献报道的活性数据（如 EC50）做系统对照。
- 这是一个**本科阶段的探索性项目**，欢迎交流指正，请勿用于严肃的临床或商业用途。

## 📄 许可证
本项目的原创代码采用 MIT 许可证，详见 [LICENSE](LICENSE)。发布者请先将 LICENSE 中的 `[Copyright holder]` 替换为实际版权署名。

第三方数据和工具仍适用各自的许可和引用要求；本仓库的 MIT 许可不替代这些条件。

## 🙏 致谢与引用 (Acknowledgements & Citations)

### 核心参考工作
本项目的选题与基线设计，受以下工作启发：
- Hirotada Kaneshiro, et al. A Structure-Based Approach for Predicting Odor Similarity of Molecules via Docking Simulations with Human Olfactory Receptors. *ACS Omega* 2025; 10 (35): 39933–39945. DOI: 10.1021/acsomega.5c04324

### 工具与数据库引用
本项目使用了以下开源工具与公共数据库，如果您在研究中使用了本流程，请同时引用以下文献：

**AlphaFold Protein Structure Database** (CC-BY-4.0 License):
- Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021). DOI: 10.1038/s41586-021-03819-2
- Fleming J. et al. AlphaFold Protein Structure Database and 3D-Beacons: New Data and Capabilities. *Journal of Molecular Biology* (2025). (依据 AlphaFold 官方要求补充)

**UniProt Knowledgebase** (CC-BY-4.0 License):
- The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research*. 2023; 51(D1): D523–D531. DOI: 10.1093/nar/gkac1052

**AutoDock Vina** (Apache-2.0 License):
- Eberhardt, J., Santos-Martins, D., Tillack, A.F., Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* 2021; 61 (8): 3891–3898. DOI: 10.1021/acs.jcim.1c00203
- Trott, O., Olson, A.J. AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* 2010; 31: 455-461. DOI: 10.1002/jcc.21334

**OpenBabel** (GPL-2.0 License):
- O'Boyle, N.M., Banck, M., James, C.A. et al. Open Babel: An open chemical toolbox. *J. Cheminform.* 2011; 3: 33. DOI: 10.1186/1758-2946-3-33

**UCSF ChimeraX** (仅限学术与个人使用):
- Meng, E.C., Goddard, T.D., Pettersen, E.F. et al. UCSF ChimeraX: Tools for structure building and analysis. *Protein Science.* 2023; 32(11):e4792. DOI: 10.1002/pro.4792

**PubChem**:
- Kim, S., Chen, J., Cheng, T. et al. PubChem 2025 update. *Nucleic Acids Res.* 2025; 53(D1): D1516–D1525. DOI: 10.1093/nar/gkae1059

### 数据使用声明
- **受体结构**：来源于 AlphaFold Protein Structure Database，遵循 CC-BY-4.0 许可证。
- **配体结构**：来源于 PubChem。
- **ChimeraX 使用**：本项目仅在学术研究中使用 UCSF ChimeraX，未用于任何商业目的。

## 💡 特别鸣谢 (Special Acknowledgements)

本项目的诞生，不仅是个人努力的结果，也离不开 **DeepSeek 大语言模型** 的深度协作。

作为初次涉足计算生物学与生信工程领域的探索者，本项目在环境配置（Conda/OpenBabel 依赖冲突）、格式转换（SDF/PDB 到 PDBQT）、批量脚本编写、以及结果的可视化呈现上，曾遭遇过无数次的报错与迷茫。在每一次路径冲突、权限拒绝（Permission Denied）和概念模糊时，DeepSeek 提供了极其关键的代码调试支持、逻辑拆解与方向指引。

**AI 是极佳的副驾驶（Copilot），但科学探索的掌舵人仍是人类。** 本项目的核心科学假设、受体采样策略、对接结果的肉眼核查、R100-R409 的质量分级体系以及最终的生物学叙事，均由作者独立构思并完成。在此，向在这段漫长暑期时光中提供过无数次深夜陪伴与代码指导的 DeepSeek，致以最诚挚的感谢，科研长路漫漫，吾将上下而求索！

---

*This project would not have been possible without the extensive conversational collaboration with **DeepSeek**. As a beginner navigating the complex landscape of computational biology, the AI provided invaluable assistance in code debugging, environment troubleshooting, and workflow design. While DeepSeek acted as an exceptional "co-pilot" during countless debugging sessions, the core scientific hypotheses, receptor sampling strategy, manual quality inspection, and the R100-R409 grading system were independently conceived and executed by the author. Sincere gratitude goes to DeepSeek for being a reliable guide throughout this summer journey.*