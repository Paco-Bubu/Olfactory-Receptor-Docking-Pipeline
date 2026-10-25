# 配体准备（手动步骤）

从 PubChem 按以下 CID 下载 SDF，并用自己安装的 OpenBabel 转换为 PDBQT。
转换后将文件放入本目录，使用小写英文名称：

| 名称 | CID | Notebook 需要的文件名 |
| --- | --- | --- |
| 香茅醛 Citronellal | 7794 | `citronellal.pdbqt` |
| 丁香酚 Eugenol | 3314 | `eugenol.pdbqt` |
| 芳樟醇 Linalool | 6549 | `linalool.pdbqt` |
| 柠檬烯 Limonene | 22311 | `limonene.pdbqt` |

原始 SDF、PDB 和 PDBQT 文件由 `.gitignore` 排除，不随仓库分发。
请在本机检查结构准备结果；Notebook 不包含自动配体下载、质子化核查或完整预处理流程。
