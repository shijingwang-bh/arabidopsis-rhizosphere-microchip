# 拟南芥根际芯片数据与代码

公开仓库：[shijingwang-bh/arabidopsis-rhizosphere-microchip](https://github.com/shijingwang-bh/arabidopsis-rhizosphere-microchip)。

本包对应 2026 年 11 月 25 日最新正文和 SI，版本为 **2026.11.25**。包含四个生长性状平均接种效应的联合 BH 校正，沿用最终图号：正文图 1–5、补充图 S1–S7（S5a、S5b 分文件）及表 S1–S5。芯片培养比较位于 S1d；RNA PCA 仅位于正文图 2c。

英文总说明见 [README.md](README.md)。正文版本及 SHA256 见 `metadata/manuscript_version.json`；逐图数据、代码及统计家族对应见 `docs/figure_data_map.csv`。

## 使用

先安装 Python 依赖，确认 R 依赖可用，然后在解压目录运行：

```sh
python -m pip install -r requirements.txt
python validate_package.py
python reproduce.py --mode all --output ../microchip_run --rscript Rscript
```

Windows 如果 Rscript 不在 PATH 中，把最后一个参数换成本机 `Rscript.exe` 的完整路径并加引号。`--output` 必须是包外的新目录或空目录。只重算用 `--mode analysis`；只从归档表重画用 `--mode figures`。运行过程不会自动安装包或改写包内数据。

## 本版统计规则

- 四个生长性状的跨重力条件平均接种效应作为一个家族联合 BH 校正。根长、成熟根直径、叶鲜重和叶面积的 q 约为 0.432、0.00408、0.00226、0.431。
- 细菌浓度及植物生长模型项报告未校正 P；效应的 95% CI 仍未作多重校正。
- 三个芯片培养性状与侧根比较仍是原来的四检验家族；图件移动不改变校正范围。
- RNA 对全部 20,449 个基因校正；GO 每种方法分别对全部 1,705 个集合校正。
- LC–MS 每个对比在全部 3,126 个特征内校正，并使用绝对效应阈值；不在展示的两种化合物内重新校正。

完整说明见 `docs/statistical_rules.md`。完整结果与出图数据分开保存，阴性结果及所有受检条目都可追溯。

## 数据起点与范围

RNA 从现有 Salmon count/TPM 矩阵开始，保留小数计数，不补造有效长度 offset 或原始定量命令。LC–MS 从供应商强度/注释表开始，其 `.xls` 实为 GB18030 制表符文本。原始 FASTQ、质谱仪采集文件和 COMSOL 模型不在这个统计复现包内；芯片和显微图版作为归档图像输入保留。

S7 提供由相同观测、CI 和完整检验家族 q 值生成的矢量重绘代码；原始矢量脚本未找到，因此同时保留提交版的精确栅格外观。该整理不改正文和 SI 文件。数据字典对独立芯片、植株、叶混样、培养皿收集及技术孔作了明确区分。

## 上传和引用

解压后将目录内容作为 GitHub 仓库文件提交，具体步骤见 `docs/UPLOAD.md`。保留 `.gitignore`、`.gitattributes`、校验清单、README 和数据目录。不要把唯一一个 ZIP 文件当成完整代码仓库。

代码沿用原有 MIT 许可；数据继续沿用原有单独使用条款，未另行指定数据许可证。作者与版本信息在 `CITATION.cff`。实际完成的验证见 `validation/VALIDATION.md`。
