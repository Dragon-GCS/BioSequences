# BioSequences

![PyPI - Downloads](https://img.shields.io/pypi/dw/Biosequences?logo=pypi&label=download) ![version](https://img.shields.io/pypi/v/biosequences?label=version&logo=pypi&color=lightgrey)
![python version](https://img.shields.io/static/v1?label=python&message=>=3.8&color=orange&logo=python) ![PyPI - License](https://img.shields.io/pypi/l/biosequences?logo=gnuprivacyguard&color=green)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/biosequences) ![GitHub last commit](https://img.shields.io/github/last-commit/Dragon-GCS/Biosequences?color=yellowgreen) ![GitHub Repo stars](https://img.shields.io/github/stars/dragon-gcs/biosequences?color=blue)

---

- [BioSequences](#biosequences)
  - [关于本项目](#关于本项目)
  - [安装](#安装)
    - [pip 安装](#pip-安装)
    - [下载源码安装](#下载源码安装)
  - [示例](#示例)
    - [加载序列信息](#加载序列信息)
    - [序列基本操作](#序列基本操作)
  - [贡献者](#贡献者)
  - [致谢](#致谢)

---

## 关于本项目

**BioSequences**是一个集合了基本的常用的生物序列分析工具的包，旨在提高日常一些基本序列分析流程的工作效率，以及为大数据分析提供一些基础支持。

完整文档请看这里[Document](https://biosequences.readthedocs.io/)。

## 安装

### pip 安装

```ps1
pip install biosequences
```

### 下载源码安装

windows下需要安装``Microsoft VC++``编译工具, Linux 需要安装gcc或其他编译工具。

```ps1
git clone https://github.com/Dragon-GCS/BioSequences.git
cd BioSequences
python -m pip install BioSequences
```

## 示例

### 加载序列信息

`bioseq`可以从标准fasta格式的文件或NCBI/Ensemble数据库读取序列信息。当`fetch`方法的参数为列表时可以批量抓取目标序列。

```python
>>> from bioseq.utils import loadFasta, fetchNCBI, fetchENS
>>> sequence1 = loadFasta("/path/to/file.fasta")
>>> bsa = fetchNCBI("NP_851335.1")
>>> actin = fetchENS("ENST00000614376")
```

### 序列基本操作

``bioseq.RNA``，``bioseq.DNA`` 和 ``bioseq.Peptide`` 都继承自 ``bioseq.Sequence``，因此三者基本操作基本一致。

* 查看序列的基本属性

    ```python
    >>> actin.GC, actin.length
    (0.5, 102)
    >>> actin.composition
    {'A': 24, 'C': 18, 'G': 33, 'T': 27}
    >>> actin.seq
    'AGAAACTTTAGCATCTGGCTAGGAGCATCTGTGGTGGCTCACCTTTCTACCTATACGTGTGAGTGGGTGACCTGAGAGGAGTACGGTGAGCATATGAGGATG'
    >>> round(bsa.weight, 1)
    69334.4
    >>> bsa.pI
    6.805
    >>> round(bsa.chargeInpH(7.4), 2)
    -13.76
    ```

* DNA序列或RNA序列可以进行转录`transcript()`，DNA序列有`translate()`方法可以翻译为RNA序列。
  还可以通过`bioseq.config.START_CODON`自定义起始密码子，以及通过修改`bioseq.config.CODON——TABLE`自定义密码子表。

    ```python
    >>> from bioseq.config import START_CODON, CODON_TABLE
    >>> actin.transcript()
    >>> START_CODON[0] = 'AGA'
    >>> actin.transcript()
    [N-RNFSIWLGASVVAHLSTYTCEWVT-C]
    >>> CODON_TABLE["AAC"] = "Y"
    >>> actin.transcript()
    [N-RYFSIWLGASVVAHLSTYTCEWVT-C]
    ```

* 两个相同类型的序列可以进行拼接

    ```python
    >>> from bioseq import DNA
    >>> dna1 = DNA("ATCG")
    >>> dna2 = DNA("GCAT")
    >>> dna1 + dna2
    "5'-ATCGGCAT-3'"
    >>> dna2 + dna1
    "5'-GCATATCG-3'"
    ```

* 通过`mutation()`方法对序列进行修改

    ```python
    >>> dna1.mutation("ATC", "GGG")
    'GGGG'
    >>> dna1.mutation(0, "AT")
    'ATGG'
    >>> dna1.mutation([0, 3], "C")
    'CTGC'
    ```

* `Sequence`用C语言实现了`Needleman-Wunsch`全局比对和`Smith-Waterman`局部比对两种基本的序列匹配算法，可以用来快速比对序列（局部比对仅返回匹配的局部序列）。

    ```python
    >>> DNA("GCATGCT").align("GATTACA")
    ('GCA-TGCT', 'G-ATTACA', -4.0)
    >>> DNA("GCATGCT").align("GATTACA", 2)
    ('AT', 'AT', 4.0)
    ```

    比对返回的前两个参数为比对后的序列，第三个参数为匹配得分，可以通过`bioseq.utils.printAlign()`来优化比对结果的显示。

    ```python
    >>> from bioseq.utils import printAlign
    >>> seq1, seq2, score = DNA("GCATGCT").align("GATTACA")
    >>> printAlign(seq1, seq2)
    1 GCA-TGCT
      ┃━┃━┃•┃•
    1 G-ATTACA
    ```

    可以通过修改`bioseq.config.AlignmentConfig`来修改匹配时的罚分，默认为`MATCH(2.0), MISMATCH(-3.0), GAP_OPEN: (-3.0), GAP_EXTEND(-3.0)`

    ```python
    >>> from bioseq.config import AlignmentConfig
    >>> AlignmentConfig.GAP_OPEN = -10
    >>> DNA("GCATGCT").align("GATTACA")
    ('GCATGCT', 'GATTACA', -6.0)
    ```

## 贡献者

<a href="https://github.com/Dragon-GCS">
<img src="https://avatars.githubusercontent.com/u/54531807?v=4" alt="@Dragon-GCS" height="40" style="border-radius: 100%; border: 2px solid">
</a>
<a href="https://github.com/laxtiz">
<img src="https://avatars.githubusercontent.com/u/3883767?v=4" alt="@laxtiz" height="40" style="border-radius: 100%; border: 2px solid">
</a>

## 致谢

* [Read the Docs](https://readthedocs.org/)
* [Img Shields](https://shields.io)
