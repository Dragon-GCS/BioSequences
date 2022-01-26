# BioSequences

![PyPI - Downloads](https://img.shields.io/pypi/dw/Biosequences?logo=pypi&label=download) ![version](https://img.shields.io/pypi/v/biosequences?label=version&logo=pypi&color=lightgrey)
![python version](https://img.shields.io/static/v1?label=python&message=>=3.8&color=orange&logo=python) ![PyPI - License](https://img.shields.io/pypi/l/biosequences?logo=gnuprivacyguard&color=green)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/biosequences) ![GitHub last commit](https://img.shields.io/github/last-commit/Dragon-GCS/Biosequences?color=yellowgreen) ![GitHub Repo stars](https://img.shields.io/github/stars/dragon-gcs/biosequences?color=blue)

---

用于分析核酸与肽段序列

## 下载源码编译

    python setup.py build_ext --inplace
    rm ./build

## pip安装

    pip install biosequences

---

# 主要功能

## bioseq.Sequence(seq="", info="")

* RNA，DNA和Peptide都基于此抽象类，因此Sequence中的属性和方法为所有序列对象公有的属性和方法。
* 相同的序列对象可以直接与同类对象或字符串进行拼接，比较。
* 所有对象都不会对seq进行检查，所以构建对象时需要主要seq中不要出现不应该出现的字符，以免发生不必要的问题

    ```python
    from bioseq import DNA, Peptide

    d1 = DNA("ATCC")
    d2 = DNA("AC")
    p1 = Peptide("MATN")

    d1  # 5'-ATCC-3'
    p1  # N-MATN-C
    d1 + d2  # 5'-ATCCAC-3'
    d2 + d1  # 5'-ACATCC-3'
    d1 + p2  # TypeError(Expected str or DNA add with DNA, not Peptide)
    d1 == d2  # False
    ```

### 属性

#### `seq`

序列信息，不可修改

#### `info`

序列的一些说明信息(可选)

#### `length`

序列的长度

#### `weight`

序列的分子量

#### `composition`

序列中各个单位的含量

### 方法

#### `align(subject, mode=1)`

    ```python
    subject(str | Sequence)：比对对象
    mode(int)：
        1 - 使用Needleman-Wunsch进行全局比对
        2 - 使用Smith-Waterman进行局部比对
    ```

#### `find(target)`

在序列中查找目标序列并返回所有匹配的起始位置

    ```python
    target(str| Sequence)：目标序列
    ```

#### `mutation(position, target)`

改变序列信息

    ```python
    position(str | int | List[int])：需要修改的单个字符或者是需要修改的字符串起始位置。
    target(str| Sequence)：目标序列
    ```

#### `toDNA()`, `toRNA()`, `toPeptide()`

将`Sequence`序列转换为对应的生物序列

## `bioseq.RNA`

用于存储RNA序列信息。

### 属性

#### `revered`

    返回序列的反向RNA序列

#### `complemented`

    返回序列的反向互补RNA序列

#### `GC`

    返回序列的GC含量

#### `orf`

    序列中的开放读码框，使用过getOrf()方法后才具有此属性

#### `peptide`

    序列转录产物，使用过tanscript()后才有此属性

### 方法

#### `revers()`

将序列自身变为其反向序列。注意：会修改序列自身

#### `complemented()`

将序列自身变为其反向互补序列。注意：会修改序列自身

#### `getOrf(topn=1, replace=False)`

获取序列上的ORF

    ```python
    topn(int)：查找序列中的所有orf，并返回最长的topn个读码框
    replace(bool)： 当multi=False时生效，是否将最长的orf替换为原序列
    ```

#### `transcript(topn=1)`

将序列翻译为肽链

    ```python
    topn(bool)：根据读码框返回长度最长的topn个翻译产物。翻译产物均为Peptide对象。
    ```

## `bioseq.DNA`

用于存储DNA序列信息。

### 方法

#### `translate()`

将DNA翻译为RNA对象并返回

#### `transcript(topn=1)`

将序列翻译为肽链

    ```python
    topn(bool)：根据读码框返回长度最长的topn个翻译产物。翻译产物均为Peptide对象。
    ```

## `bioseq.Peptide`

用于存储肽链序列信息。

### 属性

#### `pI`

基于EMBOSS数据库中氨基酸的pK值， 计算该肽链序列的等电点并返回

### 方法

#### `chargeInpH(pH: float)`

基于EMBOSS数据库中氨基酸的pK值，计算肽链在某一pH下所带的电荷量

    ```python
    pH(float): 溶液的pH值
    ```

#### `getHphob(window_size=9, show_img=True)`

基于Doolittle(1982)的氨基酸疏水性数据，计算肽链的疏水性，疏水性

    ```python
    window_size(int)：某一氨基酸的疏水性为window_size内该氨基酸位于window中心时的所有氨基酸疏水性的平均值
    show_img：绘制疏水性结果，需要matplotlib
    ```

## `bioseq.config`

可在此文件中直接修改配置数据，或通过以下函数在运行时修改部分数据

#### `setAlignPara(match = 2, mismatch = -3, gap_open = -3, gap_extend = -3)`

修改序列比对时的评分规则，需要在比对前进行设置

    ```python
    match(int) ：匹配得分(>0)
    mismath(int)：错配得分(<0)
    gap_open(int)：开口得分(<0)
    gap_extend(int)：开口延长得分(<0) 

    d1 = DNA("ATCTCGC")
    d2 = DNA("ATCCC")

    print(d1.align(d2)) #('ATCTCGC', 'ATC-C-C', 4.0)
    setAlignPara(5)
    print(d1.align(d2)) #('ATCTCGC', 'A--TCCC', -0.5)
    ```

#### `setStartCoden(coden = None)`

修改核酸序列转录时需要的起始密码子，为传入coden则将密码子初始化为*"AUG"*

    ```python
    coden(str | List(str))：密码子会在coden中寻找，如有匹配则开始进行转录

    d1 = DNA("ATCATCTCAGCATGAC")

    print(d1.transcript(filtered=False)) # []
    setStartCoden(["AUC"])
    print(d1.transcript(filtered=False)) # [N-IISA-C, N-ISA-C]
    ```

## bioseq.utils

工具

#### `printAlign(sequence1, sequence2, spacing=10, line_width=30, show_seq=True)`

在命令行中按格式输出两个比对后的序列， 可在config.SYMBOL中修改显示的符号

    ```python
    spacing(int)：序列显示间隔
    line_width(int)：每行显示的字符数
    show_sequence(bool)：是否显示序列

    d1 = DNA("ATCATCTCAGCATGAC")
    d2 = DNA("ATCATCGCATGAC")

    seq1, seq2 = d1.align(d2)
    printAlign(d1, d2)
    #    1 ATCATCTCAG CAT
    #      ┃┃┃┃┃┃•┃┃• •┃•
    #    1 ATCATCGCAT GAC
    printAlign(d1, d2, spacing=3, line_width=10, show_seq=False)
    #    1 ┃┃┃ ┃┃┃ •┃┃ •
    # 
    #   11 •┃• 
    ```

#### `loadFasta(filename)`

读取fasta文件，并返回所有读取到的(序列列表，序列名列表)**Todo：加入更多解析格式**

#### fetchNCBI(uid)

    ```python
    uid(str): NCBI中序列的唯一编号，如 NC_XXXX、NM_XXXX等，仅限于DNA(mRNA)、RNA和多肽序列，返回对应的序列对象。
    ```

> NCBI RefSeq's document: <https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole>
> some NCBI E-utilities's api: <https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/>

## Change Log

Version: **1.1.3**

fix: Wrong results from `bioseq.utils.loadFasta()`

Version: **1.1.2**

* change: `RNA.getOrf()`, `RNA.transcript()`
* add: `algoritm.pyi`

Version: **1.1.1**

* add: unittest
* fix: some bug found in unintest

Version: **1.1.0**

* add: `info`attribute for `Sequence`
* add: `toDNA()`, `toRNA()`, `toPeptide()` method for `Sequence`
* add: `utils.fetchNCBI()`
change: `utils.read_fasta()` to `utils.loadFasta` and be a generator of `Sequence`

Version: **1.0.9**

* add: add type annotations, remove `*.pyi` file
* add: `Sequence.reset_cache()` to reset some cached property, to update the value after mutation, include `weight`, `composition`, `GC(DNA, RNA)`, `orf(DNA, RNA)`, `peptide(DNA, RNA)`, `translate(DNA, RNA)`, `pI(Peptide)`, `Hphob_list(Peptide)`.
* add: warning when mutation overlaped previous mutaion
* fix: some wrong typing check in `Sequence.find()`, `Sequence.mutation()`
* remove: `return_score: bool` for `Sequence.align()`
