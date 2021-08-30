# biosequence
---

用于分析核酸与肽段序列

---
## Todo 20200830
### _Sequence
- [x] _Sequence.\__init\__(seq)
- [x] @property _Sequence.seq()
- [x] @property _Sequence.length()
- [ ] @property _Sequence.weight() 计算分子量
- [ ] _Sequence.align(subject, mode=1, boost=True, return_score=False) 序列比对
- [ ] _Sequence.analysis()  成分分析
- [ ] _Sequence.find(target) 查找序列
- [ ] _Sequence.mutation(position, target) 改变序列
- [x] _Sequence._print() 打印输出格式
- [x] _Sequence.\__str\__()
- [x] _Sequence.\__repr\__()
- [x] _Sequence.\__add\__()
- [x] _Sequence.\__radd\__()

### RNA
- [ ] @property RNA.revere 将该序列变为其反向序列并返回序列
- [ ] @property RNA.complement 将该序列变为其互补序列并返回序列
- [ ] @property RNA.GC 返回序列中GC含量,计算后保存在_GC中
- [ ] RNA.orf 序列中的开放读码框，需要先经过get_orf()计算才有此属性
- [ ] RNA.peptide 序列翻译产物，需要先经过tanscript()计算才有此属性
- [ ] RNA.reversed() 计算返回序列的反向序列
- [ ] RNA.complemented() 计算返回序列的互补序列
- [ ] RNA.getOrf(multi=False, replace=False) multi是否查找所有frame +1~+3的orf，默认值为仅查找最长的orf。 replace 当multi=False是生效，是否使用最长的orf替换原序列
- [ ] RNA.transcript(filter) filter是否仅返回最长的翻译产物。返回值为一个或多个Peptide对象。

- [ ] RNA._print() 改为5'-Seq-3'

### DNA
- [ ] DNA.translate() 返回翻译后的RNA对象
- [ ] DNA.transcript(filter) 翻译后进行转录

### Peptide
- [ ] Peptide.pl 计算肽链的等电点，第一次计算后保存在_pl中
- [ ] Peptide.Hphob 计算肽链的亲疏水性，第一次计算后保存在_Hphob中
    info from [expasy](https://web.expasy.org/protscale/)

### config
- [x] `config.setTable("conden.json")`: 使用保存在conden_table中的自定义密码子表（必须为json格式）
- [x] `config.setAlignPara(match, mismatch, gap_open, gap_extend)`: 修改比对参数

### biosequence.utils
- [ ] biosequence.utils.calculatePL（peptide）： 计算等电点
- [ ] biosequence.utils.calculateHphob（peptide）： 计算亲疏水性
- [x] biosequence.align.alignment(query, subject, mode=1, boost=True, return_score=False)
- [ ] biosequence.utils.parseFile(filename, type)：读取常见文件格式并返回对应的序列对象
- [ ] biosequence.align.printAlign(seq1, seq2, spacing=10, line_length=30)：在命令行中打印两个比对序列并显示差异
---

# 主要功能
## biosequence.Sequence

### biosequence.Sequence._Sequence

生物序列的基类，RNA，DNA，Peptide都基于此类。
序列对象可与序列对象或字符串相加，返回新的序列对象

#### 属性
##### Seq
    序列信息，不可修改（实际序列信息保存在内部属性_Seq中）
##### length
    序列的长度
##### weight
    序列的分子量（第一次计算后保存在_weight中）

#### 方法 
##### align(subject, mode=1, boost=True, return_score=False)

* subject：
    比对对象
* mode（int）：
    1 - 使用Needleman-Wunsch进行全局比对
    2 - 使用Smith-Waterman进行局部比对
    
* boost（bool）：
    是否使用c算法加速，default=True。（优化c语言算法后默认使用，并删除此参数）

##### mutation(position, target)
    改变序列

##### find(target)
    在序列中查找目标序列并返回所有匹配的起始位置

#### analysis()
    打印并返回序列中所有碱基/氨基酸的数量、含量

### biosequence.sequence.RNA
继承自biosequence.Sequence._Sequence
#### 属性
##### revere
    将该序列变为其反向序列

##### complemente
    将该序列变为其互补序列

##### orf
    序列中的开放读码框，需要先经过get_orf()计算才有此属性

##### peptide
    序列翻译产物，需要先经过tanscript()计算才有此属性


#### 方法 
##### GC()
    返回序列中GC含量

#### reversed()
    返回序列的反向序列

#### complemented()
    返回序列的互补序列

#### getOrf(multi=False, replace=False)
* multi（bool）
    是否查找所有frame +1~+3的orf，默认值为仅查找最长的orf
* replace（bool）
    当multi=False是生效，是否使用最长的orf替换原序列

### biosequence.sequence.DNA
继承自biosequence.Sequence.RNA

#### translate()
    返回翻译后的RNA对象

#### transcript(filter)
* filter（bool）    
    是否仅返回最长的翻译产物
    返回值为一个或多个Peptide对象

### biosequence.sequence.Peptide
继承自biosequence.Sequence._Sequence

#### 属性
##### pl
    肽链的等电点，第一次计算后保存在_pl中
##### Hphob
    肽链的亲疏水性，第一次计算后保存在_Hphob中
    info from [expasy](https://web.expasy.org/protscale/)

### biosequence.config
```python
from biosequence import config
# or import biosequence
# don't from biosequence.config import xxx 
```
#### config.setTable(filename)
    序列翻译前加入此代码，可以使用自定义密码子进行翻译

* filename(str)
    保存在`conden_table`文件夹中以`json`格式保存的密码子表

#### setAlignPara(match, mismatch, gap_open, gap_extend)
    序列比对前加入此代码，可以修改比对评分
* match(int) 
  匹配得分
* mismath(int)    
  错配得分
* gap_open(int)     
  开口得分
* gap_extend(int)   
  开口延长得分 
      