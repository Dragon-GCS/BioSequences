# ChangeLog

## Version: **1.1.5**

* add: `docs/` by readthedocs
* change: merge `complemented` and `complement()` to `complement`, merge `reversed` and `reverse()` to `reversed`
* change: `bioseq.config.TABLE` to `bioseq.config.CODON_TABLE`
* remove: `bioseq.config.setAlignPara`, `bioseq.config.setStartCoden`

## Version: **1.1.4**

* add: `bioseq.utils.fetchENS()`
* add: `bioseq.utils.parseFasta()`
* add: `iterator` option of `bioseq.utils.loadFasta()`
* add: multi-uids fetch support for `bioseq.utils.fetchNCBI()`

## Version: **1.1.3**

* fix: Wrong results from `bioseq.utils.loadFasta()`

## Version: **1.1.2**

* change: `RNA.getOrf()`, `RNA.transcript()`
* add: `algoritm.pyi`

## Version: **1.1.1**

* add: unittest
* fix: some bug found in unintest

## Version: **1.1.0**

* add: `info`attribute for `Sequence`
* add: `toDNA()`, `toRNA()`, `toPeptide()` method for `Sequence`
* add: `utils.fetchNCBI()`
* change: `utils.read_fasta()` to `utils.loadFasta` and be a generator of `Sequence`

## Version: **1.0.9**

* add: add type annotations, remove `*.pyi` file
* add: `Sequence.reset_cache()` to reset some cached property, to update the value after mutation, include `weight`, `composition`, `GC(DNA, RNA)`, `orf(DNA, RNA)`, `peptide(DNA, RNA)`, `translate(DNA, RNA)`, `pI(Peptide)`, `Hphob_list(Peptide)`.
* add: warning when mutation overlapped previous mutation
* fix: some wrong typing check in `Sequence.find()`, `Sequence.mutation()`
* remove: `return_score: bool` for `Sequence.align()`
