bioseq
----------------
.. automodule:: bioseq
   :members: Sequence, Peptide, RNA, RNA

.. autoclass:: bioseq.Sequence
    :members: 
        seq, align, composition, length, weight, find,
        mutation, toDNA, toRNA, toPeptide, _print
    :special-members: __init__

.. autoclass:: bioseq.Peptide
    :members: pI, chargeInpH, getHphob, _print

.. autoclass:: bioseq.RNA
    :members: 
        GC, complement, reversed, peptide, orf,
        getOrf, transcript, _print
    :undoc-members:
        peptide

.. autoclass:: bioseq.DNA
    :members: 
        GC, complement, reversed, peptide, orf,
        getOrf, translate, transcript
    :undoc-members:
        peptide

