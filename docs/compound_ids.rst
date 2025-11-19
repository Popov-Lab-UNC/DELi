================================
DELi Compound Naming Conventions
================================

Internally, DELi generates unique compound IDs for the chemicals
in the DEL using the format `<LIB_ID>-<BB1_ID>-<BB2_ID>-<BB3_ID>...`,
where `LIB_ID` is the library ID, and `BB1`, `BB2`, `BB3` are the building block IDs.

DELi requires that a DEL compound ID be generatable (and unique) from just the building blocks and library information.
Otherwise, it would not be possible to name the compound after decoding.
The best and easiest way to achieve this is to use the IDs of those components themselves. DELi does just this and uses the `-`
character to separate the components. Currently, there is no way to alter this. Therefore, a (currently) unenforced requirement
of DELi (for decoding only) is that all BB IDs and library IDs do not use the `-` character.

This constraint is not required if you are only using DELi for analysis of your own cube files.
You can use any IDs you want, as no assumptions are made by default about the IDs.

.. note::
    If you find yourself wanting support for a different compound ID format,
    raise an issue on the DELi GitHub repository to encourage us to add support for it.
