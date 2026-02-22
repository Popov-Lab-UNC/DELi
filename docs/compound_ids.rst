================================
DELi compound naming conventions
================================

Internally, DELi will generate unique compound IDs for the chemicals
in the DEL using the format <LIB_ID>-<BB1_ID>-<BB2_ID>-<BB3_ID>...
where LIB_ID is the library ID and BB1, BB2, BB3 are the building block IDs.

DELi requires that a DEL compound ID be generatable (and unique) from just the building blocks and library information,
otherwise it would not be possible to name the compound after decoding.
The best and easiest way to do this is to use the ids of those components themselves. We do just this, and use the '-'
character to separate the components. Currently, there is no way to alter this. If you find yourself wanting support for a different compound ID format,
    raise an issue on the DELi GitHub repository to encourage us to add support for it.

.. note::
    DELi does not restrict the use of the "-" character in the library and building block IDs.
    DELi **does not** implement a method to covert DEL compounds IDs back into their original library and building block
    information so it is safe to use *within DELi*.

.. note::
