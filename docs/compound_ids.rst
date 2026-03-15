.. _compound_ids:

================================
DELi compound naming conventions
================================

Making unqiue compound IDs for DELs is not nessiarily trivial. While you can just assign random numbers to each compound,
that infomration cannot be easily encoded into the DNA tags. Thus post decoding you would need to map to building block and library
information to get a compound ID, which is not ideal.

Instead DELi takes the approach of generating compound IDs from the building block and library information directly. This way, the IDs
can be generated directly during decoding.

Component ID compound IDs
-------------------------
The easiest way to do this is to just concatenate the library and building block IDs together with a separator character to make a single unquie Compound ID.
This works well since DELi already asserts that the library IDs are unique (at least within a selection) and building block IDs are unique within each cycle.
For example, if you have a compound made from library "Lib1" and building blocks "BB1", "BB2", and "BB3", the compound ID would be "Lib1-BB1-BB2-BB3".
This is easily interpretable and can help with manualy inspecting the data (one can tell if the compound has the same building block as another just by looking at the compound ID).

The one caveat is that a seperating character must be chosen to specifcy when we are changing from one component to the next. This character **MUST NOT**
be used in the library or building block IDs, otherwise the compound ID cannot be garunteed to be unique. For example if I have building blocks "BB1" and "BB1-BB2" in cycle 1 and
"BB2-BB3" and "BB3" in cycle 2, then we could have two compounds with the same compound ID "Lib1-BB1-BB2-BB3" (one made from BB1 and BB2-BB3, and one made from BB1-BB2 and BB3).
DELi defaults to using the "-" character as a separator, but this can be changed if desired (see below).
You can use the `deli data validate-ids --separator <your-separator>` command to check if a given separator is valid for your libraries and building blocks currently in your DELi data
directory. However, it is still the user's responsibility to ensure that the separator is respected in all future library and building block IDs.

Numeric in-library compound IDs
-------------------------------
DELi can also generate a numeric in-library compound ID, which is a unique integer for each compound within a given library.
This can be useful if you want to have shorter compound IDs that are still unique within the library, but do not need to be visually interpretable.
To calculate these we use the following approach: For each building block within each cycle C of size N we have building block  IDs [1, 2, ..., N].
Then for any set of building block IDs :math:`[BB_1, BB_2, \ldots, BB_K]` inside cycles
:math:`[C_1, C_2, \ldots, C_K]` with sizes :math:`[c_1, c_2, \ldots, c_K]` respectively,
we can use the formula:

.. math::

   BB_1 + (c_1 \cdot BB_2) + (c_1 \cdot c_2 \cdot BB_3) + \cdots
   + (c_1 \cdot c_2 \cdots c_{K-1} \cdot BB_K)

Or equivalently:

.. math::

   \sum_{i=1}^{K} \left( \prod_{j=1}^{i-1} c_j \right) \cdot BB_i

DELi offsets these by 1 so that 0 can represent the "null" building block (i.e. no building
block added in that cycle). The actual equation used is therefore:

.. math::

   \sum_{i=1}^{K} \left( \prod_{j=1}^{i-1} (c_j + 1) \right) \cdot BB_i

These are only unique within a given library, however, they can be combined with the library ID to make them unique across libraries. If you assign your libraries a numeric ID as well, it is also
possible to generate a single unique numeric ID for each compound across the entire DEL.
There are some downsides to this approach, the most pressing being that the IDs are now linked to the order of the building blocks within the files passed to DELi.
If the order of the files changes, so will the generated ID. This contract is also not obvious unlike a set ID column. An uninformed user may not realize the impact
modifying row order could have.

Which to use
------------
In nearly all cases, the component ID compound IDs are the best choice. They are more stable, easier to interpret, and do not have any hidden dependencies on the order of the building blocks.
The only downside to this is ensuring that the separator character is not used in the library and building block IDs. As long as this is respected, it will be garuntee that the IDs are unique
and can be used to regenerate the orginal DEL compound they encode.

Compression
^^^^^^^^^^^
At a surface glance it might seem like numeric IDs are more compact than component IDs, and that is true only when store them as integers in a non-compressed binary format.
However, even light weight compression like snappy will be able to compress the component ID strings significantly, as they are highly repetitive. In practice the component IDs
are not significantly larger than the numeric IDs when stored in a compressed format.
