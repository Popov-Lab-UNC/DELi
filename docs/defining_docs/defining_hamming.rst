Defining Hamming Codes
======================

Most modern DEL libraries utilize a hamming encoded scheme to enable fixing of single base
pair misreads during sequencing (or during PCR) of the selection. DELi support this
functionality if you libraries have it, but you need to tell DELi the how your hamming
code is set up.

DELi support but normal hamming codes AND extra parity hamming codes.
DELi does not support other error correction encoding schemes.

To keep this simple, hamming codes are stores a just two lines:
::
    hamming_order: p0,p1,p2,d3,p4,d5,d6,d7
    custom_order: p0,p1,p2,d3,p4,d5,d6,d7

This is all that is needed for DELi to understand how your code is set up,
but you will need to understand hamming codes (specifically how yours is designed) in order
to provide this file. You can read more about hamming codes `here <https://en.wikipedia.org/wiki/Hamming_code>`_

Hamming Order:
--------------
This is the "correct" order of the hamming code.
Generally hamming codes are laid out in a "power of 2" order.
So, the bit in position 1 is a parity bit, :math:`2^0=1`.
The next bit, position 2, is also a parity bit, :math:`2^1=2`
But the third bit is a data bit as there is no way for :math:`2^N` to equal 3.

Parity bits are abbreviated as ``p`` and data bits as ``d``
Thus if I had the traditional hamming(7, 4) code, ``hamming_order`` should be
``p1,p2,d3,p4,d5,d6,d7``. DELi really only needs this as a way to check that your hamming
code is something it can process, as under the hood, DELi will want to arrange your hamming
code into this order. If you cannot provide this mapping, DELi cannot do hamming decoding for you

Custom Order
------------
This is where you can provide the order of your parity bits to DELi.
Another common way for hamming codes to built is to put all the parity bits at the beginning,
and then have all the data bits follow, e.g. ``p1,p2,p4,d3,d5,d6,d7``.

The custom order allows your codes to be built this way and still be decidable by DELi.
So important parts to notice is that the number that follows a ``d`` or ``p`` bit is not
the number parity or data bit (as is often the case else where) but the *position* it has
in the list of bits for the "correct" sequence. Again, there is only 1 correct order for the
bits, but the files require you specify it as an extra check to make sure you have configure
it correctly.

Extra Parity
------------
DELi also support extra parity (sometimes called SECDED). To enable this you need to add the
``p0`` bit to the ``hamming`` and ``custom`` orders.
This would look like:
::
    hamming_order: p0,p1,p2,d3,p4,d5,d6,d7
    custom_order: p0,p1,p2,p3,d3,d5,d6,d7

This will turn on extra parity mode

.. warning::
    ``p0`` is reserved, if you use it in a non-extra parity hamming code file, it will
    cause incorrect decoding *without* raising an exception, as DELi thinks you gave it a
    valid hamming code

Example Hamming (15, 4) code
----------------------------
DELi can do a hamming code of any size, though the bigger they get the less useful they become
Here is an example file for a hamming(15, 4) code:

::
    hamming_order: p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15
    custom_order: p1,p2,d3,p4,d5,d6,d7,p8,d9,d10,d11,d12,d13,d14,d15

File Format
===========
Hamming code files do not full a common format. They are just text files, where DELi looks
for two lines: one starting with ``hamming_order:`` and the other with ``custom_order:``.
After the starting prompt, the hamming bits should follow as a comma seperated list.
The two sections should have the same number of hamming bits.
If these two lines are present only once and follow the above format, the file is valid.

Lines that start with "#" are skipped, even if they contain one of the two section starters.

DELi Data Dir
=============
Hamming code files can be saved in the DELi Data Dir sub-dir named "hamming" and end with ".txt" file extensions
