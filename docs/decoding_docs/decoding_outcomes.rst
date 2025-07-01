.. _decode-outcomes-docs:

=================
Decoding Outcomes
=================

During decoding, DELi will either fail to decode a read or it will successfully decode it into a compound ID.
While all successes are the same, DELi has various reasons why it might fail to decode a read.

Below is a list of the possible failed outcomes of decoding a read:

- ``ReadTooShort``: The read is too short to decode based on the settings.
- ``ReadTooLong``: The read is too long to decode based on the settings.
- ``LibraryLookupFailed``: The library lookup was not successful.
- ``LibraryMatchTooShort``: The library lookup resulted in a match that was too short to call.
- ``AlignmentFailed``: The alignment is not successful during calling.
- ``BuildingBlockLookupFailed``: The building block lookup was not successful.
- ``UMIMatchTooShort``: The UMI match is too short to call post decode.

The DELi Decoding report will include a summary of how often these outcomes are during decoding.
This can be useful to help debug any issues with the decoding process or to understand the quality of the reads being decoded.
