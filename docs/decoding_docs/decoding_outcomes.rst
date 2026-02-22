.. _decode-outcomes-docs:

=================
Decoding Outcomes
=================
During decoding, DELi wil either successfully decode a read or fail to decode it.
Successful decoding will result in a decoded read that includes the library, building blocks,
and UMI.
If it fails, it will return a ``FailedDecode`` object that indicates why the read could not be decoded.

Failed Decode Outcomes
----------------------
Failed decodes will include the sequence that fail and a reason for the failure in plain
text. Below are the types that can occur and a description of each.

- ``ReadTooShort``: The read is too short to decode based on the settings.
- ``ReadTooLong``: The read is too long to decode based on the settings.
- ``FailedStaticAlignment``: The static alignment failed to any matching patterns in the read.
- ``FailedLibraryBarcodeLookup``: The library barcode lookup was not successful.
- ``AlignmentFailed``: The alignment of the read to the library schema was too poor to proceed.
- ``FailedBuildingBlockCall``: The building block barcode is not present in the hashmap
- ``AmbiguousBuildingBlockBarcode``: The building block barcode matches multiple building blocks in the hashmap.
- ``ImpossibleBuildingBlockBarcode``: The decoded barcodes overlap (usually caused by many sequencing errors resulting a incorrect calls)
- ``UMIContainsAmbiguity``: The library lookup resulted in a match that was too short to call.
- ``UMIContainsINDEL``: DELi detected an INDEL in the UMI section of the read

The DELi Decoding report will include a summary of how often these outcomes are during decoding.
This can be useful to help debug any issues with the decoding process or to understand the quality
of the reads being decoded.
