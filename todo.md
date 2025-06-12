If I need to save on memery, I could convert UMI to a HexaDecimal to cut the length of UMI in half
(should cut file size in half?)


After some test, looks like I should just have the merge and counter stuff happen in Nextflow, not
in the cli.py file.

Also the pickles are too much overhead, need to just use json instead, which should be fine for
parallel because we don't care about the objects anymore, it is all file based.
Pickles seem to be added as much as 20 seconds of overhead due to long loading times.
Can just get away with jsons

Will need to think of logic for UMI verus no UMI, not that I think anyone ever does no UMI

Also notice that taking the set of a list is far faster then doing constant unions of the list

After running, the merging of the decode result is too slow when in parallel.
This is because i'm using pickle.
I need to come up with a way to avoid that.

My thoughts are to make a new class, the isn't "CalledBarcode" but just "DELCompound". This
will be what the decoder returns.

In turn, this will have functions that let you load and save raw text representations of the
compounds so we can reload them with library support later.

Also, the degen object needs to stop hashing on the compounds themselves and start
working on the ids as strings instead.
The UMI version should just be a nested dict, where the compound count is a list of umis.
This way, we can save the output as a json which will be easier to merge after the fact cause its is just raw text.
Then we can use some other util to load all the DELCompound objects if we need them (for example to generate SMILES).
This will be way better than using pickles.
