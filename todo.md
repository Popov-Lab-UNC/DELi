- [ ] work on k-mer

Need to figure out how I want to account for anchoring of reads...
maybe for now I can just skip it and move on and cycle back for speed later

Three different alignment modes
- [ ] Align to some fixed region in the barcode (if long enough) and call library then BBs
This one seems like it might be hard to generalize if a bunch of different schemas are used
- [ ] library demultiplex first then align and decode
- [ ] kmer alignment algorithm (numeric based?)

Need to make sure all library tags are at the front, need logic to reverse library barcodes if
needed to assert this logic

add headpiece part to barcode schema

think about adding mappy to make things go faster
