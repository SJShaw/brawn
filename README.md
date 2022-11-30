# brawn
A python port of [MUSCLE](https://drive5.com/muscle/downloads_v3.htm)'s
profile-profile mode for aligning sequences.

Brawn was specifically designed with repetitive insertions of query sequences
into reference alignments in mind.

Brawn takes advantage of that repetitive nature by using cache files to remove
the calculations that would need to be repeated over and over for the reference
alignments. This counters the speed loss of moving from C to python. Initially
building these cache files will take about 10 times longer than an equivalent
single run.

While results are often identical to those from MUSCLE, there will occasionally
be some minor variations. If an alignment would be very poor in the best case,
this variation will be more pronounced.

Only protein sequences are explicitly supported at the moment, full DNA/RNA
support may follow later.

## Installation
### from PyPI
```
pip install brawn
```
### from source
```
pip install .
```
### for development
```
pip install .[testing]
```

## Use
For speed, using a cache file is extremely important. They can be built and used
both programatically and the command line. They are also optional, if required.

### Basic profile-profile mode

To load an alignment from a FASTA file, use:
```
from brawn import Alignment

with open("some.fasta") as handle:
    alignment = Alignment.from_file(handle)
```

Or, to build from an existing dictionary (even if it's just one element), use:
```
name_to_sequence = {
    "A": "GT-DVG",
    "B": "GTK-VG",
}
alignment = Alignment(name_to_sequence)
```

Or, to load cached calculations directly (this is the preferable option for large alignments):
```
with open("reference.cache") as handle:
    alignment = Alignment.from_cache_file(handle)
```

As a sidenote, to build a cache for the reference alignment:
```
with open("reference.cache", "w") as handle:
    alignment.to_cache_file(handle)
```
This will add any needed computational steps prior to writing.

Once two alignments have been created/loaded, they can be combined with:
```
from brawn import combine_alignments

result = combine_alignments(alignment1, alignment2)
```

And the resulting alignment can be output to file:
```
with open("output.fasta", "w") as handle:
    result.to_file(handle=handle)

# if handle is not supplied, the output will be sent to sys.stdout
results.to_file()
```

### Simpler insertion of a sequence into an alignment

For the use case that brawn was designed for, only a single sequence is being
inserted into a reference alignment at a time. This is usually paired with
fetching a single reference sequence, and finding the matching sites.

```
# note: this still aligns with the full reference alignment
# only the reference sequence mentioned by name will be built
aligned_query, aligned_reference = get_aligned_pair(query_sequence, reference_name,
                                                    reference_alignment)

```

Or, if you prefer to have the full set of newly aligned reference sequences as a dictionary:
```
aligned_query, aligned_refs_by_name = insert_into_alignment(query_sequence, alignment)
```

### From the command line
Brawn can function as a drop-in replacement for MUSCLE's `-profile` mode, complete
command line argument conversion.

To take advantage of cached alignment calculations, first build a cache of your
alignment, with:
```
brawn input.fasta --build-cache desired_cache_path
```

From that point, you can use cached files and plain FASTA files interchangeably.
E.g.
```
brawn query.fasta --reference-alignment some_fasta_file
```
and
```
brawn query.fasta --reference-alignment some_cache_file
```
will both work as expected.

For the command line, the resulting FASTA output will to be stdout.
