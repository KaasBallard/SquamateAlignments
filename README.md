# Purpose
The purpose of this project is to soft-mask a number of snake genomes, and potentially other groups of reptiles for future ProgressiveCactus alignments.

# Installing programs
You can install the necessary programs the hard way, which I am told takes a day, or you can install in a couple of minutes using ``mamba``. You do need to make sure you have ``python 3.12`` rather than ``3.13`` since there is a package ``RepeatModeler`` relies on that doesn't have a ``3.13`` version yet. Because of this you will get a very old version of ``RepeatModeler`` from 2017 (1.0.7) when the most up to date version available normally and through ``mamba`` is 2.0.6 from late last year.

To install the correct version of python, just run the following:
```
mamba create -n RepeatMaskAnnot

mamba activate RepeatMaskAnnot

mamba install python=3.12
```

Once that is installed you can install ``RepeatModeler``:
```
mamba install repeatmodeler
```
This should install version 2.0.6. In my experience, this will also install an the second newest version of ``RepeatMasker`` version 4.1.8. There is a newer version, 4.1.9, but it came out about a week before the time of writing, meaning it will be released on ``mamba`` later.
