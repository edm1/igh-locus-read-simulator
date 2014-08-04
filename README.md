igh-locus-read-simulator
========================

Simulates single-end MiSeq IGH locus reads


### Dependancies
- Python v2.7
- Numpy
- [ART Simulation Tools](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) are required and can be downloaded using `$ bash install_ART.sh`

### Download
The repository can be downloaded using git `git clone https://github.com/edm1/igh-locus-read-simulator.git` or by following the *Download ZIP* link on the right.

### Simulate reads
The script `simulate_reads.py` is used to simulate VDJ recombination and produce MiSeq-like amplicon reads:

```
python simulate_reads.py [-h] --OutPrefix <str> [--NumReads <int>]
                         [--Vgenes <fasta>] [--Dgenes <fasta>]
                         [--Jgenes <fasta>] [--J5del <int>] [--D3del <int>]
                         [--D5del <int>] [--V3del <int>] [--VDins <int>]
                         [--DJins <int>] [--TotalReadLen <int>]
                         [--ARTbin <str>]
                         N [N ...]

```

View full list of arguments using `python 1_run_IlluminaBasecallsToFastq.py --help`

##### Required
- `--OutPrefix` - Output prefix for fasta, fastq and log files.
- `N [N ...]` - A list of integers stating what percentage of reads should be in each of the "leukemic" clones. For example, 20 15 5 would produce 3 main clones containing 20%, 15% and 5% of the reads each, the remaining 60% of the reads would be unique.

##### Optional
For a list of all optional arguments type `python simulate_reads.py --help`

#### Ouput
This script outputs:
- a fasta containing the actual clone sequences.
- a fastq containing simulated Illumina reads for each of the clones.
- a log file stating the composition of each clone.

### Notes
Deletion and insertion sizes are sampled from a Poisson distribution. More work needs to be done to see if this accurately reflects the biology but for our purposes it should do fine.
