usage: AMAS <command> [<args>]

The AMAS commands are:
  concat      Concatenate input alignments
  convert     Convert to other file format
  replicate   Create replicate data sets for phylogenetic jackknife
  split       Split alignment according to a partitions file
  summary     Write alignment summary
  remove      Remove taxa from alignment
  translate   Translate DNA alignment into protein alignment
  trim        Remove columns from alignment

Use AMAS <command> -h for help with arguments of the command of interest

positional arguments:
  command     Subcommand to run

options:
  -h, --help  show this help message and exit
usage: AMAS.py [-h] [-p CONCAT_PART] [-t CONCAT_OUT]
               [-u {fasta,phylip,nexus,phylip-int,nexus-int}]
               [-y {nexus,raxml,unspecified}] [-e] [-c CORES]
               -i IN_FILES [IN_FILES ...]
               -f {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna}

Concatenate input alignments

options:
  -h, --help            show this help message and exit
  -p, --concat-part CONCAT_PART
                        File name for the concatenated alignment partitions.
                        Default: 'partitions.txt'
  -t, --concat-out CONCAT_OUT
                        File name for the concatenated alignment. Default:
                        'concatenated.out'
  -u, --out-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        File format for the output alignment. Default: fasta
  -y, --part-format {nexus,raxml,unspecified}
                        Format of the partitions file. Default: 'unspecified'
  -e, --check-align     Check if input sequences are aligned. Default: no
                        check
  -c, --cores CORES     Number of cores used. Default: 1

required arguments:
  -i, --in-files IN_FILES [IN_FILES ...]
                        Alignment files to be taken as input. You can specify
                        multiple files using wildcards (e.g. --in-files
                        *fasta)
  -f, --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The format of input alignment
  -d, --data-type {aa,dna}
                        Type of data
usage: AMAS.py [-h] [-o SUMMARY_OUT] [-s] [-e] [-c CORES]
               -i IN_FILES [IN_FILES ...]
               -f {fasta,phylip,nexus,phylip-int,nexus-int} -d {aa,dna}

Write alignment summary

options:
  -h, --help            show this help message and exit
  -o, --summary-out SUMMARY_OUT
                        File name for the alignment summary. Default:
                        'summary.txt'
  -s, --by-taxon        In addition to alignment summary, write by
                        sequence/taxon summaries. Default: Don't write
  -e, --check-align     Check if input sequences are aligned. Default: no
                        check
  -c, --cores CORES     Number of cores used. Default: 1

required arguments:
  -i, --in-files IN_FILES [IN_FILES ...]
                        Alignment files to be taken as input. You can specify
                        multiple files using wildcards (e.g. --in-files
                        *fasta)
  -f, --in-format {fasta,phylip,nexus,phylip-int,nexus-int}
                        The format of input alignment
  -d, --data-type {aa,dna}
                        Type of data
