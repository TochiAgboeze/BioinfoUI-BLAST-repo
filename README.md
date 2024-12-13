# NCBI BLAST+ tutorial

Short introduction to using NCBI blast tools from the command line

## Why BLAST on the command line instead of using NCBI's GUI?

Sometimes, you may have to use blast on your own computer to query thousands of
sequences against a custom database of hundreds of thousands of sequences. To
do that, you will need to install Blast on your computer, format the database,
and then blast the sequences.

Below are the advantages of Command-Line BLAST:

**1. Ease of Running Multiple Queries:** with command-line BLAST, you can batch-process hundreds of queries in one script. For example, screening 100 toxin genes across 1,000 bacterial genomes becomes feasible. In GUIs, submitting queries one at a time would be tedious and time-consuming.

**2. Reproducibility and Documentation:** Command-line BLAST allows you to save the exact commands and parameters used in a script, ensuring the analysis can be replicated. In GUIs, parameters may need to be manually set each time, increasing the risk of errors or inconsistencies.

**3. Machine-Readable Results for Analysis:** Results from command-line BLAST can be output in formats like tabular (outfmt 6) or XML, which are easily parsed by downstream analysis tools. GUI results are often formatted for readability, making them less suitable for large-scale automated analyses.

**4. Custom Databases:** You can create and search against custom databases using the makeblastdb command, a feature crucial for studying specific datasets. For example, building a database of Staphylococcus aureus genomes to study MRSA strains resistant to beta-lactams.

**5. Automation:** Command-line BLAST allows automation via shell scripts. For example, you could automate BLAST queries against bacterial genomes uploaded daily to a local database.

**6. Remote Computing:** With command-line BLAST, you can leverage high-performance computing clusters to process large datasets, reducing time and computational load on your local system.


Here is a short tutorial on how to do this.

## Installing Blast+ tools

Get the compiled executables from this URL:

```
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```
e.g
```
curl -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-macosx.tar.gz

```
OR
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-macosx.tar.gz
```

Decompress the archive. For example:

```
tar xvfz ncbi-blast-2.16.0+-x64-macosx.tar.gz
```

Add the `bin` folder from the extracted archive to your path. For example, add
the following line to your `~/.bashrc` or `~/.zhrc` file:

```
export PATH="/PATH/TO/ncbi-blast-2.9.0+/bin":$PATH
```

And change the `/PATH/TO` part to the path where you have put the extracted
archive.

*conda* or *mamba* saves us from all these and we can easily use them to directly install BLAST+ and many other tools subsequently. The following command would install BLAST+ for us:
```
conda install blast
```
```
mamba install blast
```


## Example sequences to use with the tutorial

In order to test blast, you need a test fasta file. Use the following files
that come with the tutorial:

- `EcoliToxins.fasta`
- `EscherichiaDB.fasta`

## Blasting against a remote database
Instead of having to download the entirety of NR or other NCBI databases, we can BLAST against the version held on the website. This ensures we have the most up to date version but is also significantly slower and you may be denied service if you search too many query sequences in a short period of time.
We use the -remote command to do this search against the NCBI website. Lets BLAST our sequences against NR held on the NCBI website by typing:
```
blastn -query EcoliToxins.fasta -db nr/nt -remote -out result.txt -outfmt 6 -evalue 1e-30
```
**Note**: you may have to wait a long time for the remote search to finish, depending on the current server load at NCBI.


## Create blast database

The different blast tools require a formatted database to search against. In
order to create the database, we use the `makeblastdb` tool:

```
makeblastdb -in EscherichiaDB.fasta -title reference -dbtype nucl -out databases/reference
```

This will create a list of files in the `databases` folder. These are all part
of the blast database.

## Blast

We can now blast our sequences against the database. In this case, both our
query sequences and database sequences are DNA sequences, so we use the
`blastn` tool:

```
blastn -db databases/reference -query EcoliToxins.fasta -evalue 1e-3 -word_size 11 -outfmt 6 > sequences.reference
```
Often we don’t need the output alignments but we want all the details of each hit( e.g e-value, bit score, percent identity) on 1 line. Such one-line formats are also preferable for tidy data and further processing with grep or other searching commands and tools. This is achieved by changing the output type. We can do this using the -outfmt flag

You can use different output formats with the `-outmft` option:

```
  -outfmt <String>
   alignment view options:
     0 = Pairwise,
     1 = Query-anchored showing identities,
     2 = Query-anchored no identities,
     3 = Flat query-anchored showing identities,
     4 = Flat query-anchored no identities,
     5 = BLAST XML,
     6 = Tabular,
     7 = Tabular with comment lines,
     8 = Seqalign (Text ASN.1),
     9 = Seqalign (Binary ASN.1),
    10 = Comma-separated values,
    11 = BLAST archive (ASN.1),
    12 = Seqalign (JSON),
    13 = Multiple-file BLAST JSON,
    14 = Multiple-file BLAST XML2,
    15 = Single-file BLAST JSON,
    16 = Single-file BLAST XML2,
    17 = Sequence Alignment/Map (SAM),
    18 = Organism Report
```

## Modifying defaults and output types (evalue and num_alignments)
Often if we are working with many sequences, we want to make it easier to get the best results in an easy to read format. We can do this by limiting the number of results returned. Often this is performed by changing the number of alignments displayed and/or the e-value cut off.

**LIMITING THE NUMBER OF ALIGNMENTS** ensures only the most relevant results are shown, reducing noise and making it easier to identify key findings. For example, setting the limit to display only the top 10 alignments can provide a quick overview of the most relevant matches.

**E-VALUE CUTOFF**: the e-value measures the likelihood of obtaining a particular alignment by chance. A lower e-value threshold (e.g 1e-30) filters out less significant results, ensuring that only highly reliable matches are included in the output.
    - These adjustments enhance the readability and allow for a more focused analysis of sequence alignments e.g:
    
```
    blastn -query EcoliToxins.fasta -db reference -out blast_result.txt -num_alignments 1 -evalue 1e-30
```
    
  The default value (10) is quite high and not recommended as it is too permissive and will likely return false positive hits

## More options and getting help

If you need help to know the options and parameters you can pass `blastn` and
the other blast+ utilities, use the `--help` option and pipe the output into
`less`, for example:

```
blastn --help | less
```

NCBI blast tools cover more cases than DNA against DNA searches. For example,
you can search a protein database with either DNA or protein sequences. Here is
an exhaustive list of the programs that come with the blast+ distribution:

```
blast_formatter            blastn                     deltablast                 makeprofiledb              tblastn_vdb
blast_formatter_vdb        blastn_vdb                 dustmasker                 psiblast                   tblastx
blast_vdb_cmd              blastp                     get_species_taxids.sh      rpsblast                   update_blastdb.pl
blastdb_aliastool          blastx                     legacy_blast.pl            rpstblastn                 windowmasker
blastdbcheck               cleanup-blastdb-volumes.py makeblastdb                segmasker
blastdbcmd                 convert2blastmask          makembindex                tblastn
```

## Other Examples on Performing BLAST
Please replace the query and database/subject names with names of the files you have
**BLASTn (blastn -help)**

```
blastn -db nt -query nt.fasta -out result.txt
```
```
blastn -query genes.fasta -subject genomes.fasta -out results.txt
```
**Task types:** MEGABLAST vs BLASTN
```-task megablast``` is used for comparison within the same species/ highly similar sequences(it is faster)
```
blastn -query genes.fasta -subject genomes.fasta -out results.txt -task megablast
```
```-task blastn``` is for comparison of sequences from different species/somewhat similar (it is slow but finds more hits)
```
blastn -query genes.fasta -subject genomes.fasta -out results.txt -task blastn
```
**BLASTp**(blastp -help)
```
blastp -query x.faa -db dbname -out results.tsv
```
```
blastp -query x.faa -subject y.faa -out results.tsv
```
**BLASTx** (mapping nucleotides against protein database)
```
blastx -query nucl.faa -subject prot.faa -out result.txt
```

## DATABASES
**CARD**
```
https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2
```
For the CARD database when downloaded, please take note of the following while taking a look at the different files within the extracted directory:

*FASTA:*
Nucleotide and corresponding protein FASTA downloads are available as separate files for each model type.  For example, the `protein homolog model` type contains sequences of antimicrobial resistance genes that do not include mutation as a determinant of resistance- these data are appropriate for BLAST analysis of metagenomic data or searches excluding secondary screening for resistance mutations. In contrast, the `protein variant model` includes reference wild type sequences used for mapping SNPs conferring antimicrobial resistance - without secondary mutation screening, analyses using these data will include false positives for antibiotic resistant gene variants or mutants.


*INDEX FILES:*

The file `aro_index.tsv` contains a list of ARO tagging of GenBank accessions stored in CARD.


`aro_categories.tsv`: reflects AMR gene family, target drug class, and mechanism of resistance.

`aro_categories_index.tsv`: reflects AMR gene family, target drug class, and mechanism of resistance, so GenBank accessions may have more than one cross-reference.

**NCBI Database**
https://ftp.ncbi.nlm.nih.gov

