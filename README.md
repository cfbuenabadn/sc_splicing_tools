# Miscellaneous tools for alternative splicing

## GTF2psix.py

### Requirements

Requires [gtfparse](https://pypi.org/project/gtfparse/) and [tqdm](https://tqdm.github.io/), in addition to common Anaconda modules NumPy, Pandas, Argparse and gzip.

### Description

This script creates a Psix-compatible annotation of cassette exons and constitutive introns directly from a GTF file. This annotation consists of a table specifying the location (chromosome, start and end) of splice junctions. Splice junctions are annotated as supporting the inclusion of a cassette exon (\_I1 and \_I2), supporting its exclusion (\_SE), or constitutive (\_CI). You can download ready-to-use mouse (mm10) and human (hg38) annotations [here](https://github.com/lareaulab/psix/tree/master/annotation). 

The annotation is a table file with the following format:

name | intron | event  |  gene
---- | ---- | ---- | ---- 
Rpn2_1_I1 | chr2:157318043-157320120:+ | Rpn2_1 | Rpn2
Rpn2_1_I2 | chr2:157320197-157321742:+ | Rpn2_1 | Rpn2
Rpn2_1_SE | chr2:157318043-157321742:+ | Rpn2_1 | Rpn2
Rpn2_2_CI | chr2:157290713-157294152:+ | Rpn2_2 | Rpn2
... | ... | ... | ...

Each row corresponds to a splice junction or intron. The first column is a name assigned to the splice junction, that is based on the name of the gene that contains the removed intron, and whether it supports the inclusion of a cassette exon (\_I1 and \_I2), supports its exclusion (\_SE), or it is a constitutive intron (\_CI). The second column are the intron coordinates in the genome. The third column has the name of the splicing element: for example, Rpn2_1 is a cassette exon, and it is formed by three introns: Rpn2_1_I1, Rpn2_1_I2 and Rpn2_1_SE. The fourth column contains the gene name.

### Runing GTF2psix

To automatically create an annotation from a GTF file, download ```GTF2psix.py``` and run as follows

```bash
python GTF2psix.py --gtf annotation.gtf -o psix_annotation
```

### Optional arguments

```--gene_name``` specifies the tag for gene names to use in the GTF file. For example ```--gene_name gene_name``` will use the ```gene_name``` tag from the 9th column of the GTF file. The default is ```gene_id```.

```--gene_type_tag```. Some GTF files have different tags for the gene types. E.g., ```gene_type``` or ```gene_biotype```. Specify the tag with this option. Default: ```gene_type```.

```--transcript_type_tag```. Some GTF files have different tags for the transcript types. E.g., ```transcript_type``` or ```transcript_biotype```. Specify the tag with this option. Default: ```transcript_type```.

Use the argument ```--gene_type``` to limit the annotation to a specific type of genes. E.g., for an annotation of protein coding genes only, use ```--gene_type protein_coding```.

You can remove some chromosomes from the annotation using the argument ```--exclude_chromosome```. E.g., ```--exclude_chromosome chrM,chrY``` will exclude chrM and chrY from the annotation.
