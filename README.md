Installation:

# Installing 

Firstly, ensure that you have singularity installed (usually available on HPC systems, as a module). Follow the guide here: https://docs.sylabs.io/guides/3.0/user-guide/installation.html
Additionally, you will need to install conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Then, clone this repo inside your preferred location, then `cd` into it.

Install the conda environment: `fpn_env`, containing the packages required to run FindingPseudoNemo:

```
conda env create -f environment.yml
```

Activate the environment (be sure to always do this before runnning this program)

```
conda activate fpn_env
```

# Installing VEP

Go to `config.env` and set VEP_SIF (where you want the vep image to be installed) and VEP_CACHE (where vep_cache files will be stored).
Otherwise you can keep the default values set there

Install the VEP image:

```
singularity pull $VEP_SIF docker://ensemblorg/ensembl-vep
```

The command above will install the latest version (currently version 116.2). We used version version 115.2 but you are free to specify your version e.g:

```
singularity pull $VEP_SIF docker://ensemblorg/ensembl-vep:release_115.2
```

Then download the GRCh38 refseq cache/FASTA files:
 
```
singularity exec ~/Desktop/singularity_images/vep.sif INSTALL.pl \
            -c $VEP_CACHE \
            -a cf \
            -s homo_sapiens_refseq \
            -y GRCh38
```

Additionally, you will need the following plugin files (for SpliceAI, AlphaMissense and REVEL respectively):

* spliceai_scores.raw.snv.hg38.vcf.gz
* spliceai_scores.raw.indel.hg38.vcf.gz
* AlphaMissense_hg38.tsv.gz
* new_tabbed_revel_grch38.tsv.gz

Follow the guide here to set them up: https://jun2026.archive.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#alphamissense

# Installing References

Finally, download HG38, mask the homologous IKBKGP1 region and index the genome:

`./helper_scripts/install_and_set_up_hg38_index.sh`

# Running the Program:

Firstly, run the SNV detection script:

`./main_workflow.sh PATH_TO_INPUT PATH_TO_OUTPUT_FOLDER`

If you have a female WGS sample, you can also run the deletion detection module, to see if there is a potential exon 4-10 deletion, located between the mer67b repeats:

