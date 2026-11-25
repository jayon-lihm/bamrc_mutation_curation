# bamrc_mutation_curation

Screen BAM/CRAM files for possible mutations in a small genomic region of interest, using raw read counts instead of variant calls.

Use it when mutation calls are missing, or when you want to check the reads at specific positions yourself (for example, KRAS hotspots). The pipeline:

1. runs **bam-readcount** on every sample over the regions in a BED file,
2. merges each patient's read counts into one **multi-sample VCF**,
3. annotates that VCF with **SnpEff**,
4. turns the annotated VCF into a **tab-separated table** with one row per sample × alt allele × annotation.

```
samplesheet.csv ─► BAMREADCOUNT ─► BAMRC_TO_VCF ─► SNPEFF_ANNOTATION ─► SNPEFF_VCF_TO_TABLE
                   (per sample)     (per patient)    (per patient)         (per patient)
```

---

## Requirements

| Tool | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) | DSL2 |
| [Singularity / Apptainer](https://apptainer.org/) | must be on `PATH` as `singularity` |
| Python ≥ 3.6 | standard library only, no extra packages |
| bam-readcount image | e.g. built from `docker://mgibio/bam-readcount:1.0.0` |
| SnpEff 4.3t image | e.g. `docker pull lihao5657/snpeff:4.3t`; the jar must be at `/snpEff/snpEff.jar` inside the image |
| SnpEff database | e.g. `GRCh38.86`, downloaded to a local data directory |
| Reference FASTA | with a `.fai` index next to it |

### Quick setup

[`setup.sh`](setup.sh) creates a conda environment with Nextflow and Python, and pulls both images:

```bash
bash setup.sh                      # env "nf-env", images saved to ./images
bash setup.sh -n my-env -d /path/to/images
```

It skips anything that already exists and finishes by printing the `bamrc_sif` / `snpeff_sif` values to put in your config. It needs `conda` (or `mamba`) and `singularity` (or `apptainer`) on `PATH`.

To pull the images by hand instead:

```bash
singularity pull bam-readcount_1.0.0.sif docker://mgibio/bam-readcount:1.0.0
singularity pull snpeff_4.3t.sif         docker://lihao5657/snpeff:4.3t
```

Downloading the SnpEff database (example):

```bash
singularity exec --bind /path/to snpeff_4.3t.sif \
    java -jar /snpEff/snpEff.jar download -dataDir /path/to/snpeff/4.3t/data GRCh38.86
```

---

## 1. Set the paths in the config

[`nextflow.config`](nextflow.config) ships with `/path/to/...` placeholders. The pipeline will not run until you change them.

```groovy
params {
    genomes = [
        hg38: [
            fasta : '/path/to/Homo_sapiens_assembly38.fasta',   // reference FASTA (.fai next to it)
            snpeff: 'GRCh38.86',                                // SnpEff database for this genome
        ],
    ]
    bamrc_sif      = '/path/to/bam-readcount_1.0.0.sif'         // bam-readcount image
    bind_paths     = '/path/to'                                 // host folders mounted into the containers
    snpeff_sif     = '/path/to/snpeff_4.3t.sif'                 // SnpEff image
    snpeff_datadir = '/path/to/snpeff/4.3t/data'                // SnpEff data directory
}
```

| Parameter | What to set |
|---|---|
| `genomes.hg38.fasta` | Full path to the reference FASTA. A `<fasta>.fai` must sit next to it. |
| `bamrc_sif` | Full path to the bam-readcount Singularity image. |
| `snpeff_sif` | Full path to the SnpEff Singularity image. |
| `snpeff_datadir` | The SnpEff data directory, i.e. the folder that holds `GRCh38.86/`. |
| `bind_paths` | Comma-separated host folders mounted into the containers. They must cover the FASTA, the BAM/CRAM files, the SnpEff data directory and the Nextflow work directory, e.g. `'/data,/scratch'`. |

### Recommended: keep your site paths in a separate config

Leave `nextflow.config` with placeholders and put your real paths in a small override file, for example `site.config`:

```groovy
// site.config: local paths, not committed
params {
    genomes.hg38.fasta = '/data/ref/hg38/Homo_sapiens_assembly38.fasta'
    bamrc_sif          = '/data/images/bam-readcount_1.0.0.sif'
    snpeff_sif         = '/data/images/snpeff_4.3t.sif'
    snpeff_datadir     = '/data/snpeff/4.3t/data'
    bind_paths         = '/data'
}
```

Pass it with `-c`. Nextflow always reads `nextflow.config` first, and values in the `-c` file take precedence:

```bash
nextflow run main.nf -c site.config --samplesheet samplesheet.csv --outdir results_dir
```

Add the file to `.gitignore` so your paths are never committed.

### Or override on the command line

Every parameter can be passed as `--<name>`:

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --bamrc_sif /data/images/bam-readcount_1.0.0.sif \
    --snpeff_sif /data/images/snpeff_4.3t.sif \
    --snpeff_datadir /data/snpeff/4.3t/data \
    --bind_paths /data
```

(Nested values such as `genomes.hg38.fasta` are easier to set in a config file.)

### Adding another reference genome

Add an entry to `genomes`, then use its name in the samplesheet's `reference` column:

```groovy
params {
    genomes = [
        hg38: [ fasta: '/data/ref/hg38.fasta', snpeff: 'GRCh38.86' ],
        hg19: [ fasta: '/data/ref/hg19.fasta', snpeff: 'GRCh37.75' ],
    ]
    genome_aliases = [GRCh38: 'hg38', GRCh37: 'hg19']   // other names accepted (case-insensitive)
}
```

Use a region BED file with coordinates for the same genome. All samples of one patient must use the same reference.

---

## 2. Prepare the samplesheet

CSV with a header row:

| Column | Required | Description |
|---|---|---|
| `patient` | yes | Patient ID. Samples with the same ID are merged into one VCF/table. |
| `sample_name` | yes | Unique sample name. Becomes the VCF sample column. |
| `status` | yes | `tumor` or `normal` |
| `bam_file` | yes | Full path to the `.bam` or `.cram`. The index (`.bai`/`.crai`, either `x.bam.bai` or `x.bai`) must sit next to it. |
| `reference` | no | Genome name from `params.genomes` (e.g. `hg38`). If blank, `params.genome` is used. |

Example ([`example/samplesheet.csv`](example/samplesheet.csv)):

```csv
patient,sample_name,status,bam_file,reference
Pt1,Pt1_Tumor,tumor,/path/to/Pt1_Tumor_markdup_recal.cram,hg38
Pt1,Pt1_Normal,normal,/path/to/Pt1_Normal_markdup_recal.cram,hg38
```

### Samplesheet builder

Open [`samplesheet_builder.html`](samplesheet_builder.html) in a browser to build the CSV in a form. You can enter a base directory once and type only file names, the page checks required fields and file extensions, and it exports the result as `samplesheet.csv`. Nothing is uploaded: it runs entirely in the browser.

---

## 3. Choose the region

`--region` takes a BED file (`chrom  start  end`) of the positions to count. The default is [`data/KRAS_region_hg38.bed`](data/KRAS_region_hg38.bed), which covers the KRAS G12/G13, Q61 and A146 codons (hg38).

bam-readcount reads the region file as **1-based, inclusive** coordinates, unlike standard 0-based BED. Each line in the default file is one base (`start == end`).

---

## 4. Run

```bash
nextflow run main.nf \
    -c site.config \
    --samplesheet samplesheet.csv \
    --region data/KRAS_region_hg38.bed \
    --outdir /path/to/output
```

[`run_nextflow.sh`](run_nextflow.sh) is a template wrapper for this command.

On SLURM, add `-profile slurm`. Jobs go to the cluster's default partition; to choose one, add `process.queue = '<partition>'` to your site config.

Resume after a failure with `-resume`.

### All parameters

| Parameter | Default | Description |
|---|---|---|
| `samplesheet` | `./samplesheet.csv` | Input samplesheet |
| `region` | `data/KRAS_region_hg38.bed` | Regions to count |
| `outdir` | `.` (current directory) | Output root. `results/`, `work/` and `pipeline_info/` go here. |
| `genome` | `hg38` | Genome used when the samplesheet `reference` column is blank |
| `genomes` | hg38 entry | Reference FASTA and SnpEff database per genome |
| `genome_aliases` | `[GRCh38: 'hg38']` | Extra names accepted in the `reference` column |
| `bamrc_sif` | placeholder | bam-readcount image |
| `snpeff_sif` | placeholder | SnpEff image |
| `snpeff_datadir` | placeholder | SnpEff data directory |
| `snpeff_genome` | `null` | Force one SnpEff database for every patient, ignoring `genomes` |
| `bind_paths` | placeholder | Host folders mounted into containers |
| `min_base_qual` | `15` | bam-readcount `-b` (minimum base quality) |
| `max_warnings` | `1` | bam-readcount `-w` |
| `python` | `python3` | Python used for the `src/*.py` scripts |

---

## Output

```
<outdir>/
├── results/
│   ├── bamreadcount/<patient>/<sample>.bamrc.tsv          raw bam-readcount output
│   ├── vcf/<patient>.bamrc.vcf                           multi-sample VCF
│   ├── snpeff/<patient>.bamrc_ann.vcf                    SnpEff-annotated VCF
│   ├── snpeff/<patient>.bamrc.fasta                      protein sequences (SnpEff -fastaProt)
│   └── protein_coding_table/<patient>.bamrc_ann.protein_coding.tsv
├── pipeline_info/   trace.txt, report.html
└── work/            Nextflow intermediate files (safe to delete after a run)
```

### VCF ([`src/convert_bamrc_to_vcf.py`](src/convert_bamrc_to_vcf.py))

- One VCF per patient, with one sample column per sample, in samplesheet order.
- **Alt alleles** are every non-reference allele with a read count > 0 in at least one of the patient's samples. `=` and `N` are ignored.
- Sites where every sample is 100% reference are left out.
- Indels are reported from bam-readcount `-i` and normalised to VCF style. Deletions are anchored on the preceding reference base, which is why the FASTA `.fai` is needed.
- `ID` = `chrom_pos_ref_alt`. `INFO` is empty.
- `FORMAT` = `ALT:GT:DP:AD:VAF`. When there are several alt alleles, the per-allele fields are joined with `;`:

| Field | Meaning |
|---|---|
| `ALT` | alt allele(s) |
| `GT` | `1/0` if the sample has reads for that allele, `0/0` if not, `./.` if there is no data at the site |
| `DP` | total depth reported by bam-readcount |
| `AD` | `ref_count,alt_count` for that allele |
| `VAF` | `alt / (ref + alt)` |

### Table ([`src/snpeff_vcf_to_table.py`](src/snpeff_vcf_to_table.py))

One row per sample × alt allele × SnpEff annotation. The `structural_interaction_variant`, `sequence_feature` and `intragenic_variant` annotations are dropped.

| Column | Description |
|---|---|
| `patient`, `sample` | IDs |
| `chrom`, `pos`, `ref`, `alt` | variant |
| `mutation_type` | SnpEff effect, e.g. `missense_variant`, `stop_gained` |
| `gene` | gene name |
| `protein_change` | HGVS.p, e.g. `p.G12D` |
| `raw_depth` | total depth recorded by bam-readcount (`DP`) |
| `depth` | `ref_count + alt_count` |
| `ref_count`, `alt_count` | from `AD` for this alt allele |
| `vaf` | `alt_count / depth` |

> These are raw read counts with no variant calling or error model. A low-VAF alt allele can be a sequencing or alignment artefact. Check the depth, compare with the matched normal, and review candidates in IGV.
