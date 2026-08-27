# Tarefa 1 — Alinhamento de reads

## Organismo

- **Nome científico:** *Pitangus sulphuratus* (bem-te-vi)
- **Run:** `ERR17582071`
- **Tipo de sequenciamento:** WGS paired-end
- **Subset utilizado:** 5 spots

## Referência

- **Assembly:** `GCA_029286575.1_ASM2928657v1`
- **FASTA:** `GCA_029286575.1_ASM2928657v1_genomic.fna`

### Download das reads

```bash
mamba create -n sra-tools -c conda-forge -c bioconda sra-tools
mamba activate sra-tools

fastq-dump --split-files -X 5 ERR17582071

mamba deactivate
```

### Alinhamento

```bash
mamba create -n minimap2 -c conda-forge -c bioconda minimap2 samtools
mamba activate minimap2

minimap2 -d Pitangus_sulphuratus.mmi \
  GCA_029286575.1_ASM2928657v1_genomic.fna

samtools faidx GCA_029286575.1_ASM2928657v1_genomic.fna

minimap2 -ax sr \
  Pitangus_sulphuratus.mmi \
  ERR17582071_1.fastq \
  ERR17582071_2.fastq | \
samtools sort -o alinhamento_ERR17582071.bam

samtools index alinhamento_ERR17582071.bam
```

## Arquivos grandes

(https://drive.google.com/drive/folders/1pEP3cjxZ8WJUcEzd2Kppi2nfOzLeey7H?usp=sharing)
