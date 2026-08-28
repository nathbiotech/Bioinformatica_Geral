# Tarefa 1 — Do sequenciamento ao alinhamento

## Organismo

- **Nome científico:** *Mycoplasma genitalium*
- **Run:** `SRR39974648`
- **Dados:** reads paired-end

## Genoma de referência

- **Assembly:** `GCF_040556925.1_ASM4055692v1`
- **Arquivo FASTA:** `GCF_040556925.1_ASM4055692v1_genomic.fna`

### Download e conversão das reads

```bash
mamba activate sra-tools

prefetch SRR39974648

fasterq-dump --split-files SRR39974648
```

### Indexação da referência

```bash
mamba activate minimap2

minimap2 -d M_genitalium.mmi \
  GCF_040556925.1_ASM4055692v1_genomic.fna
```

### Alinhamento

```bash
minimap2 -ax sr \
  M_genitalium.mmi \
  SRR39974648_1.fastq \
  SRR39974648_2.fastq \
  > alinhamento_SRR39974648.sam
```

### Ordenação e indexação do BAM

```bash
picard SortSam \
  I=alinhamento_SRR39974648.sam \
  O=alinhamento_SRR39974648.bam \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true
```

## Arquivos restantes

[Google Drive — arquivos de entrada e saída](https://drive.google.com/drive/folders/1pEP3cjxZ8WJUcEzd2Kppi2nfOzLeey7H?usp=sharing)
