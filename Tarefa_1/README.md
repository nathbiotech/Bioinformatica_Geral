wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh

mamba create -n sra-tools -c conda-forge -c bioconda sra-tools

mamba activate sra-tools

fastq-dump --split-files -X 5 ERR17582071

mamba deactivate

mamba create -n minimap2 -c conda-forge -c bioconda minimap2 samtools

mamba activate minimap2

minimap2 -d Pitangus\_sulphuratus.mmi   
GCA\_029286575.1\_ASM2928657v1\_genomic.fna

samtools faidx GCA\_029286575.1\_ASM2928657v1\_genomic.fna

minimap2 -ax sr   
Pitangus\_sulphuratus.mmi   
ERR17582071\_1.fastq   
ERR17582071\_2.fastq |   
samtools sort -o alinhamento\_ERR17582071.bam

samtools index alinhamento\_ERR17582071.bam



\# arquivos restantes

\# https://drive.google.com/drive/folders/1pEP3cjxZ8WJUcEzd2Kppi2nfOzLeey7H?usp=sharing

