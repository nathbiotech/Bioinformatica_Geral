# Análise de scRNA-seq de quatro tecidos de camundongo

Pipeline de análise de **RNA-seq de célula única (scRNA-seq)**
utilizando o pacote **Seurat** para integrar e caracterizar células de
quatro tecidos de camundongo: olho, ovário, estômago e intestino
delgado.

O conjunto de dados utilizado é disponibilizado pela 10x Genomics:

-   [40k Mixture of Cells Dissociated from 4 Fixed Mouse
    Tissues](https://www.10xgenomics.com/datasets/40k-mixture-of-cells-dissociated-from-4-fixed-mouse-tissues-using-manual-dissociation-multiplexed-samples-4-probe-barcodes-1-standard)

## Objetivo

Processar os dados de expressão gênica de quatro tecidos em um único
objeto Seurat, realizar controle de qualidade, normalização, seleção de
genes variáveis, integração entre amostras, redução de dimensionalidade,
clusterização e anotação dos tipos celulares.

Ao final, o pipeline gera visualizações para explorar os agrupamentos e
os marcadores moleculares associados às linhagens celulares.

## Estrutura esperada dos arquivos

Disponibilize os arquivos `.h5` no diretório de trabalho do R:

``` text
.
├── eye.h5
├── ovary.h5
├── stomach.h5
├── intestine.h5
└── README.md
```

Os nomes dos arquivos devem corresponder aos utilizados no script.

## Requisitos

O código utiliza R e os seguintes pacotes:

-   [Seurat](https://satijalab.org/seurat/)
-   [ggplot2](https://ggplot2.tidyverse.org/)
-   [dplyr](https://dplyr.tidyverse.org/)

Instale os pacotes, caso ainda não estejam disponíveis:

``` r
install.packages(c("Seurat", "ggplot2", "dplyr"))
```

## Etapas da análise

1.  **Importação e criação dos objetos Seurat:** leitura das matrizes
    HDF5 de cada tecido e criação de um objeto por amostra.
2.  **Junção das amostras:** combinação dos quatro objetos em um único
    objeto Seurat, mantendo a identificação de origem de cada célula.
3.  **Controle de qualidade (QC):** cálculo da porcentagem de genes
    mitocondriais e filtragem por número de genes detectados e
    porcentagem mitocondrial.
4.  **Normalização:** aplicação de `LogNormalize`, com fator de escala
    de 10.000.
5.  **Seleção de genes variáveis:** identificação de 2.000 genes
    utilizando o método `vst`.
6.  **Escalonamento e PCA:** escalonamento dos dados e redução linear de
    dimensionalidade.
7.  **Integração:** integração das camadas com CCA (`IntegrateLayers`),
    seguida de `JoinLayers`.
8.  **Clusterização:** construção do grafo de vizinhos e identificação
    dos clusters.
9.  **UMAP:** visualização dos agrupamentos celulares e da distribuição
    por tecido.
10. **Identificação e anotação:** busca de genes marcadores e atribuição
    de identidades biológicas aos clusters.
11. **Visualizações finais:** geração de Dot Plot, Heatmap e UMAPs para
    apresentação dos resultados.

## Parâmetros principais

  Etapa               Parâmetro
  ------------------- ------------------------------------------------
  Criação do objeto   `min.cells = 3`, `min.features = 200`
  QC                  `nFeature_RNA > 200` e `nFeature_RNA < 6000`
  QC mitocondrial     `percent.mt < 15`
  Normalização        `LogNormalize`, `scale.factor = 10000`
  Genes variáveis     `selection.method = "vst"`, `nfeatures = 2000`
  Integração          `CCAIntegration`
  Clusterização       `dims = 1:15`, `resolution = 0.4`
  Marcadores          `min.pct = 0.25`, `logfc.threshold = 0.25`

Os parâmetros são os definidos no script e podem precisar de avaliação
conforme a qualidade e as características dos dados.

## Como executar

1.  Baixe o conjunto de dados no site da 10x Genomics.
2.  Disponibilize os quatro arquivos HDF5 (`eye.h5`, `ovary.h5`,
    `stomach.h5` e `intestine.h5`) no diretório de trabalho.
3.  Abra o script R no RStudio ou em outro ambiente R.
4.  Execute o pipeline em sequência.
5.  Inspecione os gráficos e as tabelas intermediárias, especialmente
    antes de confirmar a anotação dos clusters.

## Anotação dos tipos celulares

O script utiliza genes marcadores como apoio à interpretação biológica.
Entre os exemplos incluídos estão:

  Linhagem / tipo celular      Marcadores usados no script
  ---------------------------- -------------------------------------------
  Epiteliais                   `Epcam`, `Krt8`, `Krt18`, `Cdh1`, `Krt20`
  Esteroidogênicas do ovário   `Hsd3b1`, `Cyp17a1`, `Star`, `Inha`
  Retina / neuronal            `Pde6b`, `Rho`, `Sag`, `Rlbp1`
  Córnea / cristalino          `Krt12`, `Crygs`, `Pax6`
  Fibroblastos                 `Col1a1`
  Músculo liso / pericitos     `Acta2`, `Myh11`, `Rgs5`
  Células imunes               `Ptprc`, `Ccl5`, `Cd74`
  Endoteliais                  `Pecam1`, `Cdh5`
  Proliferativas               `Mki67`
  Plasmócitos                  `Jchain`
  Glia / suporte               `Plp1`
  Mesoteliais                  `Msln`
  Adipócitos                   `Adipoq`

**Importante:** os nomes atribuídos aos clusters no script são hipóteses
de anotação e devem ser validados com os marcadores, a distribuição por
tecido e os resultados de expressão. Em particular, os clusters 0 e 1
são destacados no próprio código para diagnóstico antes da anotação
definitiva. Não interprete os rótulos provisórios como resultados
biologicamente confirmados.

## Diagnóstico dos clusters 0 e 1

O script inclui uma etapa específica de inspeção antes de consolidar a
anotação:

-   Avalia a proporção de células de cada tecido em cada cluster.
-   Compara métricas de qualidade (`nFeature_RNA`, `nCount_RNA` e
    `percent.mt`) entre clusters.
-   Inspeciona os principais genes marcadores dos clusters 0 e 1.
-   Compara marcadores candidatos de epitélio, ovário, retina, córnea e
    cristalino por meio de um Dot Plot.

A interpretação sugerida no código é:

-   Se houver expressão de `Epcam`/`Krt8` e presença em vários órgãos,
    considerar uma identidade epitelial compartilhada.
-   Se houver expressão de `Hsd3b1`/`Cyp17a1`/`Star` e associação ao
    ovário, considerar células esteroidogênicas.
-   Se os indicadores apontarem baixa qualidade, considerar a
    possibilidade de células de baixa qualidade ou RNA ambiente.

Essas possibilidades devem ser verificadas nos resultados; não são
conclusões automáticas do pipeline.

## Arquivos de saída

Ao executar as etapas finais, o script salva as seguintes figuras no
diretório de trabalho:

  -----------------------------------------------------------------------
  Arquivo                             Conteúdo
  ----------------------------------- -----------------------------------
  `DotPlot_Final.png`                 Expressão dos principais marcadores
                                      por linhagem celular

  `Heatmap_Final.png`                 Assinatura transcricional dos tipos
                                      celulares anotados

  `UMAP_Anotado_Final.png`            UMAP colorido pelos tipos celulares
                                      anotados

  `UMAP_Por_Orgao_Final.png`          Distribuição dos tipos celulares
                                      separada por órgão
  -----------------------------------------------------------------------

## Observações

-   A porcentagem mitocondrial é calculada com o padrão `^mt-`, conforme
    a nomenclatura utilizada no script para genes mitocondriais de
    camundongo.
-   A integração utiliza CCA e a redução `integrated.cca` para as etapas
    de vizinhança, clusterização e UMAP.
-   O número de dimensões e a resolução de clusterização são escolhas
    analíticas e devem ser avaliadas com os gráficos e a estrutura dos
    dados.
-   A lista `novos_nomes` precisa corresponder aos IDs de cluster
    efetivamente gerados. O script apresentado contém rótulos para os
    clusters de 0 a 22.
-   A anotação celular requer validação biológica; a expressão de um
    único marcador, isoladamente, não deve ser considerada confirmação
    suficiente.

## Reprodutibilidade

O script utiliza `set.seed(42)` antes da geração do Heatmap com
subamostragem. Para facilitar a reprodutibilidade, registre também as
versões de R e dos pacotes utilizados em cada execução.

## Referência dos dados

10x Genomics. *40k Mixture of Cells Dissociated from 4 Fixed Mouse
Tissues Using Manual Dissociation, Multiplexed Samples, 4 Probe
Barcodes, 1 Standard.*

[Consultar dataset na 10x
Genomics](https://www.10xgenomics.com/datasets/40k-mixture-of-cells-dissociated-from-4-fixed-mouse-tissues-using-manual-dissociation-multiplexed-samples-4-probe-barcodes-1-standard)

------------------------------------------------------------------------

## Script completo

O código abaixo reúne as etapas de processamento, integração,
clusterização, diagnóstico, anotação e geração das figuras descritas
neste README.

``` r
library(Seurat)
library(ggplot2)
library(dplyr)

# ==============================================================================
# 1. MATRIZ DE CONTAGENS (Importação e criação do objeto Seurat)
# ==============================================================================
counts_eye       <- Read10X_h5("eye.h5")
counts_ovary     <- Read10X_h5("ovary.h5")
counts_stomach   <- Read10X_h5("stomach.h5")
counts_intestine <- Read10X_h5("intestine.h5")

seu_eye       <- CreateSeuratObject(counts = counts_eye, project = "Eye", min.cells = 3, min.features = 200)
seu_ovary     <- CreateSeuratObject(counts = counts_ovary, project = "Ovary", min.cells = 3, min.features = 200)
seu_stomach   <- CreateSeuratObject(counts = counts_stomach, project = "Stomach", min.cells = 3, min.features = 200)
seu_intestine <- CreateSeuratObject(counts = counts_intestine, project = "Small_Intestine", min.cells = 3, min.features = 200)

seurat_obj <- merge(
  x = seu_eye,
  y = list(seu_ovary, seu_stomach, seu_intestine),
  add.cell.ids = c("Eye", "Ovary", "Stomach", "Intestine"),
  project = "Multi_Tissue_scRNA"
)

# ==============================================================================
# 2. CONTROLE DE QUALIDADE (QC) E FILTRAGEM
# ==============================================================================
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Gráfico de qualidade antes da filtragem
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15
)

# ==============================================================================
# 3. NORMALIZAÇÃO
# ==============================================================================
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ==============================================================================
# 4. GENES ALTAMENTE VARIÁVEIS (Feature Selection)
# ==============================================================================
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# ==============================================================================
# 5. SCALING (Escalonamento)
# ==============================================================================
seurat_obj <- ScaleData(seurat_obj)

# ==============================================================================
# 6. PCA (Redução de dimensionalidade linear)
# ==============================================================================
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Gráfico 2D da PCA por amostra e ElbowPlot para escolher o número de dimensões (PCs)
DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident")
ElbowPlot(seurat_obj, ndims = 30)

# ==============================================================================
# OBS: INTEGRAÇÃO (Remoção de batch effects entre as amostras/tecidos)
# ==============================================================================
# No Seurat v5, como fizemos merge(), cada tecido é uma "layer".
# O comando abaixo alinha tipos celulares compartilhados (ex: células imunes, endoteliais) entre os órgãos:
seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# Junta as camadas após a integração para permitir a busca de genes marcadores
seurat_obj <- JoinLayers(seurat_obj)

# ==============================================================================
# 7. CLUSTERIZAÇÃO
# ==============================================================================
# Usamos a redução integrada ("integrated.cca") e as primeiras 15 PCs
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.cca", dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)

# ==============================================================================
# 8. UMAP (Redução de dimensionalidade não-linear)
# ==============================================================================
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.cca", dims = 1:15)

# Gráficos UMAP: um colorido por cluster e outro colorido por órgão/tecido
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")

# ==============================================================================
# 9. IDENTIFICAÇÃO DE GENES MARCADORES DE CLUSTERS
# ==============================================================================
marcadores <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Visualizando no console os top 5 genes de cada cluster para você identificar os tipos celulares
top5_tabela <- marcadores %>%
  group_by(cluster) %>%
  slice_head(n = 5)
print(top5_tabela, n = 60)

# ==============================================================================
# 10. ANOTAÇÃO DOS CLUSTERS (Identificação biológica)
# ==============================================================================
# Primeiro, verifique com FeaturePlot onde os marcadores clássicos de camundongo acendem no UMAP:
# Exemplos comuns nesses 4 tecidos:
# - Células Imunes: Ptprc (CD45), Ccl5 (Linfócitos T/NK), Cd74 (Macrófagos/Apresentadoras)
# - Músculo Liso / Pericitos: Acta2, Myh11, Rgs5
# - Células Endoteliais (vasos): Pecam1, Cdh5
# - Células Epiteliais (intestino/estômago): Epcam, Krt8, Muc2
FeaturePlot(seurat_obj, features = c("Ptprc", "Acta2", "Pecam1", "Epcam"), ncol = 2)

# Depois de olhar o FeaturePlot e a tabela 'top5_tabela', renomeie os números dos clusters!
# ATENÇÃO: Substitua os nomes abaixo de acordo com o número exato de clusters que aparecer no seu UMAP:
novos_nomes <- c(
  "0" = "Epiteliais",
  "1" = "Imunes (Linfócitos/Mieloides)",
  "2" = "Fibroblastos / Estromais",
  "3" = "Músculo Liso / Pericitos",
  "4" = "Endoteliais"
  # Continue adicionando "5" = "...", "6" = "...", etc., conforme o total de clusters gerados
)

# Aplica os nomes anotados (descomente as 3 linhas abaixo quando ajustar a lista acima):
# seurat_obj <- RenameIdents(seurat_obj, novos_nomes)
# seurat_obj$cell_type <- Idents(seurat_obj)
# DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)


# ------------------------------------------------------------------------------
# 0. VOLTAR PARA OS CLUSTERS ORIGINAIS
# Seu script anterior sobrescreveu as identidades com RenameIdents.
# A coluna seurat_clusters continua guardando os números originais (0 a 22).
# ------------------------------------------------------------------------------
Idents(seurat_obj) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# 1. DIAGNÓSTICO DOS CLUSTERS 0 E 1  (rode e LEIA antes de anotar)
# ------------------------------------------------------------------------------
# 1a. De quais órgãos vem cada cluster? (% do cluster em cada órgão)
tab_orgao <- prop.table(table(seurat_obj$seurat_clusters, seurat_obj$orig.ident), margin = 1)
print(round(tab_orgao * 100, 1))
# Se os clusters 0 e 1 vierem de 3-4 órgãos ao mesmo tempo, NÃO são
# "Esteroidogênicas do ovário" nem "Retina": são algo compartilhado
# (provavelmente epitélio) ou um cluster de baixa qualidade.

# 1b. Qualidade por cluster (cluster de baixa qualidade tem poucos genes / muito mt)
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "seurat_clusters", pt.size = 0, ncol = 3)

# 1c. Top marcadores reais dos clusters 0 e 1
marcadores %>%
  filter(cluster %in% c("0", "1")) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 15) %>%
  print(n = 30)

# 1d. DotPlot por cluster com marcadores candidatos para ovário / olho / epitélio
genes_diagnostico <- c(
  "Epcam", "Krt8", "Krt18", "Cdh1",        # epitélio geral
  "Hsd3b1", "Cyp17a1", "Star", "Inha",      # ovário (esteroidogênicas/granulosa)
  "Pde6b", "Rho", "Sag", "Rlbp1",           # retina (fotorreceptores / Müller)
  "Krt12", "Crygs", "Pax6"                  # córnea / cristalino
)
genes_diagnostico <- genes_diagnostico[genes_diagnostico %in% rownames(seurat_obj)]

DotPlot(seurat_obj, features = genes_diagnostico, group.by = "seurat_clusters") +
  RotatedAxis() +
  ggtitle("Diagnóstico - marcadores por cluster")

# ------------------------------------------------------------------------------
# 2. ANOTAÇÃO
# Ajuste "0" e "1" conforme o diagnóstico acima. Os demais clusters você já
# validou (imunes, endoteliais, músculo liso, plasmócitos, etc).
#   - Se 0 e 1 expressarem Epcam/Krt8 e vierem de vários órgãos:
#       "Epitélio / Células epiteliais"
#   - Se 0 tiver Hsd3b1/Cyp17a1/Star e vier do ovário: mantenha "Esteroidogênicas"
#   - Se forem baixa qualidade: "Baixa qualidade / RNA ambiente"
# ------------------------------------------------------------------------------
novos_nomes <- c(
  "0"  = "Esteroidogênicas (Ovário)",       # <- REVISAR
  "1"  = "Retina / Neuronal",               # <- REVISAR
  "2"  = "Fibroblastos",
  "3"  = "Epitélio (Estômago/Intestino)",
  "4"  = "Epitélio (Estômago/Intestino)",
  "5"  = "Retina / Neuronal",
  "6"  = "Epitélio (Estômago/Intestino)",
  "7"  = "Músculo Liso / Pericitos",
  "8"  = "Células Imunes",
  "9"  = "Células Imunes",
  "10" = "Epitélio (Estômago/Intestino)",
  "11" = "Células Proliferativas",
  "12" = "Endoteliais",
  "13" = "Epitélio (Estômago/Intestino)",
  "14" = "Plasmócitos (IgA/IgK)",
  "15" = "Retina / Neuronal",
  "16" = "Glia / Suporte",
  "17" = "Mesoteliais",
  "18" = "Células Proliferativas",
  "19" = "Glia / Suporte",
  "20" = "Retina / Neuronal",
  "21" = "Retina / Neuronal",
  "22" = "Adipócitos"
)

seurat_obj <- RenameIdents(seurat_obj, novos_nomes)
seurat_obj$cell_type <- Idents(seurat_obj)

# Ordem fixa das identidades (mesma ordem no DotPlot e no Heatmap)
ordem <- unique(unname(novos_nomes))
Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = ordem)
seurat_obj$cell_type <- Idents(seurat_obj)

# ------------------------------------------------------------------------------
# 3. SLIDE 10 - DOT PLOT
# Adicionei Epcam (epitélio geral), mais robusto que Krt20 sozinho.
# ------------------------------------------------------------------------------
marcadores_chave <- c(
  "Hsd3b1",  # Esteroidogênicas (Ovário)
  "Pde6b",   # Retina / Neuronal
  "Col1a1",  # Fibroblastos
  "Epcam",   # Epitélio
  "Krt20",   # Epitélio diferenciado (intestino)
  "Acta2",   # Músculo liso
  "Ptprc",   # Células imunes
  "Mki67",   # Proliferativas
  "Pecam1",  # Endoteliais
  "Jchain",  # Plasmócitos
  "Plp1",    # Glia / Suporte
  "Msln",    # Mesoteliais
  "Adipoq"   # Adipócitos
)
marcadores_chave <- marcadores_chave[marcadores_chave %in% rownames(seurat_obj)]

p_dot_final <- DotPlot(seurat_obj, features = marcadores_chave,
                       cols = c("lightgrey", "darkblue")) +
  RotatedAxis() +
  ggtitle("Dot Plot - Marcadores Principais por Linhagem")
print(p_dot_final)
ggsave("DotPlot_Final.png", plot = p_dot_final, width = 12, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. SLIDE 11 - HEATMAP
# CORREÇÃO: o ScaleData padrão só escala os genes variáveis, então Hsd3b1 e
# Pde6b sumiam do heatmap. Aqui escalamos também os marcadores escolhidos.
# ------------------------------------------------------------------------------
seurat_obj <- ScaleData(
  seurat_obj,
  features = unique(c(VariableFeatures(seurat_obj), marcadores_chave))
)

set.seed(42)
p_heat_final <- DoHeatmap(
  subset(seurat_obj, downsample = 100),
  features = marcadores_chave,
  size = 3.5
) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick3"), na.value = "white") +
  ggtitle("Heatmap - Assinatura Transcricional dos Tipos Celulares")
print(p_heat_final)
ggsave("Heatmap_Final.png", plot = p_heat_final, width = 12, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 5. SLIDE 12 - UMAP ANOTADO
# CORREÇÃO: com vários clusters agrupados na mesma identidade, o rótulo cai na
# mediana do grupo e fica no lugar errado. Usamos só a legenda (label = FALSE).
# ------------------------------------------------------------------------------
p_umap_final <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type",
                        label = FALSE, pt.size = 0.4) +
  ggtitle("UMAP - Tipos Celulares Anotados (4 Órgãos)")
print(p_umap_final)
ggsave("UMAP_Anotado_Final.png", plot = p_umap_final, width = 11, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 6. SLIDE 13 - UMAP DIVIDIDO POR ÓRGÃO
# ------------------------------------------------------------------------------
p_split_final <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type",
                         split.by = "orig.ident", ncol = 2, pt.size = 0.3) +
  ggtitle("Distribuição dos Tipos Celulares por Órgão")
print(p_split_final)
ggsave("UMAP_Por_Orgao_Final.png", plot = p_split_final, width = 13, height = 8, dpi = 300)
```
