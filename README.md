# Projeto_Aplicado_MCD

Atividade de Projeto Aplicado, desenvolvido pelo grupo My Chemical Data para o Instituto SENAI de Tecnologia Ambiental — UniSENAI.

Este projeto processa dados de identificação e abundância de compostos químicos, enriquecendo-os com informações de bases de dados públicas como PubChem e ChEBI. Inclui classificação de compostos (natural/sintético), cálculo de scores de confiança e geração de rankings.

## Como Iniciar o Projeto

1. **Baixe ou Clone o Repositório**  
   Clone este repositório para sua máquina local.

2. **Instale as Dependências**  
   Execute o comando no terminal:  
   ```
   pip install pandas requests openpyxl scikit-learn
   ```

3. **Execute o Pipeline**  
   No terminal, execute:  
   ```
   python pipeline.py
   ```

O pipeline processará os arquivos `dados_brutos/IDENTIFICACAO.xlsx` e `dados_brutos/ABUND.xlsx`, gerando saídas como `compostos_final_resultado.xlsx`, rankings e caches.

## Visão Geral do Projeto

O projeto é dividido em entregas principais, demonstradas no notebook `entregas_finais.ipynb`. Cada entrega corresponde a uma etapa do processamento de dados:

### 1. Junção de Planilhas de Identificação e Abundância
- Mescla as planilhas `IDENTIFICACAO.xlsx` e `ABUND.xlsx` pela coluna `Compound` (inner join).
- Limpa duplicatas e registros sem identificação.
- Resultado: DataFrame com compostos únicos e suas abundâncias.

### 2. Classificação: Produto Natural/Metabólito ou Composto Sintético
- Baseado na ontologia ChEBI e heurísticas por nome.
- Regras prioritárias:
  - Peptídeos detectados por abreviações de aminoácidos.
  - Termos ChEBI (ex.: "fatty acid" → Natural, "drug" → Sintético).
  - Fallback por nome ou ontologia.
- Exemplos: Sertraline (Sintético), Citric acid (Natural).

### 3. Consulta a Bases de Dados Públicas (PubChem e ChEBI)
- **PubChem PUG REST**: Propriedades moleculares (CID, fórmula, massa, IUPAC).
- **ChEBI (via PubChem)**: ID ChEBI e hierarquia taxonômica.
- Demonstração com compostos como Sertraline e Caffeine.

### 4. Obtenção da Ontologia ou Classificação Química (ChEBI)
- Hierarquia ChEBI obtida via PubChem Classification.
- Subida na árvore taxonômica até 5 níveis, excluindo nós genéricos.
- Exemplo: "secondary amino compound | organic amino compound".

### 5. Levantamento de Usos e Aplicações Conhecidas (PubChem)
- Endpoint `/compound/cid/{CID}/description/JSON` para textos descritivos.
- Inclui usos farmacológicos, propriedades biológicas, etc.

## Transformação dos Dados

### Decisão sobre Replicatas
Foi utilizada agregação por média para reduzir ruído e melhorar a confiabilidade estatística.

### APIs Utilizadas
- **PubChem**: Identificação química, propriedades moleculares, descrições e usos.
- **ChEBI**: Enriquecimento de dados com ontologia e classificação química.

### Validação
Foi realizado merge com `compostos_final.xlsx` usando Compound ID para validação dos resultados.

### Ranking
Utilizada fórmula ponderada para calcular a confiança:
```
confianca = 0.25 * media_norm + 0.20 * score_norm + 0.15 * frag_norm + 0.15 * isotope_norm + 0.15 * estabilidade + 0.10 * metadata_ok
```
Onde:
- `media_norm`: Abundância média normalizada (log-transformada).
- `score_norm`: Score de identificação normalizado.
- `frag_norm`: Fragmentation Score normalizado.
- `isotope_norm`: Isotope Similarity normalizado.
- `estabilidade`: 1 - CV normalizado (coeficiente de variação).
- `metadata_ok`: Indicador binário se há dados de PubChem ou ChEBI.

## Arquivos de Saída
O pipeline gera os seguintes arquivos:
- `pipeline_final.csv`: Dados completos enriquecidos.
- `ranking.csv`: Ranking resumido por descrição.
- `ranking_detalhado.csv`: Ranking com detalhes (confiança, média, CV, scores, etc.).
- `compostos_final_resultado.xlsx`: Planilha final no formato modelo, com colunas como ID, Metabólito/Composto, Fórmula, Massa Molecular, etc.
- `cache_pubchem.json`: Cache das consultas às APIs para evitar reprocessamento.

## Créditos
- Matheus Cruz
- Nicolas Mello
- João Paulo
