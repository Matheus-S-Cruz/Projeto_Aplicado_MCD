# Projeto_Aplicado_MCD
Atividade de Projeto Aplicado, desenvolvido pelo grupo My Chemical Data.

## Como Iniciar o Projeto:
. Baixe ou Clone o Repositório
. No terminal use o comando: pip install pandas requests openpyxl
. Inicie o projeto com o comando python: python enriquecer_pubchem.py

## Créditos:
Matheus Cruz,
Nicolas Mello,
João Paulo

## Projeto CDIA — Aula 4

### Decisão sobre replicatas
Foi utilizada agregação por média para reduzir ruído e melhorar a confiabilidade estatística.

### APIs utilizadas
- PubChem para identificação química
- (ChEBI será implementado futuramente)

### Validação
Foi realizado merge com compostos_final.xlsx usando Compound ID.

### Ranking
Utilizada fórmula:
confianca = (1 / (1 + CV)) * score * log(abundancia_media + 1)
