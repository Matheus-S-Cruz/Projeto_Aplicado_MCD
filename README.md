# QuimioAnalytics / MyChemicalData

Projeto Aplicado desenvolvido pelo grupo **My Chemical Data** para o Instituto SENAI de Tecnologia Ambiental — UniSENAI.

Aplicação web local que processa dados de **identificação** e **abundância** de compostos químicos (exportados do software Progenesis), enriquece-os com bases públicas (**PubChem**, **ChEBI**, **ClassyFire**), calcula um **score de confiança** e apresenta o resultado numa **tabela final + dashboard interativo**.

---

## Como rodar

Pré-requisito: **Python 3.10+** instalado.

**1. Clone o repositório e entre na pasta:**
```bash
git clone <url-do-repo>
cd Projeto_Aplicado_MCD
```

**2. Crie e ative um ambiente virtual (venv):**
```bash
python3 -m venv .venv
source .venv/bin/activate      # Linux/Mac
# .venv\Scripts\activate       # Windows
```

**3. Instale as dependências:**
```bash
pip install -r requirements.txt
```

**4. Rode a aplicação:**
```bash
streamlit run app.py
```

Isso abre a aplicação no navegador (normalmente `http://localhost:8501`).

> O banco de dados (`quimioanalytics.db`, SQLite) é **criado automaticamente** na primeira execução — não é preciso rodar nenhum comando de banco. Ele começa vazio; basta processar uma vez (passo abaixo) para preenchê-lo.

---

## Como usar

Na barra lateral:

1. **Envie as duas planilhas** do Progenesis — *identificação* e *abundância* — **ou** marque **"Usar planilhas de exemplo"** (usa os arquivos em `dados_brutos/`).
2. Escolha o **modo de enriquecimento**:
   - **Cache** (rápido, ideal para demonstração) — usa apenas dados já consultados antes (`cache_pubchem.json`, que já vem preenchido no repositório).
   - **APIs** (completo, mais lento) — consulta PubChem/ChEBI/ClassyFire online. **Requer internet.**
3. Clique em **Processar**.

O resultado aparece em 3 abas:
- **Tabela (Documento IST)** — a tabela final consolidada, com filtros (busca, classe, score, tag) e exportação CSV/Excel.
- **Dashboard** — indicadores e gráficos (ranking por confiança, distribuição por classe química, tags, amostra mais abundante).
- **Sobre** — descrição do sistema.

Cada processamento é salvo no **Histórico** (seletor no topo da barra lateral), permitindo revisitar análises anteriores.

---

## Estrutura do projeto

| Arquivo | Descrição |
|---|---|
| `app.py` | Interface web (Streamlit) |
| `etl.py` | Motor de ETL: leitura, enriquecimento (PubChem/ChEBI/ClassyFire), estatística, ranking e persistência |
| `models.py` | Modelagem do banco (SQLAlchemy, esquema estrela) |
| `pipeline.py` | Versão **antiga** do processamento (script único) — mantida como referência |
| `requirements.txt` | Dependências |
| `dados_brutos/` | Planilhas de exemplo (`IDENTIFICACAO.xlsx`, `ABUND.xlsx`) |
| `cache_pubchem.json` | Cache das consultas às APIs (evita reprocessar) |

Arquivos **gerados** (não versionados): `quimioanalytics.db` e `historico/`.

---

## Como funciona (resumo)

1. **Integração** — junta identificação + abundância pela coluna `Compound` (1 linha por composto).
2. **Enriquecimento** — busca fórmula, CID, ontologia (ChEBI) e classe química (ClassyFire) por composto.
3. **Estatística** — média, desvio-padrão e coeficiente de variação das abundâncias por replicata.
4. **Score de confiança** (fórmula ponderada):
   ```
   confianca = 0.25*media + 0.20*score + 0.15*fragmentacao + 0.15*isotopo + 0.15*estabilidade + 0.10*metadados
   ```
5. **Ranking e consolidação** — ordena e monta a tabela final (formato "Documento IST").

---

## Créditos
- Matheus Cruz
- Nicolas Mello
- João Paulo
