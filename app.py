"""
QuimioAnalytics — Web app local (Streamlit).

Fluxo: o usuário sobe as 2 planilhas do Progenesis (identificação + abundância)
-> o sistema roda o ETL + enriquecimento (etl.py) -> grava no banco (models.py)
-> exibe a tabela final "Documento IST" + um dashboard interativo.

Rodar:  streamlit run app.py
"""

from __future__ import annotations

import io

import pandas as pd
import plotly.express as px
import streamlit as st

import etl
import models

st.set_page_config(page_title="QuimioAnalytics", page_icon=":material/science:", layout="wide")


@st.cache_resource
def get_session():
    return models.init_db()


Session = get_session()


def tem_dados() -> bool:
    with Session() as s:
        return s.query(models.Composto).first() is not None


# =========================
# SIDEBAR — entrada e processamento
# =========================
st.sidebar.title(":material/science: QuimioAnalytics")
st.sidebar.caption("Análise e priorização de compostos químicos")

st.sidebar.header(":material/upload_file: 1. Planilhas do Progenesis")
arq_id = st.sidebar.file_uploader("Identificação (IDENTIFICACAO.xlsx)", type=["xlsx"])
arq_ab = st.sidebar.file_uploader("Abundância (ABUND.xlsx)", type=["xlsx"])
usar_exemplo = st.sidebar.checkbox("Usar planilhas de exemplo (dados_brutos/)", value=False)

st.sidebar.header(":material/tune: 2. Modo de enriquecimento")
modo = st.sidebar.radio(
    "Como buscar dados químicos?",
    ["Cache (rápido — ideal p/ demo)", "APIs (completo — mais lento)"],
    help="O modo cache usa apenas o que já foi consultado antes (cache_pubchem.json).",
)
usar_cache = modo.startswith("Cache")

if st.sidebar.button("Processar", type="primary", width="stretch", icon=":material/play_arrow:"):
    id_src = "dados_brutos/IDENTIFICACAO.xlsx" if usar_exemplo else arq_id
    ab_src = "dados_brutos/ABUND.xlsx" if usar_exemplo else arq_ab

    if id_src is None or ab_src is None:
        st.sidebar.error("Envie as duas planilhas (ou marque 'usar exemplo').")
    else:
        barra = st.sidebar.progress(0.0, text="Iniciando...")

        def _prog(msg, frac):
            barra.progress(min(frac or 0.0, 1.0), text=msg)

        try:
            resumo = etl.processar(id_src, ab_src, Session, usar_cache_apenas=usar_cache, progresso=_prog)
            st.session_state["resumo"] = resumo
            barra.empty()
            st.sidebar.success(f"Processados {resumo['compostos']} compostos!")
            st.rerun()
        except Exception as e:
            barra.empty()
            st.sidebar.error(f"Erro ao processar: {e}")

# =========================
# CORPO
# =========================
st.title("QuimioAnalytics — Painel de Compostos")

if not tem_dados():
    st.info("Envie as planilhas de identificação e abundância (ou marque *usar exemplo*) e clique em **Processar**.", icon=":material/arrow_back:")
    st.stop()

df = etl.documento_ist(Session)

if "resumo" in st.session_state:
    r = st.session_state["resumo"]
    st.caption(f"Último processamento: {r['compostos']} compostos · {r['medicoes']} medições · "
               f"{r['amostras']} amostras · modo {r['modo']}.")

aba_tabela, aba_dash, aba_sobre = st.tabs([
    ":material/table_view: Tabela (Documento IST)",
    ":material/bar_chart: Dashboard",
    ":material/info: Sobre",
])

# ---------- TABELA ----------
with aba_tabela:
    st.subheader("Tabela final — formato Documento IST")

    c1, c2, c3 = st.columns([2, 2, 1])
    busca = c1.text_input(":material/search: Buscar (composto ou descrição)")
    classes = sorted([c for c in df["Classe geral"].dropna().unique()])
    filtro_classe = c2.multiselect("Filtrar por classe", classes)
    score_min = c3.number_input("Score mínimo", value=0.0, step=5.0)

    dfx = df.copy()
    if busca:
        mask = (
            dfx["Composto"].astype(str).str.contains(busca, case=False, na=False)
            | dfx["Descrição"].astype(str).str.contains(busca, case=False, na=False)
        )
        dfx = dfx[mask]
    if filtro_classe:
        dfx = dfx[dfx["Classe geral"].isin(filtro_classe)]
    if score_min:
        dfx = dfx[dfx["Score"].fillna(0) >= score_min]

    st.caption(f"{len(dfx)} de {len(df)} compostos")
    st.dataframe(dfx, width="stretch", hide_index=True)

    csv = dfx.to_csv(index=False).encode("utf-8")
    st.download_button("Baixar CSV", csv, "documento_ist.csv", "text/csv", icon=":material/download:")
    buf = io.BytesIO()
    dfx.to_excel(buf, index=False)
    st.download_button("Baixar Excel", buf.getvalue(), "documento_ist.xlsx", icon=":material/download:")

# ---------- DASHBOARD ----------
with aba_dash:
    st.subheader("Visão geral")

    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Compostos", len(df))
    classificados = df["Classe geral"].notna() & (df["Classe geral"] != "Não classificado")
    k2.metric("Classificados", int(classificados.sum()))
    k3.metric("Score médio", f"{df['Score'].mean():.1f}")
    k4.metric("Confiança média", f"{df['Confiança'].mean():.2f}")

    g1, g2 = st.columns(2)

    # Top 15 por confiança
    top = df.nsmallest(15, "Posição")[["Composto", "Confiança"]].sort_values("Confiança")
    fig_top = px.bar(top, x="Confiança", y="Composto", orientation="h",
                     title="Top 15 compostos por confiança")
    g1.plotly_chart(fig_top, width="stretch")

    # Distribuição por classe
    cont = df["Classe geral"].fillna("(sem classe)").value_counts().reset_index()
    cont.columns = ["Classe geral", "qtd"]
    fig_classe = px.pie(cont.head(10), names="Classe geral", values="qtd",
                        title="Distribuição por classe química (top 10)")
    g2.plotly_chart(fig_classe, width="stretch")

    g3, g4 = st.columns(2)

    # Modo de aquisição
    modo_cont = df["Modo de aquisição"].fillna("Indeterminado").value_counts().reset_index()
    modo_cont.columns = ["Modo", "qtd"]
    fig_modo = px.bar(modo_cont, x="Modo", y="qtd", title="Modo de aquisição")
    g3.plotly_chart(fig_modo, width="stretch")

    # Amostra mais abundante
    if "Amostra mais abundante" in df.columns:
        am = df["Amostra mais abundante"].fillna("(sem)").value_counts().reset_index()
        am.columns = ["Amostra", "qtd"]
        fig_am = px.bar(am, x="Amostra", y="qtd", title="Compostos por amostra mais abundante")
        g4.plotly_chart(fig_am, width="stretch")

# ---------- SOBRE ----------
with aba_sobre:
    st.markdown(
        """
        ### O que é o QuimioAnalytics
        Sistema de apoio à decisão científica que processa dados de espectrometria de
        massas exportados do **Progenesis**, enriquece com bases públicas
        (**PubChem / ChEBI**), classifica os compostos e gera um **ranking de priorização**.

        **Fluxo:** upload das planilhas → ETL + enriquecimento → banco (SQLite) →
        tabela *Documento IST* + dashboard.

        **Origem dos dados:** os valores de identificação e abundância vêm do Progenesis;
        fórmula, classe química e descrição são buscadas nas APIs; o *score de confiança*,
        as estatísticas e o ranking são calculados pelo sistema.
        """
    )
