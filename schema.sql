CREATE TABLE dim_composto (
    compound_id VARCHAR PRIMARY KEY,
    nome_composto VARCHAR,
    descricao TEXT,
    pubchem_cid INT,
    chebi_id VARCHAR,
    formula VARCHAR,
    role_biologico VARCHAR,
    ontologia TEXT,
    classe_taxonomica VARCHAR
);

CREATE TABLE dim_amostra (
    id_amostra VARCHAR PRIMARY KEY,
    nome_amostra VARCHAR,
    tipo VARCHAR,
    projeto VARCHAR,
    data_coleta DATE
);

CREATE TABLE fact_abundancia (
    id_medicao SERIAL PRIMARY KEY,
    compound_id VARCHAR,
    id_amostra VARCHAR,
    abundancia_valor FLOAT,
    sinal_analitico VARCHAR,
    tipo_replicata VARCHAR,
    score_identificacao FLOAT,
    FOREIGN KEY (compound_id) REFERENCES dim_composto(compound_id),
    FOREIGN KEY (id_amostra) REFERENCES dim_amostra(id_amostra)
);