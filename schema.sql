CREATE TABLE dim_composto (
    compound_id VARCHAR(50) PRIMARY KEY,
    nome_composto VARCHAR(255) NOT NULL,
    descricao TEXT,
    pubchem_cid INT,
    chebi_id VARCHAR(50),
    formula VARCHAR(100),
    role_biologico VARCHAR(100),
    classe_taxonomica VARCHAR(100)
);

CREATE TABLE dim_amostra (
    id_amostra VARCHAR(50) PRIMARY KEY,
    nome_amostra VARCHAR(100) NOT NULL,
    tipo VARCHAR(50),
    projeto VARCHAR(100),
    data_coleta DATE
);

CREATE TABLE fact_abundancia (
    id_medicao SERIAL PRIMARY KEY,
    compound_id VARCHAR(50) NOT NULL,
    id_amostra VARCHAR(50) NOT NULL,
    abundancia_valor DOUBLE PRECISION,
    sinal_analitico VARCHAR(255),
    tipo_replicata VARCHAR(20),
    score_identificacao DOUBLE PRECISION,

    FOREIGN KEY (compound_id)
        REFERENCES dim_composto(compound_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE,

    FOREIGN KEY (id_amostra)
        REFERENCES dim_amostra(id_amostra)
        ON DELETE CASCADE
        ON UPDATE CASCADE
);

CREATE TABLE dim_ontologia (
    id SERIAL PRIMARY KEY,
    compound_id VARCHAR(50) NOT NULL,
    ontologia TEXT NOT NULL,

    FOREIGN KEY (compound_id)
        REFERENCES dim_composto(compound_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE
);

-- Índices
CREATE INDEX idx_fact_compound ON fact_abundancia(compound_id);
CREATE INDEX idx_fact_amostra ON fact_abundancia(id_amostra);
CREATE INDEX idx_ontologia_compound ON dim_ontologia(compound_id);