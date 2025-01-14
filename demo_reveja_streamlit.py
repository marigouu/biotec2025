import streamlit as st
import pandas as pd
from Bio import Entrez
import itertools

# Configuração do e-mail do Entrez
Entrez.email = "gouveiaunb@gmail.com"  # Substitua pelo seu e-mail

# Função para gerar combinações de termos
def generate_combinations(terms):
    """Gera todas as combinações possíveis de termos com operadores booleanos, agora com aspas ao redor dos termos."""
    combinations = []
    for i in range(1, len(terms) + 1):
        for subset in itertools.combinations(terms, i):
            combinations.append(' AND '.join([f'"{term.strip()}"' for term in subset]))  # Coloca aspas ao redor de cada termo
    return combinations

# Função para buscar no PubMed
def fetch_article_counts(term, year, max_results=3):
    query = f"{term} AND {year}[DP]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"]), record["IdList"]

# Função principal de Streamlit
def main():
    st.title("Revisão Bibliográfica do PubMed")

    # Entrada do usuário
    terms_input = st.text_input("Insira os termos de pesquisa separados por vírgula:", "machine learning, drug discovery, cancer")
    years = st.slider("Selecione o intervalo de anos:", 2000, 2030, (2020, 2024))
    
    # Geração de combinações
    terms = terms_input.split(",")
    combinations = generate_combinations(terms)

    results = []
    for combo in combinations:
        row = {"combination": combo}
        for year in range(years[0], years[1] + 1):
            count, _ = fetch_article_counts(combo, year)
            row[year] = count
        results.append(row)

    # Exibindo os resultados como uma tabela
    df = pd.DataFrame(results)
    st.write(df)

    # Adicionando botão de download para diferentes formatos
    st.subheader("Baixar Tabela")
    
    # Baixar como CSV
    csv_data = df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Baixar como CSV",
        data=csv_data,
        file_name="resultado_pubmed.csv",
        mime="text/csv",
    )
    
    # Baixar como TXT
    txt_data = df.to_string(index=False)
    st.download_button(
        label="Baixar como TXT",
        data=txt_data.encode('utf-8'),
        file_name="resultado_pubmed.txt",
        mime="text/plain",
    )
    
    # Baixar como XLSX
    xlsx_file = pd.ExcelWriter("resultado_pubmed.xlsx", engine='openpyxl')
    df.to_excel(xlsx_file, index=False, sheet_name="Resultados")
    xlsx_file.close()
    with open("resultado_pubmed.xlsx", "rb") as f:
        xlsx_data = f.read()
    st.download_button(
        label="Baixar como XLSX",
        data=xlsx_data,
        file_name="resultado_pubmed.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

if __name__ == "__main__":
    main()

