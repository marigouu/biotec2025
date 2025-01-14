import streamlit as st
import pandas as pd
from Bio import Entrez
import itertools

# Configuração do e-mail do Entrez
Entrez.email = "seu@email.com"  # Substitua pelo seu e-mail

# Função para gerar combinações de termos (incluindo combinações de um ou mais termos)
def generate_combinations(terms, min_terms=1):
    """Gera todas as combinações possíveis de termos com operadores booleanos."""
    combinations = []
    for i in range(min_terms, len(terms) + 1):
        for subset in itertools.combinations(terms, i):
            combinations.append(' AND '.join([f'"{term.strip()}"' for term in subset]))
    return combinations

# Função para buscar contagem de artigos no PubMed
def fetch_article_counts(term, year, max_results=3):
    query = f"{term} AND {year}[DP]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"]), record["IdList"]

# Função para buscar detalhes dos artigos no PubMed
def fetch_article_details(ids):
    handle = Entrez.esummary(db="pubmed", id=",".join(ids))
    records = Entrez.read(handle)
    handle.close()
    articles = []
    for record in records:
        pubmed_id = record["Id"]
        # Obtendo o DOI diretamente, se disponível
        doi = record.get("DOI", "N/A")
        # Formatar o DOI para ser um link clicável
        doi_link = f"[{doi}](https://doi.org/{doi})" if doi != "N/A" else "N/A"
        articles.append({
            "PubMed ID": pubmed_id,  # Exibindo apenas o PubMed ID
            "Title": record.get("Title", ""),
            "Year": record.get("PubDate", "").split(" ")[0],
            "Journal": record.get("Source", ""),
            "DOI": doi_link  # DOI como link clicável
        })
    return articles

# Função principal de Streamlit
def main():
    st.title("Revisão Bibliográfica do PubMed")

    # Entrada do usuário
    terms_input = st.text_input("Insira os termos de pesquisa separados por vírgula:", "machine learning, drug discovery, cancer")
    years = st.slider("Selecione o intervalo de anos:", 2000, 2030, (2020, 2024))
    
    # Geração de combinações de termos (para a primeira tabela)
    terms = terms_input.split(",")
    combinations_single_term_or_more = generate_combinations(terms, min_terms=1)

    results_single_term_or_more = []
    for combo in combinations_single_term_or_more:
        row = {"combination": combo}
        for year in range(years[0], years[1] + 1):
            count, ids = fetch_article_counts(combo, year)
            row[year] = count
        results_single_term_or_more.append(row)

    # Exibindo a tabela de quantidade de artigos por ano (para combinações de termos)
    st.subheader("Quantidade de Artigos por Ano (Combinações de Termos)")
    df_single_term_or_more = pd.DataFrame(results_single_term_or_more)
    st.write(df_single_term_or_more)

    # Geração de combinações de múltiplos termos para a segunda tabela (mínimo 2 termos)
    combinations_multiple_terms = generate_combinations(terms, min_terms=2)

    results_multiple_terms = []
    detailed_articles = []
    for combo in combinations_multiple_terms:
        row = {"combination": combo}
        for year in range(years[0], years[1] + 1):
            count, ids = fetch_article_counts(combo, year)
            row[year] = count
            # Buscar detalhes dos artigos mais relevantes
            if ids:
                detailed_articles.extend(fetch_article_details(ids[:3]))  # Limite de 3 artigos mais relevantes
        results_multiple_terms.append(row)

    # Exibindo a tabela de artigos mais relevantes (com combinações de termos)
    st.subheader("Artigos Mais Relevantes (Com Combinações de Termos)")
    df_relevant = pd.DataFrame(detailed_articles)
    st.write(df_relevant)

    # Adicionando botões de download
    st.subheader("Baixar Tabelas")

    # Função para baixar uma tabela como TXT
    def download_as_txt(df, filename):
        txt_data = df.to_string(index=False)
        st.download_button(
            label=f"Baixar {filename} (TXT)",
            data=txt_data.encode('utf-8'),
            file_name=f"{filename}.txt",
            mime="text/plain",
        )

    # Função para baixar uma tabela como XLSX
    def download_as_xlsx(df, filename):
        xlsx_file = pd.ExcelWriter(f"{filename}.xlsx", engine="openpyxl")
        df.to_excel(xlsx_file, index=False, sheet_name="Sheet1")
        xlsx_file.close()
        with open(f"{filename}.xlsx", "rb") as f:
            xlsx_data = f.read()
        st.download_button(
            label=f"Baixar {filename} (XLSX)",
            data=xlsx_data,
            file_name=f"{filename}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )

    # Baixar tabela de quantidade por ano (para combinações de termos)
    download_as_txt(df_single_term_or_more, "quantidade_artigos_por_ano_combinacoes")
    download_as_xlsx(df_single_term_or_more, "quantidade_artigos_por_ano_combinacoes")

    # Baixar tabela de artigos mais relevantes (para combinações de termos)
    download_as_txt(df_relevant, "artigos_mais_relevantes_combinacoes")
    download_as_xlsx(df_relevant, "artigos_mais_relevantes_combinacoes")

if __name__ == "__main__":
    main()
