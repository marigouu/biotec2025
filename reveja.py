#Este código é uma demo de uma ferramenta de revisão sistemática do Pudmed.

from Bio import Entrez
import itertools

# Configuração da API Entrez
Entrez.email = "seu@email.com"

def generate_combinations(terms):
    """Gera todas as combinações possíveis de termos com operadores booleanos."""
    combinations = []
    for i in range(1, len(terms) + 1):
        for subset in itertools.combinations(terms, i):
            combinations.append(" AND ".join(subset))
    return combinations

def fetch_article_counts_and_top_ids(term, year, max_results=3):
    """Busca contagem de artigos e os IDs mais relevantes no PubMed para um termo e ano."""
    query = f"{term} AND {year}[DP]"  # Filtra por data de publicação
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
    record = Entrez.read(handle)
    handle.close()

    count = int(record["Count"])
    ids = record["IdList"]
    return count, ids

def main():
    # Entradas: termos e intervalo de anos
    terms = ["carbon dots", "green synthesis", "polymeric nanoparticles"]
    years = list(range(2020, 2025))  # De 2020 a 2024

    # Gerar combinações de termos
    combinations = generate_combinations(terms)
    print(f"Combinações de busca: {combinations}")

    # Inicializar tabela
    results = []

    # Iterar por combinações e anos
    for combo in combinations:
        row = {"combination": combo}
        for year in years:
            count, top_ids = fetch_article_counts_and_top_ids(combo, year)
            row[year] = {"count": count, "top_ids": top_ids}
        results.append(row)

    # Exibir tabela
    print("\n=== Resultados ===")
    for row in results:
        print(f"\nCombinação: {row['combination']}")
        for year in years:
            data = row[year]
            print(f"Ano: {year} | Artigos: {data['count']} | IDs: {', '.join(data['top_ids'])}")

if __name__ == "__main__":
    main()
