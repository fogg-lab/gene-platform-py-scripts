from string import digits
from typing import TypedDict, Union
import numpy as np


class GeneInfoResult(TypedDict):
    reference_indices: list[int]
    id_type: str
    num_matches: int


def get_gene_info(
    query_genes: Union[list[str], list[int]],
    gene_reference: list[list[Union[str, int]]],
) -> GeneInfoResult:
    """Given a list of ensembl gene ids, entrez IDs, or gene symbols, find indices of the genes in the reference list.

    Args:
        query_genes (list[str]): List of ensembl gene ids, entrez IDs, or gene names.
        gene_reference (list[list[str]]): [ensembl id, entrez ID, symbol, biotype] for all genes.

    Returns:
        dict: A dictionary containing the following keys:
            - 'reference_indices' (list[int]): List of indices of the genes in the reference list, -1 if the gene is not found.
            - 'id_type' (str): id type (either 'entrezgene', 'ensembl_gene', or 'symbol').
            - 'num_matches' (int): Number of query genes found in the reference list.
            - 'duplicates' (list[int]): Indices of duplicated genes in the original query_genes list.
    """
    # Check for empty query
    if len(query_genes) == 0:
        return dict(reference_indices=[], id_type="N/A", num_matches=0)
    # Determine if Ensembl ID, Entrez ID, or HGNC symbol
    test_gene = query_genes[0]
    if isinstance(test_gene, int) or all(c in digits for c in str(test_gene)):
        id_type = "entrezgene"
    elif test_gene.lower().startswith("ens") and "0000" in test_gene:
        id_type = "ensembl_gene"
    else:
        id_type = "symbol"

    gene_record_idx = ["ensembl_gene", "entrezgene", "symbol"].index(id_type)
    formatter = int if id_type == "entrezgene" else str.lower
    ref_ids = np.array([formatter(rec[gene_record_idx]) for rec in gene_reference])
    genes = np.array([formatter(g) for g in query_genes])

    # Calculate the indices that sort genes and ref_ids
    genes_sort_idx = np.argsort(genes)
    sorted_genes = genes[genes_sort_idx]
    ref_ids_sort_idx = np.argsort(ref_ids)
    sorted_ref_ids = ref_ids[ref_ids_sort_idx]

    # Map genes to reference indices
    sorted_ref_index = 0
    gene_ref_indices = np.full(len(genes), -1, dtype=np.int32)
    for i, g in enumerate(sorted_genes):
        while sorted_ref_index < len(ref_ids) and g > sorted_ref_ids[sorted_ref_index]:
            sorted_ref_index += 1
        if sorted_ref_index == len(ref_ids):
            break
        if g == sorted_ref_ids[sorted_ref_index]:
            g_index = genes_sort_idx[i]
            ref_index = ref_ids_sort_idx[sorted_ref_index]
            gene_ref_indices[g_index] = ref_index

    return dict(
        reference_indices=gene_ref_indices.tolist(),
        id_type=id_type,
        num_matches=int(np.sum(gene_ref_indices != -1)),
    )


def get_duplicated(query_genes: Union[list[str], list[int]]) -> list[list[int]]:
    """Return indices of genes that occur more than once in query_genes, grouped by value.
    Args:
        query_genes: List of Ensembl IDs, Entrez IDs, or HGNC symbols.
    Returns:
        list[list[int]]: Groups of indices to duplicated genes in the query.
    Examples:
        >>> get_duplicated(['ab', 'op', 'ef', 'ef', 'gh', 'ab', 'op', 'cd', 'bc', 'mn'])
        [[0, 5], [2, 3], [1, 6]]
        >>> get_duplicated(['ab', 'op', 'ef'])  # No duplicates
        []
    """
    query_genes = np.array(query_genes)
    sort_idx = query_genes.argsort()
    sorted_genes = query_genes[sort_idx]
    start_idx = np.unique(sorted_genes, return_index=True)[1]
    return [
        group_idx.tolist()
        for group_idx in np.split(sort_idx, start_idx[1:])
        if len(group_idx) > 1
    ]
