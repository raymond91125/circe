"""Unit tests for the Bhattacharya 2019 innexin-expression ingest."""

from pathlib import Path

from celegans_connectome_kg.ingest.innexin_expression import read_innexin_expression

D = Path(__file__).resolve().parents[1] / "data" / "bhattacharya-2019-innexin"


def test_read_innexin_expression() -> None:
    data = read_innexin_expression(D / "innexin_expression.csv", D / "innexin_genes.csv")
    # 15 distinct innexin genes (17 fig labels collapse inx-1 a/b, inx-18 a/b)
    assert len(data.genes) == 15
    assert all(g.wbgene.startswith("WB:WBGene") and g.category == "innexin" for g in data.genes)
    # 867 verified expression entries; states constrained
    assert len(data.expressions) == 867
    assert {e.state for e in data.expressions} == {"both", "non-dauer_only", "dauer_only"}
    by = {(e.neuron_label, e.wbgene, e.isoform): e for e in data.expressions}
    # inx-1 isoforms carried as qualifiers on the same gene id
    assert ("ADA", "WB:WBGene00002123", "a") in by
    # a plasticity ('*') entry is captured (HSN inx-1a is starred in Fig 1B)
    hsn = by[("HSN", "WB:WBGene00002123", "a")]
    assert hsn.state == "non-dauer_only" and hsn.plasticity is True
