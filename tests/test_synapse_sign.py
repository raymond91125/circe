"""Unit tests for the Part C chemical-synapse sign engine (Fenyves 2020 derivation logic).

These pin the *rule* in isolation (the build integration is in test_sex_build.py). Correctness of
the rule is the real guarantee here: the aggregate coverage/E:I divergence from Fenyves on current
data is a data-coverage/vintage question decomposed in the increment-4 validation analysis, not a
rule question.
"""

from pathlib import Path

from celegans_connectome_kg.ingest.synapse_sign import (
    Receptor,
    classify_sign,
    decode_fast_transmitters,
    read_receptor_catalog,
    read_receptor_expression,
)

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data" / "cengen-ionotropic-receptors"

# A tiny hand-built catalog: one cation + one anion receptor per fast transmitter.
_CAT = {
    "WB:ach_c": Receptor("acr-x", "WB:ach_c", "acetylcholine", "cation", False),
    "WB:ach_a": Receptor("acc-x", "WB:ach_a", "acetylcholine", "anion", False),
    "WB:glu_c": Receptor("glr-x", "WB:glu_c", "glutamate", "cation", False),
    "WB:glu_a": Receptor("avr-x", "WB:glu_a", "glutamate", "anion", False),
    "WB:gab_a": Receptor("unc-49", "WB:gab_a", "GABA", "anion", False),
    "WB:pre_x": Receptor("lgc-46", "WB:pre_x", "acetylcholine", "anion", True),  # excluded
}


def test_decode_fast_transmitters() -> None:
    assert decode_fast_transmitters("a") == frozenset({"acetylcholine"})
    assert decode_fast_transmitters("l") == frozenset({"glutamate"})
    assert decode_fast_transmitters("g") == frozenset({"GABA"})
    # co-transmission keeps both fast transmitters
    assert decode_fast_transmitters("ag") == frozenset({"acetylcholine", "GABA"})
    # monoamine / betaine / serotonin / unknown letters drop out
    assert decode_fast_transmitters("as") == frozenset({"acetylcholine"})
    assert decode_fast_transmitters("lt") == frozenset({"glutamate"})
    assert decode_fast_transmitters("bl") == frozenset({"glutamate"})
    assert decode_fast_transmitters("d") == frozenset()
    assert decode_fast_transmitters("u") == frozenset()
    assert decode_fast_transmitters("") == frozenset()
    assert decode_fast_transmitters(None) == frozenset()


def test_classify_out_of_scope() -> None:
    # no fast transmitter at pre -> not an ionotropic synapse, regardless of what post expresses
    assert classify_sign(frozenset(), {"WB:ach_c"}, _CAT) == (None, "out_of_scope")


def test_classify_excitatory_cation_only() -> None:
    assert classify_sign(frozenset({"acetylcholine"}), {"WB:ach_c"}, _CAT) == (
        "excitatory",
        "excitatory",
    )


def test_classify_inhibitory_anion_only() -> None:
    assert classify_sign(frozenset({"GABA"}), {"WB:gab_a"}, _CAT) == ("inhibitory", "inhibitory")


def test_classify_complex_both_polarities() -> None:
    # post expresses both a cation and an anion ACh receptor -> complex, sign unset
    assert classify_sign(frozenset({"acetylcholine"}), {"WB:ach_c", "WB:ach_a"}, _CAT) == (
        None,
        "complex",
    )


def test_classify_unpredicted_no_match() -> None:
    # pre is GABAergic but post expresses only ACh/Glu receptors -> no match -> unpredicted
    assert classify_sign(frozenset({"GABA"}), {"WB:ach_c", "WB:glu_c"}, _CAT) == (
        None,
        "unpredicted",
    )
    # post expresses nothing -> unpredicted
    assert classify_sign(frozenset({"acetylcholine"}), set(), _CAT) == (None, "unpredicted")


def test_classify_only_matching_transmitter_counts() -> None:
    # Glu input; post co-expresses an ACh cation receptor (irrelevant) + a Glu anion receptor.
    # Only the glutamate receptor is considered -> inhibitory (not complex).
    assert classify_sign(frozenset({"glutamate"}), {"WB:ach_c", "WB:glu_a"}, _CAT) == (
        "inhibitory",
        "inhibitory",
    )


def test_classify_excluded_receptor_ignored() -> None:
    # The only matching ACh receptor post expresses is the Fenyves-excluded presynaptic lgc-46;
    # it must not be counted, leaving no match -> unpredicted.
    assert classify_sign(frozenset({"acetylcholine"}), {"WB:pre_x"}, _CAT) == (None, "unpredicted")


def test_classify_cotransmission_mixes_to_complex() -> None:
    # pre co-transmits ACh + GABA; post expresses an ACh cation + a GABA anion receptor ->
    # both polarities across the two transmitters -> complex.
    assert classify_sign(frozenset({"acetylcholine", "GABA"}), {"WB:ach_c", "WB:gab_a"}, _CAT) == (
        None,
        "complex",
    )


def test_committed_catalog_shape() -> None:
    """The committed Fenyves gene table loads: 62 distinct WBGene, acc-4/lgc-46 excluded."""
    cat = read_receptor_catalog(DATA / "ionotropic_receptor_genes.csv")
    assert len(cat) == 62  # 64 rows collapse to 62 WBGene (co-referent acr-4=des-2, acr-13=lev-8)
    assert all(r.transmitter in ("acetylcholine", "glutamate", "GABA") for r in cat.values())
    assert all(r.ion in ("cation", "anion") for r in cat.values())
    excluded = {r.symbol for r in cat.values() if r.excluded}
    assert excluded == {"acc-4", "lgc-46"}


def test_committed_expression_loads() -> None:
    expr = read_receptor_expression(DATA / "ionotropic_receptor_expression.csv")
    assert len(expr) == 290  # distinct neurons with >=1 ionotropic receptor at threshold 4
    assert all(wb.startswith("WB:WBGene") for cell in expr.values() for wb in cell)
