"""Predict the excitatory/inhibitory sign of a chemical synapse (Part C, Fenyves 2020 method).

For a chemical synapse ``pre -> post`` the predicted polarity is derived natively, the same way
CIRCE derives its neuropeptide/monoamine connectomes from expression x pairs (see
``neuropeptide.derive_neuropeptide_edges``):

    T = fast neurotransmitter(s) of ``pre``; collect the ionotropic receptor genes for T that
    ``post`` expresses (CeNGEN); the matched genes being **cation-only -> excitatory**,
    **anion-only -> inhibitory**, **mixed -> complex** (sign unset), **none -> unpredicted**
    (sign unset). A ``pre`` with no fast transmitter is out of scope (no prediction).

Inputs are the two committed CSVs of ``data/cengen-ionotropic-receptors/`` (Fenyves S1-Table
polarity + CeNGEN threshold-4 expression); the cation=excitatory / anion=inhibitory mapping is the
paper's assumption. This module is pure logic + two readers, so the rule is unit-tested in
isolation; the build (``assemble``) applies it to the hermaphrodite chemical connections. The
receptor->cation/anion calls and the connectome signed here are an independent **validation** target
(Fenyves S1 Data), so the derivation deliberately does not ingest EleganSign's output.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

#: Fenyves DOI, recorded as each predicted sign's provenance source.
SOURCE = "https://doi.org/10.1371/journal.pcbi.1007974"

#: neuron-graph NT code letter -> fast ionotropic transmitter (the three Fenyves covers). Other
#: letters (d/s/o/t monoamines, b betaine, u unknown, n none) gate no fast ionotropic receptor in
#: this layer, so they contribute no transmitter and leave the synapse out of scope.
FAST_NT: dict[str, str] = {"a": "acetylcholine", "l": "glutamate", "g": "GABA"}


@dataclass(frozen=True)
class Receptor:
    """One ionotropic postsynaptic receptor gene (Fenyves S1 Table), with its polarity."""

    symbol: str
    wbgene: str  # "WB:WBGene..."
    transmitter: str  # acetylcholine / glutamate / GABA
    ion: str  # cation / anion
    excluded: bool  # Fenyves-flagged presynaptic (acc-4, lgc-46): not a sign predictor


def read_receptor_catalog(csv_path: Path) -> dict[str, Receptor]:
    """Read ``ionotropic_receptor_genes.csv`` into {wbgene: Receptor}.

    Keyed by WBGene so the two co-referent name pairs (acr-4=des-2, acr-13=lev-8; same WBGene,
    same transmitter/ion) collapse to one entry, matching how expression is matched by WBGene.
    """
    out: dict[str, Receptor] = {}
    with Path(csv_path).open(newline="") as fh:
        for r in csv.DictReader(fh):
            wb = r["wbgene"].strip()
            out[wb] = Receptor(
                symbol=r["symbol"].strip(),
                wbgene=wb,
                transmitter=r["transmitter"].strip(),
                ion=r["ion"].strip(),
                excluded=bool(r["excluded"].strip()),
            )
    return out


def read_receptor_expression(csv_path: Path) -> dict[str, set[str]]:
    """Read ``ionotropic_receptor_expression.csv`` into {cell: {wbgene expressed}}."""
    out: dict[str, set[str]] = {}
    with Path(csv_path).open(newline="") as fh:
        for r in csv.DictReader(fh):
            out.setdefault(r["cell"].strip(), set()).add(r["wbgene"].strip())
    return out


def decode_fast_transmitters(nt_code: str | None) -> frozenset[str]:
    """Decode a neuron-graph NT code (e.g. "ag") into its fast ionotropic transmitters.

    Co-transmitters are kept (``"ag"`` -> ACh + GABA); monoamine/unknown letters drop out.
    """
    return frozenset(FAST_NT[ch] for ch in (nt_code or "") if ch in FAST_NT)


def classify_sign(
    transmitters: frozenset[str],
    expressed: set[str],
    catalog: dict[str, Receptor],
) -> tuple[str | None, str]:
    """Predict a synapse's sign from presynaptic transmitter(s) and postsynaptic expression.

    ``transmitters`` is ``decode_fast_transmitters(NT(pre))``; ``expressed`` the WBGenes ``post``
    expresses; ``catalog`` the receptor polarity table. Returns ``(sign, outcome)`` where ``sign``
    is ``"excitatory"``/``"inhibitory"`` or ``None``, and ``outcome`` is one of
    ``excitatory``/``inhibitory``/``complex``/``unpredicted``/``out_of_scope``:

    * no fast transmitter at ``pre`` -> ``(None, "out_of_scope")`` -- not an ionotropic synapse.
    * matched receptors cation-only -> ``("excitatory", "excitatory")``.
    * matched receptors anion-only -> ``("inhibitory", "inhibitory")``.
    * matched receptors of both polarities -> ``(None, "complex")``.
    * fast transmitter but no matched receptor -> ``(None, "unpredicted")``.
    """
    if not transmitters:
        return None, "out_of_scope"
    ions = {
        catalog[wb].ion
        for wb in expressed
        if wb in catalog and not catalog[wb].excluded and catalog[wb].transmitter in transmitters
    }
    if ions == {"cation"}:
        return "excitatory", "excitatory"
    if ions == {"anion"}:
        return "inhibitory", "inhibitory"
    if ions:  # both cation and anion matched
        return None, "complex"
    return None, "unpredicted"
