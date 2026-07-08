"""Phase 4 tests: load a small RDF graph into Oxigraph and run the count + sample queries."""

import pytest

from celegans_connectome_kg.build.datamodel import datamodel
from celegans_connectome_kg.export.rdf import to_turtle
from celegans_connectome_kg.verify.sparql import (
    count_summary,
    load_turtle_data,
    rows,
    sample_queries,
)


@pytest.fixture(scope="module")
def store():
    dm = datamodel()
    connectome = dm.Connectome(
        cells=[
            dm.Cell(
                id="cckg:cell/AVAL",
                name="AVAL",
                cell_type="neuron",
                anatomy="WBbt:0000001",
                sexes=["hermaphrodite", "male"],
            ),
            dm.Cell(id="cckg:cell/DA6", name="DA6", cell_type="neuron", sexes=["hermaphrodite"]),
            dm.Cell(
                id="cckg:cell/BWM-DL01",
                name="BWM-DL01",
                cell_type="muscle",
                sexes=["hermaphrodite"],
            ),
        ],
        datasets=[
            dm.Dataset(id="cckg:dataset/d1", name="D1", sex="hermaphrodite"),
            dm.Dataset(id="cckg:dataset/d2", name="D2", sex="male"),
        ],
        connections=[
            _conn(dm, "c1", "AVAL", "DA6", "chemical", 11.0, "d1"),
            _conn(dm, "c2", "AVAL", "DA6", "gap_junction", 2.0, "d1"),
            _conn(dm, "c3", "DA6", "AVAL", "gap_junction", 1.0, "d1"),
            _conn(dm, "c4", "AVAL", "DA6", "functional", 3.0, "d1"),
            _conn(dm, "c5", "AVAL", "DA6", "chemical", 5.0, "d2"),  # male-dataset edge
        ],
    )
    return load_turtle_data(to_turtle(connectome))


def _conn(dm, cid, pre, post, ctype, weight, ds):
    return dm.Connection(
        id=f"cckg:conn/{cid}",
        pre=f"cckg:cell/{pre}",
        post=f"cckg:cell/{post}",
        connection_type=ctype,
        weight=weight,
        dataset=f"cckg:dataset/{ds}",
    )


def test_count_summary(store) -> None:
    assert count_summary(store) == {
        "cells": 3,
        "cells_with_anatomy": 1,
        "datasets": 2,
        "connections": 5,
        "chemical": 2,
        "gap_junction": 2,  # both directions stored faithfully in RDF (merge is viz-only)
        "functional": 1,
    }


def test_sample_queries_present() -> None:
    queries = sample_queries()
    assert {
        "cells_by_type",
        "connections_by_type",
        "ungrounded_cells",
        "strongest_chemical_outputs",
        "cells_by_sex",
        "connections_by_sex",
        "shared_neuron_sex_partition",
    } <= set(queries)


def test_cells_by_sex_query(store) -> None:
    result = {r["sex"]: int(r["count"]) for r in rows(store, sample_queries()["cells_by_sex"])}
    assert result == {"hermaphrodite": 3, "male": 1}  # AVAL is the only both-sex cell


def test_connections_by_sex_query(store) -> None:
    result = {
        r["sex"]: int(r["count"]) for r in rows(store, sample_queries()["connections_by_sex"])
    }
    assert result == {"hermaphrodite": 4, "male": 1}


def test_shared_neuron_sex_partition(store) -> None:
    # A neuron present in both sexes (AVAL) has its edges cleanly partitioned by dataset sex.
    result = {
        r["sex"]: int(r["count"])
        for r in rows(store, sample_queries()["shared_neuron_sex_partition"])
    }
    assert result == {"hermaphrodite": 4, "male": 1}


def test_cells_by_type_query(store) -> None:
    result = {
        r["cell_type"]: int(r["count"]) for r in rows(store, sample_queries()["cells_by_type"])
    }
    assert result == {"neuron": 2, "muscle": 1}


def test_ungrounded_cells_query(store) -> None:
    names = {r["name"] for r in rows(store, sample_queries()["ungrounded_cells"])}
    assert names == {"DA6", "BWM-DL01"}  # AVAL has anatomy, so excluded


def test_strongest_chemical_outputs_query(store) -> None:
    result = rows(store, sample_queries()["strongest_chemical_outputs"])
    # AVAL->DA6 chemical summed across datasets: 11 (herm) + 5 (male) = 16
    assert result == [{"post_name": "DA6", "total_weight": "16"}]
