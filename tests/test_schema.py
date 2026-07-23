"""Phase 1: the LinkML schema is the single source of truth, so CI guards it.

These tests confirm the schema compiles, exposes the four core classes, and that the
hand-made sample instance validates against it.
"""

from pathlib import Path

import pytest
from linkml_runtime import SchemaView

SCHEMA = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "celegans_connectome_kg"
    / "schema"
    / "connectome.yaml"
)
SAMPLE = Path(__file__).resolve().parent / "data" / "sample_connectome.yaml"


@pytest.fixture(scope="module")
def view() -> SchemaView:
    return SchemaView(str(SCHEMA))


def test_schema_compiles_with_core_classes(view: SchemaView) -> None:
    classes = set(view.all_classes())
    assert {"Connectome", "Cell", "Connection", "Dataset", "Evidence"} <= classes


def test_connection_type_enum_covers_all_three(view: SchemaView) -> None:
    values = set(view.get_enum("ConnectionType").permissible_values)
    assert values == {"chemical", "gap_junction", "functional"}


def test_connectome_is_tree_root(view: SchemaView) -> None:
    assert view.get_class("Connectome").tree_root is True


def test_sex_enum_and_slots(view: SchemaView) -> None:
    # Sex enum covers both C. elegans sexes
    assert set(view.get_enum("Sex").permissible_values) == {"hermaphrodite", "male"}
    # Dataset carries the specimen sex (single-valued)
    dataset_slots = view.class_slots("Dataset")
    assert "sex" in dataset_slots
    assert view.get_slot("sex").range == "Sex"
    # Dataset carries the developmental life stage (dauer is a distinct value)
    assert "life_stage" in dataset_slots
    assert view.get_slot("life_stage").range == "LifeStage"
    life_stages = set(view.get_enum("LifeStage").permissible_values)
    assert {"L1", "L2", "dauer", "L3", "L4", "adult"} <= life_stages
    # Cell carries derived sex-presence (multivalued)
    cell_slots = view.class_slots("Cell")
    assert "sexes" in cell_slots
    sexes = view.get_slot("sexes")
    assert sexes.range == "Sex" and sexes.multivalued is True
    # Cell flags class-level placeholder endpoints (boolean)
    assert "unspecified" in cell_slots
    assert view.get_slot("unspecified").range == "boolean"
    # Reified per-sex neurotransmitter assignment (Wang 2024 male atlas)
    na_slots = view.class_slots("NeurotransmitterAssignment")
    assert {"cell", "sex", "neurotransmitter"} <= set(na_slots)
    assert view.get_class("NeurotransmitterAssignment").slot_usage["sex"].required is True


def test_sample_instance_validates() -> None:
    from linkml.validator import validate_file

    report = validate_file(str(SAMPLE), str(SCHEMA), target_class="Connectome")
    assert report.results == [], f"unexpected validation issues: {report.results}"
