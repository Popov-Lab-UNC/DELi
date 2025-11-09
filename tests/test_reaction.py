"""Tests for deli.dels.reaction module"""

import json
from pathlib import Path

import pytest
from rdkit import Chem

from deli.enumeration.reaction import (
    BBSetReactant,
    PooledReactant,
    ProductReactant,
    Reaction,
    ReactionError,
    ReactionParsingError,
    ReactionStep,
    ReactionTree,
    ReactionVial,
    StaticReactant,
)


def _mol(smiles: str):
    return Chem.MolFromSmiles(smiles)


def _smi(molobj):
    return Chem.MolToSmiles(molobj)


def test_reaction_simple_success():
    """Test making a reaction"""
    rxn = Reaction("[NH2:1].[C:2](=O)O>>[C:2](=O)[N:1]")
    # ammonia + acetic acid -> acetamide
    product = rxn.react(_mol("NC"), _mol("CC(=O)O"))
    assert _smi(product) == Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)NC"))


def test_reaction_wrong_number_of_reactants_raises():
    """Test for mismatched reactant count"""
    rxn = Reaction("[NH2:1].[C:2](=O)O>>[C:2](=O)[N:1]")
    with pytest.raises(ReactionError, match="Expected 2 reactants"):
        rxn.react(_mol("NC"))


def test_reaction_no_product_raises():
    """Test for no product on mismatched reactants"""
    rxn = Reaction("[NH2:1].[C:2](=O)O>>[C:2](=O)N")
    # pass two molecules that don't match the reactant templates
    with pytest.raises(ReactionError):
        rxn.react(_mol("CC"), _mol("CCC"))


def test_reaction_vial_and_static_reactant():
    """Test reaction vial and static reactant behavior"""
    vial = ReactionVial()
    m = _mol("CC")
    vial.add_product("foo", m)
    assert vial["foo"] is m

    static = StaticReactant(m)
    assert static.get_from_vial(vial) is m
    assert str(static).startswith("CC")


def test_bbset_and_product_reactants_get_from_vial_and_errors():
    """Test BB set and product reactant behavior and errors"""
    vial = ReactionVial()
    vial["DEL_A"] = _mol("C")
    bb = BBSetReactant("DEL_A")
    assert _smi(bb.get_from_vial(vial)) == "C"

    # missing bb set with subset id in message
    bb_subset = BBSetReactant("DEL_A:::subset1")
    with pytest.raises(ReactionError) as ei:
        bb_subset.get_from_vial(vial)
    assert "subset1" in str(ei.value)

    # product reactant missing
    prod = ProductReactant("stepX")
    with pytest.raises(ReactionError) as e:
        prod.get_from_vial(vial)
    assert "stepX" in str(e.value)

    # now add product and retrieve
    vial[prod.product_id] = _mol("CC")
    got = prod.get_from_vial(vial)
    assert _smi(got) == "CC"


def test_pooled_reactant_behavior_and_nested_error():
    """Test pooled reactant behavior and nested error"""
    vial = ReactionVial()
    vial["A"] = _mol("C")
    # pool with first missing then next available
    pool = PooledReactant([BBSetReactant("MISSING"), BBSetReactant("A")])
    got = pool.get_from_vial(vial)
    assert _smi(got) == "C"

    # pool with all missing
    pool2 = PooledReactant([BBSetReactant("X"), BBSetReactant("Y")])
    with pytest.raises(ReactionError):
        pool2.get_from_vial(vial)

    # nested pooled reactant not allowed
    with pytest.raises(ReactionParsingError):
        PooledReactant([PooledReactant([BBSetReactant("A")])])


def test_correct_number_of_reactants():
    """Test ReactionStep reactant count validation"""
    rxn1 = Reaction("[NH2:1]>>[NH2:1]")
    rxn2 = Reaction("[NH2:1].[C:2](=O)O>>[C:2](=O)N")

    # mismatched reactant count at init
    with pytest.raises(ReactionError):
        ReactionStep(
            step_name="s1",
            step_id="s1",
            reaction=rxn1,
            reactants=[StaticReactant(_mol("N")), StaticReactant(_mol("O"))],
        )

    with pytest.raises(ReactionError):
        ReactionStep(
            step_name="s2",
            step_id="s2",
            reaction=[rxn1, rxn2],
            reactants=[StaticReactant(_mol("N")), StaticReactant(_mol("CC(=O)O"))],
        )


def test_reaction_priority():
    """Test that ReactionStep picks the first matching reaction"""
    rxn1 = Reaction("[NH3:1].[CH3:2]>>[N:1][C:2]")
    rxn2 = Reaction("[NH2:1].[CH3:2]>>[N:1][C:2]")

    step = ReactionStep(
        step_name="test_step",
        step_id="test_step",
        reaction=[rxn1, rxn2],
        reactants=[BBSetReactant("TEST"), StaticReactant(_mol("CC(C=O)C(C=O)C(C=O)"))],
    )

    vial = ReactionVial()
    vial["TEST"] = _mol("N")  # ethylamine
    step.run_step(vial)
    assert _smi(vial["product_test_step"]) == _smi(_mol("NCC(C=O)C(C=O)C(C=O)"))  # should use rxn1

    vial = ReactionVial()
    vial["TEST"] = _mol("NCC")  # ethylamine
    step.run_step(vial)
    assert _smi(vial["product_test_step"]) == _smi(
        _mol("CCNCC(C=O)C(C=O)C(C=O)")
    )  # should use rxn2


DATA_DIR = Path("tests/test_data/reaction_test_data")

# Map filename -> expected exception type (or None for success)
TEST_CASES = {
    "valid_reaction.json": None,
    "valid_reaction_with_subsets_and_pools.json": None,
    "duplicate_BB_set.json": ReactionParsingError,
    "duplicated_rxn_id.json": ReactionParsingError,
    "invalid_rxn_smarts.json": Exception,  # RDKit may raise AttributeError/ValueError via Reaction
    "invalid_static_reactant.json": ReactionParsingError,
    "nonexistant_product.json": ReactionParsingError,
    "pool_nonexistant_reactant.json": ReactionParsingError,
    "reactants_nested_pool.json": ReactionParsingError,
    "reactants_type_invalid.json": ReactionParsingError,
}


@pytest.mark.parametrize("fname, expected", list(TEST_CASES.items()))
def test_reaction_tree_json_files(fname, expected):
    """Test loading ReactionTree from JSON test data files."""
    path = DATA_DIR / fname
    assert path.exists(), f"test data file {path} missing"

    raw = json.loads(path.read_text())
    possible = frozenset(raw.get("possible_bb_set_ids", []))
    reaction_steps = raw.get("reaction_steps", {})

    if expected is None:
        # should load successfully
        tree = ReactionTree.load_from_dict(reaction_steps, possible)
        assert hasattr(tree, "threads")
        assert len(tree.threads) >= 1
    else:
        with pytest.raises(expected):
            ReactionTree.load_from_dict(reaction_steps, possible)
