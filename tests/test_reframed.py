#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from reframed import *

TEST_MODEL = 'data/e_coli_core.xml.gz'
#TEST_MODEL = 'data/iML1515.xml.gz'
TEST_MODEL_COBRA = 'data/model_cobra.xml.gz'
TEST_MODEL_FBC2 = 'data/model_fbc2.xml.gz'

MIN_GROWTH = 0.1
REACTION_KO = 'R_PGI'

#set_default_solver("optlang")


class TestSBML(unittest.TestCase):
    """Test SBML import / export."""

    def test_load_cbmodel(self):
        load_cbmodel(TEST_MODEL, flavor='fbc2')

    def test_save_cbmodel_fbc2(self):
        model = load_cbmodel(TEST_MODEL)
        save_cbmodel(model, TEST_MODEL_FBC2, flavor='fbc2')

    def test_save_cbmodel_cobra(self):
        model = load_cbmodel(TEST_MODEL)
        save_cbmodel(model, TEST_MODEL_COBRA, flavor='cobra')

    def test_load_cbmodel_cobra(self):
        load_cbmodel(TEST_MODEL_COBRA, flavor='cobra')


class TestSimulation(unittest.TestCase):
    """Test SBML import / export."""

    def setUp(self):
        self.model = load_cbmodel(TEST_MODEL)
        self.obj_id = self.model.biomass_reaction

    def test_FBA(self):
        sol = FBA(self.model, constraints={REACTION_KO: 0})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_pFBA(self):
        sol = pFBA(self.model, constraints={REACTION_KO: 0})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_CAFBA(self):
        sol = CAFBA(self.model, constraints={REACTION_KO: 0})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_MOMA(self):
        sol = MOMA(self.model, constraints={REACTION_KO: 0, self.obj_id: (MIN_GROWTH, 10)})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_lMOMA(self):
        sol = lMOMA(self.model, constraints={REACTION_KO: 0, self.obj_id: (MIN_GROWTH, 10)})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_ROOM(self):
        sol = ROOM(self.model, constraints={REACTION_KO: 0, self.obj_id: (MIN_GROWTH, 10)})
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_FBAratio(self):
        self.model.add_ratio_constraint("R_PGI", "R_G6PDH2r", 1.0)
        sol = FBA(self.model)
        self.model.remove_ratio_constraint("R_PGI", "R_G6PDH2r")
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)
        self.assertAlmostEqual(sol.values["R_PGI"], sol.values["R_G6PDH2r"], 5)


class TestMedia(unittest.TestCase):

    def setUp(self):
        self.model = load_cbmodel(TEST_MODEL)
        self.obj_id = self.model.biomass_reaction

    def test_environment(self):
        env = Environment.complete(self.model, inplace=False)
        sol = FBA(self.model, constraints=env)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_minimal_medium(self):
        Environment.empty(self.model, inplace=True)
        reactions, _ = minimal_medium(self.model)
        new_env = Environment.from_reactions(reactions)
        sol = FBA(self.model, constraints=new_env)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)


class TestTransformations(unittest.TestCase):

    def setUp(self):
        self.model = load_cbmodel(TEST_MODEL)
        self.obj_id = self.model.biomass_reaction

    def test_make_irreversible(self):
        model = make_irreversible(self.model, inplace=False)
        sol = FBA(model)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_simplify(self):
        model = simplify(self.model, inplace=False)
        sol = FBA(model)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)


class TestCommunity(unittest.TestCase):

    def setUp(self):
        model1 = load_cbmodel(TEST_MODEL)
        model1.set_flux_bounds('R_ATPM', 0, 0)
        model1.id = 'm1'

        model2 = model1.copy()
        model2.id = 'm2'

#        model1.set_flux_bounds('R_NH4t', 0, 0)
#        model2.set_flux_bounds('R_GLCpts', 0, 0)

        self.community = Community('double', [model1, model2])
        self.medium = Environment.from_model(model1).get_compounds()

    def test_merge_models(self):
        self.model = self.community.merged_model

    def test_FBA(self):
        self.model = self.community.merged_model
        self.obj_id = self.model.biomass_reaction
        env = Environment.from_compounds(self.medium, fmt_func=lambda x: f"R_EX_M_{x}_e")
        env = env.apply(self.model, inplace=False)
        sol = FBA(self.model, constraints=env)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)

    def test_minimal_medium(self):
        self.model = self.community.merged_model
        self.obj_id = self.model.biomass_reaction
        Environment.empty(self.model, inplace=True)
        reactions, _ = minimal_medium(self.model)
        new_env = Environment.from_reactions(reactions)
        sol = FBA(self.model, constraints=new_env)
        self.assertGreater(sol.values[self.obj_id], MIN_GROWTH)