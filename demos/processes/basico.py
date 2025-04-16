"""
Basico/COPASI process for use with process bigraph.
"""
import os
import COPASI
from basico import (load_model, get_species, get_reactions, run_time_course)

from process_bigraph import Process, Composite, ProcessTypes
from process_bigraph.emitter import gather_emitter_results


def _set_initial_concentrations(changes, dm):
    model = dm.getModel()
    assert (isinstance(model, COPASI.CModel))

    references = COPASI.ObjectStdVector()

    for name, value in changes:
        species = model.getMetabolite(name)
        assert (isinstance(species, COPASI.CMetab))
        if species is None:
            print(f"Species {name} not found in model")
            continue
        species.setInitialConcentration(value)
        references.append(species.getInitialConcentrationReference())

    model.updateInitialValues(references)


def _get_transient_concentration(name, dm):
    model = dm.getModel()
    assert (isinstance(model, COPASI.CModel))
    species = model.getMetabolite(name)
    assert (isinstance(species, COPASI.CMetab))
    if species is None:
        print(f"Species {name} not found in model")
        return None
    return species.getConcentration()


class CopasiProcess(Process):
    """ODE component of the dfba hybrid using COPASI(basico). TODO: Generalize this to use any ode sim."""
    config_schema = {
        'model_path': 'string'
    }

    def __init__(self, config=None, core=None):
        super().__init__(config, core)
        print("Current working directory:", os.getcwd())
        print("Expected full path:", os.path.abspath(self.config['model_path']))
        print("File exists?", os.path.exists(self.config['model_path']))

        # # get the model path from the config
        # here = os.path.dirname(os.path.abspath(__file__))
        # model_path = str(os.path.join(here, '..', *self.config['model_path'].split('/')))

        # load the model
        self.model = load_model(self.config['model_path'])
        self.reaction_names = get_reactions(model=self.model).index.tolist()
        self.species_names = get_species(model=self.model).index.tolist()

    def initial_state(self):
        initial_concentrations = {
            species_name: get_species(species_name, model=self.model).initial_concentration[0]
            for species_name in self.species_names
        }

        initial_derivatives = {
            rxn_id: get_reactions(rxn_id, model=self.model).flux[0]
            for rxn_id in self.reaction_names
        }

        return {
            'species_concentrations': initial_concentrations,
            'reaction_fluxes': initial_derivatives,
            'time': 0.0
        }

    def inputs(self):
        concentrations_type = {
            name: 'float' for name in self.species_names
        }
        return {
            # 'species_concentrations': concentrations_type,
            'species_counts': {
                name: 'float' for name in self.species_names
            },
            'time': 'float'
        }

    def outputs(self):
        concentrations_type = {
            name: 'float' for name in self.species_names
        }

        reaction_fluxes_type = {
            reaction_name: 'float' for reaction_name in self.reaction_names
        }

        return {
            'species_concentrations': concentrations_type,
            'reaction_fluxes': reaction_fluxes_type,
            'time': 'float'
        }

    def update(self, inputs, interval):

        # get the new external states (18 states) from the ports
        # and set them in the model
        changes = []
        for side_id, molecules in inputs['species_counts'].items():
            for mol_id, value in molecules.items():
                if mol_id.endswith('_ext'):
                    changes.append((mol_id, value))

        _set_initial_concentrations(changes, self.copasi_model_object)

        # # spec_data_k = inputs['species_concentrations']
        # spec_data_k = inputs['species_counts']
        # for cat_id, value in spec_data_k.items():
        #     # set_type = 'concentration'
        #     set_type = "count"
        #     species_config = {
        #         'name': cat_id,
        #         'model': self.model,
        #         set_type: value
        #     }
        #     set_species(**species_config)

        # run model for "interval" length; we only want the state at the end
        tc = run_time_course(
            start_time=inputs['time'],
            duration=interval,
            update_model=True,
            model=self.model
        )

        results = {'time': interval}
        results['species_concentrations'] = {
            mol_id: float(get_species(
                name=mol_id,
                exact=True,
                model=self.model
            ).concentration[0])
            for mol_id in self.species_names
        }

        results['reaction_fluxes'] = {
            rxn_id: float(get_reactions(
                name=rxn_id,
                model=self.model
            ).flux[0])
            for rxn_id in self.reaction_names
        }

        return results


def run_basico(core):
    model_file = 'models/BIOMD0000000035_url.xml'

    spec = {
        'oscillator': {
            '_type': 'process',
            'address': 'local:basico',
            'config': {'model': model_file},
            'inputs': {
                'species_counts': ['species_counts'],
                'time': ['time'],
            },
            'outputs': {
                'species_concentrations': ['species_concentrations'],
                'reaction_fluxes': ['reaction_fluxes'],
                'time': ['time'],
            },
        },
        "emitter": {
            "_type": "step",
            "address": "local:ram-emitter",
            "config": {
                "emit": {
                    "species_concentrations": "any",
                    "global_time": "any",
                }
            },
            "inputs": {
                "species_concentrations": ['species_concentrations'],
                "global_time": ["global_time"]
            }
        }
    }

    # make the simulation
    sim = Composite({'state': spec}, core=core)

    # run the simulation
    sim.run(10)

    # get the data
    results = gather_emitter_results(sim)[("emitter",)]
    print(results)



if __name__ == '__main__':
    core = ProcessTypes()
    core.register('copasi', CopasiProcess)
    run_basico(core)
