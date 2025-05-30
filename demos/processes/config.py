
SBMLModelChangesConfig = {
    'species_changes': {
        'species_name': {
            'unit': 'maybe[string]',
            'initial_concentration': 'maybe[float]',
            'initial_particle_number': 'maybe[float]',
            'initial_expression': 'maybe[string]',
            'expression': 'maybe[string]'
        }
    },
    'global_parameter_changes': {
        'global_parameter_name': {
            'initial_value': 'maybe[float]',
            'initial_expression': 'maybe[string]',
            'expression': 'maybe[string]',
            'status': 'maybe[string]',
            'type': 'maybe[string]'  # (ie: fixed, assignment, reactions)
        }
    },
    'reaction_changes': {
        'reaction_name': {
            'parameters': {
                'reaction_parameter_name': 'maybe[int]'  # (new reaction_parameter_name value)  <-- this is done with set_reaction_parameters(name="(REACTION_NAME).REACTION_NAME_PARAM", value=VALUE)
            },
            'reaction_scheme': 'maybe[string]'   # <-- this is done like set_reaction(name = 'R1', scheme = 'S + E + F = ES')
        }
    }
}

SBMLModelConfig = {
    'model_id': 'string',
    'model_source': 'string',
    'model_language': {
        '_default': 'string',
        '_type': 'string'
    },
    'model_name': 'string',
    'model_units': 'tree[string]',
    'model_changes': SBMLModelChangesConfig  # 'tree[string]'
}

TimeCourseConfig = {
    'model': SBMLModelConfig,
    'time_config': 'tree[string]',  # ie: {start, stop, steps}
    'species_context': {
        '_default': 'concentrations',
        '_type': 'string'
    },
    'working_dir': 'string'
}