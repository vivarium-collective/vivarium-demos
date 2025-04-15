from basico import biomodels, load_model_from_string


def fetch_biomodel(model_id: str):
    # TODO: make this generalizable for those other than basico
    sbml = biomodels.get_content_for_model(model_id)
    return load_model_from_string(sbml)