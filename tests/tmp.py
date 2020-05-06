from reframed import from_cobrapy, save_cbmodel
from cobra.io import load_json_model

cwd = '/Volumes/Data/Dropbox/embl/foodtranscriptomics/models/pangenomic/'

cbmodel = load_json_model(cwd + "curated.json")
model = from_cobrapy(cbmodel)
save_cbmodel(model, cwd+"panmodel.xml", flavor="bigg")
