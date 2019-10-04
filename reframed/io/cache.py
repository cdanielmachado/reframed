from .sbml import load_cbmodel


class ModelCache:
    """ Utility to keep a cache of models with a lazy loading approach. """

    def __init__(self, ids, paths, load_args=None, post_processing=None):

        self.paths = dict(zip(ids, paths))
        self.cache = dict()
        self.load_args = load_args if load_args is not None else {}
        self.post_processing = post_processing

    def get_ids(self):
        return list(self.paths.keys())

    def get_model(self, model_id, reset_id=False):

        if model_id not in self.paths:
            raise RuntimeError("Model not in list: " + model_id)

        if model_id in self.cache:
            return self.cache[model_id]

        model = load_cbmodel(self.paths[model_id], **self.load_args)

        self.cache[model_id] = model

        if self.post_processing is not None:
            self.post_processing(model)

        if reset_id:
            model.id = model_id

        return model
