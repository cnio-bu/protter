from importlib.util import module_from_spec,spec_from_file_location
import os


def dataset_module_names():
    this_file_path = os.path.realpath(__file__)
    module_dir,this_file_name = os.path.split(this_file_path)
    return [x[:-len(".py")] for x in os.listdir(module_dir)
            if x.endswith(".py") and x != this_file_name]


def load_dataset_module(ds):
    mod_dir_path = os.path.dirname(os.path.realpath(__file__))
    mod_file_name = "{}.py".format(ds)
    mod_file_path = os.path.join(mod_dir_path,mod_file_name)

    if not os.path.isfile(mod_file_path):
        raise ModuleNotFoundError("No module found for dataset '{}'".format(ds))

    module_spec = spec_from_file_location(ds,mod_file_path)
    module = module_from_spec(module_spec)
    module_spec.loader.exec_module(module)

    return module
