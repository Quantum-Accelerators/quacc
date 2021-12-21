from monty.serialization import loadfn
import os


def load_yaml_calc(file_path):
    config = loadfn(file_path)
    if "parent" in config:
        parent_config = load_yaml_calc(
            os.path.join(os.path.dirname(file_path), config["parent"])
        )
        for k, v in parent_config.items():
            if k not in config:
                config[k] = v
            elif isinstance(v, dict):
                v_new = config.get(k, {})
                v_new.update(v)
                config[k] = v_new
    return config
