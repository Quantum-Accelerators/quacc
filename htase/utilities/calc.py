from monty.serialization import loadfn
import os


def load_yaml_calc(file_path):
    supported_args = ["inputs", "parent", "parent_magmoms", "parent_setups"]
    config = loadfn(file_path)
    for config_arg in config:
        if config_arg not in supported_args:
            raise ValueError(f"Unsupported argument: {config_arg}")
    for config_arg in ["parent", "parent_magmoms", "parent_setups"]:
        if config_arg in config:
            parent_config = load_yaml_calc(
                os.path.join(os.path.dirname(file_path), config[config_arg])
            )
            for k, v in parent_config.items():
                if k not in config:
                    config[k] = v
                elif isinstance(v, dict):
                    v_new = config.get(k, {})
                    v_new.update(v)
                    config[k] = v_new
    return config
