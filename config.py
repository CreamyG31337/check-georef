"""Settings for the tools: built-in defaults, overridden by config.toml.

config.toml lives next to the scripts (or at the path in the CHECK_GEOREF_CONFIG
environment variable) and is gitignored, so machine paths and project-specific values
stay out of the repository. config.example.toml documents every setting.
"""
import copy
import os
import tomllib

DEFAULTS = {
    "tools": {
        "las2las": r"C:\Program Files\Trimble\Trimble Business Center\Inpho_UAS\bin\las2las64.exe",
    },
    "bestfit": {"cell": 0.05, "max_samples": 40000, "search": 0.15, "shared_voxel": 0.20},
    "shift": {"suffix": " SHIFTED"},
    "control": {"files": [], "name_suffixes": []},
    "grid_compare": {
        "big_change": 1.0,
        "spot_separation": 15.0,
        "bands": [[-99, -1], [-1, -0.5], [-0.5, -0.2], [-0.2, -0.05], [-0.05, 0.05], [0.05, 0.2],
                  [0.2, 0.5], [0.5, 1], [1, 2], [2, 4], [4, 99]],
    },
}


def config_path():
    return os.environ.get("CHECK_GEOREF_CONFIG") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.toml")


def load():
    cfg = copy.deepcopy(DEFAULTS)
    path = config_path()
    if not os.path.exists(path):
        return cfg
    with open(path, "rb") as f:
        user = tomllib.load(f)
    for section, values in user.items():
        if section not in cfg or not isinstance(values, dict):
            raise ValueError(f"{path}: unknown section [{section}]")
        unknown = set(values) - set(cfg[section])
        if unknown:
            raise ValueError(f"{path}: unknown setting(s) in [{section}]: {', '.join(sorted(unknown))}")
        cfg[section].update(values)
    return cfg


CONFIG = load()
