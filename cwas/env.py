"""
Manage CWAS environment variables
"""

from collections import OrderedDict
from pathlib import Path
from typing import Optional

import dotenv


# Ref: https://python-patterns.guide/gang-of-four/singleton/
class Singleton(object):
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance


class Env(Singleton):
    def __init__(self, env_path: Path = None):
        if env_path is None:
            if hasattr(self, "path"):
                env_path = self.path
            else:
                env_path = Path.home() / ".cwas_config"

        self.set_path(env_path)
        self.env = dotenv.dotenv_values(dotenv_path=self.path)

    def set_path(self, path: Path):
        assert path is not None, "The 'path' argument cannot be None."
        old_path = getattr(self, "path", None)

        if old_path == path:
            return

        self.path = path

        if old_path is not None and old_path.exists():
            old_path.unlink()

        if not path.exists():
            path.touch()

    def get_env(self, env_key) -> Optional[str]:
        """Return None if the environment value does not exist"""
        return self.env.get(env_key)

    def set_env(self, env_key: str, env_value):
        self.env[env_key] = str(env_value).strip()

    def reset(self):
        """Make the 'env' attribute empty"""
        self.env = OrderedDict()

    def save(self):
        with self.path.open("w") as env_f:
            for k, v in self.env.items():
                print(f"{k}={v}", file=env_f)
