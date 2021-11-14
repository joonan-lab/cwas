"""
Manage CWAS environment variables
"""

from collections import OrderedDict
from pathlib import Path
from typing import Any, Optional

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
                env_path = Path.home() / ".cwas_env"

        self.path = env_path
        self.load_env_from_file()

    def get_path(self) -> Path:
        return self.path

    def set_path(self, path: Path):
        """Set the dotenv path"""
        assert path is not None, "The 'path' argument cannot be None."
        self.path = path

    def get_env(self, env_key: str) -> Optional[str]:
        """Return None if the environment value does not exist"""
        return self.env.get(env_key)

    def set_env(self, env_key: str, env_value: Any):
        """Set a new value to the environment variable (key)"""
        self.env[env_key] = str(env_value).strip()

    def reset(self):
        """Make the 'env' attribute empty"""
        self.env = OrderedDict()

    def save(self):
        """Make a new dotenv"""
        with self.path.open("w") as env_f:
            for k, v in self.env.items():
                print(f"{k}={v}", file=env_f)

    def remove_file(self):
        if self.path.exists():
            self.path.unlink()

    def load_env_from_file(self):
        self.env = dotenv.dotenv_values(dotenv_path=self.path)

    def load_env_to_os(self):
        dotenv.load_dotenv(dotenv_path=self.path)
