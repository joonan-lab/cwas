"""
Manage CWAS environment variables
"""

from pathlib import Path

import dotenv
import yaml


# Ref: https://python-patterns.guide/gang-of-four/singleton/
class Singleton(object):
    _instance = None

    def __new__(cls, *args, **kwargs):
        if Singleton._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance


class Env(Singleton):
    def __init__(self):
        self.path = Path(__file__).parent.resolve() / '.env'
        if not self.path.exists():
            self.path.touch()
        self.env = dotenv.dotenv_values(dotenv_path=self.path)

    def get_env(self, env_key):
        return self.env.get(env_key)

    def set_env(self, env_key, env_value):
        self.env[env_key] = env_value

    def save(self):
        with self.path.open('w') as env_f:
            yaml.safe_dump(self.env, env_f)
