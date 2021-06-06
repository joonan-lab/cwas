"""
Manage CWAS environment variables
"""

from pathlib import Path
import dotenv

_ENV_PATH = Path(__file__).parent.resolve() / '.env'


def init_env():
    _ENV_PATH.touch()


def env_exists() -> bool:
    return _ENV_PATH.exists()


def set_env(env_key: str, env_value: str):
    dotenv.set_key(_ENV_PATH, env_key, env_value)


def load_env() -> dict:
    return dotenv.dotenv_values(dotenv_path=_ENV_PATH)
