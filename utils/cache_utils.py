# utils/cache_utils.py
from __future__ import annotations
import json, os, time

# Functions for saving computed results in a cache file

def _load_json(path: str) -> dict:
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return {}

def _save_json(path: str, obj: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)

def get_from_pair_cache(pair_dir: str, key: str):
    cache_path = os.path.join(pair_dir, "Analysis", "cache.json")
    data = _load_json(cache_path)
    return data.get(key, None)

def put_in_pair_cache(pair_dir: str, key: str, value) -> None:
    cache_path = os.path.join(pair_dir, "Analysis", "cache.json")
    data = _load_json(cache_path)
    data[key] = value
    data["_updated"] = int(time.time())
    _save_json(cache_path, data)
