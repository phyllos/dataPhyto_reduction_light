#!/usr/bin/env python3
"""
driver_light_01.py

Run the SACLA light-data TR-SFX/phytochrome preprocessing pipeline. (?)
The structure intentionally follows the dark-data driver style: YAML config,
step-by-step logging, and per-step memory monitoring.

Usage:
    python ./driver_light_01.py --config ./inputs/config_light_upto1ps_01.yaml

Fast test from an existing C7 MATLAB/Python packed file:
    python ./driver_light_01.py --config ./inputs/config_light_from_c7.yaml
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Callable

try:
    import resource
except ImportError:  # pragma: no cover
    resource = None

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML is required. Install it first:\n"
        "  conda install pyyaml\n"
        "or\n"
        "  pip install pyyaml"
    ) from exc


# ----------------------------- small helpers -----------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run light-data preprocessing from a YAML config.")
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    return parser.parse_args()


def load_config(config_path: Path) -> dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError("Config file must be a YAML mapping/object.")
    return data


def require_key(cfg: dict[str, Any], key: str) -> Any:
    if key not in cfg:
        raise KeyError(f"Missing required config key: {key}")
    return cfg[key]


def as_path(value: Any) -> Path:
    return Path(str(value)).expanduser().resolve()


def progress(message: str) -> None:
    print(f"[progress] {message}", flush=True)


def log_banner(title: str) -> None:
    print("=" * 72, flush=True)
    print(title, flush=True)
    print("=" * 72, flush=True)


def _read_proc_status_kb(key: str) -> int | None:
    try:
        with open("/proc/self/status", "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(key):
                    return int(line.split()[1])
    except FileNotFoundError:
        return None
    return None


def _kb_to_gb(kb: int) -> str:
    return f"{kb / 1024 / 1024:.2f} GB"


def current_rss_gb() -> str:
    kb = _read_proc_status_kb("VmRSS:")
    return "N/A" if kb is None else _kb_to_gb(kb)


def process_peak_rss_gb() -> str:
    kb = _read_proc_status_kb("VmHWM:")
    if kb is not None:
        return _kb_to_gb(kb)
    if resource is None:
        return "N/A"
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = usage.ru_maxrss
    if sys.platform == "darwin":
        return f"{rss / 1024 / 1024 / 1024:.2f} GB"
    return f"{rss / 1024 / 1024:.2f} GB"


class StepMemoryMonitor:
    def __init__(self, step_name: str, sample_interval: float = 5.0, print_interval: float = 30.0):
        self.step_name = step_name
        self.sample_interval = sample_interval
        self.print_interval = print_interval
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()
        self.start_time = 0.0
        self.last_print_time = 0.0
        self.step_peak_kb = 0
        self.latest_kb = 0

    def _current_kb(self) -> int:
        kb = _read_proc_status_kb("VmRSS:")
        return 0 if kb is None else kb

    def _run(self) -> None:
        while not self._stop.is_set():
            now = time.monotonic()
            kb = self._current_kb()
            with self._lock:
                self.latest_kb = kb
                self.step_peak_kb = max(self.step_peak_kb, kb)
                if (now - self.last_print_time) >= self.print_interval:
                    elapsed = int(now - self.start_time)
                    self.last_print_time = now
                    print(
                        f"[MONITOR] {self.step_name} | elapsed={elapsed}s | "
                        f"current={_kb_to_gb(self.latest_kb)} | "
                        f"step_peak={_kb_to_gb(self.step_peak_kb)} | "
                        f"process_peak={process_peak_rss_gb()}",
                        flush=True,
                    )
            if self._stop.wait(self.sample_interval):
                break

    def start(self) -> None:
        initial = self._current_kb()
        self.start_time = time.monotonic()
        self.last_print_time = self.start_time
        self.step_peak_kb = initial
        self.latest_kb = initial
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.sample_interval + 1.0))

    def step_peak_gb(self) -> str:
        with self._lock:
            return _kb_to_gb(self.step_peak_kb)


def run_step(name: str, func: Callable[..., Any], *args: Any, enabled: bool = True) -> Any:
    if not enabled:
        print(f"[SKIP ] {name}", flush=True)
        return None

    monitor = StepMemoryMonitor(name)
    print(
        f"[START] {name} | PID={os.getpid()} | "
        f"current={current_rss_gb()} | process_peak={process_peak_rss_gb()}",
        flush=True,
    )
    monitor.start()
    try:
        result = func(*args)
    except Exception:
        monitor.stop()
        print(
            f"[FAIL ] {name} | current={current_rss_gb()} | "
            f"step_peak={monitor.step_peak_gb()} | process_peak={process_peak_rss_gb()}",
            flush=True,
        )
        traceback.print_exc()
        raise
    finally:
        gc.collect()

    monitor.stop()
    print(
        f"[ END ] {name} | current={current_rss_gb()} | "
        f"step_peak={monitor.step_peak_gb()} | process_peak={process_peak_rss_gb()}",
        flush=True,
    )
    return result


# ------------------------------ config logic -----------------------------

def validate_config(raw: dict[str, Any], config_path: Path) -> dict[str, Any]:
    project_root = as_path(raw.get("project_root", config_path.parent.parent))
    functions_dir = as_path(raw.get("functions_dir", project_root / "Functions_python"))
    output = as_path(require_key(raw, "output"))
    output.mkdir(parents=True, exist_ok=True)

    start_from = str(raw.get("start_from", "C1")).upper()
    if start_from not in {"C1", "C8"}:
        raise ValueError("start_from must be either C1 or C8")

    cfg: dict[str, Any] = {
        "project_root": project_root,
        "functions_dir": functions_dir,
        "output": output,
        "start_from": start_from,
        "data_name": str(raw.get("data_name", "upto1ps")),
        "a": float(raw.get("a", 54.22)),
        "b": float(raw.get("b", 115.78)),
        "c": float(raw.get("c", 117.08)),
        "wavelength": float(raw.get("wavelength", 1.77)),
        "ss": int(raw.get("ss", 155)),
        "tmin": float(raw.get("tmin", -84.0)),
        "tmax": float(raw.get("tmax", 550.0)),
        "seed": None if raw.get("seed", None) is None else int(raw.get("seed")),
        "remove_negative_pixels": bool(raw.get("remove_negative_pixels", False)),
        "make_drl_mask": bool(raw.get("make_drl_mask", False)),
        "use_drl_mask": bool(raw.get("use_drl_mask", False)),
        "generate_hkl_avg": bool(raw.get("generate_hkl_avg", True)),
    }

    if not functions_dir.is_dir():
        raise FileNotFoundError(f"functions_dir not found: {functions_dir}")

    if start_from == "C1":
        for key in ["stream", "params", "delays"]:
            p = as_path(require_key(raw, key))
            if not p.is_file():
                raise FileNotFoundError(f"{key} file not found: {p}")
            cfg[key] = p
    else:
        c7_file = as_path(require_key(raw, "c7_file"))
        if not c7_file.is_file():
            raise FileNotFoundError(f"c7_file not found: {c7_file}")
        cfg["c7_file"] = c7_file

    return cfg


def add_functions_to_path(functions_dir: Path) -> None:
    sys.path.insert(0, str(functions_dir))


def stage_file(output: Path, data_name: str, suffix: str) -> str:
    return str(output / f"dataPhyto_{data_name}_{suffix}.hdf5")


def run_pipeline(cfg: dict[str, Any]) -> None:
    add_functions_to_path(cfg["functions_dir"])

    from C1_read_partialator_params_light import read_partialator_params_light
    from C2_read_stream_files_light import read_stream_files_light
    from C3_grab_hklI_sacla2021_light import grab_hklI_sacla2021_light
    from C4_find_hkl_max_light import find_hkl_max_light
    from C5_find_hkl_redundancy_light import find_hkl_redundancy_light
    from C6_pack_phyto_allInfo_light import pack_phyto_allInfo_light
    from C7_pack_phyto_final_light import pack_phyto_final_light
    from C8_generate_uniform_delay_sacla2021_hdf5 import generate_uniform_delay_sacla2021_hdf5
    from C9_pack_phyto_light_lpf import pack_phyto_light_lpf
    from C10_generate_maskDRL_light import generate_maskDRL_light
    from C11_pack_phyto_light_lpf_drl import pack_phyto_light_lpf_drl
    from C12_pack_phyto_light_lpf_drl_scl import pack_phyto_light_lpf_drl_scl
    from C13_pack_phyto_light_lpf_drl_scl_bst import pack_phyto_light_lpf_drl_scl_bst

    output: Path = cfg["output"]
    data_name = cfg["data_name"]

    log_banner("SACLA light-data preprocessing started")
    for key in ["project_root", "functions_dir", "output", "start_from", "data_name"]:
        print(f"{key}: {cfg[key]}", flush=True)
    print(f"a,b,c = {cfg['a']}, {cfg['b']}, {cfg['c']} A", flush=True)
    print(f"wavelength = {cfg['wavelength']} A", flush=True)
    print(f"uniform delay: ss={cfg['ss']}, tmin={cfg['tmin']}, tmax={cfg['tmax']}, seed={cfg['seed']}", flush=True)
    print(f"Initial current RSS = {current_rss_gb()}", flush=True)
    print(f"Initial process peak RSS = {process_peak_rss_gb()}", flush=True)

    record: dict[str, str] = {}

    if cfg["start_from"] == "C1":
        c1_file = run_step(
            "C1 read_partialator_params_light",
            read_partialator_params_light,
            str(cfg["params"]),
            str(output),
            data_name,
        )
        record["C1_file"] = str(c1_file)

        c2_file = run_step(
            "C2 read_stream_files_light",
            read_stream_files_light,
            str(cfg["stream"]),
            data_name,
            str(output),
        )
        record["C2_file"] = str(c2_file)

        hklI_folder = run_step(
            "C3 grab_hklI_sacla2021_light",
            grab_hklI_sacla2021_light,
            str(cfg["stream"]),
            data_name,
            str(output),
        )
        record["C3_hklI_folder"] = str(hklI_folder)

        hkl_max_file = run_step(
            "C4 find_hkl_max_light",
            find_hkl_max_light,
            str(hklI_folder),
            data_name,
            str(output),
        )
        record["C4_file"] = str(hkl_max_file)

        redundancy_file = run_step(
            "C5 find_hkl_redundancy_light",
            find_hkl_redundancy_light,
            str(hklI_folder),
            str(hkl_max_file),
            data_name,
            str(output),
        )
        record["C5_file"] = str(redundancy_file)

        allinfo_h5 = run_step(
            "C6 pack_phyto_allInfo_light",
            pack_phyto_allInfo_light,
            str(c2_file),
            str(c1_file),
            str(cfg["delays"]),
            data_name,
            str(output),
        )
        record["C6_file"] = str(allinfo_h5)
        allinfo_dat = str(Path(allinfo_h5).with_suffix(".dat"))

        c7_file = run_step(
            "C7 pack_phyto_final_light",
            pack_phyto_final_light,
            allinfo_dat,
            str(redundancy_file),
            str(hklI_folder),
            data_name,
            str(output),
            bool(cfg["remove_negative_pixels"]),
        )
        record["C7_file"] = str(c7_file)
    else:
        c7_file = str(cfg["c7_file"])
        record["C7_file"] = c7_file
        print(f"[INFO ] start_from=C8, using existing C7 file: {c7_file}", flush=True)

    c8_file = run_step(
        "C8 generate_uniform_delay_sacla2021_hdf5",
        generate_uniform_delay_sacla2021_hdf5,
        str(c7_file),
        stage_file(output, data_name, "C8_unifdelay"),
        int(cfg["ss"]),
        float(cfg["tmin"]),
        float(cfg["tmax"]),
        cfg["seed"],
    )
    record["C8_file"] = str(c8_file)

    c9_file = run_step(
        "C9 pack_phyto_light_lpf",
        pack_phyto_light_lpf,
        str(c8_file),
        stage_file(output, data_name, "C9_LPF"),
        float(cfg["a"]),
        float(cfg["b"]),
        float(cfg["c"]),
        float(cfg["wavelength"]),
    )
    record["C9_file"] = str(c9_file)

    mask_file = None
    if cfg["make_drl_mask"] or cfg["use_drl_mask"]:
        mask_file = run_step(
            "C10 generate_maskDRL_light",
            generate_maskDRL_light,
            str(c9_file),
            str(output / f"maskDRL_{data_name}.hdf5"),
            float(cfg["a"]),
            float(cfg["b"]),
            float(cfg["c"]),
        )
        record["C10_mask_file"] = str(mask_file)

    c11_file = run_step(
        "C11 pack_phyto_light_lpf_drl",
        pack_phyto_light_lpf_drl,
        str(c9_file),
        stage_file(output, data_name, "C11_LPF_DRL"),
        str(mask_file) if cfg["use_drl_mask"] and mask_file is not None else None,
    )
    record["C11_file"] = str(c11_file)

    c12_file = run_step(
        "C12 pack_phyto_light_lpf_drl_scl",
        pack_phyto_light_lpf_drl_scl,
        str(c11_file),
        stage_file(output, data_name, "C12_LPF_DRL_SCL"),
        float(cfg["a"]),
        float(cfg["b"]),
        float(cfg["c"]),
        bool(cfg["generate_hkl_avg"]),
    )
    record["C12_file"] = str(c12_file)

    c13_file = run_step(
        "C13 pack_phyto_light_lpf_drl_scl_bst",
        pack_phyto_light_lpf_drl_scl_bst,
        str(c12_file),
        stage_file(output, data_name, "C13_LPF_DRL_SCL_BST"),
    )
    record["C13_file"] = str(c13_file)

    log_banner("Pipeline complete")
    for k, v in record.items():
        print(f"{k}: {v}", flush=True)
    print(f"Final current RSS = {current_rss_gb()}", flush=True)
    print(f"Final process peak RSS = {process_peak_rss_gb()}", flush=True)


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        raise SystemExit(f"Config file not found: {config_path}")
    raw = load_config(config_path)
    cfg = validate_config(raw, config_path)
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
