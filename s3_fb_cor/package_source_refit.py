"""Package paired S3 source/A/B/q candidates without activating them."""
from pathlib import Path
import argparse
import shutil
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from s3_fb_cor.source_refit import read, write, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", type=Path, required=True)
    args = parser.parse_args()
    task = args.task.resolve()
    bundle = task / "source_calibration"
    bundle.mkdir(exist_ok=True)
    source_files = {}
    for period, source_run, name in [
        ("high1", "0460", "cali_si_r460.txt"),
        ("low", "0457", "cali_si_r457_low_current.txt"),
        ("high2", "0457", "cali_si_r457_high_current.txt")]:
        source = read(ROOT / f"cali_si/parameters/source_{source_run}_peak_centers_20261102/candidate.json")
        original = ROOT / "pars" / name
        rows = np.loadtxt(original, skiprows=1)
        keys = (10000 * rows[:, 0] + 100 * rows[:, 1] + rows[:, 2]).astype(int)
        before = rows.copy()
        changed = []
        for c in source["channels"].values():
            if c["detector"] not in ["S3 ring", "S3 sector"]:
                continue
            key = c["low_beam_key"] if period == "low" else c["key"]
            index = np.flatnonzero(keys == key)
            assert len(index) == 1 and c["window_checked"]
            rows[index[0], 3:6] = [c["c0"], c["c1"], c["c2"]]
            # The final field is a legacy source-fit chi value. Preserve it as
            # historical metadata; this new center-line fit has a different test.
            changed.append(int(index[0]))
        assert len(changed) == 56
        unchanged = sorted(set(range(len(rows))) - set(changed))
        assert np.array_equal(before[unchanged], rows[unchanged])
        output = bundle / name
        header = "cid sid ch p0 p1 p2 legacy_source_chi_not_refitted"
        np.savetxt(output, rows, header=header, comments="", fmt=["%d", "%d", "%d", "%.17g", "%.17g", "%.17g", "%.17g"])
        assert np.array_equal(np.loadtxt(output, skiprows=1), rows)
        source_files[period] = {"file": name, "sha256": sha(output), "baseline_sha256": sha(original),
            "source_run": source_run, "s3_channels_replaced": 56, "other_channels_preserved": len(unchanged)}
    package = read(task / "coefficients.json")
    package["required_source_calibration"] = source_files
    package["source_calibration_scope"] = "Only S3 source coefficients updated; other channels preserved. Use the paired source_calibration/ files together."
    write(task / "coefficients.json", package)
    fb = ROOT / "s3_fb_cor/parameters/source_recalibrated_20261102"
    att = ROOT / "s3_energy_attenuation/parameters/source_recalibrated_20261102"
    for destination in [fb, att]:
        destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(task / "coefficients.json", fb / "coefficients.json")
    shutil.copy2(task / "coefficients.csv", fb / "coefficients.csv")
    shutil.copytree(bundle, fb / "source_calibration", dirs_exist_ok=True)
    shutil.copy2(task / "position_factors.csv", att / "ring_position_factors.csv")
    shutil.copy2(task / "run_summary.csv", att / "run_summary.csv")
    write(att / "complete_package.json", {"path": "../../../s3_fb_cor/parameters/source_recalibrated_20261102/coefficients.json",
        "sha256": sha(fb / "coefficients.json"), "formal_replacement": False})
    readme = """# 新源刻度配套 S3 A/B/q 候选（2026-11-02）

这是 0244 的完整候选，尚未切换 `pars/active.json`。全部173轮；高流强 ring1–8，低流强 ring1–24。高流强外环系数仅用于完整结构，禁止作为生产范围启用。

`coefficients.json` 与 `source_calibration/` 内三个基础刻度文件必须配套；只更新S3源系数，Si array保持现用系数。文本最后一列保留历史源拟合χ值，不是本次两点一次刻度的χ值。

E0=hardware_gain*(c0+c1*ADC)，Er=q[matched_ring]*(Ar+Br*E0r)，Es=q[matched_ring]*(As+Bs*E0s)。不再叠加旧正背面和衰减链。q作用于匹配ring的两面，不强迫不同ring等能。

全量系数和采纳状态见CSV；重拟未通过者保留等价换算。新基础刻度不独立决定高能绝对能标，参考轮0472/0616/0680的现用物理能标固定。数据为E盘保存的全173轮分段ADC配对样本，本次未扫描全量C2。

程序：`s3_fb_cor/recalibrate_source.cpp`、`source_refit.py`、`s3_energy_attenuation/refit_position.py`；完整日志/方案/核验/拟合在 `codex/analysis/s3_source_calibration_refit_20261102/`；阅读版 `jupyter/codex_0244_150sm_s3_new_source_calibration_abq.ipynb`。
"""
    (fb / "README.md").write_text(readme)
    (att / "README.md").write_text("# 新源刻度下的逐 ring 衰减候选\n\n完整 A/B/q 与三个配套源刻度文件位于 `s3_fb_cor/parameters/source_recalibrated_20261102/`。这里维护逐 ring q、重拟候选、留出检查和选择状态。`complete_package.json`记录完整包及SHA；未启用正式参数。\n")
    files = list(fb.rglob("*")) + list(att.rglob("*"))
    write(task / "package_manifest.json", {"formal_replacement": False,
        "files": {str(p.relative_to(ROOT)): sha(p) for p in files if p.is_file()}, "source_files": source_files})
    print(f"Packaged {len(package['runs'])} runs; active parameters unchanged")


if __name__ == "__main__":
    main()
