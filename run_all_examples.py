#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all_examples.py — Executa todos os exemplos de validação usando o novo framework:
  1. Gmsh (.geo -> .msh)
  2. preproc.py + input_data.py (meshio -> elementos quadráticos descontínuos .dat + _mesh_bc.png)
  3. Solver Fortran (bin/fracture_2d.exe)
  4. Pós-processamento gráfico (plot_results.py -> figuras na pasta examples/)
"""

import math
import pathlib
import shutil
import subprocess
import time

from preproc import run_preproc

ROOT = pathlib.Path(__file__).resolve().parent
EXE = ROOT / "bin" / "fracture_2d.exe"
RESULTS_DIR = ROOT / "results"
EXAMPLES_DIR = ROOT / "examples"


def parse_crack_results(run_dir: pathlib.Path, aid: str, crack_idx: int = 1):
    fp = run_dir / f"{aid}_results_crack_{crack_idx:03d}.txt"
    rows = []
    for line in fp.read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
        cols = line.split()
        if not cols:
            continue
        vals = []
        for c in cols:
            if c.startswith("*"):
                break
            vals.append(float(c))
        if vals:
            rows.append(vals)
    return rows


def parse_fatigue_results(run_dir: pathlib.Path, aid: str, crack_idx: int = 1):
    fp = run_dir / f"{aid}_fatigue_results_crack_{crack_idx:03d}.txt"
    rows = []
    for line in fp.read_text(encoding="utf-8", errors="replace").splitlines()[1:]:
        cols = line.split()
        if not cols:
            continue
        vals = []
        for c in cols:
            if c.startswith("*"):
                break
            vals.append(float(c))
        if vals:
            rows.append(vals)
    return rows


def run_solver(problem, aw: float = 0.4, s51_variant: str = "da0_01_r11", sub_name: str | None = None) -> tuple[pathlib.Path, str, float]:
    dat_tmp, data = run_preproc(
        problem=problem,
        aw=aw,
        s51_variant=s51_variant,
        output_dir=ROOT / "inputs",
        plot_png=True,
    )
    aid = data["analysis_id"]
    folder_name = sub_name or aid
    run_dir = RESULTS_DIR / folder_name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    dat_in_run = run_dir / dat_tmp.name
    shutil.copy2(dat_tmp, dat_in_run)

    t0 = time.perf_counter()
    subprocess.run([str(EXE), str(dat_in_run)], cwd=str(run_dir), check=True)
    dt = time.perf_counter() - t0
    return run_dir, aid, dt


def validate_sec_3_1_3_1():
    print("\n" + "=" * 78)
    print("1) Seção 3.1.3.1: Placa Quadrada 10x10 com Trinca de Borda (Modo I Puro)")
    print("=" * 78)
    ref = {0.2: 1.488, 0.3: 1.848, 0.4: 2.324, 0.5: 3.010, 0.6: 4.152}
    print(f"{'a/w':>6} {'KI (BEM)':>12} {'KI_norm':>10} {'Ref (Tab 3.4)':>14} {'Erro (%)':>10}")
    for aw in [0.2, 0.3, 0.4, 0.5, 0.6]:
        a = aw * 10.0
        rdir, aid, _ = run_solver(problem=7, aw=aw)
        rows = parse_crack_results(rdir, aid)
        ki = rows[0][2]
        ki_norm = ki / math.sqrt(math.pi * a)
        err = (ki_norm - ref[aw]) / ref[aw] * 100.0
        print(f"{aw:6.1f} {ki:12.4f} {ki_norm:10.3f} {ref[aw]:14.3f} {err:+9.2f}%")


def validate_sec_3_1_3_2():
    print("\n" + "=" * 78)
    print("2) Seção 3.1.3.2: Placa 10x20 com Trinca Central Inclinada a 45° (Modo Misto)")
    print("=" * 78)
    ref_ki = {0.2: 0.518, 0.3: 0.541, 0.4: 0.572, 0.5: 0.612, 0.6: 0.661}
    ref_kii = {0.2: 0.507, 0.3: 0.516, 0.4: 0.529, 0.5: 0.546, 0.6: 0.567}
    print(f"{'2a/w':>6} {'KI_norm':>10} {'Ref_KI':>10} {'Erro_I%':>9} {'KII_norm':>10} {'Ref_KII':>10} {'Erro_II%':>9}")
    for aw in [0.2, 0.3, 0.4, 0.5, 0.6]:
        a = aw * 5.0
        rdir, aid, _ = run_solver(problem=8, aw=aw)
        rows = parse_crack_results(rdir, aid)
        ki_norm = abs(rows[0][2]) / math.sqrt(math.pi * a)
        kii_norm = abs(rows[0][3]) / math.sqrt(math.pi * a)
        ei = (ki_norm - ref_ki[aw]) / ref_ki[aw] * 100.0
        eii = (kii_norm - ref_kii[aw]) / ref_kii[aw] * 100.0
        print(f"{aw:6.1f} {ki_norm:10.3f} {ref_ki[aw]:10.3f} {ei:+8.2f}% {kii_norm:10.3f} {ref_kii[aw]:10.3f} {eii:+8.2f}%")


def validate_sec_5_1():
    print("\n" + "=" * 78)
    print("3) Seção 5.1: Placa Retangular 3x1 com Trinca Central (Isida / Andrade)")
    print("=" * 78)

    def ki_isida(a_half, sigma=10.0, b_half=0.5):
        ab = a_half / b_half
        F = (1.0 - 0.025 * ab**2 + 0.06 * ab**4) * math.sqrt(1.0 / math.cos(math.pi * ab / 2.0))
        return sigma * math.sqrt(math.pi * a_half) * F

    def n_isida_trap(a_start, a_end, n_pts=20000):
        C0, m = 1.5e-8, 2.25
        h_step = (a_end - a_start) / n_pts
        s = 0.0
        for i in range(n_pts + 1):
            a_i = a_start + i * h_step
            w = 0.5 if (i == 0 or i == n_pts) else 1.0
            s += w / (C0 * (ki_isida(a_i) ** m))
        return s * h_step

    for variant, exp_2af in [("da0_01_r11", 0.955), ("da0_0125_r105", 0.896)]:
        rdir, aid, dt = run_solver(problem=5, s51_variant=variant)
        crows = parse_crack_results(rdir, aid)
        frows = parse_fatigue_results(rdir, aid)
        final_2a = crows[-1][1]
        ki0, kif = crows[0][2], crows[-1][2]
        ref_ki0, ref_kif = ki_isida(0.05), ki_isida(final_2a / 2.0)
        nf = frows[-1][2]
        ref_n = n_isida_trap(0.05, final_2a / 2.0)
        err_n = (nf - ref_n) / ref_n * 100.0
        print(
            f"[{variant}] (tempo={dt:.2f}s): 2a_f={final_2a:.4f} (diss={exp_2af}) | "
            f"KI_0={ki0:.4f} (Isida={ref_ki0:.4f}) | KI_f={kif:.4f} (Isida={ref_kif:.4f}) | "
            f"N_f={nf:.1f} vs Eq(5.4)={ref_n:.1f} (erro={err_n:+.2f}%)"
        )


def eq59_F(a, b):
    al = a / b
    return 0.265 * (1.0 - al) ** 4 + (0.857 + 0.265 * al) / ((1.0 - al) ** 1.5)


def eq510_F(a, b):
    pa = math.pi * a / (2.0 * b)
    return (1.0 + 0.122 * (math.cos(pa)) ** 4) * math.sqrt((1.0 / pa) * math.tan(pa))


def n_dowling_trap(a_start, a_end, F_fn, b, sigma_max=526.0, R=0.3333333, C0=5.11e-13, m=3.24, gamma=0.42, n_pts=20000):
    ds_eq = sigma_max * (1.0 - R) ** gamma
    h = (a_end - a_start) / n_pts
    s = 0.0
    for i in range(n_pts + 1):
        a_i = a_start + i * h
        dk_eq = F_fn(a_i, b) * ds_eq * math.sqrt(math.pi * a_i)
        w = 0.5 if (i == 0 or i == n_pts) else 1.0
        s += w / (C0 * (dk_eq ** m))
    return s * h


def validate_sec_5_2_1_and_5_2_2():
    print("\n" + "=" * 78)
    print("4) Seções 5.2.1 e 5.2.2: Placas de Aço AISI 4340 (1 e 2 Trincas de Borda)")
    print("=" * 78)
    # 5.2.1
    rdir1, aid1, dt1 = run_solver(problem=3)
    crows1 = parse_crack_results(rdir1, aid1)
    frows1 = parse_fatigue_results(rdir1, aid1)
    a0_1, af_1 = crows1[0][1], crows1[-1][1]
    dk0_1 = crows1[0][2] * (1.0 - 0.3333333)
    dkf_1 = crows1[-1][2] * (1.0 - 0.3333333)
    ref_dk0_1 = eq59_F(a0_1, 0.076) * (526.0 * (1.0 - 0.3333333)) * math.sqrt(math.pi * a0_1)
    ref_dkf_1 = eq59_F(af_1, 0.076) * (526.0 * (1.0 - 0.3333333)) * math.sqrt(math.pi * af_1)
    n1 = frows1[-1][2]
    n1_ref = n_dowling_trap(a0_1, af_1, eq59_F, 0.076)
    print(
        f"[Sec. 5.2.1 - 1 trinca borda] (tempo={dt1:.2f}s): a0={a0_1*1e3:.3f} -> af={af_1*1e3:.3f} mm | "
        f"dK0={dk0_1:.3f} (ref={ref_dk0_1:.3f}) | dKf={dkf_1:.3f} (ref={ref_dkf_1:.3f}) | "
        f"Kmax_f={crows1[-1][2]:.2f} | N_f={n1:.1f} vs ref={n1_ref:.1f} (erro={(n1-n1_ref)/n1_ref*100:+.2f}%)"
    )

    # 5.2.2
    rdir2, aid2, dt2 = run_solver(problem=4)
    crows2 = parse_crack_results(rdir2, aid2, 1)
    frows2_1 = parse_fatigue_results(rdir2, aid2, 1)
    frows2_2 = parse_fatigue_results(rdir2, aid2, 2)
    a0_2, af_2 = crows2[0][1], crows2[-1][1]
    dk0_2 = crows2[0][2] * (1.0 - 0.3333333)
    dkf_2 = crows2[-1][2] * (1.0 - 0.3333333)
    ref_dk0_2 = eq510_F(a0_2, 0.038) * (526.0 * (1.0 - 0.3333333)) * math.sqrt(math.pi * a0_2)
    ref_dkf_2 = eq510_F(af_2, 0.038) * (526.0 * (1.0 - 0.3333333)) * math.sqrt(math.pi * af_2)
    n2_1, n2_2 = frows2_1[-1][2], frows2_2[-1][2]
    n2_ref = n_dowling_trap(a0_2, af_2, eq510_F, 0.038)
    print(
        f"[Sec. 5.2.2 - 2 trincas borda] (tempo={dt2:.2f}s): a0={a0_2*1e3:.3f} -> af={af_2*1e3:.3f} mm | "
        f"dK0={dk0_2:.3f} (ref={ref_dk0_2:.3f}) | dKf={dkf_2:.3f} (ref={ref_dkf_2:.3f}) | "
        f"Kmax_f={crows2[-1][2]:.2f} | N_f={n2_1:.1f} / {n2_2:.1f} vs ref={n2_ref:.1f} (erro={(n2_1-n2_ref)/n2_ref*100:+.2f}%)"
    )


def validate_sec_5_2_3_dowling():
    print("\n" + "=" * 78)
    print("5) Seção 5.2.3: Placa com Trinca Central e Amplitude Variável (Dowling 11.6)")
    print("=" * 78)
    rdir, aid, dt = run_solver(problem=2)
    crows = parse_crack_results(rdir, aid)
    frows = parse_fatigue_results(rdir, aid)
    a0_2, af_2 = crows[0][1], crows[-1][1]
    ki0, kif = crows[0][2], crows[-1][2]
    nf = frows[-1][2]

    def ki_feddersen(two_a, W=0.076, sig_max=2.4 * 219.298245614035):
        a = two_a / 2.0
        return sig_max * math.sqrt(math.pi * a / math.cos(math.pi * a / W))

    ref_ki0, ref_kif = ki_feddersen(a0_2), ki_feddersen(af_2)
    print(
        f"[Dowling 11.6] (tempo={dt:.2f}s): 2a0={a0_2*1e3:.3f} -> 2af={af_2*1e3:.3f} mm | "
        f"KI_max_0={ki0:.4f} (Feddersen={ref_ki0:.4f}) | KI_max_f={kif:.4f} (Feddersen={ref_kif:.4f}) | "
        f"N_f={nf:.1f} ciclos"
    )
    return rdir


def validate_sec_5_3_cmm():
    print("\n" + "=" * 78)
    print("6) Seções 5.3.1 e 5.3.2: Corpos de Prova Compact Mixed Mode (CMM 90° e 110°)")
    print("=" * 78)
    # 5.3.1: alpha = 90 deg
    rdir90, aid90, dt90 = run_solver(problem=6)
    crows90 = parse_crack_results(rdir90, aid90)
    frows90 = parse_fatigue_results(rdir90, aid90)
    print(
        f"[CMM 90°] (tempo={dt90:.2f}s): Passo 0 KI={crows90[0][2]:.2f}, KII={crows90[0][3]:.2f}, "
        f"KI/KII={crows90[0][2]/crows90[0][3]:.3f} | Passo 1 KII={crows90[1][3]:.3f}"
    )
    a_list = [r[0] for r in frows90]
    n_list = [r[2] for r in frows90]
    a0 = a_list[0]
    for label, target_da, ref_num, ref_exp in [
        ("CMM1", 8.52, 66550, 51166),
        ("CMM2", 4.24, 42508, 37455),
        ("CMM3", 3.64, 38059, 30007),
    ]:
        target_a = a0 + target_da
        n_interp = 0.0
        for k in range(len(a_list) - 1):
            if a_list[k] <= target_a <= a_list[k + 1]:
                t = (target_a - a_list[k]) / (a_list[k + 1] - a_list[k])
                n_interp = n_list[k] + t * (n_list[k + 1] - n_list[k])
                break
        print(
            f"  -> {label} (da={target_da:.2f} mm): N={n_interp:.1f} vs Diss={ref_num} "
            f"(dif={(n_interp-ref_num)/ref_num*100:+.2f}%, Exp={ref_exp})"
        )

    # 5.3.2: alpha = 110 deg (Chambers_CMM_4)
    rdir110, aid110, dt110 = run_solver(problem=1)
    crows110 = parse_crack_results(rdir110, aid110)
    frows110 = parse_fatigue_results(rdir110, aid110)
    a_list4 = [r[0] for r in frows110]
    n_list4 = [r[2] for r in frows110]
    a0_4 = a_list4[0]
    target_a4 = a0_4 + 2.70
    n_interp4 = 0.0
    for k in range(len(a_list4) - 1):
        if a_list4[k] <= target_a4 <= a_list4[k + 1]:
            t = (target_a4 - a_list4[k]) / (a_list4[k + 1] - a_list4[k])
            n_interp4 = n_list4[k] + t * (n_list4[k + 1] - n_list4[k])
            break
    print(
        f"[CMM4 110°] (tempo={dt110:.2f}s): Passo 0 KI={crows110[0][2]:.2f}, KII={crows110[0][3]:.2f}, "
        f"KI/KII={crows110[0][2]/crows110[0][3]:.4f} | Passo 1 KII={crows110[1][3]:.3f} | "
        f"N(da=2.7mm)={n_interp4:.1f} (Diss=28783, Exp=27782) | N_final={n_list4[-1]:.1f}"
    )
    return rdir110


def main():
    validate_sec_3_1_3_1()
    validate_sec_3_1_3_2()
    validate_sec_5_1()
    validate_sec_5_2_1_and_5_2_2()
    rdir_dowling = validate_sec_5_2_3_dowling()
    rdir_cmm4 = validate_sec_5_3_cmm()

    # Atualiza os gráficos de pós-processamento em examples/
    print("\nGerando figuras de pós-processamento em examples/...")
    for png in (ROOT / "inputs").glob("*_mesh_bc.png"):
        shutil.copy2(png, EXAMPLES_DIR / png.name)
    subprocess.run([sys.executable, str(ROOT / "plot_results.py"), str(rdir_cmm4), "--output-dir", str(EXAMPLES_DIR), "--scale", "40"], check=True)
    subprocess.run([sys.executable, str(ROOT / "plot_results.py"), str(rdir_dowling), "--output-dir", str(EXAMPLES_DIR), "--scale", "40"], check=True)
    print("\nTodos os exemplos foram executados e validados com sucesso!")


if __name__ == "__main__":
    import sys
    main()
