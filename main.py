import solver
import short
import units

W = 52  # terminal width

def header():
    print("=" * W)
    print("POLYTROPIC TOV SOLVER".center(W))
    print("=" * W)

def section(title):
    print("\n" + f"  {title}")
    print("  " + "-" * (W - 2))

def row(label, value, unit=""):
    print(f"  {label:<20} {value:>20}  {unit}")

def divider():
    print("=" * W)


if __name__ == "__main__":

    header()

    # ── inputs ────────────────────────────────────────────────────────────
    section("EOS PARAMETERS")
    K     = float(input("  Polytropic Constant K  (GEOM) : "))
    Gamma = float(input("  Polytropic Index Gamma        : "))
    rho_c = float(input("  Central Density rho_c  (GEOM) : "))

    section("SOLVER SELECTION")
    print("  1. Direct  (integrate P directly)")
    print("  2. Log     (integrate log P — stable near surface)")
    print("  3. Enthalpy (integrate enthalpy — stable near surface)")
    SOLVER_TYPE = int(input("\n  Choice: "))
    flag = SOLVER_TYPE
    method_map = {1: 'direct', 2: 'log', 3: 'enthalpy'}
    methods = method_map[SOLVER_TYPE]


    # ── solve ─────────────────────────────────────────────────────────────
    print(f"\n  Running {methods}...\n")

    R, M, sol = solver.solver(rho_c, K, Gamma, method = methods)


    # ── results ───────────────────────────────────────────────────────────
    section("RESULTS — GEOMETRIZED UNITS  (G = c = 1)")
    row("Radius  R",  f"{R:.9f}",  "[M_sun]")
    row("Mass    M",  f"{M:.9f}",  "[M_sun]")

    section("RESULTS — CGS UNITS")
    row("Radius  R",  f"{R * units.LENGTH_CGS / 1e5:.9f}",  "[km]")
    row("Mass    M",  f"{M * units.M_sun:.9e}",             "[g]")
    row("Mass    M",  f"{M:.9f}",                           "[M_sun]")

    # ── short file ────────────────────────────────────────────────────────
    section("CONSTRUCTING GR1D SHORT FILE")
    fname   = input("  Output filename        : ")
    short.create(sol, K, Gamma, fname, flag)

    divider()
    print("  File saved to:".ljust(22) + fname)
    divider()
    print()
