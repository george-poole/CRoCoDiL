from lucifex.sim import Simulation, run_from_cli, np_postprocess
from co2_pkg.sim import solutal_rectangle


def posthook(
    sim: Simulation,
    delete_h5_xdmf: bool | None = None,
) -> None:
    if delete_h5_xdmf is None:
        return    
    np_postprocess(sim, delete_h5_xdmf=delete_h5_xdmf)


def main():
    run_from_cli(solutal_rectangle, posthook)


if __name__ == "__main__":
    main()

