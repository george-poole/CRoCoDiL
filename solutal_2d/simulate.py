from lucifex.sim import integrate_from_cli, postprocess_grids

from co2_dissolution_pkg.sim import solutal_rectangle


def main():
    simulation = integrate_from_cli(solutal_rectangle)
    postprocess_grids(
        simulation, 
        delete_h5_xdmf=True,
    )


if __name__ == "__main__":
    main()

