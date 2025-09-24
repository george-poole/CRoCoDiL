from lucifex.sim import integrate_from_cli, write_structured

from co2_dissolution_pkg.sim import solutal_2d


def main():
    simulation = integrate_from_cli(solutal_2d)
    write_structured(
        simulation, 
        delete_h5_xdmf=True,
    )


if __name__ == "__main__":
    main()

