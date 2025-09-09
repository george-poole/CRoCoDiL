from lucifex.sim import integrate_from_cli, write_structured

from co2_dissolution_pkg.simulation import solutal_convective_dissolution_2d


def main():
    simulation = integrate_from_cli(solutal_convective_dissolution_2d)
    write_structured(
        simulation, 
        delete_h5_xdmf=True,
    )


if __name__ == "__main__":
    main()

