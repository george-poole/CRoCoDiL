from lucifex.sim import integrate_from_cli, postprocess_structured

from co2_dissolution_pkg.simulation import solutal_convective_dissolution_2d


def main():
    simulation = integrate_from_cli(solutal_convective_dissolution_2d)
    postprocess_structured(
        simulation, 
        exclude=['u'],
    )


if __name__ == "__main__":
    main()

