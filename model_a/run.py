from lucifex.sim import Simulation, run_from_cli, xdmf_to_npz
from model.dns import dns_model_a


def posthook(
    sim: Simulation,
    delete_xdmf: bool | None = None,
) -> None:
    if delete_xdmf is None:
        return    
    xdmf_to_npz(sim, delete_xdmf=delete_xdmf)


def main():
    run_from_cli(dns_model_a, posthook)


if __name__ == "__main__":
    main()

