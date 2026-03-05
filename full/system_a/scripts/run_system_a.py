import os
import signal
from lucifex.sim import Simulation, run_from_cli, xdmf_to_npz
from crocodil.dns.system_a import dns_system_a


def posthook(
    sim: Simulation,
    delete_xdmf: bool | None = None,
) -> None:
    print(f"PID: {os.getpid()} finished")
    if delete_xdmf is not None: 
        xdmf_to_npz(sim, delete_xdmf=delete_xdmf)
    print(f"PID: {os.getpid()} postprocessed")
    os.kill(os.getpid(), signal.SIGTERM)


def main():
    run_from_cli(dns_system_a, posthook)


if __name__ == "__main__":
    main()

