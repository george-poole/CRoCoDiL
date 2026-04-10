import os
import signal
import datetime
from typing import Any

from lucifex.fem import Constant
from lucifex.fdm import ConstantSeries, FE
from lucifex.sim import Simulation, Stopper, run_from_cli, xdmf_to_npz, stopper
from crocodil.dns.system_a import dns_system_a


def prehook(_: Simulation) -> None:
    print(f"PID: {os.getpid()} started at {str(datetime.datetime.now())}", flush=True)


def posthook(
    sim: Simulation,
    delete_xdmf: bool | None = None,
    kill_pid: bool = False,
) -> None:
    print(f"PID: {os.getpid()} finished at {str(datetime.datetime.now())}", flush=True)
    if delete_xdmf is not None: 
        xdmf_to_npz(sim, delete_xdmf=delete_xdmf)
        print(f"PID: {os.getpid()} postprocessed at {str(datetime.datetime.now())}", flush=True)
    if kill_pid:
        print(f"PID: {os.getpid()} to be killed", flush=True)
        os.kill(os.getpid(), signal.SIGTERM)


def is_greater_than(
    u: float | Constant | ConstantSeries,
    u_thresh: float,
) -> bool:
    return float(FE(u)) > u_thresh


def c_gt_stopper(
    sim: Simulation,
    *args,
) -> Stopper:
    cMin: Constant = sim['cMin']

    # cMin: float

    cMin
    
    return stopper(cMin, is_greater_than)()


c_gt_stopper_ = stopper('u', is_greater_than)


def main():
    run_from_cli(dns_system_a, prehook=prehook, posthook=posthook)


if __name__ == "__main__":
    main()

