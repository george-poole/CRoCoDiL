from memory_profiler import profile
from dolfinx.fem import FunctionSpace
from lucifex.mesh import rectangle_mesh
from lucifex.fem import Function


@profile
def test_dolfinx_memory(
    Nx: int, 
    Ny: int,
    n_alloc: int,
) -> None:
    mesh = rectangle_mesh(1.0, 1.0, Nx, Ny)
    fs = FunctionSpace(mesh, ('P', 1))

    f_list = []
    for _ in range(n_alloc):
        f = Function((mesh, 'P', 1))
        f_list.append(f)

if __name__ == '__main__':
    test_dolfinx_memory(160, 160, 500)
