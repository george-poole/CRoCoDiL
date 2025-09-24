from dolfinx.mesh import Mesh

from lucifex.mesh import rectangle_mesh, mesh_boundary, MeshBoundary


def create_rectangle_domain(
    Lx: float,
    Ly: float,
    Nx: int,
    Ny: int,
    cell: str,
    name: str = 'LxLy',
    clockwise_names: tuple[str, str, str, str] = ('upper', 'right', 'lower', 'left'),
) -> tuple[Mesh, MeshBoundary]:
    
    mesh = rectangle_mesh(Lx, Ly, Nx, Ny, cell, name=name)
    boundary = mesh_boundary(
        mesh,
        {
            clockwise_names[0]: lambda x: x[1] - Ly,
            clockwise_names[1]: lambda x: x[0] - Lx,
            clockwise_names[2]: lambda x: x[1],
            clockwise_names[3]: lambda x: x[0],
        },
    )

    return mesh, boundary