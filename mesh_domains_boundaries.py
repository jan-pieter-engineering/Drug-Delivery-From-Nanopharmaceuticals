# -----------------------------------------------------------------
#  M E S H ,  D O M A I N S  &  B O U N D A R I E S
# -----------------------------------------------------------------
from fenics import *

def define_mesh(nx, r_i, r_o, element_degree, element_family):
    """ Mesh is defined by number of knots, left & right boundary coordinate, element family & degree. """

    # Create mesh and define function space
    mesh = IntervalMesh(nx, r_i, r_o)
    space_dim = mesh.geometry().dim()
    V = FunctionSpace(mesh, element_family, element_degree)


    return mesh, space_dim, V



def define_boundaries(mesh, space_dim):
    """ Boundaries are defined. """

    # Initialize mesh function for boundary domains
    boundaries = MeshFunction("size_t", mesh, space_dim - 1)
    boundaries.set_all(0)

    # Define new measures associated with the exterior boundaries
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)


    return ds



def define_r(element_degree):
    """ r defines 1d-cyclinder coordinates. """

    return Expression('x[0]', degree=element_degree)



def collect_mesh_parameters(r_i, r_o, nx, element_degree, element_family, boundary_tol):

    _mesh, _space_dim, _V = define_mesh(nx, r_i, r_o, element_degree, element_family)
    _ds = define_boundaries(_mesh, _space_dim)
    _r = define_r(element_degree)
    _n_facet = FacetNormal(_mesh)  # define normal vector


    return _mesh, _V, _ds, _r, _n_facet