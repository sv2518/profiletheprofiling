from firedrake import *
from firedrake.petsc import PETSc
import time
import sys


commit = sys.argv[1]
path_to_dir = sys.argv[2]
run = sys.argv[3]
problem = sys.argv[4]

  
def problem1():
    """ Solve equations for lowest order RT discretisation of Poisson problem with Hybridization and CG+jacobi on the Trace system"""
    mesh = UnitSquareMesh(20, 20)
    U = FunctionSpace(mesh, "RT", 1)
    V = FunctionSpace(mesh, "DG", 0)
    W = U * V
    sigma, u = TrialFunctions(W)
    tau, v = TestFunctions(W)
    n = FacetNormal(mesh)

    # Define the source function
    x, y = SpatialCoordinate(mesh)
    f = Function(V)
    f.interpolate(10*exp(-(pow(x - 0.5, 2) + pow(y - 0.5, 2)) / 0.02))

    # Define the variational forms
    a = (inner(sigma, tau) + inner(u, div(tau)) + inner(div(sigma), v)) * dx
    L = -inner(f, v) * dx + Constant(0.0) * dot(conj(tau), n) * (ds(3) + ds(4))

    # Compare hybridized solution with non-hybridized
    w = Function(W)
    bc1 = DirichletBC(W[0], as_vector([0.0, -sin(5*x)]), 1)
    bc2 = DirichletBC(W[0], as_vector([0.0, sin(5*y)]), 2)
    bcs = [bc1, bc2]

    params = {'mat_type': 'matfree',
                        'ksp_type': 'preonly',
                        'pc_type': 'python',
                        'pc_python_type': 'firedrake.HybridizationPC',
                        'hybridization': {'ksp_type': 'cg',
                                        'pc_type': 'jacobi',
                                        'ksp_rtol': 1e-8}}

    solve(a == L, w, bcs=bcs, solver_parameters=params)


def problem2():
    """ Solve equations for lowest order RT discretisation of Poisson problem with Hybridization and CG+GTMG on the Trace system.
        This will produce more local kernels than problem1. """
    nlevels = 2
    n=2
    s=1
    base = SquareMesh(n, n, s, quadrilateral=True)
    basemh = MeshHierarchy(base, nlevels)
    mh = ExtrudedMeshHierarchy(basemh, s, base_layer=n)
    mesh = mh[-1]
    x = SpatialCoordinate(mesh)

    def get_p1_space():
        return FunctionSpace(mesh, "CG", 1)

    def get_p1_prb_bcs():
        return [DirichletBC(get_p1_space(), zero(), "on_boundary"),
                DirichletBC(get_p1_space(), zero(), "top"),
                DirichletBC(get_p1_space(), zero(), "bottom")]

    def p1_callback():
        P1 = get_p1_space()
        p = TrialFunction(P1)
        q = TestFunction(P1)
        return inner(grad(p), grad(q))*dx

    p = 0
    RT = FiniteElement("RTCF", quadrilateral, p+1)
    DG_v = FiniteElement("DG", interval, p)
    DG_h = FiniteElement("DQ", quadrilateral, p)
    CG = FiniteElement("CG", interval, p+1)
    HDiv_ele = EnrichedElement(HDiv(TensorProductElement(RT, DG_v)),
                            HDiv(TensorProductElement(DG_h, CG)))
    U = FunctionSpace(mesh, HDiv_ele)
    V = FunctionSpace(mesh, "DQ", p)
    W = U * V
    sigma, u = TrialFunctions(W)
    tau, v = TestFunctions(W)

    L=1
    exact = 100*x[0]*(L-x[0])*x[1]*(L-x[1])*x[2]*(L-x[2])
    f = -div(grad(exact))

    a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
    L = inner(f, v)*dx

    w = Function(W)
    
    mgmatfree_mtf = {'snes_type': 'ksponly',
                    'ksp_type': 'preonly',
                    'mat_type': 'matfree',
                    'pc_type': 'mg',
                    'pc_mg_type': 'full',
                    'ksp_rtol': 1e-8,
                    'mg_coarse_ksp_type': 'preonly',
                    'mg_coarse_pc_type': 'python',
                    'mg_coarse_pc_python_type': 'firedrake.AssembledPC',
                    'mg_coarse_assembled_pc_type': 'lu',
                    'mg_coarse_assembled_pc_factor_mat_solver_type': 'superlu_dist',
                    'mg_levels_ksp_type': 'chebyshev',
                    'mg_levels_ksp_norm_type': 'unpreconditioned',
                    'mg_levels_ksp_max_it': 3,
                    'mg_levels_pc_type': 'none',
                    'ksp_norm_type': 'preconditioned'
                    }

    params = {'mat_type': 'matfree',
                        'ksp_type': 'preonly',
                        'pc_type': 'python',
                        'pc_python_type': 'firedrake.HybridizationPC',
                        'hybridization': {'ksp_type': 'cg',
                                            'pc_type': 'python',
                                            'ksp_rtol': 1e-5,
                                            'ksp_atol': 1e-8,
                                            'mat_type': 'matfree',
                                            'localsolve': {'ksp_type': 'preonly',
                                                            'pc_type': 'fieldsplit',
                                                            'pc_fieldsplit_type': 'schur'},
                                            'ksp_norm_type': 'preconditioned',
                                            'pc_python_type': 'firedrake.GTMGPC',
                                            'gt': {'mg_levels': {'ksp_type': 'chebyshev',
                                                                 'pc_type': 'none',
                                                                 'ksp_max_it': 3},
                                                    'mg_coarse': mgmatfree_mtf,
                                                    'mat_type':'matfree'
                                                    }
                                            }
                        }
    appctx = {'get_coarse_operator': p1_callback,
                'get_coarse_space': get_p1_space,
                'coarse_space_bcs': get_p1_prb_bcs()}


    if ("hybridization" in params.keys() 
        and "localsolve" in params["hybridization"].keys() 
        and "mat_type" in params["hybridization"]["localsolve"].keys()):
        if params["hybridization"]["localsolve"]["mat_type"] == "matfree":
            fcp = {"slate_compiler": {"replace_mul": True}}
            appctx["form_compiler_parameters"] = fcp
    else:
        fcp = {}

    solve(a == L, w, solver_parameters=params, appctx=appctx, form_compiler_parameters=fcp)

def timing(path_to_dir, problem):
    if problem == "problem1":
        problem = problem1
    elif problem == "problem2":
        problem = problem2
    start_time = time.time()
    problem()
    timed = (time.time() - start_time)
    descr = "a+"
    with open(path_to_dir, descr) as myfile:
        myfile.write(str(timed)+",")
    

log = PETSc.Log.isActive()

timing(f"{path_to_dir}time_log{log}_{run}_{commit}.txt", problem)
