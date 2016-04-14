from dolfin import *

mesh = UnitIntervalMesh(50)
V = FunctionSpace(mesh, "CG", 12)

alpha = 1
gamma = 1.2
dt = 0.005
k = Constant(dt)
t = dt
T = 0.5
g = Expression("1",
             gamma = gamma, t = 0)

f = Expression("1")


u0 = project(f, V)


bc = DirichletBC(V, 0.0, "x[0] < DOLFIN_EPS && on_boundary")

u = TrialFunction(V)
v = TestFunction(V)
#F = ((u-u0)/k*v + inner(u.dx(0), v))*dx
a = u*v*dx + k/2*inner(u.dx(0), v)*dx
L = u0*v*dx - k/2*inner(u0.dx(0), v)*dx
u = Function(V)


while (t <= T):
    solve(a == L, u, bc)

    u0.assign(u)
    t += dt
    p = plot(u0)


interactive()
