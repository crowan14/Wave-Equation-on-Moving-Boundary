# Wave-Equation-on-Moving-Boundary

The simplest form of the wave equation is

$$ \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2}, \quad u(a,t)=u(b,t)=0, \quad u(x,0) = f(x), \quad \frac{\partial u}{\partial t}(x,0) = 0 $$

This could model lateral vibrations in a tensioned string, or the axial vibrations of an elastic bar. Constructing an approximate solution to this equation is a canonical problem in numerical PDE's. First, we discretize the solution in space with a set of basis functions modulated by time-dependent coefficients:

$$ u(x,t) = \sum_i a_i(t) w_i(x) $$

For Galerkin methods, the weak form of the PDE is then obtained by integrating against the same spatial functions which are used to discretize the solution. This reads


$$ \int_a^b \frac{\partial^2 u}{\partial t^2} w_j dx = \int_a^b \frac{\partial^2 u}{\partial x^2} w_j dx $$

This discretization of the displacement is substituted, integration by parts is used to transfer a spatial derivative from the solution to the test function, and we obtain the following system of ordinary differential equations:

$$ M_{ij} \ddot a_j + K_{ij} a_j = 0 $$


