# Wave-Equation-on-Moving-Boundary

The simplest form of the wave equation is

$$ \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2}, \quad u(a,t)=u(b,t)=0, \quad u(x,0) = f(x), \quad \frac{\partial u}{\partial t}(x,0) = 0 $$

This could model lateral vibrations in a tensioned string, or the axial vibrations of an elastic bar. Constructing an approximate solution to this equation is a canonical problem in numerical PDE's. First, we discretize the solution in space with a set of basis functions modulated by time-dependent coefficients:

$$ u(x,t) = \sum_i a_i(t) w_i(x) $$

For Galerkin methods, the weak form of the PDE is then obtained by integrating against the same spatial functions which are used to discretize the solution. This reads


$$ \int_a^b \frac{\partial^2 u}{\partial t^2} w_j dx = \int_a^b \frac{\partial^2 u}{\partial x^2} w_j dx $$

This discretization of the displacement is substituted, integration by parts is used to transfer a spatial derivative from the solution to the test function, and we obtain the following system of ordinary differential equations:

$$ M_{ij} \ddot a_j + K_{ij} a_j = 0 $$

Once one has gone through this process many times, it is very formulaic. So far so good. We can now ask a seemingly innocent question: what if the domain on which the PDE is solved changes with time? For example, imagine the vibrations of a string where a clamp (which zeros the displacement) is slid along the string's length. This effectively changes the length of the domain while the string is vibrating. This is not a traditional problem in PDE's, and is a bit tough to think about. This is a first pass at making sense of this problem. It is interesting that a seemingly valid physical scenario seems to so clearly fall outside the purview of traditional PDEs. 

The moving domain problem can be written as

$$ \frac{\partial^2 u}{\partial t^2} = \frac{\partial^2 u}{\partial x^2}, \quad u(a(t),t)=u(L,t)=0, \quad u(x,0) = f(x), \quad \frac{\partial u}{\partial t}(x,0) = 0 $$


