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

We will assume that $a(t)\geq0$ which means that $L$ is the total length of the string and $a(t)$ determines the position of the clamp enforcing the zero displacement boundary condition. For simplicitly, we imagine expanding the domain from some $a_0=a(0)$ to a final position $a_f=a(T)<a_0$ where $T$ is the total elapsed time. Decreasing the coordinate specifying the position of the clamp corresponds to increasing the size of the domain. 

Our way of simulating this problem is a (potentially sketchy) appeal to intuition rather than the result of surveying the literature, which proved to be quite mathematical! We will simply treat the part of the string hidden by the clamp as constrained to zero displacement until it is exposed as a result of the clamp's motion. This constraint will be built into the displacement discretization. As such, the displacement will be discretized with

$$ u(x,t) = \sum_i a_i(t) w_i(x) H(x-a(t)) $$

where $H(x-a(t))$ is a Heaviside step function is zero for $x<a(t)$ and 1 for $x\geq a(t) $. This discretization ensures that only the part of the string exposed by the moving clamp has a non-zero contribution from the shape functions. Note that this means the shape functions discretizing the displacement change in time, which is a departure from usual discretization techniques. No problem! We can naively soldier on and see what happens. We have to be a bit more careful in constructing the weak form because the space and time components of the solution are no longer decoupled. For ease of notation, we will introduce the following convention:

$$ W_i(x,t) = w_i(x) H(x-a(t)) $$

The weak form of the wave equation with this discretization can be written as 


$$ \int_0^L \frac{\partial^2 u}{\partial t^2} W_j dx + \int_0^L \frac{\partial u}{\partial x}\frac{\partial W_j}{\partial x} dx = 0 $$ 

