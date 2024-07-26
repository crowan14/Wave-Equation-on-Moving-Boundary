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

We will assume that $a(t)\geq0$ which means that $L$ is the total length of the string and $a(t)$ determines the position of the clamp enforcing the zero displacement boundary condition. For simplicitly, we imagine expanding the domain from some $a_0=a(0)$ to a final position $a_f=a(T)\leq a_0$ where $T$ is the total elapsed time. Decreasing the coordinate specifying the position of the clamp corresponds to increasing the size of the domain. 

Our way of simulating this problem is a (potentially sketchy) appeal to intuition rather than the result of surveying the literature, which proved to be quite mathematical! We will simply treat the part of the string hidden by the clamp as constrained to zero displacement until it is exposed as a result of the clamp's motion. This constraint will be built into the displacement discretization. As such, the displacement will be discretized with

$$ u(x,t) = \sum_i a_i(t) w_i(x) H(x-a(t)) $$

where $H(x-a(t))$ is a Heaviside step function is zero for $x\leq a(t)$ and 1 for $x\geq a(t) $. This discretization ensures that only the part of the string exposed by the moving clamp has a non-zero contribution from the shape functions. Note that this means the shape functions discretizing the displacement change in time, which is a departure from usual discretization techniques. No problem! We can naively soldier on and see what happens. We have to be a bit more careful in constructing the weak form because the space and time components of the solution are no longer decoupled. For ease of notation, we will introduce the following convention:

$$ W_i(x,t) = w_i(x) H(x-a(t)) $$

The weak form of the wave equation with this discretization can be written as 

$$ \int_0^L \frac{\partial^2 u}{\partial t^2} W_j dx + \int_0^L \frac{\partial u}{\partial x}\frac{\partial W_j}{\partial x} dx = 0 $$ 

Note that we have to be more careful taking spatial derivatives of the shape functions $W_j(x,t)$ because their spatial dependence comes both from the underlying basis set $w_i(x)$ and the Heaviside step function. This also suggests that we need to regularize the step functions otherwise there would be infinite derivatives at the position of the step. The step functions can be smoothly approximated using hyperbolic tangents. Plugging in the displacement discretization, we have

$$ \sum_i \ddot a_i \Big( \int_0^L W_i(x,t) W_j(x,t) dx\Big) + \sum_i a_i \Big( \int_0^L\frac{\partial W_i}{\partial x} \frac{\partial W_j}{\partial x} dx \Big) = 0 $$

This shows that the mass and stiffness matrices are time-dependent quantities as a result of the time-varying step functions multiplying basis functions. The governing system of ODE's is then

$$ M_{ij}(t) \ddot a_j(t) + K_{ij}(t) a_j(t) = 0 $$

We can time integrate this system with a forward Euler scheme. This reads

$$ M_{ij}(t) \Big( \frac{a_j(t+1) - 2a_j(t) + a_j(t-1)}{\Delta t^2}\Big) + K_{ij} a_j(t) = 0 $$

The time updating scheme for the displacement degrees of freedom is then

$$ a(t+1) = -\Delta t^2 M^{-1}(t) K(t) a(t) + 2a(t) - a(t-1) $$

The stiffness and mass matrices are evaluated at the current time. This method is implemented in the attached code. Note that when reconstructing the solution the degrees of freedom need to be matched with their shape functions at a correspond time. This is because the basis changes with time. It is not clear whether this is "allowed." The results of this method indicate that either that it is difficult to have intuition for a problem of this sort, or that the method is extremely unstable. One thing that we can check is the total energy of the solution with time. Given that we are modeling free vibrations, there does not seem to be any energy input into the system even though the domain expands. Thus, the total energy of the solution should remain constant. The total energy at a given time is the sum of the strain energy and then kinetic energy, namely


$$ U(t) = \int_0^L \frac{1}{2}\Big( \frac{\partial u}{\partial t} \Big)^2 + \frac{1}{2} \Big( \frac{\partial u}{\partial x} \Big)^2 dx = T + V $$


The strain energy for the discrete solution can be written simply as

$$ V = \frac{1}{2} a_i K_{ij} a_j $$

