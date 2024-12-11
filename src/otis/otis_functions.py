#============================================================#
#         Welcome to the Orbiting TImescale for Sinking library
#============================================================
#         O  T  I  S
#         /\     /\     Created by Alice Damiano
#        /  \___/  \     contact: alice.damiano@inaf.it
#       (    o o    )
#       /     ^     \
#===============================================================
        
#===============================================================


from otis.otis_classes import *
from otis.otis_utils import *





def maxwellian_distribution(x,r, Halo_Type=NFWHalo):
    """Calculate the value of the maxwellian distribution given a dispersion"""
    sigma_function_vectorized = np.vectorize(Halo_Type.sigmav_nfw)
    sigma = sigma_function_vectorized(r)

    return ((np.sqrt(1/(2*np.pi*(sigma**2))**3)*x*x* np.exp(-(x**2) / (2 * sigma**2))))

def interpolation(radius, vel, Halo_Type, velocity_distribution_function=maxwellian_distribution):

    # Define the integrand for integration
    def integrand(v, r):
        return velocity_distribution_function(v, r, Halo_Type)

    # Perform the integration for a given distance r
    def compute_integral(v, r):
        integral, error = quad(integrand, 0.001, v, args=(r,))
        return integral

    # Vectorize the compute_integral function for grid processing
    compute_integral_vect = np.vectorize(compute_integral)

    # Generate the grid of v and r
    V, R = np.meshgrid(vel, radius, indexing='ij')  # Ensure consistent indexing
    print("Computing integral on the grid...")

    # Compute Z over the grid
    Z = compute_integral_vect(V, R)
    print("Integral computed.")

    # Create the interpolator
    print("Interpolating...")
    interp = RegularGridInterpolator((vel, radius), Z)
    print("Interpolation done.")
    return interp


################## heart of the code: produce the analytical times

############ you can provide an eccentricity and the code translate it into a velocity value!!!!!!!! - IMPROVEMENT

def analytical_times(BH,interp, Halo_Type=NFWHalo, bmax = 0., bmin = 0.,  interpolating_technique='Radau', t_span = (0,14), Nstep = 1000):
    


    def v_circ(r):
        res = Grav*Halo_Type.mass_nfw(r)/r               ############ The halo type shoud have the same functions 
        return np.sqrt(res)

    
    
    def acceleration_friction_noncircular(r,v, bmin, bmax):

        ####### Make this in input..
        
        bmin = Grav*BH.mass/((v**2)+((2/3)*(Halo_Type.sigmav_nfw(r)**2)))
        bmax = Halo_Type.r_s
        res = -16 * np.pi*np.pi*(np.log(bmax/bmin))*Grav*Grav*(BH.mass)*Halo_Type.rho_nfw(r)*interp((v,r))/(v*v)
        return res
    

    

    # Equations of motion for the system
    def equations(t, cond):
        [x0,y0, vx0, vy0 ]=cond
        r0 = np.longdouble(np.sqrt((x0**2)+(y0**2)))
        v = np.longdouble(np.sqrt((vx0**2)+(vy0**2)))
        dvx = Halo_Type.nfw_acc(r0)*(x0/r0)+acceleration_friction_noncircular(r0, v, bmin,bmax )*(vx0/v)
        dvy = Halo_Type.nfw_acc(r0)*(y0/r0)+acceleration_friction_noncircular(r0, v, bmin, bmax)*(vy0/v)
      
        return [np.longdouble(vx0), np.longdouble(vy0), dvx, dvy]

    
    #######Initial conditions
    


      
    t_eval = np.linspace(t_span[0], t_span[1], Nstep)  # Time points where solution is evaluated
    sol = solve_ivp(equations, t_span, [BH.initial_position[0], BH.initial_position[1], BH.initial_velocity[0], BH.initial_velocity[1]], t_eval=t_eval, method=interpolating_technique)
    
    
    # Solve the equations of motion
    #sol = solve_ivp(equations, t_span, [20, 0, 0, vtheta],t_eval=t_eval, method='RK45', min_step = 0.1)
    
    # Extract the solution
    x = sol.y[0]
    y = sol.y[1]
    vx = sol.y[2]
    vy= sol.y[3]
    v = np.sqrt(vx**2+vy**2)
    #ax = (nfw_acc(r0)*(x0/r0))+acceleration_friction_noncircular(r0, b_max, b_min, v)*(vx0/(v))
    #ay =(nfw_acc(r0)*(y0/r0))+acceleration_friction_noncircular(r0, b_max, b_min, v)*(vy0/(v))
    
    print(sol.message)
    return(sol.t, x, y, vx, vy)
