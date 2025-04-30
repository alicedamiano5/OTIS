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





def maxwellian_distribution(x, r, Halo_Type, bulge=None, Bulge_Type=None):
    """Calculate the value of the Maxwellian distribution given a dispersion."""
    
    # Compute the halo velocity dispersion
    sigma_halo = np.vectorize(Halo_Type.sigmav_nfw)(r)

    if bulge and Bulge_Type:
        # Compute the bulge velocity dispersion
        sigma_bulge = np.vectorize(Bulge_Type.sigma_bulge)(r)

        # Correctly compute the total velocity dispersion
        sigma_total = np.sqrt((sigma_halo**2 + sigma_bulge**2))
    else:
        sigma_total = sigma_halo
    ############ this is making the assumption that 
    
    # Maxwellian Distribution Formula
    return (np.sqrt(1 / (2 * np.pi * (sigma_total ** 2)) ** 3) *
            x ** 2 * np.exp(-(x ** 2) / (2 * sigma_total ** 2)))

def interpolation(radius, vel, Halo_Type, velocity_distribution_function=maxwellian_distribution, bulge=False, Bulge_Type=None):

    # Define the integrand for integration
    def integrand(v, r):
        return velocity_distribution_function(v, r, Halo_Type, bulge, Bulge_Type)

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

def analytical_times(BH,interp, Halo_Type=NFWHalo, bmax = 0., bmin = 0.,  interpolating_technique='Radau', t_span = (0,14), Nstep = 1000, r_limit= 0.01, v_limit= 0.01, DF = False, bulge = False, Bulge_Type=None):
    


    def v_circ(r):
        res = Grav*Halo_Type.mass_nfw(r)/r               ############ The halo type shoud have the same functions 
        return np.sqrt(res)
    
    if bulge == True: 
        def sigma_tot(r): 
            sigma = np.sqrt((Halo_Type.sigmav_nfw(r)**2)+(Bulge_Type.sigma_bulge(r)**2))
            return sigma
        
        def acc_tot(r):
            acc = Bulge_Type.bulge_acc(r)+Halo_Type.nfw_acc(r)
            return acc
        def rho_tot(r):
            rho = Halo_Type.rho_nfw(r) +  Bulge_Type.rho_bulge(r)
            return rho
    else:
        def sigma_tot(r):
            sigma = Halo_Type.sigmav_nfw(r)
            return sigma
        def acc_tot(r):
            acc = Halo_Type.nfw_acc(r)
            return acc
        def rho_tot(r):
            rho = Halo_Type.rho_nfw(r) 
            return rho
        
        

    
    
    def acceleration_friction_noncircular(r,v, bmin, bmax):

        ####### Make this in input..
        
        bmin = Grav*BH.mass/((v**2)+((2/3)*(sigma_tot(r)**2)))
        bmax = Halo_Type.r_s
        res = -16 * np.pi*np.pi*(np.log(bmax/bmin))*Grav*Grav*(BH.mass)*Halo_Type.rho_nfw(r)*interp((v,r))/(v*v)
        return res
    

    

    # Equations of motion for the system
    def equations(t, cond):
        [x0,y0, vx0, vy0 ]=cond
        r0 = np.longdouble(np.sqrt((x0**2)+(y0**2)))
        v = np.longdouble(np.sqrt((vx0**2)+(vy0**2)))
        print(r0, v)
        if DF == False: 
             dvx = Halo_Type.nfw_acc(r0)*(x0/r0)
             dvy = Halo_Type.nfw_acc(r0)*(y0/r0)
        else:
            dvx = acc_tot(r0)*(x0/r0)+acceleration_friction_noncircular(r0, v, bmin,bmax )*(vx0/v)
            dvy = acc_tot(r0)*(y0/r0)+acceleration_friction_noncircular(r0, v, bmin, bmax)*(vy0/v)
      
        return [np.longdouble(vx0), np.longdouble(vy0), dvx, dvy]

    
    #######Initial conditions

    def stop_condition(t, y):
        """Stops integration when r < r_limit or v < v_limit"""
        # Extract state variables
        r = y[0]**2 + y[1]**2
        v = y[2]**2 + y[3]**2

        terminal_condition = 0

        if r<r_limit or v < v_limit:
            terminal_condition = 1
    
        return terminal_condition  # Stop when either reaches its limit

    
    def stop_condition(t, y):
        """Stops integration when r < r_limit or v < v_limit"""
        # Extract state variables
        r = np.sqrt(y[0]**2 + y[1]**2)
        v = np.sqrt(y[2]**2 + y[3]**2)
    
        return min(r - r_limit, v - v_limit)
    # Tell `solve_ivp` to stop when the event is triggered
    stop_condition.terminal = True  
    stop_condition.direction = -1  # Detect when values decrease past limits




      
    t_eval = np.linspace(t_span[0], t_span[1], Nstep)  # Time points where solution is evaluated
    
    #sol = solve_ivp(equations, t_span, [BH.initial_position[0], BH.initial_position[1], BH.initial_velocity[0], BH.initial_velocity[1]], t_eval=t_eval, method=interpolating_technique)
    #rk45 = RK45(fun, t0, y0, t_bound=1000, max_step=0.1)



    # Tell `solve_ivp` to stop when the event is triggered
    #stop_condition.terminal = True  
    #stop_condition.direction = -1  # Detect when values decrease past limits
    
    # Solve the ODEs
    sol = solve_ivp(
        equations, 
        t_span, 
        [BH.initial_position[0], BH.initial_position[1], BH.initial_velocity[0], BH.initial_velocity[1]], 
        t_eval=t_eval, 
        method=interpolating_technique, 
        events=stop_condition  # Pass the stopping function
    )

    
        
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