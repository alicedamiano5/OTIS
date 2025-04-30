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

from otis.otis_utils import *

##################### since initially suited to work on Gadget code units, this is a conversion function 

########Initialize function, all that you should initialize at the beginning
Grav = 43079.93181818183



def convert_in_gyr(time_code):
    conversion_fac = (u_l/u_vel)
    time_in_s = time_code *conversion_fac
    time_in_gyr = time_in_s/(3600*24*365*1e9)
    return time_in_gyr


######### NFW halo parameters
class NFWHalo:

    def __init__(self, r_vir, conc, m_vir):
        self.r_vir = r_vir
        self.conc = conc
        self.m_vir = m_vir
        self.r_s = r_vir / conc  # Scale radius
        

    def rho_s(self):
        """Calculate the scale density."""
        numerator = self.m_vir
        denominator = 4 * np.pi * (self.r_s ** 3) * (
            np.log(1 + self.conc) - (self.conc / (1 + self.conc)))
        return numerator / denominator

    def mass_nfw(self, r):
        """Calculate the mass enclosed within radius r."""
        x = r / self.r_s
        res = (
            4
            * np.pi
            * self.rho_s()
            * (self.r_s ** 3)
            * (np.log(1 + x) - (x / (1 + x)))
        )
        return res

    def sigmav_nfw(self, r):
        """Calculate the velocity dispersion at radius r."""
        # Define the function to integrate
        g = lambda x: (
            (np.log(1 + x / self.r_s) / ((1 + x / self.r_s) ** 2 * (x / self.r_s) ** 3))
            - (1 / (((1 + x / self.r_s) ** 3) * x ** 2 / (self.r_s ** 2))))

        # Prefactor for the calculation
        prefac = 4 * np.pi * Grav * self.rho_s() * self.r_s

        # Perform the integration
        lower_limit = r
        result, error = quad(g, lower_limit, np.inf)

        # Adjust the result
        result *= (r / self.r_s) * (1 + (r / self.r_s)) ** 2

        # Final calculation
        sigmas = result * prefac

        return np.sqrt(sigmas)

    
    def rho_nfw(self, r):
        """Calculate the NFW density at radius r."""
        x = r / self.r_s
        return self.rho_s() / (x * (1 + x) ** 2)

    def nfw_potential(self, r):
        return -Grav * 4*np.pi*(self.r_s**3)*self.rho_s() * np.log(1 + r / self.r_s) / (r)
    
    # Derivative of the potential (force)
    def nfw_acc(self,r):
        x=r/self.r_s
        first = self.m_vir/(np.log(1+self.conc)-(self.conc/(1+self.conc)))
        second=(1/(r))*(-((np.log(1+x))/r)+(1/((r+self.r_s))))
        return np.longdouble(first*second*Grav)
    
    def v_circ(self,r):
        res = Grav*self.mass_nfw(r)/r
        return np.sqrt(res)


class BH:
    def __init__(self, mass, initial_position, initial_velocity):
        self.mass     = mass
        self.initial_position = initial_position
        self.initial_velocity = initial_velocity


class Bulge:
    def __init__(self, r_bulge, m_bulge):
        self.r_bulge = r_bulge
        self.m_bulge = m_bulge
    
    def rho_bulge(self,r):

        dens = (self.m_bulge/(2*np.pi*r*((self.r_bulge+r)**3)))*self.r_bulge
        return dens

    def mass_bulge(self,r):

        mass = self.m_bulge*r*r/((r+self.r_bulge)**2)
        return mass
       
    def sigma_bulge(self,r):

        prefac = -Grav*self.m_bulge*r*((self.r_bulge+r)**3)
        one = (np.log(abs(r))-np.log(abs(r+self.r_bulge)))/(self.r_bulge**5)
        two = 1/((self.r_bulge**4)*(r+self.r_bulge))
        three = 1/(2*(self.r_bulge**3)*((r+self.r_bulge)**2))
        four=1/(3*self.r_bulge*self.r_bulge*((r+self.r_bulge)**3))
        five=1/(4*self.r_bulge*((r+self.r_bulge)**4))
        sigma = (one + two + three+four+five)*prefac
        return np.sqrt(sigma)
    
    def v_circ_bulge(self,r):
        res = Grav*self.mass_bulge(r)/r
        return np.sqrt(res)
    
    def bulge_potential(self, r): 
        res = -Grav*self.m_bulge/(r+self.r_bulge)
        return res


    def bulge_acc(self, r):
        res = -Grav*self.m_bulge/((r+self.r_bulge)*(r+self.r_bulge))
        return res

