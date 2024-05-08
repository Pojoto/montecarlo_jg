
import numpy as np

    
class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        """Constructor for BitString objects

        Parameters
        ----------
        N   : number
            Length of bitstring
        """
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __str__(self):
        """
        Print dunder method to display bitstring
        """
        string = ""
        for zero in self.config:
            string += str(zero)
        return string

    def __eq__(self, other):        
        """
        Check if two BitString objects are equal

        Parameters
        ----------
        other   : Bitstring
            other Bitstring
        Returns
        -------
        equality  : boolean
            Boolean value showing equality or not
        """
        self_dec = self.int()
        other_dec = other.int()

        return self_dec == other_dec

    
    def __len__(self):
        """
        Calculate the length of a BitString object

        Returns
        -------
        length  : number
            Length of the bitstring
        """
        return len(self.config)

    def on(self):
        """
        Calculate the number of bits of value 1 in the bitstring

        Returns
        -------
        ons  : number
            Number of 1 bits
        """
        total = 0
        for bit in self.config:
            if bit == 1:
                total += 1
        return total
    
    def off(self):
        """
        Calculate the number of bits of value 1 in the bitstring

        Returns
        -------
        ons  : number
            Number of 1 bits
        """
        total = 0
        for bit in self.config:
            if bit == 0:
                total += 1
        return total
    
    def flip_site(self,i):
        """
        Flip a bit at a given index

        Parameters
        ----------
        i   : number
            index to flip
        """
        self.config[i] = (self.config[i] + 1) % 2
            
    
    def int(self):
        """
        Calculate the actual integer value of the current bitstring

        Returns
        ----------
        int   : number
            Integer value
        """
        temp = self.config[::-1]
        total = 0
        for i, bit in enumerate(temp):
            total += 2**i * bit

        return total
 

    def set_config(self, s:list[int]):
        """
        Set the bitstring by inputting a list

        Parameters
        ----------
        s   : list of ints
            List of bit values to configure
        """
        self.config = s;
        
    def set_int_config(self, dec:int):
        """
        Set the bitstring to match a given integer

        Parameters
        ----------
        dec   : integer
            Integer for the bitstring to match value of
        """
        bits = []
        for bit in self.config:
            remainder = dec % 2
            bits.insert(0, remainder)
            dec = dec // 2
        self.set_config(bits)




class IsingHamiltonian:
    """
    Class that represents an Ising Hamiltonian and provides various operations along with it.
    """
    def __init__(self, G: list[list[int]], mus: list[int]):
        """
        Constructor for IsingHamiltonians
        
        Parameters
        ----------
        G    : Graph
            input graph defining the Hamiltonian
        mus   : Bitstring
            input configuration
        """
        self.G = G

        self.bs = BitString(len(mus))
        self.bs.set_config(mus)
    
    def energy(self):
        """Compute energy of configuration, `bs`

            .. math::
                E = \\left<\\hat{H}\\right>

        Returns
        -------
        energy  : float
            Energy of the input configuration
        """
        sum = 0.0

        for i, node in enumerate(self.G):
            for edge_tuple in node:
                spin1 = 2 * self.bs.config[i] - 1
                spin2 = 2 * self.bs.config[edge_tuple[0]] - 1
                weight = edge_tuple[1]

                sum += spin1 * spin2 * weight
        
        return sum / 2
    


    def compute_average_values(self, T: float):
        """
        Compute the average value of Energy, Magnetization, 
        Heat Capacity, and Magnetic Susceptibility 

            .. math::
                E = \\left<\\hat{H}\\right>

        Parameters
        ----------
        T    : float
            temperature of the system
        Returns
        -------
        energy  : float
        magnetization  : float
        heat capacity  : float
        magnetic susceptibility  : float
        """

        Z = 0

        for i in range(2 ** len(self.bs.config)):
            self.bs.set_int_config(i)
            e = self.energy()
            Z += np.exp(-(1/T) * e)

        
        E = 0
        M = 0
        E_squared = 0
        M_squared = 0

        for i in range(2 ** len(self.bs.config)):
            self.bs.set_int_config(i)
            e = self.energy()

            this_P = np.exp(-(1/T) * e) / Z 
            this_M = self.bs.on() - self.bs.off()

            E += e * this_P

            M += this_M * this_P

            E_squared += e * e * this_P
            M_squared += this_M * this_M  * this_P


        HC = (E_squared - E ** 2) * (T ** -2)
        MS = (M_squared - M ** 2) * (T ** -1)
        

        return E, M, HC, MS






