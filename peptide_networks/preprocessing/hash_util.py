from pandas import read_csv
from numpy import float64

class HashString(object):
    AMINO_ACIDS = {'A':0, 'R': 1, 'N':2,'D':3,'C':4, 'Q':5,'E':6,'G':7, 'H':8, 'I':9, 'L':10, 'K': 11, 'M': 12, 'F':13,'P':14, 'S':15, 'T':16, 'W':17,'Y':18,'V':19}
    MOD = 2**64-1
    BASE = 20
    def __init__(self, string_to_hash=""):
        """
        Class to compute rolling hash values for peptide sequences. 
        The hash function can be rolled over longer sequences by repeatedly
            using the pop_front and insert methods.
        The size and hash (probably) uniquely identify the peptide.
        Use the get-methods to fetch attributes. 

        Parameters
        ----------
        string_to_hash : TYPE, optional
            Optional starting peptide sequence. The default is "".

        Returns
        -------
        None.

        """
        string_to_hash = string_to_hash.upper()
        self.hash_value = 0
        self.first_index = 0
        self.size = 0
        self.charstring = []
        
        while(self.size < len(string_to_hash)):
            self.insert(string_to_hash[self.size])

    def _po(self,a,b):
        """
        Fast way of computing large powers.

        Parameters
        ----------
        a : float
            Base.
        b : int
            Exponent.

        Returns
        -------
        float
            a**b.

        """
        if b == 0:
            return 1
        c = self._po(a,b // 2)
        if b % 2 == 1:
            return (((c*c)%self.MOD)*a) % self.MOD
        else:
            return (c*c)%self.MOD

    def __eq__(self,other):
        return (self.size == other.size) and (self.hash_value == other.hash_value)
    
    def __str__(self):
        return ''.join(self.charstring) + " : " + str(self.hash_value)

    def insert(self,char):
        """
        Inserts a character at the end of the sequence.

        Parameters
        ----------
        char : char
            Character to insert.

        Returns
        -------
        None.

        """
        char = char.upper()
        self.hash_value *= self.BASE
        self.hash_value += self.AMINO_ACIDS[char]
        self.hash_value %= self.MOD
        self.charstring.append(char)
        self.size += 1

    def pop_front(self):
        """
        Removes the first character in the string.

        Raises
        ------
        IndexError
            When the hash string already is of length 0.

        Returns
        -------
        char
            The removed character.
        """
        
        
        if self.size == 0:
            raise IndexError('Unable to pop HashString of length 0')
            
        self.hash_value -= self.AMINO_ACIDS[self.charstring[self.first_index]]*self._po(self.BASE, self.size - 1)
        while(self.hash_value < 0):
            self.hash_value += self.MOD
        self.hash_value %= self.MOD
        self.first_index += 1
        self.size -= 1
        return self.charstring[self.first_index-1]

    def pop_back(self):
        """
        Removes last character in sequence.

        Raises
        ------
        IndexError
            When sequence already is of length 0.

        Returns
        -------
        char
            The removed character.

        """
        if self.size == 0:
            raise IndexError('Unable to pop HashString of length 0')
            
        self.hash_value -= (self.AMINO_ACIDS[self.charstring[self.size-1]])
        self.hash_value = self.hash_value // self.BASE 
        while(self.hash_value < 0):
            self.hash_value += self.MOD 
        self.hash_value %= self.MOD 
        self.size -= 1
        return self.charstring.pop()

    def getString(self):
        return ''.join(self.charstring)
    
    def getHash(self):
        return self.hash_value
    
    def getSize(self):
        return self.size

class PeptideSequence(HashString):
    def __init__(self, sequence):
        """
        Utility class to easily fetch and hash peptide sequences. 
        When querying the clean dataframes, use the loc method as follows:
            peptide_sequence = PeptideSequence('GATCA')
            info = clean_df.loc[peptide_sequence.loc(),:]

        Parameters
        ----------
        sequence : string
            Peptide sequence to be hashed.

        Returns
        -------
        None.
        """
        super().__init__(sequence)
    
    def loc(self):
        return (self.getSize(), self.getHash())
    
    
def import_clean_data(path):
    return read_csv(path, index_col=[0,1], header=0)