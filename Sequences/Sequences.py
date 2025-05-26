class SequenceUser:
    """Base class for biological sequence representation entered by the user.

        Attributes:
            _sequence (str): The stored sequence in uppercase
        """
    def __init__(self,sequence : str):
        """Initialize the sequence.

                Args:
                    sequence (str): The input biological sequence
                """
        self._sequence : str = sequence.upper()

    def __str__(self):
        return self._sequence

    def __repr__(self):
        return str(self)

    def seq(self):
        """Get the stored sequence.

                       Returns:
                           str: The stored sequence
                       """
        return self._sequence

class ProteinSequence(SequenceUser):
    """Class representing a protein sequence with FASTA identifier.

        Attributes:
            identifier (str): FASTA header/identifier line
            data (str): The protein sequence data
        """
    def __init__(self, identifier: str, data: str):
        self.identifier = identifier
        self.data = data

    def __str__(self):
        return f"{self.identifier}\n{self.data}"

def fasta_to_pro_seq(fasta_file):
    """Read a protein sequence from a FASTA file.

       Args:
           fasta_file (str): Path to the FASTA file

       Returns:
           ProteinSequence: Protein sequence object or None if reading fails
       """
    try:
        with open(fasta_file, "r") as file:
            identifier = file.readline()
            content = file.read()
            if identifier and content:
                return ProteinSequence(identifier, content)
            else:
                print("missing identifier or _sequence.")
                return None
    except FileNotFoundError:
        print("File Not Found")
        return None
    except PermissionError:
        print("Permission Denied")
        return None

def convert_to_sequence(seq_str):
    """Convert input string to SequenceUser object consistently"""
    return SequenceUser(seq_str.upper())
