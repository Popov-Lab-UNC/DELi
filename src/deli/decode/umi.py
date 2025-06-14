"""umi functions and classes"""

nuc_to_binary = {"A": "00", "T": "01", "G": "10", "C": "11"}


class UMI:
    """Call class for UMIs"""

    def __init__(self, umi_tag: str):
        """
        Initialize a UMICall

        Parameters
        ----------
        umi_tag: str
            the DNA sequence of the umi tag
        """
        self.umi_tag = umi_tag

    def __str__(self):
        """Return the DNA sequence of the UMI as a string"""
        return self.umi_tag

    def __len__(self):
        """Return the length of the UMI tag"""
        return len(self.umi_tag)

    def to_ascii_code(self) -> str:
        """
        Convert the UMI tag to a printable ASCII string representation

        This will covert every 3 byte nucleotide section of the UMI tag into a single byte.
        This can help reduce the size of the UMI tag when storing in raw files.

        Notes
        -----
        If the UMI tag is not a multiple of 3, it will be padded with
        '0's to make it a multiple of 3.

        Returns
        -------
        str
            the printable ASCII string representation of the UMI tag
        """
        bit_str = "".join([nuc_to_binary[nuc] for nuc in self.umi_tag])
        code = ""
        for i in range(0, len(bit_str), 6):
            _sub_bit_str = bit_str[i : i + 6].ljust(6, "0")
            code += chr(
                int(_sub_bit_str, 2) + 35
            )  # ASCII offset for printable characters skipping '"'
        return code
