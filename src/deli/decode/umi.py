"""umi functions and classes"""

nuc_to_binary = {"A": "00", "T": "01", "G": "10", "C": "11"}

COMPRESSION_BIT_SIZE = 6
COMPRESSION_ASCII_OFFSET = 35


class UMI:
    """Call class for UMIs"""

    def __init__(self, umi_tag: str):
        """
        Initialize a UMICall

        Parameters
        ----------
        umi_tag: str
            the DNA observed_seq of the umi tag
        """
        self.umi_tag = umi_tag

    def __str__(self):
        """Return the DNA observed_seq of the UMI as a string"""
        return self.umi_tag

    def __len__(self):
        """Return the length of the UMI tag"""
        return len(self.umi_tag)

    def to_ascii_code(self) -> str:
        """
        Convert the UMI tag to a printable ASCII string representation

        Each nucleotide is converted to 2 bits (A=00, G=01, ...), then groups of 6 bits
        are mapped to printable ASCII characters, starting with 35 ('#').
        This can help reduce the size of the UMI tag when storing in raw files
        by achieving 3:1 compression.

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
        for i in range(0, len(bit_str), COMPRESSION_BIT_SIZE):
            _sub_bit_str = bit_str[i : i + COMPRESSION_BIT_SIZE].ljust(COMPRESSION_BIT_SIZE, "0")
            code += chr(
                int(_sub_bit_str, 2) + COMPRESSION_ASCII_OFFSET
            )  # ASCII offset for printable characters skipping '"'
        return code
