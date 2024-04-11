# SequenceParser.py
from abc import ABC, abstractmethod
from collections import Counter
import logging
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
import os
# Import all definitions from iupac_codes.py
from iupac_codes import *

class SequenceParser(ABC):
    """An abstract base class for parsing sequence files and analyzing sequence data."""

    def __init__(self, filename):
        """
        Initialize the SequenceParser with the filename.

        Parameters:
        - filename (str): The path to the sequence file.
        """
        self.filename = filename # Store the filename
        self.sequences = []  # Initialize an empty list to store sequences

        # Configure logging
        self._configure_logging()

    def _configure_logging(self):
        """Configure logging settings."""
        log_file = os.path.splitext(os.path.basename(self.filename))[0] + '_parser.log'
        logging.basicConfig(filename=log_file, level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')

    @abstractmethod
    def read_sequences(self):
        """
        Read sequences from the file.

        This method must be implemented by subclasses to handle specific file formats.
        """
        pass  # Placeholder method to be implemented by subclasses

    def _validate_file(self):
        """Validate the input file."""
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"File not found: {self.filename}")
        if os.path.getsize(self.filename) == 0:
            raise ValueError("File is empty")

    def _validate_sequences(self):
        """Validate the parsed sequences."""
        if not self.sequences:
            raise ValueError("No sequences found")
        if any(not seq for seq in self.sequences):
            raise ValueError("Empty sequences found")
    
    def get_sequences(self):
        """
        Get the list of sequences.

        Returns:
        - list: A list of sequences.
        """
        return self.sequences  # Return the list of sequences

    def get_sequence_lengths(self):
        """
        Get the lengths of the sequences.

        Returns:
        - list: A list of sequence lengths.
        """
        return [len(seq) for seq in self.sequences]  # Return a list of sequence lengths
    
    def get_sequence_count(self):
        """
        Get the number of sequences in the file.

        Returns:
        - int: The number of sequences.
        """
        # Validate sequences
        self._validate_sequences()

        # Return the number of sequences
        return len(self.sequences)

    def get_nucleotide_distribution(self):
        """
        Get the distribution of nucleotides in the sequences.

        Returns:
        - dict: A dictionary where keys are nucleotides (A, T, G, C) and values are lists of nucleotide counts at each position.
        """
        # Validate sequences
        self._validate_sequences()

        # Initialize a dictionary to store nucleotide counts
        nucleotide_distribution = {base: [] for base in Standard_Nucleotides}
        # Loop over each sequence and count occurrences of each nucleotide in the sequence
        for seq in self.sequences:
            counts = Counter(seq)
            # Loop over each nucleotide and append count to the respective nucleotide's list
            for base in Standard_Nucleotides:
                nucleotide_distribution[base].append(counts[base])
        # Return the nucleotide distribution dictionary
        return nucleotide_distribution

    def plot_nucleotide_distribution(self):
        """
        Plot the distribution of nucleotides in the sequences.
        """
        # Validate sequences
        self._validate_sequences()

        # Get nucleotide distribution dictionary
        nucleotide_distribution = self.get_nucleotide_distribution()
        if not nucleotide_distribution:
            logging.warning(
                "No sequences to analyze for nucleotide distribution.")
            print("No sequences to analyze for nucleotide distribution.")
            return
        # Create array of sequence indices
        positions = np.arange(len(self.sequences))
        width = 0.2
        # Loop over each nucleotide and plot bar chart for each nucleotide's count at each sequence index
        for base in nucleotide_distribution:
            plt.bar(positions + width * (ord(base) - ord('A')),
                    nucleotide_distribution[base], width, label=base)

        plt.xlabel('Sequence Index')
        plt.ylabel('Count')
        plt.title('Nucleotide Distribution')
        plt.legend()
        plt.show()

    def compute_gc_content(self):
        """
        Compute the GC content for each sequence.

        Returns:
        - list: A list of GC percentages for each sequence.
        """
        # Validate sequences
        self._validate_sequences()

        # Compute GC content for each sequence using Biopython's GC function
        gc_percentages = [GC(seq) for seq in self.sequences]
        return gc_percentages  # Return list of GC percentages
    
    def find_gc_by_position(self):
        """
        Find the GC ratio at each position in the reads and plot the results.

        Returns:
        - tuple: A tuple containing two lists: positions and corresponding GC content.
        """
        # Compute the GC content for each sequence
        gc_content = self.compute_gc_content()

        # Transpose the list of GC content to get a list of lists where each sublist contains GC content at the same position
        gc_content_by_position = list(zip(*gc_content))

        # Compute the average GC content at each position
        avg_gc_content = [sum(gc) / len(gc) for gc in gc_content_by_position]

        # Plot the GC content by position
        plt.plot(range(len(avg_gc_content)), avg_gc_content)
        plt.xlabel('Position in Read')
        plt.ylabel('GC Content')
        plt.title('GC Content by Position')
        plt.show()

        return range(len(avg_gc_content)), avg_gc_content


class FastaParser(SequenceParser):
    """A class for parsing FASTA sequence files."""


    def __init__(self, filename):
        """
        Initialize the FastaParser with the filename.

        Parameters:
        - filename (str): The path to the FASTA sequence file.
        """
        super().__init__(filename)

    def read_sequences(self):
        """
        Read sequences from the FASTA file.

        This method reads sequences from the FASTA file using Biopython's SeqIO module.
        """
        # Validate the file
        self._validate_file()

        # Use Biopython's SeqIO to parse the FASTA file
        with open(self.filename, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                self.sequences.append(str(record.seq))

    def read_sequences(self):
        """Read sequences from the FASTA file and store them in a dictionary."""
        # Validate the input file
        self._validate_file()

        # Parse the FASTA file and store records in a dictionary
        self.record_dict = SeqIO.to_dict(SeqIO.parse(self.filename, "fasta"))

        # Update the list of sequences
        self.sequences = list(self.record_dict.values())

    def get_sequence_by_header(self, header):
        """
        Get the sequence corresponding to a given header.

        Parameters:
        - header (str): The header of the sequence.

        Returns:
        - str: The sequence corresponding to the header.
        """
        if self.record_dict is None:
            raise ValueError("Sequences have not been read from file yet")

        if header not in self.record_dict:
            raise KeyError(f"Header '{header}' not found in the FASTA file")

        return str(self.record_dict[header].seq)


class FastaParser(SequenceParser):
    """A class for parsing FASTA sequence files."""

    def read_sequences(self):
        """Read sequences from the FASTA file."""
        # Validate the file
        self._validate_file()

        # Reset sequences list
        self.sequences = []

        try:
            # Read sequences from the FASTA file
            with open(self.filename, 'r') as file:
                for record in SeqIO.parse(file, "fasta"):
                    self.sequences.append(str(record.seq))
        except FileNotFoundError:
            logging.error(f"File not found: {self.filename}")
            raise
        except Exception as e:
            logging.error(f"Error reading sequences from FASTA file: {e}")
            raise

    def get_headers(self):
        """Get the headers of the sequences."""
        try:
            # Read headers from the FASTA file
            with open(self.filename, 'r') as file:
                headers = [record.id for record in SeqIO.parse(file, "fasta")]
            return headers
        except FileNotFoundError:
            logging.error(f"File not found: {self.filename}")
            raise
        except Exception as e:
            logging.error(f"Error reading headers from FASTA file: {e}")
            raise

    def get_sequence_by_header(self, header):
        """
        Get the sequence corresponding to the given header.

        Parameters:
        - header (str): The header of the sequence to retrieve.

        Returns:
        - str: The sequence corresponding to the header.
        """
        # Validate header
        if not header:
            raise ValueError("Header is empty")

        try:
            # Read sequences and headers from the FASTA file
            with open(self.filename, 'r') as file:
                for record in SeqIO.parse(file, "fasta"):
                    if record.id == header:
                        return str(record.seq)
            # If header not found
            raise ValueError(f"Header '{header}' not found in FASTA file")
        except FileNotFoundError:
            logging.error(f"File not found: {self.filename}")
            raise
        except Exception as e:
            logging.error(f"Error reading sequences from FASTA file: {e}")
            raise

    def plot_sequence_length_distribution(self):
        """
        Plot the distribution of sequence lengths in the FASTA file.
        """
        # Validate sequences
        self._validate_sequences()

        # Get sequence lengths
        sequence_lengths = self.get_sequence_lengths()

        # Plot histogram of sequence lengths
        plt.hist(sequence_lengths, bins=20, color='skyblue', edgecolor='black')
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.show()

    def compute_sequence_statistics(self):
        """
        Compute statistics of sequences in the FASTA file.

        Returns:
        - dict: A dictionary containing sequence statistics (e.g., min length, max length, average length).
        """
        # Validate sequences
        self._validate_sequences()

        # Get sequence lengths
        sequence_lengths = self.get_sequence_lengths()

        # Compute statistics
        stats = {
            'min_length': min(sequence_lengths),
            'max_length': max(sequence_lengths),
            'average_length': sum(sequence_lengths) / len(sequence_lengths)
        }
        return stats
