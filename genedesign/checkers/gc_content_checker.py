class GCContentChecker:
    """
    Checks if the GC content of a given DNA sequence falls within an acceptable range.
    Standard synthetic biology constraints usually prefer 40% to 60% GC content 
    to ensure stable DNA synthesis and sequencing.
    """
    def __init__(self, min_gc=0.40, max_gc=0.60):
        # Set our acceptable biological thresholds
        self.min_gc = min_gc
        self.max_gc = max_gc

    def run(self, sequence):
        """
        Calculates the GC content of the sequence.
        
        Args:
            sequence (str or list): The DNA sequence as a string or a list of codons.
            
        Returns:
            passes_check (bool): True if GC content is within range, False otherwise.
            gc_content (float): The actual calculated GC percentage.
        """
        # If the sequence was passed as a list of codons, join it into a single string
        if isinstance(sequence, list):
            sequence = "".join(sequence)
            
        # Handle the edge case of an empty sequence to avoid dividing by zero!
        if len(sequence) == 0:
            return False, 0.0

        # Standardize to uppercase just in case
        sequence = sequence.upper()
        
        # Count the Gs and Cs
        g_count = sequence.count('G')
        c_count = sequence.count('C')
        
        # Calculate the percentage
        gc_content = (g_count + c_count) / len(sequence)
        
        # Check if it passes our defined thresholds
        passes_check = self.min_gc <= gc_content <= self.max_gc
        
        return passes_check, gc_content