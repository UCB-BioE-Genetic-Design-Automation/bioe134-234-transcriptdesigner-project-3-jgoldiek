import random
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker  
from genedesign.checkers.gc_content_checker import GCContentChecker
from genedesign.models.transcript import Transcript
from genedesign.seq_utils.hairpin_counter import hairpin_counter

# Importing the hairpin checker directly as a function
from genedesign.checkers.hairpin_checker import hairpin_checker 


class TranscriptDesigner:
    def __init__(self):
        # 1. The original required variables (Don't change these names!)
        self.aminoAcidToCodon = {}
        self.rbsChooser = None

        # 2. Your new checker tools
        self.codon_checker = CodonChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.gc_checker = GCContentChecker()
    
    def initiate(self):
        """
        Signals readiness to the benchmarker and sets up the required data.
        """
        # 1. Populate the dictionary the original code expects
        self.aminoAcidToCodon = {
            'A': ['GCG', 'GCC', 'GCA', 'GCT'],
            'C': ['TGC', 'TGT'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGC', 'GGT', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
            'M': ['ATG'],
            'N': ['AAC', 'AAT'],
            'P': ['CCG', 'CCA', 'CCT', 'CCC'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'S': ['AGC', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT'],
            'T': ['ACC', 'ACG', 'ACT', 'ACA'],
            'V': ['GTG', 'GTT', 'GTC', 'GTA'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC'],
            '*': ['TAA', 'TGA', 'TAG']
        }
        
        # 2. Initiate all the class-based checkers so they load their rules!
        self.codon_checker.initiate()
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()
        #self.gc_checker.initiate()
        
        # (Note: hairpin_checker is a function, so it doesn't need to be initiated)
        
        # 3. If you have an RBS chooser, initiate it if it exists
        if self.rbsChooser is not None:
            self.rbsChooser.initiate()

    def naive_translate(self, protein_sequence):
        cds = []
        for i, aa in enumerate(protein_sequence):
            options = self.aminoAcidToCodon.get(aa.upper(), ['NNN'])
            safe = [c for c in options if c not in self.codon_checker.rare_codons]
            pool = safe if safe else options
            
            best_codon = pool[0]
            best_count = float('inf')
            for candidate in pool:
                test = cds + [candidate]
                local_seq = ''.join(test)
                start = max(0, len(local_seq) - 50)
                count, _ = hairpin_counter(local_seq[start:], 3, 4, 9)
                if count < best_count:
                    best_count = count
                    best_codon = candidate
            cds.append(best_codon)
        return cds

    def mutate_sequence(self, current_cds, protein_sequence, mutation_rate=0.05):
        """Randomly swaps a few codons for synonymous ones to fix biological errors."""
        new_cds = current_cds.copy()
        
        num_mutations = max(1, int(len(protein_sequence) * mutation_rate))

        for _ in range(num_mutations):
            idx = random.randint(0, len(protein_sequence) - 1)
            aa = protein_sequence[idx].upper()

            possible_codons = self.aminoAcidToCodon.get(aa, [current_cds[idx]])

            if len(possible_codons) > 1:
                choices = [c for c in possible_codons if c != current_cds[idx]]
                if choices:
                    new_cds[idx] = random.choice(choices)

        return new_cds
    
    def smart_restart(self, safe_mutation_map, peptide_len):
        from genedesign.seq_utils.reverse_complement import reverse_complement
        codons = []
        for i in range(peptide_len):
            options = safe_mutation_map[i].copy()
            random.shuffle(options)
            # Avoid codons that are reverse complements of anything in the last 15 codons (45bp window)
            window = codons[max(0, i-15):i]
            bad = {reverse_complement(c) for c in window}
            preferred = [c for c in options if c not in bad]
            codons.append(preferred[0] if preferred else options[0])
        return codons

    def run(self, peptide: str, ignores: set):
        """
        Translates the peptide to DNA and optimizes it using an Ultra-Fast 
        Error-State Random Walk with Short-Circuiting and Lazy Evaluation.
        """
        try:
            import random
            import re
            from collections import Counter
            
            clean_regex = re.compile(r'[^ATCGatcg]')
            peptide_len = len(peptide)
            
            # Pre-compute safe mutations to avoid lists in the loop
            safe_mutation_map = {}
            for i, aa in enumerate(peptide):
                aa_upper = aa.upper()
                possible_codons = self.aminoAcidToCodon.get(aa_upper, [])
                # WITH THIS:
                safe = [c for c in possible_codons 
                        if c not in self.codon_checker.rare_codons
                        and self.codon_checker.codon_frequencies.get(c, 0.01) > 0.0]
                safe_mutation_map[i] = safe if safe else possible_codons

            # ==========================================
            # 1. THE SHORT-CIRCUIT EVALUATOR
            # ==========================================
            def get_errors_and_diags(test_codons, best_total_errors=None):
                cds_string = ''.join(test_codons)
                
                # Run Fast Checkers First
                c_res = self.codon_checker.run(test_codons[:-1])
                f_res = self.forbidden_checker.run(cds_string)
                p_res = self.promoter_checker.run(cds_string)
                g_res = self.gc_checker.run(cds_string)
                
                p_c = c_res[0] if isinstance(c_res, tuple) else c_res
                diversity = c_res[1] if isinstance(c_res, tuple) and len(c_res) > 1 else 0.0
                rare_count = c_res[2] if isinstance(c_res, tuple) and len(c_res) > 2 else 0
                cai = c_res[3] if isinstance(c_res, tuple) and len(c_res) > 3 else 0.0
                
                p_f = f_res[0] if isinstance(f_res, tuple) else f_res
                p_p = p_res[0] if isinstance(p_res, tuple) else p_res
                p_g = g_res[0] if isinstance(g_res, tuple) else g_res
                
                fast_errors = 0
                if not p_f: fast_errors += 1
                if not p_p: fast_errors += 1
                if not p_g: fast_errors += 1
                if rare_count > 3: fast_errors += 1
                if diversity < 0.5: fast_errors += 1
                if cai < 0.2: fast_errors += 1
                
                # THE SHORT CIRCUIT: If fast errors alone make this worse than our best, 
                # ABORT instantly. Do not run the hairpin checker!
                if best_total_errors is not None and fast_errors > best_total_errors:
                    return fast_errors + 999, cai, False, None
                
                # Run Slow Checker (Only if the fast checkers passed the threshold)
                h_res = hairpin_checker(cds_string)
                p_h = h_res[0] if isinstance(h_res, tuple) else h_res
                
                slow_errors = 0
                if not p_h: slow_errors += 1
                
                total_errors = fast_errors + slow_errors
                passed_all = (total_errors == 0 and p_c)
                
                return total_errors, cai, passed_all, (c_res, f_res, p_res, g_res, h_res)

            # ==========================================
            # 2. THE LAZY TARGETER
            # ==========================================
            def get_bad_indices(test_codons, diags):
                # We only run this expensive Regex parsing when we actually KEEP a mutation
                c_res, f_res, p_res, g_res, h_res = diags
                cds_string = ''.join(test_codons)
                
                p_c = c_res[0] if isinstance(c_res, tuple) else c_res
                diversity = c_res[1] if isinstance(c_res, tuple) and len(c_res) > 1 else 0.0
                rare_count = c_res[2] if isinstance(c_res, tuple) and len(c_res) > 2 else 0
                
                p_f = f_res[0] if isinstance(f_res, tuple) else f_res
                p_p = p_res[0] if isinstance(p_res, tuple) else p_res
                p_h = h_res[0] if isinstance(h_res, tuple) else h_res
                
                bad_indices = set()
                def extract(data):
                    bad = []
                    if isinstance(data, str):
                        for line in data.split('\n'):
                            if not line.strip(): continue
                            seq_str = line.split(':')[1] if ':' in line else line
                            clean_seq = clean_regex.sub('', seq_str).upper()
                            if len(clean_seq) >= 4:
                                idx = cds_string.find(clean_seq)
                                if idx != -1:
                                    bad.extend(list(range(idx // 3, (idx + len(clean_seq)) // 3 + 1)))
                    return bad

                
                '''if not p_h:
                    chunk_size = 50
                    overlap = 25
                    cds_string = ''.join(test_codons)
                    for i in range(0, len(cds_string) - chunk_size + 1, overlap):
                        chunk = cds_string[i:i + chunk_size]
                        count, _ = hairpin_counter(chunk, 3, 4, 9)
                        if count > 1:
                            start_codon = i // 3
                            end_codon = min((i + chunk_size) // 3 + 1, peptide_len)
                            bad_indices.update(range(start_codon, end_codon))
                            break  # Target the first failing chunk; subsequent mutations fix the rest'''
                if not p_h: bad_indices.update(extract(h_res[1] if isinstance(h_res, tuple) else h_res))
                if not p_f: bad_indices.update(extract(f_res[1] if isinstance(f_res, tuple) else f_res))
                if not p_p: bad_indices.update(extract(p_res[1] if isinstance(p_res, tuple) else p_res))
                
                if not bad_indices:
                    if rare_count > 3:
                        for i, codon in enumerate(test_codons[:-1]):
                            if codon in self.codon_checker.rare_codons:
                                bad_indices.add(i)
                    elif diversity < 0.5:
                        counts = Counter(test_codons[:-1])
                        most_common = counts.most_common(1)[0][0]
                        for i, codon in enumerate(test_codons[:-1]):
                            if codon == most_common:
                                bad_indices.add(i)
                                
                return list(set([i for i in bad_indices if i < peptide_len]))

            # ==========================================
            # 3. THE OPTIMIZATION LOOP
            # ==========================================
            best_codons = self.naive_translate(peptide)
            best_codons.append("TAA")
            
            best_errors, best_cai, passed, diags = get_errors_and_diags(best_codons)
            if passed:
                selected_rbs = self.rbsChooser.run(''.join(best_codons), ignores) if self.rbsChooser else None
                return Transcript(selected_rbs, peptide, best_codons)
                
            best_bad_indices = get_bad_indices(best_codons, diags)
            
            plateau_count = 0
            PLATEAU_LIMIT = 150

            for attempt in range(1000):
                test_codons = best_codons.copy()
                # ADD THESE TWO LINES BACK:
                h_res = diags[4]
                hairpin_failed = not (h_res[0] if isinstance(h_res, tuple) else h_res)

                num_mutations = 1
                if best_bad_indices:
                    if hairpin_failed:
                        max_burst = min(5, len(best_bad_indices))
                        num_mutations = random.randint(2, max_burst)
                    else:
                        max_burst = min(3, len(best_bad_indices))
                        num_mutations = random.randint(1, max_burst)
                    target_indices = random.sample(best_bad_indices, num_mutations)
                else:
                    target_indices = [random.randint(0, peptide_len - 1)]

                valid_mutation_made = False
                for idx in target_indices:
                    available_moves = [c for c in safe_mutation_map[idx] if c != best_codons[idx]]
                    # REPLACE:
                    if available_moves:
                        if hairpin_failed:
                            from genedesign.seq_utils.reverse_complement import reverse_complement
                            window = test_codons[max(0, idx-15):idx] + test_codons[idx+1:min(peptide_len, idx+16)]
                            bad = {reverse_complement(c) for c in window}
                            preferred = [c for c in available_moves if c not in bad]
                            test_codons[idx] = random.choice(preferred) if preferred else random.choice(available_moves)
                        else:
                            test_codons[idx] = random.choice(available_moves)
                        valid_mutation_made = True
                if not valid_mutation_made:
                    continue 
                
                # Fast Evaluation (Keep the rest of your loop exactly the same below this!)
                test_errors, test_cai, test_passed, test_diags = get_errors_and_diags(test_codons, best_errors)
                
                if test_passed:
                    selected_rbs = self.rbsChooser.run(''.join(test_codons), ignores) if self.rbsChooser else None
                    return Transcript(selected_rbs, peptide, test_codons)
                    
                # Determine if we keep the mutation
                keep_mutation = False
                if test_errors < best_errors:
                    keep_mutation = True
                elif test_errors == best_errors:
                    if best_errors > 0:
                        keep_mutation = True
                    elif test_cai >= best_cai:
                        keep_mutation = True
                        
                # Only if we KEEP the mutation do we overwrite best state and run the Regex Lazy Targeter
                if keep_mutation:
                    plateau_count = 0
                    best_codons = test_codons
                    best_errors = test_errors
                    best_cai = test_cai
                    diags = test_diags
                    best_bad_indices = get_bad_indices(best_codons, diags)
                  
                else:
                    plateau_count += 1

                if plateau_count >= PLATEAU_LIMIT:
                    best_codons = self.smart_restart(safe_mutation_map, peptide_len)  # ← changed
                    best_codons.append("TAA")
                    best_errors, best_cai, passed, diags = get_errors_and_diags(best_codons)
                    plateau_count = 0
                    if passed:
                        selected_rbs = self.rbsChooser.run(''.join(best_codons), ignores) if self.rbsChooser else None
                        return Transcript(selected_rbs, peptide, best_codons)
                    best_bad_indices = get_bad_indices(best_codons, diags)

            # If it fails, unpack the diags and crash cleanly
            c_res, f_res, p_res, g_res, h_res = diags
            print(f"\n--- DIAGNOSTIC FOR {peptide[:10]} ---")
            print(f"Codon Output: {c_res}")
            print(f"Forbidden Output: {f_res}")
            print(f"Promoter Output: {p_res}")
            print(f"GC Output: {g_res}")
            print(f"Hairpin Output: {h_res}")
            raise ValueError(f"Failed on {peptide[:10]}. See diagnostic above.")

        except Exception as e:
            import traceback
            print("\n=== CRASH TRACEBACK ===")
            traceback.print_exc()
            raise e


'''
    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Codon table with highest CAI codon for each amino acid (for E. coli)
        self.aminoAcidToCodon = {
            'A': "GCG", 'C': "TGC", 'D': "GAT", 'E': "GAA", 'F': "TTC",
            'G': "GGT", 'H': "CAC", 'I': "ATC", 'K': "AAA", 'L': "CTG",
            'M': "ATG", 'N': "AAC", 'P': "CCG", 'Q': "CAG", 'R': "CGT",
            'S': "TCT", 'T': "ACC", 'V': "GTT", 'W': "TGG", 'Y': "TAC"
        }

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)


'''
