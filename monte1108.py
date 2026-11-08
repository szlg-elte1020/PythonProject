# whats monte carlo project?
# a computational technique that uses repeated random sampling to
# calculate probabilities and observe trends in complex, unpredictable systems.
# Instead of solving an equation mathematically,
# a Monte Carlo algorithm "rolls the dice" thousands of times in a computer simulation
# to see what outcomes happen most frequently.

import random  # a standard library #module for generating pseudo-random numbers
# in this case to simulate real-world genetic variation
# Provided code: codon table and translator
# codon_table is a python #dictionary where key are the codons
# and values are amino acids the codons encode
CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}
BASES = "ACGT"
TEST_SEQUENCE = (
    "ATGGCCTTTAAAGGCTGGCTCGAGAAACCC"
    "GGGTACCATCGTAGCGATCAAATTGTGTAG"
)

#leftover bases at the end (<3) are ignored when translating DNA into AA
def translate(seq):
    protein = []  # initialize empty list to store sequence
    for i in range(0, len(seq) - 2, 3):  # DNA string starts at index 0 moving with 3 steps
        # stops at len(seq) and ignores partial codons
        amino_acid = CODON_TABLE[seq[i:i + 3]]  # extracts bases 0,1,2 =slice
        if amino_acid == "Stop":
            break  # checks if its a stop codon
        protein.append(amino_acid)  # appends valid amino acid names to the protein list
    return protein  # returns the array of translated amino acids

# Role 1: Mutator # mutator functions handle the altering DNA sequences
# but python strings are immutable thus string slicing is needed to construct a new modified string
def substitute(seq):  # whats len and getitem?
    # swaps out one base at random position for a different one
    pos = random.randrange(len(seq))  # picks a valid index from 0 up to len(seq)-1
    new_base = random.choice([b for b in BASES if b != seq[pos]])  # using list comprehension
    # to filter bases so only the different ones from seq[pos] are kept
    return seq[:pos] + new_base + seq[pos + 1:], pos  # first we slice up to index pos
    # + we take the new letter #+from index after pos to the end of sequence
    # we glue them togethet -> only replacing one letter

# mimicks natural transition/transversion bias (in 2:1 ratio)
def realistic_substitute(seq): #accepts DNA string as input
    pos = random.randrange(len(seq)) #selects a random integer index
    old_base = seq[pos] #stoes the original base at index pos
    transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}  # natural transition pairs
    #defines swapping within same chemical class (purine with purine, pyrimidine <-> pyrimidine)
    transversions = {  # transversion alternatives #swapping across chemical classes
        'A': ['C', 'T'],
        'G': ['C', 'T'],
        'C': ['A', 'G'],
        'T': ['A', 'G']
    } #purine with pyrimidine giving 2 possible outcomes for each
    # 66% chance of transition, 34% chance of transversion
    if random.random() < 0.66:  #geenrates random float between 0.0 -1.0
        new_base = transitions[old_base] #assigns direct transition partner
    else:
        new_base = random.choice(transversions[old_base]) #randonly picks one of the 2 options
    return seq[:pos] + new_base + seq[pos + 1:], pos #mut. seq and mut. index
#reconstructs string using strict slicing returning 2 element tuple

def insert(seq):  # unlike substitute which skips the character to replace it
    # insert adds an extra base making the seq. one character longer
    pos = random.randrange(len(seq) + 1)  # -1 would be until stop signal
    # while +1 is at the end of sequence, since there are always n+1 possible insertions spots
    # python thus "chooses" spots equally
    new_base = random.choice(BASES)
    return seq[:pos] + new_base + seq[pos:], pos
    # so insert doesn't lose letters, while substitution skips the original charachter

def delete(seq):  # we remove one base at a random position
    pos = random.randrange(len(seq))  # picks a random valid index
    return seq[:pos] + seq[pos + 1:], pos  # again index [pos] isnt included
    # the left and right pieces are glued directly onto each other

def random_mutation(seq):  # consider it like a master function
    # thus we dont need to manually pick which functions are needed
    kind = random.choice(["substitution", "insertion", "deletion"])  # kind means equal probability to choose from
    if kind == "substitution":  # checks this condition first, if its true it runs
        new_seq, pos = substitute(seq)
    elif kind == "insertion":  # if the *if* condition was false
        new_seq, pos = insert(seq)
    else:  # if none of the conditions were true, the this is default fallback
        new_seq, pos = delete(seq)
    return new_seq, kind, pos  # returns 3-item tuple back
    # kind means here a string stating which type of mutation occured
    # pos is the index integer where the mutation occured

# AA biochemical properties
AA_PROPERTIES = {  #dictionary as a lookup table
    # non-polar and hydrophobic
    "Phe": "nonpolar", "Leu": "nonpolar", "Ile": "nonpolar", "Val": "nonpolar",
    "Met": "nonpolar", "Trp": "nonpolar", "Ala": "nonpolar", "Pro": "nonpolar",
    "Ser": "polar", "Thr": "polar", "Tyr": "polar", "Asn": "polar", "Gln": "polar", "Cys": "polar",  # polar uncharged
    "Lys": "basic", "Arg": "basic", "His": "basic",  # basic has positive charge
    "Asp": "acidic", "Glu": "acidic",  # acidic has negative charge
    "Gly": "special"  # glycine
}
def classify_missense_impact(orig_aa, mut_aa):  # determines if missense mutation is convervative or radical
    prop_orig = AA_PROPERTIES.get(orig_aa, "unknown") #.get() prevents key errors
    prop_mut = AA_PROPERTIES.get(mut_aa, "unknown") #we look up biochemical category
    if prop_orig == prop_mut: #are they the same chemical class?
        return "conservative_missense" #no functional impact
    return "radical_missense" #properties differ protein is destabilized or altered

# Role 2: Classifier #did mutation alter the protein?
# it compares original with the mutated protein list and checks the mutation type
def classify(original_protein, mutated_protein, kind):
    # kind indicates the mutations #def classifies
    if original_protein == mutated_protein:
        return "silent"  # substitution is equal to original list
        # this is when the codons encode the same amino acid
    if kind in ("insertion", "deletion"):  # if added a nucleotide base
        return "frameshift"  # this one is simpler, no need to compare protein string
        # it returns frameshift if we add or delete a base #no protein list needed
    elif len(mutated_protein) < len(original_protein):
        return "nonsense"  # translation stops earlier than usual
    else:
        return "missense"

# Role 3: Simulator #is a program that can model processes repeatedly under conditions
# here its random sampling to estimate probabilities because a single seq testing doesnt reveal a full picture
def run_simulation(seq, trials, mode="any"):
    # this automated steps completes four steps over many (1000) iterations/trials
    counts = {"silent": 0, "missense": 0, "nonsense": 0, "frameshift": 0}
    # counts a dictionary score keeper starting at zero
    original_protein = translate(seq)  # creating reference point
    # the original sequence never changes during simulation
    # translates the original unmutated DNA into a protein list once before testing mutations
    for _ in range(trials):  # loop variable: we don't track iteration number
        if mode == "substitution":  # on each iteration which mode to run
            mut_seq, pos = substitute(seq)  # only tests base swaps
            kind = "substitution"
        else:
            mut_seq, kind, pos = random_mutation(seq)  # randomly pick between modes
        mutated_protein = translate(mut_seq)  # translate into amino acid list
        effect = classify(original_protein, mutated_protein, kind)
        # to determine which category returned it compares new to original list
        counts[effect] += 1  # adds score to dictionary
    return counts

# Role 4: Reporter #taking raw numerical data into human-readable output
def print_report(counts):  # converts dictionary to a clean terminal table
    total = sum(counts.values())  # calculates total completed trials/iterations
    if total == 0:  # to avoid ZeroDivisionError trials/=0
        print("No trials executed.")
        return  # to protect the program against crashes
    print(f"{'Effect':<12} | {'Count':<6} | {'Percentage':<10}")
    # we use f-string alignment modifiers to create columns 12-6-10 characters wide
    print("-" * 35)  # repeats hyphen cha 35 times drawing a horizontal line under heaeder
    for category, count in counts.items():
        pct = (count / total) * 100  # calculates each category's percentage
        print(f"{category.capitalize():<12} | {count:<6} | {pct:>6.2f}%")
        # first letter capitalized, right aligns 6 cha, 2f rounds floating-point numbers to 2 decimals

def plot_results(counts, title="Effects of random mutations"):
    import matplotlib.pyplot as plt  # imports plotting module locally inside the function
    # under standard alias plt
    # its only loaded to memory when the graph is actually generated
    plt.bar(list(counts.keys()), list(counts.values()))
    # converts keys to list to be labels along X-axis and extracts numbers for Y-axis
    plt.title(title)
    plt.xlabel("Effect on protein")
    plt.ylabel("Number of mutations")
    plt.show()  # rendering figure #user has to close the pop-up

# prints a single mutation step by terminal color highlight
def demo_one_mutation_colored(seq): #define function by single string input
    RED = "\033[91m" #ANSI escape codes
    GREEN = "\033[92m"
    ORANGE = "\033[38;5;208m"
    BLUE = "\033[34m"
    PINK="\033[95m"
    RESET = "\033[0m" #resets to default formatting
    orig_protein = translate(seq)
    mut_seq, kind, pos = random_mutation(seq) #unpacks tuple into variables
    mut_protein = translate(mut_seq)
    effect = classify(orig_protein, mut_protein, kind) #categorizes outcome
    if kind == "substitution":  # highlights mutated cha in DNA
        highlighted_mut = mut_seq[:pos] + PINK + mut_seq[pos] + RESET + mut_seq[pos + 1:]
    elif kind == "insertion":
        highlighted_mut = mut_seq[:pos] + GREEN + mut_seq[pos] + RESET + mut_seq[pos + 1:]
    else:  # deletion repr. for removed base
        highlighted_mut = mut_seq[:pos] + RED + "[*DELETED*]" + RESET + mut_seq[pos:]
    print("\n" + "=" * 45)
    print(f"{BLUE} SINGLE MUTATION DEMO{RESET}")
    print("=" * 45)
    print(f"Original DNA:\n{seq}")
    print(f"Mutated DNA:\n{highlighted_mut}")
    print(f"Mutation Event:\n{ORANGE}{kind.upper()}{RESET} at index {pos}")
    print(f"Original Protein: {'-'.join(orig_protein)}") #hyphen-separated strings
    print(f"Mutated Protein: {'-'.join(mut_protein)}") #like Met-Ala-Phe
    print(f"Outcome Effect: {PINK if effect in ('nonsense', 'frameshift') else GREEN}{effect.upper()}{RESET}")
    print("=" * 45 + "\n") #color codes the outcome category

# side by side chart for comparison
import matplotlib.pyplot as plt #importing visualization library
def plot_comparative_results(sub_counts, all_counts): #2 dict. arguments
    categories = ["silent", "missense", "nonsense", "frameshift"] #order on x axis
    sub_values = [sub_counts.get(c, 0) for c in categories] #extracts numerical counts
    all_values = [all_counts.get(c, 0) for c in categories] #for substitution trials in order
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True) #both charts share the same vertical scale
    # left substitutions only
    ax1.bar(categories, sub_values, color=['#2ecc71', '#3498db', '#e74c3c', '#95a5a6'])
    ax1.set_title("Substitutions Only (Point Mutations)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Count (Trials)")
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    # right every mutation
    ax2.bar(categories, all_values, color=['#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
    ax2.set_title("All Mutations (Sub + Insertion + Deletion)", fontsize=12, fontweight='bold')
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    plt.suptitle("Monte Carlo Simulation: Impact of DNA Mutations", fontsize=14, fontweight='bold')
    plt.tight_layout() #prevents subplots margins overlaps
    plt.show() #opens figure

def demo_one_mutation(seq):
    orig_protein = translate(seq)  # translate to base list of amino acids
    mut_seq, kind, pos = random_mutation(seq)  # unpacks tuple
    # new string-type of mutation-position
    # puts random operation in random position in sequence
    mut_protein = translate(mut_seq)
    effect = classify(orig_protein, mut_protein, kind)  # save results = into effect variable
    print("=== Single Mutation Demo ===")
    print(f"Original Sequence: {seq}")
    print(f"Original Protein: {orig_protein}")
    print(f"Mutation Applied: {kind} at index {pos}")
    print(f"Mutated Sequence: {mut_seq}")
    print(f"Mutated Protein: {mut_protein}")
    print(f"Effect Category: {effect}")
    print("======\n")

# Main program - play button, launches the program
if __name__ == "__main__":  # executed directly?prevents autorunning
    # executing indirectly means importing from a file as a module into another script
    # doesn't run as main script
    demo_one_mutation_colored(TEST_SEQUENCE)  # uses pre-defined test_seq
    print("=== Substitutions only ===")  # section divider
    sub_counts = run_simulation(TEST_SEQUENCE, 1000, mode="substitution")
    print_report(sub_counts)
    print()  # blank line for readibility
    print("=== All mutation types ===")
    all_counts = run_simulation(TEST_SEQUENCE, 1000, mode="any")
    print_report(all_counts)
    plot_results(sub_counts, "Effects of 1,000 random substitutions") #single plot for point mutations
    # substitutions is separate from random table to isolate effects of a single-letter point mutation
    # no frameshifts vs random table has frameshift
    plot_comparative_results(sub_counts, all_counts) #side by side comparison

# substitutions is separate from random table to isolate effects of a single-letter point mutation
# no frameshifts vs random table has frameshift
# this project is useful because bioinformaticians use randomized computational simulations
# to predict evolution or disease while #learning to build end-to-end data processing pipeline