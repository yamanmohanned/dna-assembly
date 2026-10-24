# Glossary

Generated from `docs/glossary.yaml`; do not edit by hand.

## accession and version

**In plain words.** The catalogue number of a sequence in the library. The number after the dot says which edition it is, so the same number always gives the same letters.

**Precisely.** An identifier such as NC_001422.1. The part before the dot names the record, the part after it is the version; a changed sequence gets a new version number.

## ALNS (adaptive large neighbourhood search)

**In plain words.** Repeatedly smash part of the solution (destroy) and rebuild it (repair), keep the result if it is good, and learn which smashing and rebuilding tools work best.

**Precisely.** Destroy operators: random, worst-connection, segment. Repair operators: greedy insertion and regret-k insertion. Simulated-annealing acceptance with a calibrated start temperature.

## ambiguity code (IUPAC)

**In plain words.** A letter that means "not sure which one", for example N for "any of the four". Real data sometimes contains it where the machine could not decide.

**Precisely.** An IUPAC nucleotide symbol other than A, C, G, T: R, Y, S, W, K, M, B, D, H, V or N. This project rejects genomes that contain any (decision D-016) rather than guessing.

## base (nucleotide letter)

**In plain words.** One of the four letters, A, C, G or T, that spell out DNA.

**Precisely.** A single symbol of the alphabet {A, C, G, T}. In this project a sequence is a string over this alphabet.

## bootstrap confidence interval

**In plain words.** Re-draw the results at random many times and see how much the answer wobbles; the middle 95% of the wobble is the interval.

**Precisely.** Percentile interval of the median paired difference over 10,000 resamples of the read sets.

## cache

**In plain words.** A saved copy kept close by, so we do not have to download it again and can work without internet.

**Precisely.** Here, the folder data/raw/ holding the FASTA file, its NCBI metadata and SHA256SUMS for each genome. Loading reads only the cache unless a download is explicitly requested.

## candidate list

**In plain words.** A short list of the most promising partners for each piece, so the search tries the good moves first instead of guessing among all of them.

**Precisely.** For every read the k reads with the largest overlap into it (predecessors) and out of it (successors). In tabu_cl a sampled insert move puts a read next to one of its candidates with probability candidate_fraction (decision D-045).

## checksum (SHA-256)

**In plain words.** A fingerprint of a file. If even one letter changes, the fingerprint changes completely, so we can tell that a file is untouched.

**Precisely.** The 256-bit output of the SHA-256 hash function, written as 64 hexadecimal digits. The same bytes always give the same value.

## contig

**In plain words.** A stretch of the genome that the reads cover without a gap, so it can be put together in one piece.

**Precisely.** A true contig is a maximal group of reads chained by true overlaps of at least theta bases, counted from the known read positions. An assembler's contigs are judged against these.

## coverage and depth

**In plain words.** How many layers of strips lie on top of each spot of the long paper, on average.

**Precisely.** With N reads of length l from a genome of length G, the nominal coverage is c = N * l / G. The depth at a position is the number of reads that cover it; its average over the genome is exactly c when every read lies inside the genome.

## data card

**In plain words.** A one-page fact sheet about a data set, like the label on a food package.

**Precisely.** A JSON and Markdown summary generated from the cached files: accession, version, organism, length, base counts, GC content and checksums.

## delta evaluation

**In plain words.** After a small change we do not recount everything; we only look at the few joins that changed.

**Precisely.** For a swap or insert move the change of F is the sum of the new edge weights minus the old ones, computed in constant time and tested against a full recomputation.

## DNA

**In plain words.** The instruction book of a living thing, written as a very long chain of four kinds of building blocks.

**Precisely.** A molecule made of two strands of nucleotides. Each nucleotide carries one of four bases: adenine (A), cytosine (C), guanine (G) or thymine (T).

## dummy node

**In plain words.** An extra imaginary read, a start-and-finish point, that turns the chain of reads into a closed loop.

**Precisely.** Node 0 with c(0, j) = len_j and c(i, 0) = 0. It is what makes the cycle cost equal to L.

## edge effect

**In plain words.** Near the two ends of the genome fewer strips can cover a spot, because every strip has to fit inside the genome.

**Precisely.** For a position x closer than l - 1 to an end, the probability q(x) that a random read covers it is smaller than l / (G - l + 1), so uncovered positions are more likely than exp(-c) says.

## effect size (A12, Cliff's delta)

**In plain words.** How big a difference is, not only whether it exists: the chance that a random result of method A beats a random result of method B.

**Precisely.** Vargha-Delaney A12 = P(A > B) + 0.5 P(A = B) (0.5 means no effect); Cliff's delta = 2 * A12 - 1.

## Entrez

**In plain words.** The website service that lets a program ask the NCBI library for a sequence.

**Precisely.** The NCBI search and retrieval system. Programs use its E-utilities (for example esummary and efetch) through Biopython's Bio.Entrez module.

## evaluation (budget)

**In plain words.** One look at how good a candidate change would be. A fair race gives every solver the same time and also reports how many looks it used.

**Precisely.** One delta of a move, one insertion-gain computation, or one full recomputation. Wall-clock time is the main budget; evaluations are reported as well.

## FASTA

**In plain words.** A very simple text format for DNA: one line that starts with ">" and names the sequence, followed by lines of letters.

**Precisely.** A text format with a header line beginning with ">" followed by one or more lines of sequence. Parsing rules used here are in genomeasm/io.py and docs/DATA_LOADING.md.

## GC content

**In plain words.** Out of all the letters, the share that are G or C.

**Precisely.** (number of G + number of C) / sequence length, a number between 0 and 1, usually shown as a percentage.

## genome

**In plain words.** The complete text of the instruction book of one organism or virus.

**Precisely.** The full DNA sequence of an organism or virus. Here, one string over {A, C, G, T} of length G, taken from an NCBI record.

## genome fraction, identity, N50 and misjoin

**In plain words.** How much of the real genome the assembly rebuilt, how correct its letters are, how long the typical piece is, and whether a piece jumps between distant places.

**Precisely.** Defined in genomeasm/experiments/bio.py: windows of 100 letters are located on the reference by shared 15-letter words; genome fraction counts positions covered by windows with identity >= 0.9; N50 is the contig length such that longer contigs hold half of the assembled letters.

## greedy and nearest neighbour

**In plain words.** Greedy: glue the best-matching pair first, then the next best, and so on. Nearest neighbour: start somewhere and always step to the best-matching unused read.

**Precisely.** Greedy edge merging accepts overlaps in decreasing order when the first read has no successor yet, the second no predecessor, and no cycle is closed. Both are deterministic.

## ground truth

**In plain words.** The answer key. Because we made the reads ourselves, we know exactly where each one came from and which mistakes were put in.

**Precisely.** For every simulated read, its true start position, the number of reference bases it covers and its counts of substitutions, insertions and deletions.

## Hamming distance and error-tolerant overlap

**In plain words.** The Hamming distance counts the places where two equally long texts differ. An error-tolerant overlap accepts a matching edge that has a few typing mistakes in it.

**Precisely.** The largest k >= tau whose suffix and prefix differ in at most floor(tolerance * k) places; its score is k itself. Insertions and deletions are not handled.

## heuristic and exact solver

**In plain words.** An exact solver is guaranteed to find the very best order but only for a handful of reads. A heuristic is a clever shortcut that finds a very good order quickly, without a guarantee.

**Precisely.** Brute force (n <= 9) and Held-Karp dynamic programming (n <= 15) are exact. Greedy, tabu search, ALNS and LKH-3 are heuristics; the exact solvers are used to check them on small instances.

## Holm correction

**In plain words.** When you run several tests, some will look significant by luck. This correction makes the bar higher so that the whole group is judged fairly.

**Precisely.** The step-down Holm-Bonferroni adjustment of p-values, applied within each family of comparisons.

## hyper-parameter tuning (random search)

**In plain words.** Every method has dials. We try many random dial settings on practice cases and keep the best, then use that setting on the real test cases, which the tuning never saw.

**Precisely.** 24 configurations per algorithm (trial 0 is the default), evaluated on 18 tuning read sets from a region of phiX174 that the evaluation does not use; the configuration with the smallest mean regret to the best known is kept.

## indel (insertion and deletion)

**In plain words.** A typing mistake where a letter is skipped (deletion) or an extra letter slips in (insertion).

**Precisely.** Deletion: a reference base is dropped with probability p_del. Insertion: after each reference base, one random base is inserted with probability p_ins.

## Lander-Waterman model

**In plain words.** A classic piece of mathematics that predicts how many spots stay bare, and how many separate pieces appear, when random strips are thrown onto a long paper.

**Precisely.** Reads start at uniformly random positions. It predicts depth ~ Poisson(c), uncovered fraction ~ exp(-c) and about N * exp(-c * (1 - theta / l)) contigs, ignoring the genome ends. Original source: Lander and Waterman (1988); the citation is TODO-VERIFY until it is checked.

## layout and consensus

**In plain words.** Layout: slide each read along so it sits on top of the one before it. Consensus: where several reads lie on top of each other, let them vote on each letter and keep the most common one.

**Precisely.** offset(pi_(k+1)) = offset(pi_k) + len(pi_k) - ov(pi_k, pi_(k+1)); a contig ends where an overlap is 0. Each column takes the most frequent letter, ties going to the first of A, C, G, T.

## LKH-3

**In plain words.** A famous program written by a researcher for travelling-salesman problems. We use it as the strong reference to compare our solvers with.

**Precisely.** Helsgaun's Lin-Kernighan-Helsgaun heuristic, version 3.0.6 here, run through Ubuntu on a TSPLIB file; every tour it returns is re-checked against the cost matrix.

## minimum overlap (tau)

**In plain words.** The shortest matching edge we trust. Very short matches happen by luck between reads that have nothing to do with each other.

**Precisely.** The parameter tau in the overlap definition. The project picks the smallest tau for which fewer than one chance overlap is expected among n(n-1) ordered pairs in a random-sequence model, namely below n(n-1) * (4/3) * 4^(-tau).

## Monte Carlo simulation

**In plain words.** Finding out how something behaves by running it many times with random dice and looking at all the results.

**Precisely.** Estimating an expected value by averaging over R independent random runs. The uncertainty of the average is its standard error.

## multi-armed bandit and policy

**In plain words.** Several slot machines with unknown payouts: pull one, see the reward, and decide which to pull next. Here each machine is a pair of destroy and repair tools. The rule that chooses is the policy.

**Precisely.** Policies: random, round robin, roulette wheel (classical ALNS weights), epsilon-greedy and UCB1 (both on a sliding window of the last 100 rewards, because usefulness changes over time).

## NCBI

**In plain words.** A large public library in the United States that stores DNA sequences for scientists.

**Precisely.** The National Center for Biotechnology Information. It hosts GenBank, RefSeq and the Entrez search and download system.

## objective (fitness) function F and merged length L

**In plain words.** The score of an order: how many letters we save by gluing neighbouring reads. Bigger is better. Equivalently, how short the glued text becomes.

**Precisely.** F(pi) = sum of ov(pi_k, pi_(k+1)) over k = 1..n-1, to be maximised. L(pi) = (sum of all read lengths) - F(pi), the length of the merged sequence, to be minimised.

## overlap (suffix-prefix)

**In plain words.** How far the end of one read matches the beginning of another, like the matching edge of two jigsaw pieces. It depends on the direction: A then B is not the same as B then A.

**Precisely.** ov(i, j) is the largest k with tau <= k < min(len_i, len_j) such that the last k letters of read i equal the first k letters of read j; 0 if there is none.

## parser

**In plain words.** A small program that reads a file and checks that it follows the rules of its format.

**Precisely.** A function that turns text into structured data and rejects text that violates the format. genomeasm.io.parse_fasta is the FASTA parser.

## permutation (solution)

**In plain words.** A list that puts every read in some order, each read used exactly once.

**Precisely.** A bijection pi of the n reads. The search space has n! elements; a solution is valid exactly when it is a permutation.

## plateau

**In plain words.** A stretch of the search in which nothing gets better any more.

**Precisely.** A phase of an ALNS run in which no move improves the current solution, so every reward is zero and every arm looks equally good to a bandit.

## Poisson distribution

**In plain words.** The pattern of counts you get when things land at random, like how many raindrops hit each paving stone.

**Precisely.** P(k) = exp(-c) * c^k / k! for k = 0, 1, 2, ... Lander-Waterman approximates the depth at a position by a Poisson distribution with mean c.

## post-hoc variant

**In plain words.** A version of a method that was added after seeing the first results. Its results are reported, but they count as exploration and never as a confirmation.

**Precisely.** An algorithm or analysis not fixed in the configuration before the first main run: tabu_cl (D-045) and alns_ucb1_rt (D-046). Its tests form a separate family and are not merged into the pre-registered one.

## provenance

**In plain words.** The label on a result that says exactly how and when it was made.

**Precisely.** The record stored with every result: configuration hash, git commit and whether the code was modified, random seed, library versions, machine information and a timestamp.

## random seed

**In plain words.** A starting number that makes a computer's dice rolls repeat exactly, so a random experiment can be run again and give identical results.

**Precisely.** An integer given to numpy's SeedSequence. The same seed and configuration always give the same output; independent streams are spawned from it for read positions and for errors.

## read

**In plain words.** A short piece of DNA text that the sequencing machine reads in one go, like one strip torn from a very long roll of paper.

**Precisely.** A string over {A, C, G, T} obtained from one position of the genome, possibly with errors. Here the typical length is l = 100.

## read set and paired design

**In plain words.** One simulated batch of reads is one test case. Every solver gets exactly the same batch, so we compare them fairly, case by case.

**Precisely.** A read set is an instance setting plus a seed. Solvers are compared paired on the same read set; the read set is the unit of replication (30 per instance setting).

## RefSeq

**In plain words.** A shelf in the NCBI library with carefully checked reference copies.

**Precisely.** NCBI's curated collection of reference sequences. Accessions starting with NC_ are RefSeq records of complete molecules.

## relative gap to the best known

**In plain words.** How far a solver's score is below the best score anybody reached on the same test case, in percent. Zero means as good as the best.

**Precisely.** gap = 100 * (best - F) / best, where best is the largest fitness any solver reached on that read set. The gap to LKH-3 uses the LKH-3 fitness instead.

## simulated annealing

**In plain words.** Accept a worse solution sometimes, often at the beginning and almost never at the end, like a metal cooling slowly.

**Precisely.** A worse solution (change delta < 0) is accepted with probability exp(delta / T); T falls geometrically from T0 to T0 * 1e-3 as the budget is used.

## standard error

**In plain words.** How much an average would wobble if the whole experiment were repeated. Smaller means more trust in the average.

**Precisely.** The sample standard deviation across R independent runs divided by the square root of R. Checks in this project pass when the observed value is within a stated number of standard errors of the expected value.

## substitution error

**In plain words.** A typing mistake where one letter is swapped for a different letter.

**Precisely.** A kept base replaced by one of the three other bases, each with probability 1/3. The substitution rate is the probability of this for each kept base.

## substring-free reads

**In plain words.** A set of reads where no read is a copy of another or sits completely inside another, because such reads add no new letters.

**Precisely.** Duplicates are removed (the first copy is kept) and every read that is a substring of a longer read is removed. The numbers removed are always recorded.

## tabu search

**In plain words.** Keep taking the best small change, even when it makes things worse, but remember what you just moved and forbid moving it again for a while, so you do not walk in circles.

**Precisely.** Local search over swap and insert moves with a tabu list (tenure), an aspiration criterion (a tabu move is allowed if it gives a new best) and restarts from the shaken best solution.

## tie-breaking

**In plain words.** What a chooser does when two options look exactly equally good: always take the first one, or flip a coin.

**Precisely.** UCB1 as first implemented took the lowest arm index among equal indices. With all rewards equal it then played arm 0 for ever (a defect, decision D-046); ucb1_rt breaks ties uniformly at random.

## travelling salesman problem (TSP, ATSP)

**In plain words.** Find the cheapest round trip through a list of towns that visits each town once. In the asymmetric version, the price from A to B can differ from the price from B to A.

**Precisely.** Find a Hamiltonian cycle of minimum total cost in a complete directed graph. Here town j after town i costs len_j - ov(i, j), plus a dummy node so that the cycle through all reads costs exactly L(pi).

## true ordering

**In plain words.** The order the reads really have along the genome, which we know because we simulated them.

**Precisely.** The reads sorted by their true start position. Its fitness F is the reference line against which the solvers are compared.

## TSPLIB file

**In plain words.** A standard text file format for travelling-salesman problems that many solvers, such as LKH-3, can read.

**Precisely.** A file with TYPE: ATSP, EDGE_WEIGHT_TYPE: EXPLICIT, EDGE_WEIGHT_FORMAT: FULL_MATRIX and the integer cost matrix in the EDGE_WEIGHT_SECTION.

## UCB1 (upper confidence bound)

**In plain words.** Pick the tool whose average result plus a bonus for being tried rarely is highest: optimism in the face of uncertainty.

**Precisely.** Index = mean reward + sqrt(2 ln N / n), with n the arm's rewards in its window and N the total. Citation Auer, Cesa-Bianchi and Fischer (2002) is TODO-VERIFY.

## uncovered fraction

**In plain words.** The share of spots on the long paper that no strip happened to cover.

**Precisely.** The fraction of genome positions with depth 0. Lander-Waterman predicts exp(-c); the exact finite-genome model of this project also accounts for the genome ends.

## Wilcoxon signed-rank test

**In plain words.** A test that asks, for pairs of results on the same case, whether one method wins more or by more than luck would allow.

**Precisely.** A non-parametric paired test on the differences, two-sided here; zero differences are discarded.

## winner's curse

**In plain words.** If you pick the best of many noisy trials, the best one usually looks better than it really is.

**Precisely.** The selection bias of choosing the configuration with the smallest mean regret among 24 random trials on a small tuning set: its tuning regret is optimistic. All algorithms are tuned with the same number of trials, so the bias is similar for all.
