import numpy as np

from Bio.Seq import Seq, MutableSeq
from Bio.SeqIO import SeqRecord, write

from tqdm import trange


DNA_ALPHABET = ('A', 'T', 'G', 'C')

rng = np.random.default_rng(2024)


def mutate_seq(seq, prob, rng):
    seq_len = len(seq)
    seq_list = list(seq)
    mut_poses = rng.binomial(1, prob, size=seq_len)
    for j in range(seq_len):
        if mut_poses[j]:
            seq_list[j] = rng.choice(DNA_ALPHABET)

    return ''.join(seq_list)


def write_test(genome, reads, test_id):
    genome = SeqRecord(
        genome,
        id=f'genome_{test_id}',
        description=f'{len(genome)}bp'
    )
    reads = [
        SeqRecord(seq, id=str(i), description='')
        for i, seq in enumerate(reads)
    ]

    with open(f'genome_{test_id}.fasta', 'w') as f_out:
        write(genome, f_out, 'fasta')

    with open(f'reads_{test_id}.fasta', 'w') as f_out:
        write(reads, f_out, 'fasta')


# Genome length: 12bp
# Reads: all 3-mers exactly once
# No errors

genome0 = Seq(''.join(rng.choice(DNA_ALPHABET, 12)))
reads0 = [genome0[i:i + 3] for i in range(len(genome0) - 2)]
rng.shuffle(reads0)
write_test(genome0, reads0, 0)


# Genome length: 300bp
# Reads: Random position and length (20 - 40bp)
# Coverage: 10
# No errors

genome1 = Seq(''.join(rng.choice(DNA_ALPHABET, 300)))
reads1 = []
for i in range(100):
    read_len = rng.integers(20, 40, endpoint=True)
    start = rng.integers(301 - read_len)
    reads1.append(genome1[start:start + read_len])

write_test(genome1, reads1, 1)


# Genome length: 1000bp
# Reads: Random position and length (40 - 60bp)
# Coverage: 20
# Errors: 0.1%

genome2 = Seq(''.join(rng.choice(DNA_ALPHABET, 1000)))
reads2 = []

for i in range(400):
    read_len = rng.integers(40, 60, endpoint=True)
    start = rng.integers(1001 - read_len)
    read = MutableSeq(genome2[start:start + read_len])

    error_poses = rng.binomial(1, 0.001, size=read_len)
    for j in range(read_len):
        if error_poses[j]:
            read[j] = rng.choice(DNA_ALPHABET)

    reads2.append(read)

write_test(genome2, reads2, 2)


# Genome length: 100_000bp
# Reads: Random position and length (150 - 250bp)
# Coverage: 40
# Errors: 0.1%

genome3 = Seq(''.join(rng.choice(DNA_ALPHABET, 100_000)))
reads3 = []

for i in range(20_000):
    read_len = rng.integers(150, 250, endpoint=True)
    start = rng.integers(100_001 - read_len)
    read = MutableSeq(genome3[start:start + read_len])

    error_poses = rng.binomial(1, 0.001, size=read_len)
    for j in range(read_len):
        if error_poses[j]:
            read[j] = rng.choice(DNA_ALPHABET)

    reads3.append(read)

write_test(genome3, reads3, 3)


# Genome: 1 chromosome 1Mbp with different number and lengths of
#                           repeating regions
#         2 plasmids 10Kbp  with some common regions and genes
#                           in one and two copies respectively
#                           0.1% mutations between copies of same plasmid
#         2 plasmids 2Kbp   in one and three copies respectively
#                           0.05% mutations between copies of same plasmid
# Reads: Random position
#        length 200bp
# Coverage:
# Errors: 0.05%

genome4 = []

chromosome = ''.join(rng.choice(DNA_ALPHABET, 1_000_000))
for rep_length, rep_num in [
    (5, 200),
    (600, 3),
    (20000, 2),
]:
    rep_seq = ''.join(rng.choice(DNA_ALPHABET, rep_length))
    for _ in range(rep_num):
        rep_start = rng.integers(1_000_001 - rep_length)
        chromosome = chromosome[:rep_start] + rep_seq + chromosome[rep_start + rep_length:]

genome4.append(Seq(chromosome))

plasmid1 = ''.join(rng.choice(DNA_ALPHABET, 10_000))
plasmid2 = ''.join(rng.choice(DNA_ALPHABET, 10_000))

for rep_length, rep_num1, rep_num2 in [
    (20, 5, 15),
    (500, 1, 1),
    (2000, 1, 1),
]:
    rep_seq = ''.join(rng.choice(DNA_ALPHABET, rep_length))
    for _ in range(rep_num1):
        rep_start = rng.integers(10_001 - rep_length)
        plasmid1 = plasmid1[:rep_start] + rep_seq + plasmid1[rep_start + rep_length:]
    for _ in range(rep_num2):
        rep_start = rng.integers(10_001 - rep_length)
        plasmid2 = plasmid2[:rep_start] + rep_seq + plasmid2[rep_start + rep_length:]

genome4.extend((
    Seq(plasmid1),
    Seq(plasmid2),
    Seq(mutate_seq(plasmid2, 0.001, rng)),
))

plasmid3 = ''.join(rng.choice(DNA_ALPHABET, 2_000))
plasmid4 = ''.join(rng.choice(DNA_ALPHABET, 2_000))

genome4.extend((
    Seq(plasmid3),
    Seq(plasmid4),
    Seq(mutate_seq(plasmid4, 0.0005, rng)),
    Seq(mutate_seq(plasmid4, 0.0005, rng)),
))

genome4_lengths = np.array([len(seq) for seq in genome4])
choice_probs = genome4_lengths / genome4_lengths.sum()

reads4 = []

for i in trange(200_000):
    seq_n = rng.choice(len(genome4), p=choice_probs)
    read_source = genome4[seq_n]
    read_len = 200
    if seq_n == 0:
        start = rng.integers(genome4_lengths[seq_n] - read_len + 1)
        read = MutableSeq(read_source[start:start + read_len])
    else:
        start = rng.integers(genome4_lengths[seq_n])
        read_source = str(read_source) * 2 # plasmids are circular
        read = MutableSeq(read_source[start:start + read_len])

    error_poses = rng.binomial(1, 0.0005, size=read_len)
    for j in range(read_len):
        if error_poses[j]:
            read[j] = rng.choice(DNA_ALPHABET)

    reads4.append(read)

genome4_ids = [
    'chromosome',
    'plasmid1',
    'plasmid2_1',
    'plasmid2_2',
    'plasmid3',
    'plasmid4_1',
    'plasmid4_2',
    'plasmid4_3',
]

genome4 = [
    SeqRecord(seq, id=id_, description=f'{seq_len}bp')
    for seq, id_, seq_len in zip(genome4, genome4_ids, genome4_lengths)
]
reads4 = [
    SeqRecord(seq, id=str(i), description='')
    for i, seq in enumerate(reads4)
]

with open(f'genome_4.fasta', 'w') as f_out:
    write(genome4, f_out, 'fasta')

with open(f'reads_4.fasta', 'w') as f_out:
    write(reads4, f_out, 'fasta')
