from schreiersims import *
from random import Random
import re

# Three permutations that generate the Rubik's Cube.
# Shamelessly stolen from
# https://github.com/runjak/2020-06-06.enthusiasticon/blob/master/src/permutations.ts#L58

X = [
   9, 10, 11, 12, 13, 14, 15, 16, 17,
  45, 46, 47, 48, 49, 50, 51, 52, 53,
  24, 21, 18, 25, 22, 19, 26, 23, 20,
   8,  7,  6,  5,  4,  3,  2,  1,  0,
  38, 41, 44, 37, 40, 43, 36, 39, 42,
  35, 34, 33, 32, 31, 30, 29, 28, 27,
]
Y = [
  20, 23, 26, 19, 22, 25, 18, 21, 24,
  11, 14, 17, 10, 13, 16,  9, 12, 15,
  47, 50, 53, 46, 49, 52, 45, 48, 51,
  33, 30, 27, 34, 31, 28, 35, 32, 29,
   2,  5,  8,  1,  4,  7,  0,  3,  6,
  38, 41, 44, 37, 40, 43, 36, 39, 42,
]
R = [
   0,  1, 11,  3,  4, 14,  6,  7, 17,
   9, 10, 47, 12, 13, 50, 15, 16, 53,
  24, 21, 18, 25, 22, 19, 26, 23, 20,
   8, 28, 29,  5, 31, 32,  2, 34, 35,
  36, 37, 38, 39, 40, 41, 42, 43, 44,
  45, 46, 33, 48, 49, 30, 51, 52, 27,
]

center_cubelet_faces = list(range(4, 54, 9))


def fmt_large_num(x):
    return re.subn(r'(?<=\d)(?=(\d{3})+$)', ',', str(x))[0]


def add_rubiks_gens(grp):
    for gen in [X, Y, R]:
        grp.add_gen(gen)


def print_all_stats(grp):
    print_sgs_stats(grp)
    print_performance_stats(grp)


def print_sgs_stats(grp):
    print(f"  group order = {fmt_large_num(grp.order())}")
    print(f"  strong generating set base = {grp.base()}")
    print(f"  strong generating set size = {len(grp.generators())}")


def print_performance_stats(grp):
    print(f"  took {fmt_large_num(grp.cfg.stats.products)} group products")
    if grp.cfg.stats.rounds:
        print(f"  took {grp.cfg.stats.rounds} sifting rounds")


print("""
We're going to build a strong generating set for the Rubik's Cube group,
including reorientations of the whole cube. This is a permutation group of the
9 * 6 = 54 visible cublet faces.

It is generated using the operations X, Y which allow us to arbitrarily
reorient the whole cube without turning any sides.
""")

print('X =', fmt_perm(X))
print('Y =', fmt_perm(Y))

print("""
And the operation R, which turns a single side. By reorienting the cube we can
use this to turn any side and thus get all possible ways to permute the cube.
""")
print('R =', fmt_perm(R))

print("""
------------------------------------------------------------------------------
First we're going to use the deterministic Schreier-Sims algorithm without
limiting the depth of Schreier trees.
""")

grp = Group(Config(monte_carlo=False, schreier_tree='deep'))
add_rubiks_gens(grp)
print_all_stats(grp)

known_group_order = grp.order()

print("""
------------------------------------------------------------------------------
Next we're going to use shallow Schreier trees, which require more work to
build, but ensure that sifting can be done efficiently. That usually pays off.
""")

grp = Group(Config(monte_carlo=False))
add_rubiks_gens(grp)
print_all_stats(grp)

print("""
------------------------------------------------------------------------------
Using the Monte Carlo Random Schreier-Sims algorithm is a lot more efficient.
""")

grp = Group()
add_rubiks_gens(grp)
grp.build()
print_all_stats(grp)

print("""
------------------------------------------------------------------------------
It doesn't guarantee that the result is correct though, as the algorithm
terminates after successfully sifting k random elements. If we could generate
uniform samples of our group, the chance of failure would be limited by 2^-k.
As we're using a heuristic to generate random elements, we do not even have
that guarantee. In practice though, the chance of failure is often smaller.

To demonstrate failure, we set k (called exit_rounds in the code) to a small
value like 1 and repeatedly build a strong generating set until we get the
wrong group order.
""")

while True:
    # The default config fixes the random seed, so we need to set rng
    grp = Group(Config(exit_rounds=1))
    add_rubiks_gens(grp)
    grp.build()
    if grp.order() != known_group_order:
        break

print_all_stats(grp)

print("""
------------------------------------------------------------------------------
If we know the order of the group though, we can turn it into a Las Vegas
algorithm. That is even more efficient, as it terminates as soon as the strong
generating set is complete and thus needs fewer sifting rounds.
""")

grp = Group()
add_rubiks_gens(grp)
grp.build(known_group_order)
print_all_stats(grp)

print("""
------------------------------------------------------------------------------
If we don't know the order of the group we can turn the Monte Carlo algorithm
into a Las Vegas algorithm by performing a verification step. If the
verification fails, it will provide us with a group element that doesn't sift
through the incomplete strong generating set. In that case we add the sifting
residue and perform some more Random Schreier-Sims rounds.

The simplest way to verify a strong generating set is to generate all Schreier
generators and sift them. It's not very efficient though, as that's basically
what the deterministic Schreier-Sims algorithm does, so we lose the advantage
of the faster Random Schreier-Sims approach.
""")

grp = Group(Config(verify='schreier'))
add_rubiks_gens(grp)
failures = grp.build_verified()

print(f'  verification failures = {failures}')
print_all_stats(grp)

print("""
------------------------------------------------------------------------------
Luckily there are more efficient tests for a strong generating set. One that is
often used in practice is Sims's Verify routine. It's quite a bit more
complicated than the Random Schreier-Sims algorithm and requires among other
things base changes (which we'll cover next). It's also quite a bit more
expensive than just the Random Schreier-Sims algorithm, but often still a lot
better than the deterministic Schreier-Sims algorithm.
""")

grp = Group()
add_rubiks_gens(grp)
failures = grp.build_verified()

print(f'  verification failures = {failures}')
print_all_stats(grp)

print("""
------------------------------------------------------------------------------
To illustrate that the verification does indeed work, we can again set
exit_rounds to a small value of 1.
""")

while True:
    # The default config fixes the random seed, so we need to set rng
    grp = Group(Config(exit_rounds=1))
    add_rubiks_gens(grp)
    failures = grp.build_verified()
    if failures:
        break

print(f'  verification failures = {failures}')
print_all_stats(grp)

saved_grp = grp

print("""
------------------------------------------------------------------------------
If you happen to be a Rubik's Cube enthusiast, you might already have noticed,
that the printed group order is larger than the roughly 43 quintillion that is
often quoted.

This is because we're also counting different orientations of the cube. If we
don't want to count them, we can stabilize the center cubelet faces of each
cube face, thereby fixing them in place.

One way to do this is to specify the center cublet faces as a prefix of the
base before running the Schreier-Sims algorithm and then take the first
subgroup in the stabilizer chain that fixes these points.

Here we're using the known order of the group including cube reorientations to
get the faster Las Vegas variant.
""")

grp = Group(base=center_cubelet_faces)
print(f'  base perfix = {grp.base()}')
add_rubiks_gens(grp)
grp.build(known_group_order)

print_sgs_stats(grp)

print()
print(f'subgroup stabilizing the cube orientation:')
print_sgs_stats(grp.stabilizer_chain()[len(center_cubelet_faces)])

print_performance_stats(grp)

print("""
------------------------------------------------------------------------------
If we already generated a strong generating set but with a different base, it's
often possible to perform a base change faster than generating a new strong
generating set from scratch.

This uses a heuristic combination of conjugating the strong generating set and
partially rebuilding the strong generating set using the Las Vegas (as we know
the group orders) Random Schreier-Sims algorithm.
""")

grp = saved_grp
grp.cfg.stats = Stats()  # Reset stats
print(f'  old base = {grp.base()}')
grp.change_base(center_cubelet_faces)
print_sgs_stats(grp)
print()
print(f'subgroup stabilizing the cube orientation:')
print_sgs_stats(grp.stabilizer_chain()[len(center_cubelet_faces)])
print_performance_stats(grp)
