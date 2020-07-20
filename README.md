# Experiments with Permutation Group Algorithms

I have several planned projects that require me to deal with permutation
groups. Computing with permutation groups often requires building a strong
generating set. This can be done using the Schreier-Sims algorithm. As there
are a lot of different variants and trade-offs when implementing the
Schreier-Sims algorithm, I decided to write a small Python prototype that
allows me to quickly explore these to get a better understanding of the
variants described in the literature.

Optimized implementations of this will hopefully appear in a Rust crate in the
future.

My goal when documenting this was to make sure that it's not too hard to follow
the implementation together with the existing literature. The documentation
isn't a self contained explanation of the algorithms though.

There is also a commented example script in `example_rubiks.py`. Its output is
also at the end of this document.

## Disclaimer

This was written to explore different variants and trade-offs and gather some
performance statistics about them. It is not optimized for runtime performance
at all though! For larger groups, computing group products will dominate the
runtime, so I could simply count these, without having to optimize anything to
benchmark different variants.

I do not recommend using this code for anything but exploring implementation
details of these algorithms. Well tuned implementations of permutation group
algorithms can be found in the open source GAP computer algebra system
(gap-system.org). SageMath also contains a Python interface to GAP.

## Example Output

This is a transcript of running `example_rubiks.py`:

We're going to build a strong generating set for the Rubik's Cube group,
including reorientations of the whole cube. This is a permutation group of the
9 * 6 = 54 visible cublet faces.

It is generated using the operations X, Y which allow us to arbitrarily
reorient the whole cube without turning any sides.

    X = (0 9 45 35)(1 10 46 34)(2 11 47 33)(3 12 48 32)(4 13 49 31)(5 14 50 30)(6 15 51 29)(7 16 52 28)(8 17 53 27)(18 24 26 20)(19 21 25 23)(36 38 44 42)(37 41 43 39)
    Y = (0 20 53 42)(1 23 52 39)(2 26 51 36)(3 19 50 43)(4 22 49 40)(5 25 48 37)(6 18 47 44)(7 21 46 41)(8 24 45 38)(9 11 17 15)(10 14 16 12)(27 33 35 29)(28 30 34 32)

And the operation R, which turns a single side. By reorienting the cube we can
use this to turn any side and thus get all possible ways to permute the cube.

    R = (2 11 47 33)(5 14 50 30)(8 17 53 27)(18 24 26 20)(19 21 25 23)

------------------------------------------------------------------------------
First we're going to use the deterministic Schreier-Sims algorithm without
limiting the depth of Schreier trees.

    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 6, 5, 8, 7, 15, 12, 17, 14, 16, 1, 23, 25, 4, 34, 32, 26, 13, 3]
    strong generating set size = 47
    took 238,782 group products

------------------------------------------------------------------------------
Next we're going to use shallow Schreier trees, which require more work to
build, but ensure that sifting can be done efficiently. That usually pays off.

    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 6, 5, 8, 12, 7, 14, 16, 15, 17, 1, 23, 25, 4, 34, 26, 32, 13, 3]
    strong generating set size = 43
    took 202,928 group products

------------------------------------------------------------------------------
Using the Monte Carlo Random Schreier-Sims algorithm is a lot more efficient.

    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 1, 3, 4, 5, 6, 8, 7, 12, 15, 13, 14, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 27
    took 1,077 group products
    took 35 sifting rounds

------------------------------------------------------------------------------
It doesn't guarantee that the result is correct though, as the algorithm
terminates after successfully sifting k random elements. If we could generate
uniform samples of our group, the chance of failure would be limited by 2^-k.
As we're using a heuristic to generate random elements, we do not even have
that guarantee. In practice though, the chance of failure is often smaller.

To demonstrate failure, we set k (called exit_rounds in the code) to a small
value like 1 and repeatedly build a strong generating set until we get the
wrong group order.

    group order = 173,008,013,097,959,424,000
    strong generating set base = [0, 2, 1, 3, 5, 4, 6, 8, 7, 12, 13, 14, 15, 17, 23, 16, 25, 32, 26]
    strong generating set size = 34
    took 789 group products
    took 32 sifting rounds

------------------------------------------------------------------------------
If we know the order of the group though, we can turn it into a Las Vegas
algorithm. That is even more efficient, as it terminates as soon as the strong
generating set is complete and thus needs fewer sifting rounds.

    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 1, 3, 4, 5, 6, 7, 8, 12, 13, 14, 16, 15, 17, 25, 23, 26, 32, 34]
    strong generating set size = 36
    took 859 group products
    took 34 sifting rounds

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

    verification failures = 0
    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 1, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 23, 26, 25, 32, 34]
    strong generating set size = 31
    took 451,849 group products
    took 40 sifting rounds

------------------------------------------------------------------------------
Luckily there are more efficient tests for a strong generating set. One that is
often used in practice is Sims's Verify routine. It's quite a bit more
complicated than the Random Schreier-Sims algorithm and requires among other
things base changes (which we'll cover next). It's also quite a bit more
expensive than just the Random Schreier-Sims algorithm, but often still a lot
better than the deterministic Schreier-Sims algorithm.

    verification failures = 0
    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 1, 3, 4, 5, 8, 6, 7, 12, 13, 14, 15, 16, 23, 17, 25, 32, 26, 34]
    strong generating set size = 35
    took 46,475 group products
    took 43 sifting rounds

------------------------------------------------------------------------------
To illustrate that the verification does indeed work, we can again set
exit_rounds to a small value of 1.

    verification failures = 1
    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [0, 2, 1, 3, 4, 5, 6, 8, 7, 12, 13, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 32
    took 41,279 group products
    took 30 sifting rounds

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

    base perfix = [4, 13, 22, 31, 40, 49]
    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [4, 13, 22, 31, 40, 49, 2, 0, 1, 3, 5, 6, 7, 8, 12, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 31

subgroup stabilizing the cube orientation:

    group order = 43,252,003,274,489,856,000
    strong generating set base = [2, 0, 1, 3, 5, 6, 7, 8, 12, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 28
    took 204 group products
    took 28 sifting rounds

------------------------------------------------------------------------------
If we already generated a strong generating set but with a different base, it's
often possible to perform a base change faster than generating a new strong
generating set from scratch.

This uses a heuristic combination of conjugating the strong generating set and
partially rebuilding the strong generating set using the Las Vegas (as we know
the group orders) Random Schreier-Sims algorithm.

    old base = [0, 2, 1, 3, 4, 5, 6, 8, 7, 12, 13, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    group order = 1,038,048,078,587,756,544,000
    strong generating set base = [4, 13, 22, 31, 40, 49, 0, 2, 1, 3, 5, 6, 8, 7, 12, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 26

subgroup stabilizing the cube orientation:

    group order = 43,252,003,274,489,856,000
    strong generating set base = [0, 2, 1, 3, 5, 6, 8, 7, 12, 14, 15, 16, 17, 23, 25, 26, 32, 34]
    strong generating set size = 24
    took 60 group products
    took 16 sifting rounds
