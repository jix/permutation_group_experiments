"""
This implements a few variants of the Schreier-Sims algorithm for generating
strong generating sets of permutation groups.

This was written to explore different variants and trade-offs and gather some
performance statistics about them. It is not optimized for runtime performance
at all though! For larger groups, computing group products will dominate the
runtime, so I could simply count these, without having to optimize anything to
benchmark different variants.

I do not recommend using this code for anything but exploring implementation
details of these algorithms. Well tuned implementations of permutation group
algorithms can be found in the open source GAP computer algebra system
(gap-system.org). SageMath also contains a Python interface to GAP.
"""


from dataclasses import dataclass, field
from typing import Literal
from random import Random
import copy
import re
import math
import sys

if sys.version_info < (3, 7):
    raise RuntimeError('python before version 3.7 is not supported')


def check_perm(p):
    """Raise an exception if p is not a permutation.
    """
    if set(p) != set(range(len(p))):
        raise ValueError('not a permutation')


def fmt_inv_pair_helper(p):
    """Handle permutation/inverse paris.

    For a pair of a permutation and its inverse, return just the permutation.
    Return other values unchanged.
    """
    if isinstance(p, tuple) and isinstance(p[0], (tuple, list)):
        check_perm(p[1])
        if not is_id_perm(mult_perm(p[0], p[1])):
            raise ValueError('expected a permutation or an inverse pair')
        p = p[0]
    return p


def fmt_perm(p):
    """Display a permutation as a product of cycles.
    """
    p = fmt_inv_pair_helper(p)
    check_perm(p)
    seen = set()
    out = []
    for i in range(len(p)):
        if i in seen:
            continue
        if p[i] == i:
            continue
        cycle = [i]
        j = p[i]
        while j != i:
            seen.add(j)
            cycle.append(j)
            j = p[j]
        out.append('(%s)' % ' '.join(map(str, cycle)))

    if not out:
        return '()'
    return ''.join(out)


def fmt_perm_list(ps):
    """Display a list of permutations, each as a product of cycles.
    """
    return ', '.join(map(fmt_perm, ps))


ELEMENT_SEP_RE = r' *[, ] *'
CYCLE_RE = rf'\(( *\d+({ELEMENT_SEP_RE}\d+)* *)?\) *'
PERM_LIST_TOKEN_RE = rf'({CYCLE_RE})+|, *'


def parse_perm(s, n=0):
    """Parse a permutation given as a product of cycles.
    """
    cycles = []
    stripped = re.subn(r'\s', ' ', s)[0].strip()
    for match in re.finditer(CYCLE_RE + r'|.', stripped):
        cycle = match.group().strip()
        if len(cycle) == 1:
            raise ValueError(f"could not parse permutation {s!r}")
        cycle = cycle[1:-1].strip()
        if not cycle:
            continue
        cycle = list(map(int, re.split(ELEMENT_SEP_RE, cycle)))
        if cycle:
            cycles.append(cycle)
    n = max(n, max(map(max, cycles), default=-1) + 1)
    out = id_perm(n)
    for cycle in cycles:
        out = mult_perm(out, cycle_perm(n, cycle))
    return out


def parse_perm_list(s, n=0):
    """Parse a list of permutations, each given as a product of cycles.
    """
    stripped = re.subn(r'\s', ' ', s)[0].strip()
    perms = []
    for match in re.finditer(PERM_LIST_TOKEN_RE + r'|.', stripped):
        token = match.group().strip()
        if token == ',':
            continue
        perm = parse_perm(token, n)
        n = max(n, len(perm))
        perms.append(perm)
    return [perm + list(range(len(perm), n)) for perm in perms]


def to_gap(p):
    """Shift permutation to act on 1..n for exporting to GAP.
    """
    return [0, *(1 + i for i in p)]


def from_gap(p):
    """Shift permutation fixing 0 to act on 0..n-1 for importing from GAP.
    """
    if p[0] != 0:
        raise ValueError("permutation doesn't fix 0")
    return [i - 1 for i in p[1:]]


def gap_fmt_perm(p):
    """Display a permutation as a product of cycles in GAP's syntax.
    """
    return fmt_perm(to_gap(fmt_inv_pair_helper(p))).replace(' ', ',')


def gap_fmt_perm_list(ps):
    """Display a permutation as a product of cycles in GAP's syntax.
    """
    return ', '.join(map(gap_fmt_perm, ps))


def gap_parse_perm(s, n=0):
    """Equivalent of parse_perm using GAP's syntax.
    """
    return from_gap(parse_perm(s, n + 1))


def gap_parse_perm_list(s, n=0):
    """Equivalent of parse_perm_list using GAP's syntax.
    """
    return list(map(from_gap, parse_perm_list(s, n + 1)))


def id_perm(n):
    """Return identity permutation on 0..n-1.
    """
    return list(range(n))


def is_id_perm(p):
    """Return whether p is the identity permutation.
    """
    return p is None or all(i == j for i, j in enumerate(p))


def cycle_perm(n, cycle):
    """Return a given cycle acting on 0..n-1.
    """
    p = id_perm(n)
    if not cycle:
        return p
    for i, j in zip(cycle, cycle[1:]):
        p[i] = j
    p[cycle[-1]] = cycle[0]
    return p


def mult_perm(p, q):
    """Multiply two permutations.

    If exactly one of the permutations is None, returns the other.
    """
    if p is None:
        return q
    elif q is None:
        return p
    return [q[i] for i in p]


def mult_perms(xs):
    """Multiply a list of permutations.
    """
    w = None
    for x in xs:
        w = mult_perm(w, x)
    return w


def inv_perm(p):
    """Inverse of a permutation.
    """
    if p is None:
        return None
    inv = [None] * len(p)
    for i, j in enumerate(p):
        inv[j] = i
    return inv


def exp_perm(p, n):
    """Take a permutation to the nth power.
    """
    if p is None:
        return None
    if n == 0:
        return id_perm(len(p))
    elif n < 0:
        return exp_perm(inv_perm(p), -n)

    q = exp_perm(p, n >> 1)
    if n & 1:
        q = mult_perm(p, q)

    return q


@dataclass
class Stats:
    products: int = 0
    rounds: int = 0


@dataclass
class Config:
    # Strategy to build shallow Schreier trees.
    #
    # 'gap':   Use the algorithm implemented by GAP and described in Seress's
    #          "Permutation Group Algorithms" 4.4.3.
    # 'halve': Variant of 'gap' which adds a generator for the deepest point of
    #          the orbit instead of the first exceeding the threshold.
    # 'deep':  Don't build shallow Schreier trees.
    schreier_tree: Literal['gap', 'halve', 'deep'] = 'gap'

    # Apply generators from both sides to sift generators further down the
    # stabilizer chain.
    deep_sift_gens: bool = False

    # Use the randomized Monte Carlo Schreier-Sims algorithm.
    monte_carlo: bool = True

    # Exit the Monte Carlo algorithm after the given number of rounds without
    # progress.
    exit_rounds: int = 10

    # For Monte Carlo, whenever a sift residue is added to some stabilizer
    # subgroup, also add this many random schreier generators.
    #
    # Not very effecitve as currently implemented.
    random_schreier_gens: int = 0

    # Which method should be used by the verify method to check whether a
    # strong generating set is complete.
    #
    # 'schreier': Sift all Schreier generators. Slow but simple.
    # 'sims':     Use Sims's verify routine described by Seress in "Permutation
    #             Group Algorithms" 8.2.
    verify: Literal['schreier', 'sims'] = 'sims'

    # Parameters for black box random element generation. See GAP's
    # documentation for ProductReplacer.
    rng_accus: int = 5
    rng_extra_slots: int = 5
    rng_scramble: int = 30
    rng_scramble_factor: int = 4

    stats: Stats = field(default_factory=lambda: Stats())
    rng: Random = field(default_factory=lambda: Random())


class NotInOrbit(ValueError):
    pass


class IncompleteStrongGeneratingSet(ValueError):
    def __init__(self, message, witness=None):
        super().__init__(message)
        self.witness = witness


class GroupRng:
    """Generate approximately uniform random elements of a group.

    This uses the "rattle" algorithm as described in GAP's documentation.
    """
    def __init__(self, cfg=None):
        self.cfg = cfg or Config()
        self.reservoir = []
        self.accus = []
        self.accu = 0

    def add_gen(self, gen):
        """Add a generator to the group.
        """
        if not self.reservoir:
            self.reservoir = [id_perm(len(gen))] * self.cfg.rng_extra_slots
            self.accus = [id_perm(len(gen))] * self.cfg.rng_accus
        self.reservoir.append(gen)
        self.new_gens = True

    def sample(self):
        """Sample a random element of the group.
        """
        if self.new_gens:
            self.new_gens = False
            self.scramble()

        return self.stir()

    def stir(self):
        """Perform a random replacement step.
        """

        i = self.cfg.rng.randrange(1, len(self.reservoir))
        j = self.cfg.rng.randrange(1, len(self.reservoir))

        p = self.reservoir[i]
        if self.cfg.rng.randrange(2):
            p = inv_perm(p)

        self.reservoir[0] = c = mult_perm(self.reservoir[0], p)
        self.cfg.stats.products += 1

        if self.cfg.rng.randrange(2):
            c = inv_perm(c)

        self.reservoir[j] = q = mult_perm(self.reservoir[j], c)
        self.cfg.stats.products += 1

        if self.cfg.rng.randrange(2):
            q = inv_perm(q)

        self.accu = (self.accu + 1) % len(self.accus)
        self.accus[self.accu] = r = mult_perm(self.accus[self.accu], q)
        self.cfg.stats.products += 1
        return r

    def scramble(self):
        """Perform a scramble of the reservoir"""

        gen_size = len(self.reservoir) - self.cfg.rng_extra_slots

        steps = max(
            self.cfg.rng_scramble, self.cfg.rng_scramble_factor * gen_size)

        for _ in range(steps):
            self.stir()

    def apply_map(self, f):
        """Apply a function to every permutation of the rng state.

        Mostly useful if f is a group morphism.
        """

        self.reservoir = [f(g) for g in self.reservoir]
        self.accus = [f(g) for g in self.accus]


class Group:
    """A permutation group.

    The group is represented using a strong generating set built using a
    configurable variant of the Schreier-Sims algorithm.

    The base parameter specifies a prefix of the base that will be used for the
    strong generating set. It will be extended automatically as necessary.
    """
    def __init__(self, cfg=None, base=()):
        self.cfg = cfg or Config()
        self.gens = []
        self.basepoint = None

        self.tree_gens = []
        self.tree_expand = []
        self.tree = {}
        self.tree_depth = 0

        self.stab = None
        self.rng = None

        if base:
            self.basepoint = base[0]
            self.tree = {self.basepoint: None}
            self.stab = Group(cfg, base=base[1:])

    def copy(self):
        """Returns a shallow copy of the group sharing the stabilizer chain.
        """
        group = copy.copy(self)
        group.gens = list(self.gens)
        group.rng = copy.deepcopy(self.rng)

        return group

    def deep_copy(self):
        """Return a deep copy of the group, not sharing the stabilizer chain.
        """
        group = self.copy()
        if group.stab:
            group.stab = self.stab.deep_copy()
        return group

    def stabilizer_chain(self, prefix=None):
        """Return the chain of stabilizer groups as a list.
        """
        if prefix is None:
            prefix = []
        prefix.append(self)
        if self.stab is None:
            return prefix
        return self.stab.stabilizer_chain(prefix)

    def generators(self):
        """Return a set of generators for this group.

        When cfg.monte_carlo is False or after successful verification this is
        a strong generating set.
        """
        if self.stab is None:
            return self.gens
        else:
            return self.stab.generators() + self.gens

    def base(self):
        """Return the base of the strong generating set for this group.
        """
        return [stab.basepoint for stab in self.stabilizer_chain()[:-1]]

    def remove_redundand_basepoints(self):
        """Remove redundant points from the base.
        """
        if not self.stab:
            return

        self.stab.remove_redundand_basepoints()

        if len(self.tree) == 1:
            stab = self.stab
            rng = self.rng
            self.__dict__.clear()
            self.__dict__.update(stab.__dict__)
            self.gens = list(self.gens)
            self.rng = rng

    def add_gen(self, gen):
        """Add a generator to this group.
        """
        gen = self.sift_gen(gen)

        if not is_id_perm(gen):
            if self.cfg.monte_carlo:
                if self.rng is None:
                    self.rng = GroupRng(self.cfg)
                    for g, g_inv in self.generators():
                        self.rng.add_gen(g)
                self.rng.add_gen(gen)

            self.add_nonmember_gen(gen)

    def build(self, known_order=None):
        """Run the Monte Carlo algorithm until the exit condition is met.

        If known_order is given, it is used as exit condition. Otherwise it
        exits after cfg.exit_rounds many rounds without progress.
        """
        if not self.cfg.monte_carlo:
            return
        if known_order is None:
            stationary_rounds = 0
            while stationary_rounds < self.cfg.exit_rounds:
                stationary_rounds += 1
                if self.monte_carlo_round():
                    stationary_rounds = 0
        else:
            while self.order() < known_order:
                self.monte_carlo_round()

    def monte_carlo_round(self):
        """One iteration of the randomized Monte Carlo Schreier-Sims algorithm.
        """
        if not self.cfg.monte_carlo:
            raise RuntimeError('cfg.monte_carlo is False')
        self.cfg.stats.rounds += 1
        p = self.rng.sample()
        p = self.sift_gen(p)
        p_is_missing = not is_id_perm(p)
        if p_is_missing:
            self.add_nonmember_gen(p)
        return p_is_missing

    def sift_gen(self, gen):
        """Perform sifting for adding a new generator.
        """
        if self.cfg.deep_sift_gens:
            return self.deep_sift(gen)
        else:
            return self.sift(gen)

    def sift(self, p, trace=None, depth=math.inf):
        """Perform sifting on p.

        Applies generators to p to fix p[basepoint] = basepoint and if
        successful recurses in stabilizer subgroup.

        When trace is a list, the applied generators are appended to it where
        each generator is represented as (i, p) where i is the index into
        generators() and p is True if the generator is inverted.
        """
        if not depth:
            return None

        if self.basepoint is None or p is None:
            return p

        try:
            p_stab = self.move_to_basepoint(p[self.basepoint], p, trace)
        except NotInOrbit:
            return p

        return self.stab.sift(p_stab, trace, depth=depth - 1)

    def deep_sift(self, p):
        """Perform deep sifting on p.

        Applies generators to both sides of p to fix p[basepoint] = basepoint
        and if successful recurses in stabilizer subgroup.
        """

        if self.basepoint is None:
            return p

        try:
            p_stab = self.move_to_basepoint(p[self.basepoint], p)
        except NotInOrbit:
            for a in sorted(self.tree.keys()):
                if p[a] in self.tree:
                    p_tmp = inv_perm(self.move_to_basepoint(a, inv_perm(p)))
                    p_stab = self.move_to_basepoint(p_tmp[self.basepoint], p)
                    return self.stab.deep_sift(p_stab)
            return p

        return self.stab.deep_sift(p_stab)

    def order(self):
        """Compute the order of the group.
        """
        if self.basepoint is None:
            return 1
        return len(self.tree) * self.stab.order()

    def move_to_basepoint(self, a, p=None, trace=None):
        """Apply a product q of tree_gens with q[a] == basepoints to p.
        """
        if a not in self.tree:
            raise NotInOrbit("a not in basepoints's orbit")

        steps = []

        while a != self.basepoint:
            edge_i, edge_pol = self.tree[a]
            edge_gen = self.tree_gens[edge_i][edge_pol]
            if trace is not None:
                trace.extend(self.tree_expand[edge_i][edge_pol])
            a = edge_gen[a]
            p = mult_perm(p, edge_gen)
            self.cfg.stats.products += 1

        return p

    def add_nonmember_gen(self, gen, las_vegas=False):
        """Add a generator that is not a member of this group yet.
        """
        if self.basepoint is None:
            for i, j in enumerate(gen):
                if i != j:
                    self.basepoint = i
                    break
            if self.basepoint is None:
                return

            self.stab = Group(self.cfg)

        if gen[self.basepoint] == self.basepoint:
            # we can add this generator directly to the stabilizer subgroup
            self.stab.add_nonmember_gen(gen, las_vegas)
        else:
            self.gens.append((gen, inv_perm(gen)))

        self.rebuild_schreier_tree()

        if self.cfg.monte_carlo or las_vegas:
            self.add_random_schreier_gens()
        else:
            self.add_all_schreier_gens()

    def add_all_schreier_gens(self):
        """Add all Schreier generators to the stabilizer subgroup.
        """
        for gen, inv_gen in self.generators():
            for a in sorted(self.tree.keys()):
                p = inv_perm(self.move_to_basepoint(a, inv_gen))
                schreier_gen = self.move_to_basepoint(p[self.basepoint], p)
                self.stab.add_gen(schreier_gen)

    def add_random_schreier_gens(self):
        """Add some random Schreier generators to the stabilizer subgroup.

        Used for the monte carlo variant if cfg.random_schreier_gens is
        nonzero.
        """
        gens = self.generators()
        orbit = sorted(self.tree.keys())
        for i in range(self.cfg.random_schreier_gens):
            gen, inv_gen = self.cfg.rng.choice(gens)
            a = self.cfg.rng.choice(orbit)
            p = inv_perm(self.move_to_basepoint(a, inv_gen))
            schreier_gen = self.move_to_basepoint(p[self.basepoint], p)

            # We can't use add_gen recursively for the monte carlo variant
            # so we manually sift and check for identity
            schreier_gen = self.sift_gen(schreier_gen)

            if not is_id_perm(schreier_gen):
                self.add_nonmember_gen(schreier_gen)

    def rebuild_schreier_tree(self):
        """Rebuild the Schreier tree for the orbit containing basepoint.
        """

        self.tree_gens = []
        self.tree_expand = []

        for i, gen in enumerate(self.generators()):
            self.tree_gens.append(gen)
            self.tree_expand.append(([(i, False)], [(i, True)]))

        while True:
            a = self.schreier_tree_bfs()
            if a is None:
                break
            # we add a new (tree only) generator to reach the element `a`
            # that was too far from the root in a single step.
            trace = []
            new_gen = self.move_to_basepoint(a, trace=trace)
            assert new_gen is not None
            inv_trace = [(i, not p) for i, p in reversed(trace)]
            self.tree_gens.append((inv_perm(new_gen), new_gen))
            self.tree_expand.append((inv_trace, trace))

    def schreier_tree_bfs(self):
        """Attempt to build a Schreier tree using a bfs search.

        Depending on cfg.schreier_tree it checks whether the tree is too deep
        (and might abort early) and returns a point of the orbit that is too
        deep (but reachable in the tree built so far).

        """
        mode = self.cfg.schreier_tree

        self.tree = {self.basepoint: None}
        queue = [(self.basepoint, 0)]

        while queue:
            a, depth = queue.pop(0)

            for i, gen_pair in enumerate(self.tree_gens):
                for pol, gen in enumerate(gen_pair):
                    b = gen[a]
                    if b in self.tree:
                        continue

                    self.tree[b] = (i, 1 - pol)
                    queue.append((b, depth + 1))

                    if mode == 'gap' and depth + 1 >= 2 * len(self.tree_gens):
                        # We want a shallow tree, so we break the outer loop to
                        # add a new generator to reach `a` in a single step.
                        # This is what GAP does, see also Seress's "Permutation
                        # Group Algorithms" 4.4.3
                        return b

        if mode == 'halve' and depth >= 2 * len(self.tree_gens):
            return a

        self.tree_depth = depth

    def split_generators(self):
        """Split any generator with a non-prime-power sized basepoint orbit.

        Was only used for debugging.

        For every maximal prime power p^i of the orbit size n, add the
        generator g^(n/p^i).
        """
        for stab in reversed(self.stabilizer_chain()):
            stab.split_generators_step()

    def split_generators_step(self):
        """One step of split_generators."""
        new_gens = []
        for g, g_inv in self.gens:
            orbit_size = 0
            a = g[self.basepoint]
            while a != self.basepoint:
                orbit_size += 1
                a = g[a]
            orbit_factors = {}
            p = 2
            n = orbit_size
            while p*p <= n:
                if n % p == 0:
                    n //= p
                    orbit_factors.setdefault(p, 1)
                    orbit_factors[p] *= p
                else:
                    p = (p + 1) | 1
            orbit_factors.setdefault(n, 1)
            orbit_factors[n] *= n

            for pn in sorted(orbit_factors.values()):
                n = orbit_size // pn
                new_g = exp_perm(g, n)
                new_gens.append((new_g, inv_perm(new_g)))

        self.gens = new_gens
        self.rebuild_schreier_tree()

    def remove_redundand_generators(self):
        """Removes redundant generators of the strong generating set.

        For each group in the stabilizer chain, it checks each generator in
        turn and removes it if it's generated by the stabilizer subgroup
        together with the other generators.

        Can raise a IncompleteStrongGeneratingSet exception, but is not
        guaranteed to do so for an incomplete strong generating set.
        """
        for stab in reversed(self.stabilizer_chain()):
            stab.remove_redundand_generators_step()

    def remove_redundand_generators_step(self):
        """Removes redundant top-level generators of the strong generating set.

        It checks each top-level generator in turn and removes it if it's
        generated by the stabilizer subgroup together with the other top-level
        generators.

        Can raise a IncompleteStrongGeneratingSet exception, but is not
        guaranteed to do so for an incomplete strong generating set.
        """
        if self.basepoint is None:
            return

        require_rebuild = False

        i = 0
        while i < len(self.gens):
            gen_i, gen_i_inv = self.gens.pop(i)
            self.rebuild_schreier_tree()
            require_rebuild = False
            if gen_i[self.basepoint] in self.tree:
                witness = self.sift(gen_i)
                if is_id_perm(witness):
                    continue
                else:
                    raise IncompleteStrongGeneratingSet(
                        'incomplete strong generating set detected '
                        'while removing redundant generators',
                        witness=witness)
            self.gens[i:i] = [(gen_i, gen_i_inv)]
            require_rebuild = True
            i += 1

        if require_rebuild:
            self.rebuild_schreier_tree()

    def apply_map(self, f_perm, f_base):
        """Apply f_perm to every generator and f_base to every base point.

        The resulting state is only valid if f_perm is a group morphism and if
        `f_base(b g) = f_base(b) f_perm(g)` for all group elements g and all
        base points b.
        """

        if self.basepoint is None:
            return

        self.stab.apply_map(f_perm, f_base)

        def f_gen(g):
            g = f_perm(g)
            return (g, inv_perm(g))

        self.gens = [f_gen(g) for g, g_inv in self.gens]
        self.basepoint = f_base(self.basepoint)

        gens = self.generators()

        extra_tree_gens = self.tree_gens[len(gens):]
        self.tree_gens = gens + [f_gen(g) for g, g_inv in extra_tree_gens]

        self.tree = {f_base(b): step for b, step in self.tree.items()}

    def conjugate(self, w):
        """Conjugate this group.
        """
        if w is None:
            return

        w_inv = inv_perm(w)

        def f_perm(g):
            self.cfg.stats.products += 2
            return mult_perms((w_inv, g, w))

        def f_base(b):
            return w[b]

        self.apply_map(f_perm, f_base)

    def insert_toplevel_block_stabilizer(self):
        """Insert a block stabilizer subgroup directly below this group.

        If this group is generated by the stabilizer subgroup and more than a
        single generator, this will insert a subgroup having the current
        stabilizer as proper subgroup such that this group is generated by the
        new subgroup and a _single_ generator.

        This is done by finding a block system such that the block containing
        the first base point is invariant under the action of all but one
        generator.

        Then the set on which this group acts is extended by points
        representing the orbit of that block, so that stabilizing the point
        corresponding to the block containing our basepoint stabilizes that
        block.

        This is used for SGS verification methods that can only (or are more
        efficient when they) verify adding a single generator to an already
        verified subgroup.

        See also Seress's "Permutation Group Algorithms" 8.2. "general case"

        This construction may fail if the strong generating set is incomplete.
        In that case an IncompleteStrongGeneratingSet exception containing a
        missing group element as witness is raised.
        """

        if len(self.gens) < 2:
            return False

        gens = self.generators()

        block_stab = Group(cfg=self.cfg)

        block_stab.stab = self.stab
        block_stab.gens = self.gens[:-1]
        block_stab.basepoint = self.basepoint
        block_stab.rebuild_schreier_tree()

        self.stab = block_stab
        self.gens = [self.gens[-1]]

        n = len(self.gens[0][0])

        repr_map = [None] * n

        base_block = list(self.stab.tree)
        base_repr = min(base_block)
        for point in base_block:
            repr_map[point] = base_repr

        queue = [base_block]
        block_paths = {base_repr: []}

        while queue:
            block_a = queue.pop()
            repr_a = min(block_a)
            for g_pair in gens:
                for g in g_pair:
                    block_b = list(g[a] for a in block_a)
                    repr_b = min(block_b)
                    if repr_map[repr_b] is None:
                        queue.append(block_b)
                        block_paths[repr_b] = block_paths[repr_a] + [g]
                    for b in block_b:
                        if repr_map[b] not in (None, repr_b):
                            repr_c = repr_map[b]

                            # find v in stab so that path_b'[b] v = basept
                            # find w in stab so that path_c'[b] w = basept
                            # v' path_b path_c' w is a witness

                            path_b = block_paths[repr_a] + [g]
                            path_c = block_paths[repr_c]

                            self.cfg.stats.products += (
                                len(path_b) + len(path_c) - 1)
                            path_b = mult_perms(path_b)
                            path_c = mult_perms(path_c)

                            part_b = self.move_to_basepoint(
                                b, inv_perm(path_b))
                            part_c = self.move_to_basepoint(
                                b, inv_perm(path_c))

                            witness = self.stab.sift(
                                mult_perm(inv_perm(part_b), part_c))

                            assert not is_id_perm(witness)

                            raise IncompleteStrongGeneratingSet(
                                'incomplete strong generating set detected'
                                ' while checking blocks', witness=witness)

                        repr_map[b] = repr_b

        blocks = sorted(set(repr_map) - set([None]))

        if len(blocks) == 1:
            residue = self.stab.sift(self.gens[0][0])
            if not is_id_perm(residue):
                raise IncompleteStrongGeneratingSet(
                    'incomplete strong generating set detected'
                    ' while checking blocks', witness=residue)

            self.stab.rng = self.rng
            self.__dict__ = self.stab.__dict__

            return self.insert_toplevel_block_stabilizer()

        block_names = {block: n + i for i, block in enumerate(blocks)}

        def f_perm(g):
            return g + [block_names[repr_map[g[block]]] for block in blocks]

        self.apply_map(f_perm, lambda b: b)

        self.basepoint = block_names[repr_map[self.basepoint]]
        self.rebuild_schreier_tree()

        return True

    def sample(self, depth=math.inf):
        """Return a uniform random element or coset representative.

        The depth parameter specifies for which subgroup on the chain a coset
        representative should be sampled. If it isn't given or points past the
        last non-trivial subgroup this returns a uniformly sampled group
        element.

        This is only correct if the strong generating set is complete.
        """
        if not depth or self.basepoint is None:
            return None

        if self.stab is None:
            stab_sample = None
        else:
            stab_sample = self.stab.sample(depth - 1)

        target = self.cfg.rng.choice(sorted(self.tree))

        traversal_sample = self.move_to_basepoint(target)

        self.cfg.stats.products += 1

        return mult_perm(stab_sample, inv_perm(traversal_sample))

    def orbit(self, a):
        """Compute the orbit of a point.

        Returns a dictionary with the orbit points as key and a word that moves
        a to that orbit point as value.
        """
        if a in self.tree:
            return sorted(self.tree)

        orbit = {a: []}
        queue = [a]

        gens = self.generators()

        while queue:
            a = queue.pop()
            for i, gen in enumerate(gens):
                for p, g in enumerate(gen):
                    b = g[a]
                    if b not in orbit:
                        orbit[b] = orbit[a] + [(i, bool(p))]
                        queue.append(b)

        return orbit

    def stabilized_depth(self, a):
        """Return the depth of the first subgroup that does not move a.
        """
        if self.basepoint is None:
            return 0

        stab_depth = self.stab.stabilized_depth(a)

        if stab_depth == 0:
            for g, g_inv in self.gens:
                if g[a] != a:
                    break
            else:  # if no generator moved a
                return 0

        return 1 + stab_depth

    def change_base(self, base):
        """Change the base to start with the given prefix.

        Often this is much faster than generating a new strong generating set
        with the given base prefix. Sometimes it may be slower though, more
        tuning of the involved heuristics might help here.
        """

        # Inserting redundand base points is free, so we can remove them here
        # to get better estimates
        self.remove_redundand_basepoints()

        # How many levels are we allowed to regenerate before building the
        # remaining chain from scratch?
        budget = len(self.base()) + len(base)

        # Perform a base change that permutes the new base using group elements
        # to save work
        w = self.change_base_conjugated(base, budget)

        # Undo the permutation of the new base by conjugating the group
        self.conjugate(inv_perm(w))

    def change_base_conjugated(self, base, budget):
        """Change the base to start with a permutation of the given prefix.

        Returns the permutation applied to the base perfix. This permutation
        will be an element of the group.

        Used to implement change_base.
        """

        # TODO it would be nice to refactor this into smaller parts that can be
        # used independently

        # Nothing to do
        if not base:
            return

        # Trivial group, rest are redundant base points
        if self.basepoint is None:
            group = Group(self.cfg, base)
            group.rng = self.rng
            self.__dict__ = group.__dict__
            return

        # Can we permute the base so that the first basepoint matches?
        if base[0] in self.tree:
            w = self.move_to_basepoint(base[0])
            if w:
                base = [w[b] for b in base]
            assert base[0] == self.basepoint

            old_gens = self.generators()

            w_stab = self.stab.change_base_conjugated(base[1:], budget)

            if self.generators() != old_gens:
                self.rebuild_schreier_tree()

            w = mult_perm(w, w_stab)
            if w:
                self.cfg.stats.products += 1
            return w

        # Is the first base point of the old base redundant?
        # TODO is this still necessary? do we avoid this in all cases?
        if len(self.tree) == 1:
            group = self.stab.copy()
            group.rng = self.rng
            self.__dict__ = group.__dict__
            return self.change_base_conjugated(base)

        # If the first base point of the new base is redundant we can insert it
        # for free
        if self.stabilized_depth(base[0]) == 0:
            group = Group(self.cfg, [base[0]])
            group.rng, self.rng = self.rng, None
            self.__dict__, group.__dict__ = group.__dict__, self.__dict__
            self.stab = group
            # no need to rebuild schreier tree, this is a redundant basepoint
            return self.stab.change_base_conjugated(base[1:], budget)

        # Otherwise we need to regenerate at least some levels of the
        # stabilizer chain

        # We use the known order to be able to do this as a Las Vegas algorithm
        target_order = self.order()

        # We can reduce the number of levels to regenerate, by permuting the
        # base. We're going to permute it such that the first basepoint of the
        # new base is the first one that is stabilized in the old stabilizer
        # chain.
        orbit = self.orbit(base[0])

        depth, point, word = min(
            (self.stabilized_depth(a), a, w) for a, w in orbit.items())

        gens = self.generators()

        w_adjust = mult_perms(gens[i][p] for i, p in word)
        if word:
            self.cfg.stats.products += len(word) - 1
            base = [w_adjust[b] for b in base]

        assert base[0] == point

        # When regenerating all levels up to the point where our new first
        # basepoint is stabilized, we can also process all following basepoints
        # that would also be stabilized up to that point.
        cut = 1

        while cut < len(base) and self.stabilized_depth(base[cut]) <= depth:
            cut += 1

        full_depth = len(self.stabilizer_chain())

        # keep order, but remove duplicates
        new_base = list(dict.fromkeys(base[:cut] + self.base()[:depth]))
        sift_depth = len(new_base)

        # We pay for the levels that we're regenerating that are not part of
        # our target prefix
        budget -= sift_depth - cut

        # If we're out of budget, we regenerate all remaining levels
        if budget < 0:
            # It seems that we can reduce the number of random samples we need
            # by keeping the old base as suffix and sifting the old strong
            # generating set first
            new_base = list(dict.fromkeys(base + self.base()))
            group = Group(self.cfg, new_base)

            # Las Vegas sample, sift and add if missing loop:
            while group.order() < target_order:
                if gens:
                    g = gens.pop()[0]
                else:
                    g = self.sample()
                self.cfg.stats.rounds += 1
                residue = group.sift(g)
                if not is_id_perm(residue):
                    group.add_nonmember_gen(residue, las_vegas=True)

            group.rng = self.rng
            self.__dict__ = group.__dict__

            self.stabilizer_chain()[len(base)].remove_redundand_basepoints()

            return w_adjust

        # If we still have budget we only generate some levels and reuse the
        # suffix of the stabilizer chain that we don't regenerate
        group = Group(self.cfg, new_base)

        new_chain_part = group.stabilizer_chain()[:-1]
        new_chain_part[-1].stab = self.stabilizer_chain()[depth]

        for subgroup in reversed(new_chain_part):
            subgroup.rebuild_schreier_tree()

        # We start by sifting the old strong generators from the part that
        # we're regenerating, again this seems to be quite effective in
        # reducing the number of random samples needed.
        gens = []
        for subgroup in self.stabilizer_chain()[:depth]:
            gens.extend(subgroup.gens)

        # Las Vegas sample, sift and add if missing loop:
        while group.order() < target_order:
            if gens:
                g = gens.pop(0)[0]
            else:
                g = self.sample(depth)
            self.cfg.stats.rounds += 1
            residue = group.sift(g, depth=sift_depth)
            if not is_id_perm(residue):
                group.add_nonmember_gen(residue, las_vegas=True)

        group.rng = self.rng
        self.__dict__ = group.__dict__

        self.stabilizer_chain()[cut].remove_redundand_basepoints()

        # This recursive call will quickly shift through the prefix we already
        # handled and continues after that
        w = self.change_base_conjugated(base, budget)

        w = mult_perm(w_adjust, w)
        self.cfg.stats.products += 1

        return w

    def verify(self):
        """Perform a verification for every step of the stabilizer chain.

        Uses the verification method specified in the Config.
        """
        if self.cfg.verify == 'schreier':
            for stab in reversed(self.stabilizer_chain()):
                self.schreier_verify_step()
        elif self.cfg.verify == 'sims':
            group = self.deep_copy()
            group.remove_redundand_generators()
            for stab in reversed(group.stabilizer_chain()):
                group.simss_verify_step()
        else:
            raise ValueError(f'unknown verification method {self.cfg.verify}')

    def build_verified(self, known_order=None):
        """Run the Monte Carlo algorithm until verification succeeds.

        If known_order is passed, no verification is necessary and thus no
        verification is performed.

        Returns the number of verification failures.
        """
        if known_order is not None:
            self.build(known_order)
            return 0

        failures = 0
        self.build()
        while True:
            try:
                self.verify()
            except IncompleteStrongGeneratingSet as e:
                failures += 1
                self.add_nonmember_gen(e.witness)
                self.build()
            else:
                break
        return failures

    def schreier_verify_step(self):
        """Slow verification by sifting all Schreier generators.

        Raises an IncompleteStrongGeneratingSet exception if verification
        fails.
        """
        if self.basepoint is None:
            return True

        for gen, inv_gen in self.generators():
            for a in sorted(self.tree.keys()):
                p = inv_perm(self.move_to_basepoint(a, inv_gen))
                schreier_gen = self.move_to_basepoint(p[self.basepoint], p)
                residue = self.stab.sift(schreier_gen)
                if not is_id_perm(residue):
                    raise IncompleteStrongGeneratingSet(
                        "incomplete strong generating set detected"
                        " while sifting Schreier generators",
                        witness=residue)

    def simss_verify_step(self):
        """Verifies the strong generating set using Sims's verify routine.

        Assumes stab has a strong generating set, i.e. this needs to be called
        for every subgroup on the stabilizer chain to verify a strong
        generating set.

        This follows the description given by Seress in "Permutation Group
        Algorithms" 8.2.

        Raises an IncompleteStrongGeneratingSet exception if verification
        fails.
        """
        if not self.gens:
            return

        if len(self.gens) > 1:
            # handles only a single additional generator, so we insert block
            # stabilizers, see insert_toplevel_block_stabilizer
            group = self.deep_copy()
            try:
                if group.insert_toplevel_block_stabilizer():
                    group.stab.simss_verify_step()
                group.simss_verify_step()
            except IncompleteStrongGeneratingSet as e:
                e.witness = e.witness[:len(self.gens[0][0])]
                raise e
            return

        stab = self.stab.deep_copy()

        z, z_inv = self.gens[0]
        gamma = z_inv[self.basepoint]

        stab.change_base([gamma])

        stab_gamma = stab.stab

        gens = [g for g, g_inv in stab.generators()] + [z_inv]

        z_idx = len(gens) - 1
        stab_idx = z_idx - len(stab.gens)

        words = {self.basepoint: []}
        sequence = [self.basepoint]
        gen_pos = [1] * z_idx + [0]

        gen_idx = z_idx

        while gen_idx < len(gens):
            if gen_pos[gen_idx] == len(sequence):
                gen_idx += 1
                continue

            point = sequence[gen_pos[gen_idx]]
            gen_pos[gen_idx] += 1
            gen = gens[gen_idx]
            next_point = gen[point]

            if next_point not in words:
                sequence.append(next_point)
                words[next_point] = words[point] + [gen_idx]
                gen_idx = 0

        deltas = []
        u = {}
        y = {}
        v = {}
        gammas = []

        for point in sequence:
            word = words[point]
            if not word or word[-1] >= stab_idx:
                gamma_ij = point
                if not word or word[-1] == z_idx:
                    # point is also delta_i
                    deltas.append(point)

                gammas.append(gamma_ij)

                stab_orbit_repr_prefix = max(
                    (i + 1 for i, idx in enumerate(word) if idx == z_idx),
                    default=0)

                prefix = word[:stab_orbit_repr_prefix]
                suffix = word[stab_orbit_repr_prefix:]

                u_delta_i = mult_perms(gens[i] for i in prefix)
                y_gamma_ij = mult_perms(gens[i] for i in suffix)

                u_gamma_ij = mult_perm(u_delta_i, y_gamma_ij)
                self.cfg.stats.products += max(0, len(word) - 1)

                y[gamma_ij] = y_gamma_ij or id_perm(len(z))
                u[gamma_ij] = u_gamma_ij or id_perm(len(z))

        for gamma_ij in gammas:
            z_word = words[z[gamma_ij]]

            stab_gamma_orbit_repr_prefix = max(
                (i + 1 for i, idx in enumerate(z_word) if idx >= stab_idx),
                default=0)

            prefix = z_word[:stab_gamma_orbit_repr_prefix]
            suffix = z_word[stab_gamma_orbit_repr_prefix:]

            stab_gamma_orbit_repr = self.basepoint
            for idx in prefix:
                stab_gamma_orbit_repr = gens[idx][stab_gamma_orbit_repr]

            v_gamma_ij_z = mult_perm(
                u[stab_gamma_orbit_repr],
                mult_perms(gens[i] for i in suffix))

            self.cfg.stats.products += len(suffix)

            v[z[gamma_ij]] = v_gamma_ij_z or id_perm(len(z))

        for gamma_ij in gammas:
            assert v[z[gamma_ij]][self.basepoint] == z[gamma_ij]

        # because the following two assertions hold, we can skip checking
        # delta_0 when checking (c)
        assert deltas[0] == self.basepoint
        assert u[self.basepoint] == id_perm(len(z))

        # because the following two assertions hold, we automatically check (a)
        # while checking (c)
        assert deltas[1] == gamma
        assert u[gamma] == z_inv

        # check (c), we can skip delta_0
        for delta_i in deltas[1:]:
            stab.change_base([delta_i])

            u_delta_i = u[delta_i]
            u_delta_i_inv = inv_perm(u_delta_i)

            for g, g_inv in stab.stab.generators():
                g_u_delta_i_inv = mult_perms((
                    u_delta_i, g, u_delta_i_inv))
                self.cfg.stats.products += 2

                residue = stab.sift(g_u_delta_i_inv)

                if not is_id_perm(residue):
                    raise IncompleteStrongGeneratingSet(
                        "incomplete strong generating set detected"
                        " during Sims's verification step",
                        witness=residue)

        # check (b)
        for gamma_ij in gammas:
            b = mult_perms((u[gamma_ij], z, inv_perm(v[z[gamma_ij]])))
            residue = stab.sift(b)
            if not is_id_perm(residue):
                raise IncompleteStrongGeneratingSet(
                    "incomplete strong generating set detected"
                    " during Sims's verification step",
                    witness=residue)
