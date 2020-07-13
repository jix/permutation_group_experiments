from dataclasses import dataclass, field
from typing import Literal
from random import Random
import re


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
    return all(i == j for i, j in enumerate(p))


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
    # 'gap':   Use the algorithm implemented by GAP and described in Seress'
    #          "Permutation Group Algorithms" 4.4.3.
    # 'halve': Variant of 'gap' which adds a generator for the deepest point of
    #          the orbit instead of the first exceeding the threshold.
    # 'deep':  Don't build shallow Schreier trees.
    schreier_tree: Literal['gap', 'halve', 'deep'] = 'gap'

    # Don't keep a list of generators for each stabilizer, instead use the
    # strong generating set to generate Schreier trees and Schreier generators.
    #
    # This is very inefficient when using the deterministic Schreier-Sims
    # algorithm.
    incomplete_gens: bool = True

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

    # Parameters for black box random element generation. See GAP's
    # documentation for ProductReplacer.
    rng_accus: int = 5
    rng_extra_slots: int = 5
    rng_scramble: int = 30
    rng_scramble_factor: int = 4

    stats: Stats = field(default_factory=lambda: Stats())
    rng: Random = field(default_factory=lambda: Random(0))


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


class Group:
    """A permutation group.

    The group is represented using a strong generating set built using a
    configurable variant of the Schreier-Sims algorithm.
    """
    def __init__(self, cfg=None):
        self.cfg = cfg or Config()
        self.gens = []
        self.gens_complete = True
        self.basepoint = None

        self.tree_gens = []
        self.tree_expand = []
        self.tree = {}
        self.tree_depth = 0

        self.stab = None
        self.rng = None
        self.blocks = None

    def stabilizer_chain(self, prefix=None):
        """Return the chain of stabilizer groups as a list.
        """
        if prefix is None:
            prefix = []
        prefix.append(self)
        if self.stab is None:
            return prefix
        return self.stab.stabilizer_chain(prefix)

    def strong_generators(self):
        """Return a strong generating set for this group.
        """
        strong_gens = []
        seen = set()
        for group in reversed(self.stabilizer_chain()):
            for gen in group.gens:
                key = tuple(gen[0])
                if key in seen:
                    continue
                seen.add(key)
                strong_gens.append(gen)
        return strong_gens

    def generators(self):
        """Return a set of generators for this group.
        """
        if self.gens_complete or self.stab is None:
            return self.gens
        else:
            return self.stab.generators() + self.gens

    def add_gen(self, gen):
        """Add a generator to this group.
        """
        gen = self.sift_gen(gen)

        if not is_id_perm(gen):
            if self.cfg.monte_carlo:
                if self.rng is None:
                    self.rng = GroupRng(self.cfg)
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

    def sift(self, p, trace=None):
        """Perform sifting on p.

        Applies generators to p to fix p[basepoint] = basepoint and if
        successful recurses in stabilizer subgroup.

        When trace is a list, the applied generators are appended to it.
        """
        if self.basepoint is None or p is None:
            return p

        try:
            p_stab = self.move_to_basepoint(p[self.basepoint], p, trace)
        except NotInOrbit:
            return p

        return self.stab.sift(p_stab, trace)

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

    def add_nonmember_gen(self, gen):
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

        self.gens.append((gen, inv_perm(gen)))

        if gen[self.basepoint] == self.basepoint:
            # we can add this generator directly to the stabilizer subgroup
            self.stab.add_nonmember_gen(gen)

            if self.cfg.incomplete_gens:
                self.gens.pop()
                self.gens_complete = False

        self.rebuild_schreier_tree()

        if self.cfg.monte_carlo:
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

    def verify_all_schreir_gens(self):
        """Slow verification by sifting all schreier generators.

        This is performed recursively down the stabilizer chain.
        Returns True if the strong generating set is complete.
        """
        if self.basepoint is None:
            return True

        for gen, inv_gen in self.generators():
            for a in sorted(self.tree.keys()):
                p = inv_perm(self.move_to_basepoint(a, inv_gen))
                schreier_gen = self.move_to_basepoint(p[self.basepoint], p)
                if not is_id_perm(self.stab.sift(schreier_gen)):
                    return False

        return self.stab.verify_all_schreir_gens()

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

        for gen in self.generators():
            self.tree_gens.append(gen)
            self.tree_expand.append(([gen[0]], [gen[1]]))

        while True:
            a = self.schreier_tree_bfs()
            if a is None:
                break
            # we add a new (tree only) generator to reach the element `a`
            # that was too far from the root in a single step.
            trace = []
            new_gen = self.move_to_basepoint(a, trace=trace)
            assert new_gen is not None
            inv_trace = [inv_perm(g) for g in trace[::-1]]
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
                        # This is what GAP does, see also Seress' "Permutation
                        # Group Algorithms" 4.4.3
                        return b

        if mode == 'halve' and depth >= 2 * len(self.tree_gens):
            return a

        self.tree_depth = depth

    def split_generators(self):
        """Split any generator with a non-prime-power sized basepoint orbit.

        Was only used for debugging build_block_stabilizer_chain.

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

    def build_block_stabilizer_chain(self):
        """Insert block stabilizer into the chain.

        This inserts groups into the stabilizer chain, so that each group is
        generated by the following group and a single additional generator.
        This is realized by first testing for and removing redundand
        generators. If multiple generators g_1, ..., g_k remain, insert the
        groups <H, g_1>, <H, g_1, g_2>, ... <H, g_1, ..., g_{k-1}>. Each of
        these groups stabilizes a block of the following group. The stabilized
        block is the orbit of basepoint in the smaller group.

        See also Seress' "Permutation Group Algorithms" 8.2. "general case"

        This construction may fail if the strong generating set is incomplete.
        In that case an IncompleteStrongGeneratingSet exception containing a
        missing group element as witness is raised.

        The subgroups are not inserted as objects, but are defined implicitly
        by taking a prefix of the `gens` field. The corresponding blocks are
        stored in the `blocks` field.
        """
        for stab in reversed(self.stabilizer_chain()):
            stab.build_block_stabilizer_chain_step()

    def build_block_stabilizer_chain_step(self):
        """Insert block stabilizers for one step in the stabilizer chain.

        See build_block_stabilizer_chain for a description.
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

        self.blocks = [None]

        for i in range(1, len(self.gens)):
            stab_gens = self.stab.generators() + self.gens[:i]

            queue = [self.basepoint]
            block = set(queue)

            while queue:
                a = queue.pop()
                for g_pair in stab_gens:
                    for g in g_pair:
                        b = g[a]
                        if b not in block:
                            block.add(b)
                            queue.append(b)

            block = frozenset(block)

            point_to_block = {a: block for a in block}

            queue = [block]
            blocks = set(queue)
            block_paths = {block: []}

            gens = stab_gens + [self.gens[i]]

            while queue:
                block_a = queue.pop()
                for g_pair in gens:
                    for g in g_pair:
                        block_b = frozenset(g[a] for a in block_a)
                        for b in block_b:
                            point_to_block.setdefault(b, block_b)
                        for b in block_b:
                            if point_to_block[b] != block_b:
                                block_c = point_to_block[b]

                                # find v in stab so that path_b'[b] v = basept
                                # find w in stab so that path_c'[b] w = basept
                                # v' path_b path_c' w is a witness

                                path_b = block_paths[block_a] + [g]
                                path_c = block_paths[block_c]

                                self.cfg.stats.products += (
                                    len(path_b) + len(path_c) - 1)
                                path_b = mult_perms(path_b)
                                path_c = mult_perms(path_c)

                                saved_gens = self.gens
                                self.gens = self.gens[:i]
                                self.rebuild_schreier_tree()

                                part_b = self.move_to_basepoint(
                                    b, inv_perm(path_b))
                                part_c = self.move_to_basepoint(
                                    b, inv_perm(path_c))

                                witness = self.sift(
                                    mult_perm(inv_perm(part_b), part_c))

                                self.gens = saved_gens
                                self.rebuild_schreier_tree()

                                assert not is_id_perm(witness)

                                raise IncompleteStrongGeneratingSet(
                                    'incomplete strong generating set detected'
                                    ' while checking blocks', witness=witness)

                        if block_b not in blocks:
                            blocks.add(block_b)
                            queue.append(block_b)
                            block_paths[block_b] = block_paths[block_a] + [g]

            self.blocks.append((blocks, point_to_block))
