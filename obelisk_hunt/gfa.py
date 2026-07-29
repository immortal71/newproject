"""Minimal GFA1 parser for assembly-graph circularity evidence.

We only need one signal out of the graph: does a segment have a link back to
itself (either orientation)? That is the de Bruijn-graph analogue of the
sequence-level terminal redundancy in circularity.py, and it is the same
"contig links to its own start" check used by Recycler/Unicycler to call
circular plasmids. It is independent evidence: a graph self-loop can exist
even when the sequence-level overlap is shorter than our detection window
(e.g. a bare k-1 overlap), and a sequence-level overlap can exist without a
matching graph edge on partially-scaffolded assemblies.

Only S (segment) and L (link) lines are parsed; everything else is ignored.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class SelfLoop:
    segment: str
    from_orient: str
    to_orient: str
    cigar: str


@dataclass
class GfaGraph:
    segment_lengths: dict = field(default_factory=dict)
    self_loops: dict = field(default_factory=dict)  # segment -> list[SelfLoop]

    def has_self_loop(self, segment: str) -> bool:
        return bool(self.self_loops.get(segment))


def parse_gfa(path: str) -> GfaGraph:
    graph = GfaGraph()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            tag = fields[0]

            if tag == "S":
                name = fields[1]
                seq = fields[2]
                graph.segment_lengths[name] = len(seq) if seq != "*" else None

            elif tag == "L":
                # L  From  FromOrient  To  ToOrient  CIGAR
                from_seg, from_orient, to_seg, to_orient, cigar = fields[1:6]
                if from_seg == to_seg:
                    graph.self_loops.setdefault(from_seg, []).append(
                        SelfLoop(from_seg, from_orient, to_orient, cigar)
                    )

    return graph
