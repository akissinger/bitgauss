//! Graph realization over GF(2): rewriting a matrix so that every column has
//! Hamming weight at most 2.
//!
//! A matrix with this property and the same rowspace as a given matrix exists
//! if and only if the binary matroid represented by the columns of the matrix
//! is *graphic*; in that case the rewritten matrix can be taken to be the
//! vertex-edge incidence matrix of a realizing graph (with one redundant row
//! per connected component removed). Graphicness is decided, and the graph
//! constructed, with the algorithm of
//!
//! > R. Bixby and D. Wagner, *An almost linear-time algorithm for graph
//! > realization*, Mathematics of Operations Research 13 (1988), pp. 99-123.
//!
//! This module is a direct port of a Java implementation of the Bixby-Wagner
//! algorithm by Takahiro Ohto (Copyright (C) 2001-2003).  It is
//! dual licensed under 3-clause BSD (following the original) and Apache 2.0
//! (following the rest of the BitGauss project).
//!
//! The port was produced by Claude Fable 5 and maintains the same control flow
//! and data structures, so it can be checked against the original line by line,
//! see `docs/graphic/GRP.java`. The original 18 unit tests were also produced by
//! the LLM and checked by a human.

use crate::bitmatrix::BitMatrix;

impl BitMatrix {
    /// Returns a matrix `g` with the same rowspace as `self` in which every
    /// column has Hamming weight at most 2, or `None` if no such matrix
    /// exists.
    ///
    /// Such a matrix exists precisely when the binary matroid represented by
    /// the columns of `self` is graphic, which is decided with the
    /// Bixby-Wagner graph-realization algorithm. On success `g` has full row
    /// rank, i.e. `self.rank()` rows.
    ///
    /// Equivalent to `self.graphic_form_with_options(false, false).0`; see
    /// [`graphic_form_with_options`](Self::graphic_form_with_options) for
    /// partial realizations and the basis-change matrix.
    pub fn graphic_form(&self) -> Option<BitMatrix> {
        self.graphic_form_with_options(false, false).0
    }

    /// Generalization of [`graphic_form`](Self::graphic_form) that can
    /// optionally allow partial solutions and compute the basis-change matrix.
    ///
    /// This returns a triple `(g, b, h)` where `g` and `b` are optional matrices
    /// and `h` is a vector of column indices.
    ///
    /// If `partial == false`, then this method returns the graphic form of the
    /// matrix as `g`, or `None` if no such form exists.
    ///
    /// If `partial == true`, then a matrix is always returned for `g`, but some
    /// columns may have Hamming weight greater than 2. The indices of these
    /// columns are returned as `h`.
    ///
    /// When `basis_change == true` and `g != None`, an invertible matrix `b`
    /// is returned such that `g == b * self`.
    pub fn graphic_form_with_options(
        &self,
        partial: bool,
        basis_change: bool,
    ) -> (Option<BitMatrix>, Option<BitMatrix>, Vec<usize>) {
        let cols = self.cols();
        let mut rref = self.clone();
        let mut row_ops = basis_change.then(|| BitMatrix::identity(self.rows()));
        let pcols = match row_ops.as_mut() {
            Some(ops) => rref.gauss_helper(true, 1, ops),
            None => rref.gauss_helper(true, 1, &mut ()),
        };
        if pcols.is_empty() {
            let b = basis_change.then(|| BitMatrix::zeros(0, self.rows()));
            return (Some(BitMatrix::zeros(0, cols)), b, Vec::new());
        }

        let mut circuits = fundamental_circuits(&rref, &pcols);
        let mut skipped: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut td = loop {
            match try_realize(pcols.len(), &circuits) {
                Ok(td) => break td,
                Err(_) if !partial => return (None, None, Vec::new()),
                // drop the offending circuit and rerun: a failed insertion
                // leaves the t-decomposition mid-mutation, and removing a
                // circuit can also split a block, so the cheapest sound
                // recovery is to redo the (near-linear) realization
                Err(ci) => skipped.push(circuits.remove(ci)),
            }
        };
        let mut result = assemble_incidence(&mut td, &pcols, cols);

        // Any basis of the rowspace forces each skipped column to be the XOR
        // of the spanning-forest columns in its fundamental circuit.
        for (j, circuit) in &skipped {
            for &i in circuit {
                for row in 0..result.rows() {
                    if result.bit(row, pcols[i]) {
                        let b = result.bit(row, *j);
                        result.set_bit(row, *j, !b);
                    }
                }
            }
        }

        let mut skipped_cols: Vec<usize> = skipped.into_iter().map(|(j, _)| j).collect();
        skipped_cols.sort_unstable();
        let b = row_ops.map(|ops| self::basis_change(&result, &ops, &pcols));
        (Some(result), b, skipped_cols)
    }
}

/// Expresses each row of `realized` as a combination of the rows of the
/// original matrix, returning `b` with `realized == b * original`.
///
/// The rows of `realized` lie in the rowspace, and the RREF rows form a basis
/// in which each basis vector is 1 on its own pivot column and 0 on the other
/// pivot columns, so row `r` of `realized` is the XOR of the RREF rows `i`
/// with `realized[(r, pcols[i])]` set. `row_ops` is the row-operation record
/// of the elimination (the same operations applied to an identity matrix,
/// so RREF row `i` equals `row_ops` row `i` times the original matrix).
fn basis_change(realized: &BitMatrix, row_ops: &BitMatrix, pcols: &[usize]) -> BitMatrix {
    let mut b = BitMatrix::zeros(realized.rows(), row_ops.cols());
    for r in 0..realized.rows() {
        for (i, &p) in pcols.iter().enumerate() {
            if realized.bit(r, p) {
                b.add_bits_to_row(row_ops.row(i), r);
            }
        }
    }
    b
}

/// Reads the fundamental circuits off a reduced row echelon form: the pivot
/// columns are a basis ("tree edges", numbered by their row `0..rank`), and
/// each nonzero non-pivot column is a "cotree edge" whose circuit is the set
/// of tree edges it depends on. Returns `(column, circuit)` pairs.
fn fundamental_circuits(rref: &BitMatrix, pcols: &[usize]) -> Vec<(usize, Vec<usize>)> {
    let mut is_pivot = vec![false; rref.cols()];
    for &p in pcols {
        is_pivot[p] = true;
    }

    let mut circuits: Vec<(usize, Vec<usize>)> = Vec::new();
    for (j, &pivot) in is_pivot.iter().enumerate() {
        if pivot {
            continue;
        }
        let circuit: Vec<usize> = (0..pcols.len()).filter(|&i| rref.bit(i, j)).collect();
        if !circuit.is_empty() {
            circuits.push((j, circuit));
        }
    }
    circuits
}

/// Runs the Bixby-Wagner realization on all circuits, block by block. Returns
/// the resulting t-decomposition, or the index (into `circuits`) of the first
/// circuit that could not be realized.
fn try_realize(rank: usize, circuits: &[(usize, Vec<usize>)]) -> Result<TDecomp, usize> {
    let blocks = decompose_blocks(rank, circuits);
    let mut td = TDecomp::new(rank);
    for block in &blocks {
        let block_circuits: Vec<(usize, &[usize])> = block
            .iter()
            .map(|&ci| (circuits[ci].0, circuits[ci].1.as_slice()))
            .collect();
        if let Err(pos) = td.realization(&block_circuits) {
            return Err(block[pos]);
        }
    }
    Ok(td)
}

/// Builds the output matrix from a completed realization: extracts the
/// endpoints of the edge realizing each matrix column and assembles the
/// full-row-rank incidence matrix.
fn assemble_incidence(td: &mut TDecomp, pcols: &[usize], cols: usize) -> BitMatrix {
    // Endpoints (vertex pair) of the edge realizing each matrix column.
    // Zero columns are realized by self-loops and stay `None` (all-zero).
    let mut endpoints: Vec<Option<(usize, usize)>> = vec![None; cols];
    let mut num_vertices = td.extract_endpoints(pcols, &mut endpoints);

    // Tree edges that occur in no circuit are isolated bridges: realize
    // each as an edge between two fresh vertices of its own.
    for i in 0..pcols.len() {
        if td.tree_edge[i].is_none() {
            endpoints[pcols[i]] = Some((num_vertices, num_vertices + 1));
            num_vertices += 2;
        }
    }

    incidence_with_full_row_rank(num_vertices, &endpoints)
}

/// Splits the circuits into 1-connected blocks (port of Java `MatrixDivide`):
/// two circuits belong to the same block when they are linked by a chain of
/// shared tree edges. Returns the blocks as lists of circuit indices, each in
/// the breadth-first order expected by [`TDecomp::realization`].
fn decompose_blocks(rank: usize, circuits: &[(usize, Vec<usize>)]) -> Vec<Vec<usize>> {
    let mut tree_to_circuits: Vec<Vec<usize>> = vec![Vec::new(); rank];
    for (ci, (_, circuit)) in circuits.iter().enumerate() {
        for &t in circuit {
            tree_to_circuits[t].push(ci);
        }
    }

    let mut tree_seen = vec![false; rank];
    let mut circuit_seen = vec![false; circuits.len()];
    let mut blocks = Vec::new();

    for start in 0..circuits.len() {
        if circuit_seen[start] {
            continue;
        }
        circuit_seen[start] = true;
        // `block` doubles as the BFS queue: elements before `pos` are processed
        let mut block = vec![start];
        let mut pos = 0;
        while pos < block.len() {
            let ci = block[pos];
            pos += 1;
            for &t in &circuits[ci].1 {
                if tree_seen[t] {
                    continue;
                }
                tree_seen[t] = true;
                for &cj in &tree_to_circuits[t] {
                    if !circuit_seen[cj] {
                        circuit_seen[cj] = true;
                        block.push(cj);
                    }
                }
            }
        }
        blocks.push(block);
    }
    blocks
}

/// Builds the incidence matrix of the realized graph from the per-column
/// endpoint pairs, dropping one (redundant) vertex row per connected component
/// so that the result has full row rank.
fn incidence_with_full_row_rank(
    num_vertices: usize,
    endpoints: &[Option<(usize, usize)>],
) -> BitMatrix {
    fn find(parent: &mut [usize], mut v: usize) -> usize {
        while parent[v] != v {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        v
    }

    // union-find over vertices to identify the connected components
    let mut parent: Vec<usize> = (0..num_vertices).collect();
    for &(u, v) in endpoints.iter().flatten() {
        let (ru, rv) = (find(&mut parent, u), find(&mut parent, v));
        if ru != rv {
            parent[ru] = rv;
        }
    }

    // drop each component's union-find root vertex; the incidence rows of one
    // component sum to zero, so the remaining rows span the same space
    let mut row_of: Vec<Option<usize>> = vec![None; num_vertices];
    let mut rows = 0;
    for (v, row) in row_of.iter_mut().enumerate() {
        if find(&mut parent, v) != v {
            *row = Some(rows);
            rows += 1;
        }
    }

    let mut result = BitMatrix::zeros(rows, endpoints.len());
    for (j, ep) in endpoints.iter().enumerate() {
        if let Some((u, v)) = *ep {
            if let Some(ru) = row_of[u] {
                result.set_bit(ru, j, true);
            }
            if let Some(rv) = row_of[v] {
                result.set_bit(rv, j, true);
            }
        }
    }
    result
}

type EdgeId = usize;
type NodeId = usize;
type MemberId = usize;
type MarkerId = usize;

/// Marker orientation constants (Java `Marker.head` / `Marker.tail`)
const HEAD: bool = false;
const TAIL: bool = true;

/// The three kinds of member graphs in a t-decomposition (Java
/// `Member.polygon` / `Member.bond` / `Member.prime`)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum GraphType {
    Polygon,
    Bond,
    Prime,
}

/// An edge of the graph under construction (Java class `Edge`). Real edges
/// realize a tree edge or cotree edge of the matroid; virtual marker edges
/// (with `marker` set, `name == -1` in Java) glue members together.
#[derive(Clone, Copy, Debug, Default)]
struct EdgeData {
    /// owning member, before union-find resolution (Java `memberName`)
    member: Option<MemberId>,
    head: Option<NodeId>,
    tail: Option<NodeId>,
    /// the marker this edge is half of, if it is a virtual edge (Java `markerName`)
    marker: Option<MarkerId>,
    path_edge_flag: bool,
}

/// A vertex of the graph under construction (Java class `Node`). In a polygon
/// member, `head`/`tail` link the node into the circular edge list.
#[derive(Clone, Debug)]
struct NodeData {
    uf_parent: NodeId,
    uf_size: usize,
    /// final vertex number, assigned during endpoint extraction (Java `name`)
    name: Option<usize>,
    /// scratch counter used when typing prime members (Java `degree`)
    degree: usize,
    /// the edge this node is the head of (Java `head`)
    head: Option<EdgeId>,
    /// the edge this node is the tail of (Java `tail`)
    tail: Option<EdgeId>,
    /// scratch list of incident path edges used when typing prime members
    /// (Java `endEdgeArrayList`)
    end_edges: Vec<EdgeId>,
}

/// A member of the t-decomposition (Java class `Member`)
#[derive(Clone, Debug)]
struct MemberData {
    uf_parent: MemberId,
    uf_size: usize,
    /// scratch depth used by `make_depth_partition` (Java `depth`, -1 = unset)
    depth: Option<usize>,
    graph_type: GraphType,
    /// number of edges, maintained for bond members only (Java `edgeSize`)
    edge_size: isize,
    /// the marker joining this member to its parent member (Java `parentMarker`)
    parent_marker: Option<MarkerId>,
    /// path edges of the circuit currently being added (Java `pathEdge`)
    path_edges: Vec<EdgeId>,
    /// child markers typed 2 or 3 for the current circuit (Java `type2or3`)
    type2or3: Vec<EdgeId>,
    /// child markers typed 4 for the current circuit (Java `type4`)
    type4: Vec<EdgeId>,
    /// scratch flag used by `find_root_pair` (Java `findRootFlag`)
    find_root_flag: bool,
    /// parent marker edge of the root member, which has no parent marker
    /// (Java `rootParentMarkerEdge`)
    root_parent_marker_edge: Option<EdgeId>,
}

/// A virtual edge pair gluing a child member to its parent member (Java class
/// `Marker`)
#[derive(Clone, Copy, Debug)]
struct MarkerData {
    /// the virtual edge lying in the parent member
    parent: EdgeId,
    /// the virtual edge lying in the child member
    child: EdgeId,
    /// false: head joins head; true: head joins tail (Java `orient`)
    orient: bool,
    child_orient: bool,
}

/// The t-decomposition (Java class `TDecomposition`), with all `Edge`, `Node`,
/// `Member` and `Marker` objects arena-allocated and referenced by index.
struct TDecomp {
    edges: Vec<EdgeData>,
    nodes: Vec<NodeData>,
    members: Vec<MemberData>,
    markers: Vec<MarkerData>,
    /// edge realizing tree edge `i` (pivot column `pcols[i]`), once added
    /// (Java `treeEdge`, 1-based there)
    tree_edge: Vec<Option<EdgeId>>,
    /// edges realizing cotree columns, in insertion order (Java `coTreeEdge`)
    cotree_edges: Vec<EdgeId>,
    /// original matrix column of each cotree edge, parallel to `cotree_edges`
    cotree_cols: Vec<usize>,
    /// members located by rules 1-3 while typing the current circuit (Java `k`)
    k: [Option<MemberId>; 2],
    /// nodes located by rules 1-3 while typing the current circuit (Java `u`)
    u: [Option<NodeId>; 2],
}

impl TDecomp {
    fn new(rank: usize) -> Self {
        TDecomp {
            edges: Vec::new(),
            nodes: Vec::new(),
            members: Vec::new(),
            markers: Vec::new(),
            tree_edge: vec![None; rank],
            cotree_edges: Vec::new(),
            cotree_cols: Vec::new(),
            k: [None, None],
            u: [None, None],
        }
    }

    fn new_edge(&mut self) -> EdgeId {
        self.edges.push(EdgeData::default());
        self.edges.len() - 1
    }

    fn new_node(&mut self) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(NodeData {
            uf_parent: id,
            uf_size: 1,
            name: None,
            degree: 0,
            head: None,
            tail: None,
            end_edges: Vec::new(),
        });
        id
    }

    /// Java `new Marker()`: allocates the parent/child virtual edge pair
    fn new_marker(&mut self) -> MarkerId {
        let parent = self.new_edge();
        let child = self.new_edge();
        let id = self.markers.len();
        self.edges[parent].marker = Some(id);
        self.edges[child].marker = Some(id);
        self.markers.push(MarkerData {
            parent,
            child,
            orient: false,
            child_orient: false,
        });
        id
    }

    /// Java `new Member()` (the empty-polygon constructor), also used as the
    /// common setup of the other Member constructors
    fn new_member(&mut self, graph_type: GraphType) -> MemberId {
        let id = self.members.len();
        self.members.push(MemberData {
            uf_parent: id,
            uf_size: 1,
            depth: None,
            graph_type,
            edge_size: 0,
            parent_marker: None,
            path_edges: Vec::new(),
            type2or3: Vec::new(),
            type4: Vec::new(),
            find_root_flag: false,
            root_parent_marker_edge: None,
        });
        id
    }

    // ------------------------------------------------------------------
    // union-find (Java `Node.find/union`, `Member.find/union`)
    // ------------------------------------------------------------------

    /// Java `Node.find()`
    fn node_find(&mut self, n: NodeId) -> NodeId {
        let mut root = n;
        while self.nodes[root].uf_parent != root {
            root = self.nodes[root].uf_parent;
        }
        let mut cur = n;
        while cur != root {
            let next = self.nodes[cur].uf_parent;
            self.nodes[cur].uf_parent = root;
            cur = next;
        }
        root
    }

    /// Java `Node.union(Node)`
    fn node_union(&mut self, a: NodeId, b: NodeId) {
        let r1 = self.node_find(a);
        let r2 = self.node_find(b);
        if self.nodes[r1].uf_size > self.nodes[r2].uf_size {
            self.nodes[r2].uf_parent = r1;
            self.nodes[r1].uf_size += self.nodes[r2].uf_size;
        } else {
            self.nodes[r1].uf_parent = r2;
            self.nodes[r2].uf_size += self.nodes[r1].uf_size;
        }
    }

    /// Java `Member.find()`
    fn member_find(&mut self, m: MemberId) -> MemberId {
        let mut root = m;
        while self.members[root].uf_parent != root {
            root = self.members[root].uf_parent;
        }
        let mut cur = m;
        while cur != root {
            let next = self.members[cur].uf_parent;
            self.members[cur].uf_parent = root;
            cur = next;
        }
        root
    }

    /// Java `Member.union(Member, byte)`: merges two members; the union-find
    /// loser's `parent_marker` is discarded so the merged member keeps the
    /// parent marker of `b`'s side (the parent in every call site)
    fn member_union(&mut self, a: MemberId, b: MemberId, gtype: GraphType) {
        let m1 = self.member_find(a);
        let m2 = self.member_find(b);
        if self.members[m1].uf_size > self.members[m2].uf_size {
            self.members[m2].uf_parent = m1;
            self.members[m1].uf_size += self.members[m2].uf_size;
            self.members[m1].parent_marker = self.members[m2].parent_marker;
            match gtype {
                GraphType::Prime => self.members[m1].graph_type = GraphType::Prime,
                GraphType::Bond => self.members[m1].edge_size += self.members[m2].edge_size - 2,
                GraphType::Polygon => {}
            }
        } else {
            self.members[m1].uf_parent = m2;
            self.members[m2].uf_size += self.members[m1].uf_size;
            match gtype {
                GraphType::Prime => self.members[m2].graph_type = GraphType::Prime,
                GraphType::Bond => self.members[m2].edge_size += self.members[m1].edge_size - 2,
                GraphType::Polygon => {}
            }
        }
    }

    // ------------------------------------------------------------------
    // basic accessors
    // ------------------------------------------------------------------

    /// Java `Edge.getMember()`
    fn edge_member(&mut self, e: EdgeId) -> MemberId {
        let m = self.edges[e].member.expect("edge without member");
        self.member_find(m)
    }

    /// Java `Member.getParent()`
    fn member_parent(&mut self, m: MemberId) -> Option<MemberId> {
        let mk = self.members[m].parent_marker?;
        let pe = self.markers[mk].parent;
        Some(self.edge_member(pe))
    }

    /// Java `Member.parentMarkerEdge()`: the virtual edge representing the
    /// parent marker inside this member. `None` for a root member whose
    /// `root_parent_marker_edge` was never set (where Java returns null).
    fn parent_marker_edge(&self, m: MemberId) -> Option<EdgeId> {
        match self.members[m].parent_marker {
            Some(mk) => Some(self.markers[mk].child),
            None => self.members[m].root_parent_marker_edge,
        }
    }

    /// Java `Member.addPathEdge(Edge)`
    fn add_path_edge(&mut self, m: MemberId, e: EdgeId) {
        self.members[m].path_edges.push(e);
        self.edges[e].path_edge_flag = true;
    }

    /// Java `Marker.setChildOrient(boolean)`
    fn set_child_orient(&mut self, mk: MarkerId, b: bool) {
        self.markers[mk].child_orient = b;
    }

    /// Java `Marker.setParentOrient(boolean)`
    fn set_parent_orient(&mut self, mk: MarkerId, parent_orient: bool) {
        self.markers[mk].orient = self.markers[mk].child_orient ^ parent_orient;
    }

    /// Java `Marker.union(byte)`: contracts the marker, merging the child
    /// member into the parent and identifying the virtual edges' endpoints
    fn marker_union(&mut self, mk: MarkerId, gtype: GraphType) {
        let MarkerData {
            parent,
            child,
            orient,
            ..
        } = self.markers[mk];
        let cm = self.edges[child]
            .member
            .expect("marker child without member");
        let pm = self.edges[parent]
            .member
            .expect("marker parent without member");
        self.member_union(cm, pm, gtype);
        let ch = self.edges[child].head.unwrap();
        let ct = self.edges[child].tail.unwrap();
        let ph = self.edges[parent].head.unwrap();
        let pt = self.edges[parent].tail.unwrap();
        if orient {
            self.node_union(ch, pt);
            self.node_union(ct, ph);
        } else {
            self.node_union(ch, ph);
            self.node_union(ct, pt);
        }
    }

    // ------------------------------------------------------------------
    // edge and node operations (Java `Edge` / `Node` methods)
    // ------------------------------------------------------------------

    /// Java `Edge.nextEdge()`: the following edge in the circular polygon list
    fn next_edge(&self, e: EdgeId) -> EdgeId {
        self.nodes[self.edges[e].tail.unwrap()].tail.unwrap()
    }

    /// Java `Edge.preEdge()`: the preceding edge in the circular polygon list
    fn pre_edge(&self, e: EdgeId) -> EdgeId {
        self.nodes[self.edges[e].head.unwrap()].head.unwrap()
    }

    /// Java `Edge.swap(Edge)`: exchanges the positions of two edges of a polygon
    fn swap_edges(&mut self, a: EdgeId, b: EdgeId) {
        let a_h = self.edges[a].head.unwrap();
        let a_t = self.edges[a].tail.unwrap();
        let b_h = self.edges[b].head.unwrap();
        let b_t = self.edges[b].tail.unwrap();
        self.nodes[a_h].tail = Some(b);
        self.nodes[a_t].head = Some(b);
        self.nodes[b_h].tail = Some(a);
        self.nodes[b_t].head = Some(a);
        self.edges[b].head = Some(a_h);
        self.edges[b].tail = Some(a_t);
        self.edges[a].head = Some(b_h);
        self.edges[a].tail = Some(b_t);
    }

    /// Java `Edge.isEnd(Node, Node)`: whether the edge runs between `n1` and
    /// `n2` (which must be union-find roots)
    fn edge_is_end(&mut self, e: EdgeId, n1: NodeId, n2: NodeId) -> bool {
        let t = self.node_find(self.edges[e].tail.unwrap());
        let h = self.node_find(self.edges[e].head.unwrap());
        (t == n1 && h == n2) || (t == n2 && h == n1)
    }

    /// Java `Edge.isSameEnd(Edge)`: whether two edges are parallel
    fn edges_same_end(&mut self, a: EdgeId, b: EdgeId) -> bool {
        let a_t = self.node_find(self.edges[a].tail.unwrap());
        let a_h = self.node_find(self.edges[a].head.unwrap());
        let b_t = self.node_find(self.edges[b].tail.unwrap());
        let b_h = self.node_find(self.edges[b].head.unwrap());
        (a_t == b_t && a_h == b_h) || (a_t == b_h && a_h == b_t)
    }

    /// Java `Edge.setNextEdge(Edge)`
    fn set_next_edge(&mut self, e: EdgeId, target: EdgeId) {
        let t = self.edges[e].tail.unwrap();
        self.nodes[t].tail = Some(target);
        self.edges[target].head = Some(t);
    }

    /// Java `Edge.setNextEdge(Node, Edge)`
    fn set_next_edge_via(&mut self, e: EdgeId, target_node: NodeId, target_edge: EdgeId) {
        self.nodes[target_node].head = Some(e);
        self.nodes[target_node].tail = Some(target_edge);
        self.edges[target_edge].head = Some(target_node);
        self.edges[e].tail = Some(target_node);
    }

    /// Java `Edge.setPreEdge(Edge)`
    fn set_pre_edge(&mut self, e: EdgeId, target: EdgeId) {
        let h = self.edges[e].head.unwrap();
        self.nodes[h].head = Some(target);
        self.edges[target].tail = Some(h);
    }

    /// Java `Edge.setPreEdge(Node, Edge)`
    fn set_pre_edge_via(&mut self, e: EdgeId, target_node: NodeId, target_edge: EdgeId) {
        self.nodes[target_node].tail = Some(e);
        self.nodes[target_node].head = Some(target_edge);
        self.edges[target_edge].tail = Some(target_node);
        self.edges[e].head = Some(target_node);
    }

    /// Java `Node.isJoined(Node)`: whether two polygon nodes are adjacent
    fn nodes_joined(&self, n: NodeId, m: NodeId) -> bool {
        let h = self.nodes[n].head.unwrap();
        let t = self.nodes[n].tail.unwrap();
        self.edges[h].head == Some(m) || self.edges[t].tail == Some(m)
    }

    /// Java `Node.jointEdge(Node)`: an edge between two adjacent polygon nodes
    fn joint_edge(&self, n: NodeId, m: NodeId) -> Option<EdgeId> {
        let t = self.nodes[n].tail.unwrap();
        if self.edges[t].tail == Some(m) {
            return Some(t);
        }
        let h = self.nodes[n].head.unwrap();
        if self.edges[h].head == Some(m) {
            return Some(h);
        }
        None
    }

    /// Java `Node.isAnotherEnd(Node)`: follows the path of end edges starting
    /// at `start` and tests whether it terminates at `test` (all nodes are
    /// union-find roots)
    fn is_another_end(&mut self, start: NodeId, test: NodeId) -> bool {
        let mut wkn = start;
        let mut wke = self.nodes[wkn].end_edges[0];
        let h = self.node_find(self.edges[wke].head.unwrap());
        wkn = if h == wkn {
            self.node_find(self.edges[wke].tail.unwrap())
        } else {
            h
        };
        while self.nodes[wkn].end_edges.len() == 2 {
            wke = if self.nodes[wkn].end_edges[0] == wke {
                self.nodes[wkn].end_edges[1]
            } else {
                self.nodes[wkn].end_edges[0]
            };
            let h = self.node_find(self.edges[wke].head.unwrap());
            wkn = if h == wkn {
                self.node_find(self.edges[wke].tail.unwrap())
            } else {
                h
            };
        }
        wkn == test
    }

    // ------------------------------------------------------------------
    // member constructors (Java `Member` constructors)
    // ------------------------------------------------------------------

    /// Java `Member(Edge, ArrayList, Edge[], ArrayList)`: builds a new polygon
    /// member containing `marker_edge` plus a fresh edge for each element of
    /// `c` (elements `< rank` are tree edges, the rest is the new cotree edge)
    fn new_polygon_member(&mut self, marker_edge: EdgeId, c: &[usize]) -> MemberId {
        let member = self.new_member(GraphType::Polygon);
        self.edges[marker_edge].member = Some(member);
        self.members[member].parent_marker = self.edges[marker_edge].marker;
        let mut n = self.new_node();
        self.edges[marker_edge].head = Some(n);
        self.nodes[n].tail = Some(marker_edge);
        for &j in c {
            let wn = n;
            let e = self.new_edge();
            self.edges[e].member = Some(member);
            if j < self.tree_edge.len() {
                self.tree_edge[j] = Some(e);
            } else {
                self.cotree_edges.push(e);
            }
            n = self.new_node();
            self.edges[e].head = Some(n);
            self.edges[e].tail = Some(wn);
            self.nodes[wn].head = Some(e);
            self.nodes[n].tail = Some(e);
        }
        self.edges[marker_edge].tail = Some(n);
        self.nodes[n].head = Some(marker_edge);
        member
    }

    /// Java `Member(Edge, Edge, Edge)`: new 3-edge bond member with parent
    /// marker edge `f`
    fn new_bond3_member(&mut self, f: EdgeId, f1: EdgeId, f2: EdgeId) -> MemberId {
        let head_to_tail = if self.edges[f1].head.is_some() && self.edges[f2].head.is_some() {
            let h1 = self.node_find(self.edges[f1].head.unwrap());
            let t2 = self.node_find(self.edges[f2].tail.unwrap());
            h1 == t2
        } else {
            false
        };
        let n1 = self.new_node();
        let n2 = self.new_node();
        let member = self.new_member(GraphType::Bond);
        self.members[member].edge_size = 3;
        self.edges[f].member = Some(member);
        self.edges[f1].member = Some(member);
        self.edges[f2].member = Some(member);
        self.edges[f].head = Some(n1);
        self.edges[f].tail = Some(n2);
        self.edges[f1].head = Some(n1);
        self.edges[f1].tail = Some(n2);
        if head_to_tail {
            self.edges[f2].head = Some(n2);
            self.edges[f2].tail = Some(n1);
        } else {
            self.edges[f2].head = Some(n1);
            self.edges[f2].tail = Some(n2);
        }
        self.members[member].parent_marker = self.edges[f].marker;
        member
    }

    /// Java `Member(Edge, Edge)`: new 2-edge bond member rooted at `f1`
    fn new_bond2_member(&mut self, f1: EdgeId, f2: EdgeId) -> MemberId {
        let n1 = self.new_node();
        let n2 = self.new_node();
        let member = self.new_member(GraphType::Bond);
        self.members[member].edge_size = 2;
        self.edges[f1].member = Some(member);
        self.edges[f2].member = Some(member);
        self.edges[f1].head = Some(n1);
        self.edges[f2].head = Some(n1);
        self.edges[f1].tail = Some(n2);
        self.edges[f2].tail = Some(n2);
        self.members[member].root_parent_marker_edge = Some(f1);
        member
    }

    // ------------------------------------------------------------------
    // polygon surgery (Java `Member.divide/squeeze/rootSqueeze/makeNewPolygon`)
    // ------------------------------------------------------------------

    /// Java `Member.isParentMarker(Edge, Edge)`: whether the parent marker of
    /// `t`'s member lies strictly between `s` and `t` in polygon order
    fn is_parent_marker_between(&mut self, s: EdgeId, t: EdgeId) -> bool {
        let m = self.edge_member(t);
        let mk = self.members[m]
            .parent_marker
            .expect("polygon member without parent marker");
        let pm = self.markers[mk].child;
        let mut wke = self.next_edge(s);
        while wke != t {
            if wke == pm {
                return true;
            }
            wke = self.next_edge(wke);
        }
        false
    }

    /// Java `Member.makeNewPolygon(...)`: splits this polygon between nodes
    /// `h` and `t`, moving one side into `new_member` with `edge_in_new`
    /// closing the new polygon and `edge_in_old` closing the old one
    #[allow(clippy::too_many_arguments)]
    fn make_new_polygon(
        &mut self,
        old_member: MemberId,
        h: NodeId,
        t: NodeId,
        edge_in_new: EdgeId,
        edge_in_old: EdgeId,
        new_member: MemberId,
        change_end: bool,
    ) {
        let new_node1 = self.new_node();
        let new_node2 = self.new_node();
        let h1 = self.nodes[h].tail.unwrap();
        let t2 = self.nodes[h].head.unwrap();
        let t1 = self.nodes[t].head.unwrap();
        let h2 = self.nodes[t].tail.unwrap();
        let mut wke = h1;
        while wke != h2 {
            self.edges[wke].member = Some(new_member);
            wke = self.next_edge(wke);
        }
        if change_end {
            self.set_next_edge(t1, edge_in_new);
            self.set_pre_edge(h1, edge_in_new);
            self.set_next_edge_via(t2, new_node1, edge_in_old);
            self.set_pre_edge_via(h2, new_node2, edge_in_old);
        } else {
            self.set_next_edge_via(t1, new_node1, edge_in_new);
            self.set_pre_edge_via(h1, new_node2, edge_in_new);
            self.set_next_edge(t2, edge_in_old);
            self.set_pre_edge(h2, edge_in_old);
        }
        self.edges[edge_in_new].member = Some(new_member);
        self.edges[edge_in_old].member = Some(old_member);
    }

    /// Java `Member.divide(Node, Node, Edge)`: splits the polygon `m` at nodes
    /// `s` and `t` and hangs a new bond containing `f` between the two halves
    fn divide(&mut self, m: MemberId, s: NodeId, t: NodeId, f: EdgeId) {
        let new_member = self.new_member(GraphType::Polygon);
        let m1 = self.new_marker();
        let m2 = self.new_marker();
        // Java `isParentMarker(Node s, Node t)` forwards to `(s.head, t.tail)`
        let s_head = self.nodes[s].head.unwrap();
        let t_tail = self.nodes[t].tail.unwrap();
        if self.is_parent_marker_between(s_head, t_tail) {
            let (m2_parent, m2_child) = (self.markers[m2].parent, self.markers[m2].child);
            let (m1_parent, m1_child) = (self.markers[m1].parent, self.markers[m1].child);
            self.make_new_polygon(m, s, t, m2_parent, m1_child, new_member, true);
            self.members[new_member].parent_marker = self.members[m].parent_marker;
            self.members[m].parent_marker = Some(m1);
            self.new_bond3_member(m2_child, f, m1_parent);
        } else {
            let (m2_parent, m2_child) = (self.markers[m2].parent, self.markers[m2].child);
            let (m1_parent, m1_child) = (self.markers[m1].parent, self.markers[m1].child);
            self.make_new_polygon(m, s, t, m2_child, m1_parent, new_member, true);
            self.members[new_member].parent_marker = Some(m2);
            self.new_bond3_member(m1_child, f, m2_parent);
        }
    }

    /// Java `Member.squeeze(Node)`: shrinks the polygon `m` so that `n` is
    /// adjacent to the parent marker, splitting off the rest as new members
    fn squeeze_node(&mut self, m: MemberId, n: NodeId) {
        let p_marker = self.parent_marker_edge(m);
        if self.members[m].graph_type != GraphType::Polygon {
            return;
        }
        let p_marker = p_marker.expect("polygon member without parent marker edge");
        let tail_adjacent = self.edges[p_marker].tail == Some(n)
            || self.edges[self.next_edge(p_marker)].tail == Some(n);
        if !tail_adjacent {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let pm_tail = self.edges[p_marker].tail.unwrap();
            self.make_new_polygon(m, pm_tail, n, mk_child, mk_parent, new_member, false);
            self.members[new_member].parent_marker = Some(new_marker);
        }
        let head_adjacent = self.edges[p_marker].head == Some(n)
            || self.edges[self.pre_edge(p_marker)].head == Some(n);
        if !head_adjacent {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let pm_head = self.edges[p_marker].head.unwrap();
            self.make_new_polygon(m, pm_head, n, mk_parent, mk_child, new_member, true);
            self.members[new_member].parent_marker = self.members[m].parent_marker;
            self.members[m].parent_marker = Some(new_marker);
        }
    }

    /// Java `Member.rootSqueeze(Node, Edge)`: like `squeeze_node` but for the
    /// hypopath root member, squeezing between `n` and the edge `e1`
    fn root_squeeze_node(&mut self, m: MemberId, n: NodeId, e1: EdgeId) {
        if self.members[m].graph_type != GraphType::Polygon {
            return;
        }
        // Java `isParentMarker(Node s, Edge t)` forwards to `(s.head, t)`
        let n_head = self.nodes[n].head.unwrap();
        let is_parent = self.is_parent_marker_between(n_head, e1);
        let n_tail = self.nodes[n].tail.unwrap();
        if !(n_tail == e1 || self.next_edge(n_tail) == e1) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let e1_head = self.edges[e1].head.unwrap();
            if is_parent {
                self.make_new_polygon(m, n, e1_head, mk_parent, mk_child, new_member, false);
                self.members[new_member].parent_marker = self.members[m].parent_marker;
                self.members[m].parent_marker = Some(new_marker);
            } else {
                self.make_new_polygon(m, n, e1_head, mk_child, mk_parent, new_member, false);
                self.members[new_member].parent_marker = Some(new_marker);
            }
        }
        let n_head = self.nodes[n].head.unwrap();
        if !(n_head == e1 || self.pre_edge(n_head) == e1) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let e1_tail = self.edges[e1].tail.unwrap();
            if is_parent {
                self.make_new_polygon(m, n, e1_tail, mk_parent, mk_child, new_member, true);
                self.members[new_member].parent_marker = self.members[m].parent_marker;
                self.members[m].parent_marker = Some(new_marker);
            } else {
                self.make_new_polygon(m, n, e1_tail, mk_child, mk_parent, new_member, true);
                self.members[new_member].parent_marker = Some(new_marker);
            }
        }
    }

    /// Java `Member.squeeze(Edge)`: shrinks the polygon `m` so that the child
    /// marker `c_marker` is adjacent to the parent marker
    fn squeeze_edge(&mut self, m: MemberId, c_marker: EdgeId) {
        let p_marker = self.parent_marker_edge(m);
        if self.members[m].graph_type != GraphType::Polygon {
            return;
        }
        let p_marker = p_marker.expect("polygon member without parent marker edge");
        let next = self.next_edge(p_marker);
        if !(next == c_marker || self.next_edge(next) == c_marker) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let pm_tail = self.edges[p_marker].tail.unwrap();
            let cm_head = self.edges[c_marker].head.unwrap();
            self.make_new_polygon(m, pm_tail, cm_head, mk_child, mk_parent, new_member, false);
            self.members[new_member].parent_marker = Some(new_marker);
        }
        let pre = self.pre_edge(p_marker);
        if !(pre == c_marker || self.pre_edge(pre) == c_marker) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let pm_head = self.edges[p_marker].head.unwrap();
            let cm_tail = self.edges[c_marker].tail.unwrap();
            self.make_new_polygon(m, pm_head, cm_tail, mk_parent, mk_child, new_member, true);
            self.members[new_member].parent_marker = self.members[m].parent_marker;
            self.members[m].parent_marker = Some(new_marker);
        }
    }

    /// Java `Member.rootSqueeze(Edge, Edge)`: shrinks the hypopath root
    /// polygon so that the child markers `e1` and `e2` become adjacent
    fn root_squeeze_edges(&mut self, m: MemberId, e1: EdgeId, e2: EdgeId) {
        if self.members[m].graph_type != GraphType::Polygon {
            return;
        }
        let (mut wke1, mut wke2) = (e1, e2);
        while wke1 != e2 && wke2 != e1 {
            wke1 = self.next_edge(wke1);
            wke2 = self.next_edge(wke2);
        }
        let (s, t) = if wke1 == e2 { (e1, e2) } else { (e2, e1) };
        let is_parent = self.is_parent_marker_between(s, t);
        let s_next = self.next_edge(s);
        if !(s_next == t || self.next_edge(s_next) == t) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let s_tail = self.edges[s].tail.unwrap();
            let t_head = self.edges[t].head.unwrap();
            if is_parent {
                self.make_new_polygon(m, s_tail, t_head, mk_parent, mk_child, new_member, false);
                self.members[new_member].parent_marker = self.members[m].parent_marker;
                self.members[m].parent_marker = Some(new_marker);
            } else {
                self.make_new_polygon(m, s_tail, t_head, mk_child, mk_parent, new_member, false);
                self.members[new_member].parent_marker = Some(new_marker);
            }
        }
        let s_pre = self.pre_edge(s);
        if !(s_pre == t || self.pre_edge(s_pre) == t) {
            let new_member = self.new_member(GraphType::Polygon);
            let new_marker = self.new_marker();
            let (mk_parent, mk_child) = (
                self.markers[new_marker].parent,
                self.markers[new_marker].child,
            );
            let s_head = self.edges[s].head.unwrap();
            let t_tail = self.edges[t].tail.unwrap();
            if is_parent {
                self.make_new_polygon(m, s_head, t_tail, mk_parent, mk_child, new_member, true);
                self.members[new_member].parent_marker = self.members[m].parent_marker;
                self.members[m].parent_marker = Some(new_marker);
            } else {
                self.make_new_polygon(m, s_head, t_tail, mk_child, mk_parent, new_member, true);
                self.members[new_member].parent_marker = Some(new_marker);
            }
        }
    }
}

impl TDecomp {
    // ------------------------------------------------------------------
    // the Bixby-Wagner incremental realization (Java `TDecomposition`)
    // ------------------------------------------------------------------

    /// Java `TDecomposition.realization(ArrayList)`: realizes one block of
    /// circuits, given as `(original column, fundamental circuit)` pairs.
    /// On failure returns the position within `block` of the offending
    /// circuit; the t-decomposition is then left mid-mutation and must not be
    /// used for further insertions.
    fn realization(&mut self, block: &[(usize, &[usize])]) -> Result<(), usize> {
        let (first_col, first_circuit) = block[0];
        let first_edge = first_circuit[0];
        // seed bond: a dummy edge (Java `treeEdge[0]`) parallel to the first
        // tree edge of the block
        let dummy = self.new_edge();
        let fe = self.new_edge();
        self.tree_edge[first_edge] = Some(fe);
        self.new_bond2_member(dummy, fe);
        // Java ignores the result of the first addCircuit: with only the
        // seed bond present it cannot fail
        let _ = self.add_circuit(first_col, first_circuit);
        for (pos, &(col, circuit)) in block.iter().enumerate().skip(1) {
            if !self.add_circuit(col, circuit) {
                return Err(pos);
            }
        }
        Ok(())
    }

    /// Java `TDecomposition.addCircuit(int[])`
    fn add_circuit(&mut self, col: usize, s: &[usize]) -> bool {
        self.k = [None, None];
        self.u = [None, None];
        let mut p: Vec<usize> = Vec::new();
        let mut c: Vec<usize> = Vec::new();
        for &si in s {
            match self.tree_edge[si] {
                Some(e) => {
                    p.push(si);
                    self.edges[e].path_edge_flag = true;
                }
                None => c.push(si),
            }
        }
        // the new cotree edge, numbered past the tree edges as in Java
        c.push(self.tree_edge.len() + self.cotree_cols.len());
        self.cotree_cols.push(col);
        self.hypopath(&p) && self.update(&c)
    }

    /// Java `TDecomposition.hypopath(ArrayList)`: checks that the tree edges
    /// `p` can be rearranged into a (hypo)path of the graph built so far,
    /// typing every member holding path edges bottom-up
    fn hypopath(&mut self, p: &[usize]) -> bool {
        let mut path_members: Vec<MemberId> = Vec::new();
        for &pi in p {
            let e = self.tree_edge[pi].unwrap();
            let wkm = self.edge_member(e);
            if self.members[wkm].path_edges.is_empty() {
                path_members.push(wkm);
            }
            self.add_path_edge(wkm, e);
        }

        let root = self.find_root_list(&path_members);
        let mut partition = self.make_depth_partition(&path_members, root);

        let mut depth = partition.len();
        while depth > 1 {
            depth -= 1;
            while let Some(wkm) = partition[depth].pop() {
                let type2or3_size = self.members[wkm].type2or3.len();
                let type4_size = self.members[wkm].type4.len();
                if type4_size > 1 || type2or3_size > 2 || (type4_size == 1 && type2or3_size > 0) {
                    return false;
                }
                if !self.typing(wkm) {
                    return false;
                }
                self.clear_path_state(wkm);
            }
        }
        let type2or3_size = self.members[root].type2or3.len();
        let type4_size = self.members[root].type4.len();
        if type4_size > 1 || type2or3_size > 2 || (type4_size == 1 && type2or3_size > 0) {
            return false;
        }
        if !self.is_path(root) {
            return false;
        }
        self.clear_path_state(root);
        self.rule5();
        true
    }

    /// Clears the per-circuit scratch state of a member after typing (the
    /// trailing loops of Java `hypopath`)
    fn clear_path_state(&mut self, m: MemberId) {
        let path_edges = std::mem::take(&mut self.members[m].path_edges);
        for e in path_edges {
            self.edges[e].path_edge_flag = false;
        }
        self.members[m].type2or3.clear();
        self.members[m].type4.clear();
    }

    /// Java `TDecomposition.findRoot(ArrayList)`: nearest common ancestor of a
    /// list of members, computed by pairwise reduction
    fn find_root_list(&mut self, p: &[MemberId]) -> MemberId {
        let mut queue: std::collections::VecDeque<MemberId> = p.iter().copied().collect();
        while queue.len() > 1 {
            let m0 = queue.pop_front().unwrap();
            let m1 = queue.pop_front().unwrap();
            let root = self.find_root_pair(m0, m1);
            queue.push_back(root);
        }
        queue
            .pop_front()
            .expect("hypopath must contain at least one member")
    }

    /// Java `TDecomposition.findRoot(Member, Member)`: nearest common ancestor
    /// of two members, found by alternately walking up from both and flagging
    fn find_root_pair(&mut self, m1: MemberId, m2: MemberId) -> MemberId {
        let mut wk = [Some(m1), Some(m2)];
        let mut visited: Vec<MemberId> = Vec::new();
        let mut result: Option<MemberId> = None;

        'phase1: while wk[0].is_some() && wk[1].is_some() {
            for slot in &mut wk {
                let m = slot.unwrap();
                if self.members[m].find_root_flag {
                    result = Some(m);
                    break 'phase1;
                }
                self.members[m].find_root_flag = true;
                visited.push(m);
                *slot = self.member_parent(m);
            }
        }
        if result.is_none() {
            'phase2: for slot in &mut wk {
                while let Some(m) = *slot {
                    if self.members[m].find_root_flag {
                        result = Some(m);
                        break 'phase2;
                    }
                    self.members[m].find_root_flag = true;
                    visited.push(m);
                    *slot = self.member_parent(m);
                }
            }
        }
        for m in visited {
            self.members[m].find_root_flag = false;
        }
        let root = result.expect("members of one block must share a root");
        self.member_find(root)
    }

    /// Java `TDecomposition.makeDepthPartition(ArrayList, Member)`: groups the
    /// path members (and their ancestors up to `root`) by depth below `root`
    fn make_depth_partition(&mut self, path: &[MemberId], root: MemberId) -> Vec<Vec<MemberId>> {
        let mut partition: Vec<Vec<MemberId>> = Vec::new();
        let mut visited: Vec<MemberId> = Vec::new();
        let mut copy: Vec<MemberId> = path.to_vec();
        let mut depth_stack: Vec<MemberId> = Vec::new();

        self.members[root].depth = Some(0);
        visited.push(root);
        partition.push(vec![root]);
        let mut max_depth = 0;

        while let Some(mut wk) = copy.pop() {
            while self.members[wk].depth.is_none() {
                depth_stack.push(wk);
                wk = self
                    .member_parent(wk)
                    .expect("path member must have an ancestor below the root");
            }
            let mut depth = self.members[wk].depth.unwrap();
            while let Some(w) = depth_stack.pop() {
                depth += 1;
                if depth > max_depth {
                    partition.push(Vec::new());
                }
                self.members[w].depth = Some(depth);
                visited.push(w);
                partition[depth].push(w);
            }
            if depth > max_depth {
                max_depth = depth;
            }
        }
        for m in visited {
            self.members[m].depth = None;
        }
        partition
    }

    // ------------------------------------------------------------------
    // rules for locating the attachment members/nodes of the new edge
    // ------------------------------------------------------------------

    /// Java `TDecomposition.rule1`
    fn rule1(&mut self, m: MemberId, n: NodeId) -> bool {
        if self.u[0].is_none() {
            self.k[0] = Some(m);
            self.u[0] = Some(n);
        } else if self.u[1].is_none() {
            self.k[1] = Some(m);
            self.u[1] = Some(n);
        } else {
            return false;
        }
        true
    }

    /// Java `TDecomposition.rule2`
    fn rule2(&mut self, m: MemberId, n1: NodeId, n2: NodeId) -> bool {
        if self.u[0].is_some() {
            false
        } else {
            self.k[0] = Some(m);
            self.k[1] = Some(m);
            self.u[0] = Some(n1);
            self.u[1] = Some(n2);
            true
        }
    }

    /// Java `TDecomposition.rule3`: walks the marker chain from `k[0]` up to
    /// `m`, following the node `n` through each marker identification
    fn rule3(&mut self, m: MemberId, n: NodeId) -> bool {
        if self.u[1].is_some() {
            return false;
        }
        let mut test_members: Vec<MemberId> = Vec::new();
        let mut wkm = self.k[0].expect("rule3 requires k[0]");
        let mut wkn = n;
        while wkm != m {
            test_members.push(wkm);
            wkm = self.member_parent(wkm).expect("rule3 walk must reach m");
        }
        while let Some(w) = test_members.pop() {
            wkm = w;
            let mk = self.members[wkm].parent_marker.unwrap();
            let parent_edge = self.markers[mk].parent;
            let child_edge = self.parent_marker_edge(wkm).unwrap();
            let ph = self.node_find(self.edges[parent_edge].head.unwrap());
            if ph == wkn {
                wkn = if self.markers[mk].orient {
                    self.node_find(self.edges[child_edge].tail.unwrap())
                } else {
                    self.node_find(self.edges[child_edge].head.unwrap())
                };
                continue;
            }
            let pt = self.node_find(self.edges[parent_edge].tail.unwrap());
            if pt == wkn {
                wkn = if self.markers[mk].orient {
                    self.node_find(self.edges[child_edge].head.unwrap())
                } else {
                    self.node_find(self.edges[child_edge].tail.unwrap())
                };
            } else {
                self.k[1] = self.member_parent(wkm);
                self.u[1] = Some(wkn);
                return true;
            }
        }
        self.k[1] = Some(wkm);
        self.u[1] = Some(wkn);
        true
    }

    /// Java `TDecomposition.rule5`: moves the attachment pair `k`/`u` off a
    /// parent marker (or onto an adjacent bond) where possible
    fn rule5(&mut self) {
        if self.k[0] != self.k[1] {
            return;
        }
        let Some(k0) = self.k[0] else { return };
        let u0 = self.u[0].expect("rule5 requires u[0]");
        let u1 = self.u[1].expect("rule5 requires u[1]");
        let pme = self
            .parent_marker_edge(k0)
            .expect("rule5 member without parent marker edge");
        if self.edge_is_end(pme, u0, u1) && self.members[k0].graph_type != GraphType::Bond {
            let mk = self.members[k0].parent_marker.unwrap();
            let parent_edge = self.markers[mk].parent;
            self.u[0] = Some(self.node_find(self.edges[parent_edge].tail.unwrap()));
            self.u[1] = Some(self.node_find(self.edges[parent_edge].head.unwrap()));
            self.k[0] = self.member_parent(k0);
            self.k[1] = self.k[0];
        } else if self.members[k0].graph_type == GraphType::Polygon && self.nodes_joined(u0, u1) {
            let joint = self.joint_edge(u0, u1).unwrap();
            // Java: `name == -1` (a virtual marker edge) whose child member is a bond
            if let Some(mk) = self.edges[joint].marker {
                let child = self.markers[mk].child;
                let child_member = self.edge_member(child);
                if self.members[child_member].graph_type == GraphType::Bond {
                    self.u[0] = Some(self.node_find(self.edges[child].head.unwrap()));
                    self.u[1] = Some(self.node_find(self.edges[child].tail.unwrap()));
                    self.k[0] = Some(child_member);
                    self.k[1] = self.k[0];
                }
            }
        }
    }

    /// Java `TDecomposition.update(ArrayList)`: inserts the new cotree edge
    /// (and any tree edges not yet present, as a polygon hanging off it)
    /// between the attachment nodes `u[0]`/`u[1]` found by the rules
    fn update(&mut self, c: &[usize]) -> bool {
        let f: EdgeId;
        if c.len() == 1 {
            f = self.new_edge();
            self.cotree_edges.push(f);
        } else {
            let mk = self.new_marker();
            f = self.markers[mk].parent;
            let child = self.markers[mk].child;
            self.new_polygon_member(child, c);
        }
        let k0 = self.k[0].expect("rules must set k[0] before update");
        let k1 = self.k[1].expect("rules must set k[1] before update");
        if k0 == k1 {
            if self.members[k0].graph_type != GraphType::Polygon {
                if self.members[k0].graph_type == GraphType::Bond {
                    self.members[k0].edge_size += 1;
                    let wke = self
                        .parent_marker_edge(k0)
                        .expect("bond without parent marker edge");
                    let h = self.node_find(self.edges[wke].head.unwrap());
                    let t = self.node_find(self.edges[wke].tail.unwrap());
                    self.u[0] = Some(h);
                    self.u[1] = Some(t);
                }
                // f.setEndNode(u[0], u[1]) — (head, tail)
                self.edges[f].head = self.u[0];
                self.edges[f].tail = self.u[1];
                self.edges[f].member = Some(k0);
            } else {
                let u0 = self.u[0].unwrap();
                let u1 = self.u[1].unwrap();
                if self.nodes_joined(u0, u1) {
                    // replace the edge joining u[0] and u[1] by a marker whose
                    // child is a new bond of that edge and f
                    let f1 = self.joint_edge(u0, u1).unwrap();
                    let m2 = self.new_marker();
                    let m2_parent = self.markers[m2].parent;
                    let m2_child = self.markers[m2].child;
                    self.edges[m2_parent].member = Some(k0);
                    let f1_head = self.edges[f1].head.unwrap();
                    let f1_tail = self.edges[f1].tail.unwrap();
                    self.nodes[f1_head].tail = Some(m2_parent);
                    self.nodes[f1_tail].head = Some(m2_parent);
                    self.edges[m2_parent].tail = Some(f1_tail);
                    self.edges[m2_parent].head = Some(f1_head);
                    self.new_bond3_member(m2_child, f, f1);
                } else {
                    self.divide(k0, u1, u0, f);
                }
            }
        } else {
            let first_root = self.find_root_pair(k0, k1);
            let mut r1: Vec<MemberId> = Vec::new();
            let mut r2: Vec<MemberId> = Vec::new();
            let mut wkm = k0;
            while wkm != first_root {
                r1.push(wkm);
                wkm = self.member_parent(wkm).expect("walk to root");
            }
            wkm = k1;
            while wkm != first_root {
                r2.push(wkm);
                wkm = self.member_parent(wkm).expect("walk to root");
            }
            let mut root = first_root;
            let wkm_a0 = *r1.last().expect("k[0] must lie strictly below the root");
            let mut e1 = self.markers[self.members[wkm_a0].parent_marker.unwrap()].parent;
            if root != k0 && root != k1 {
                let wkm_a1 = *r2.last().unwrap();
                let mut e2 = self.markers[self.members[wkm_a1].parent_marker.unwrap()].parent;
                if self.members[root].graph_type == GraphType::Prime && self.edges_same_end(e1, e2)
                {
                    // e1 and e2 are parallel in a prime member: bundle them
                    // into a new bond hanging off a fresh marker
                    let m1 = self.new_marker();
                    let f1 = self.markers[m1].parent;
                    self.edges[f1].member = Some(root);
                    let e1_tail = self.node_find(self.edges[e1].tail.unwrap());
                    let e1_head = self.node_find(self.edges[e1].head.unwrap());
                    self.edges[f1].tail = Some(e1_tail);
                    self.edges[f1].head = Some(e1_head);
                    let m1_child = self.markers[m1].child;
                    root = self.new_bond3_member(m1_child, e1, e2);
                    if self.members[wkm_a0].graph_type == GraphType::Bond {
                        let mk = self.members[wkm_a0].parent_marker.unwrap();
                        self.marker_union(mk, GraphType::Bond);
                        r1.pop();
                        root = self.member_find(root);
                        let last = *r1.last().unwrap();
                        e1 = self.markers[self.members[last].parent_marker.unwrap()].parent;
                    }
                    if self.members[wkm_a1].graph_type == GraphType::Bond {
                        let mk = self.members[wkm_a1].parent_marker.unwrap();
                        self.marker_union(mk, GraphType::Bond);
                        r2.pop();
                        root = self.member_find(root);
                        let last = *r2.last().unwrap();
                        e2 = self.markers[self.members[last].parent_marker.unwrap()].parent;
                    }
                }
                if self.members[root].graph_type == GraphType::Bond
                    && self.members[root].edge_size > 3
                {
                    // split e1 and e2 off a large bond into a 3-edge bond
                    self.members[root].edge_size -= 1;
                    let m1 = self.new_marker();
                    let f1 = self.markers[m1].parent;
                    self.edges[f1].member = Some(root);
                    let e1_tail = self.node_find(self.edges[e1].tail.unwrap());
                    let e1_head = self.node_find(self.edges[e1].head.unwrap());
                    self.edges[f1].tail = Some(e1_tail);
                    self.edges[f1].head = Some(e1_head);
                    let m1_child = self.markers[m1].child;
                    root = self.new_bond3_member(m1_child, e1, e2);
                }
                self.root_squeeze_edges(root, e1, e2);
            }
            let root = self.find_root_pair(k0, k1);
            let mut marker_stack: Vec<MarkerId> = Vec::new();
            let mut wkm = k0;
            while wkm != root {
                marker_stack.push(self.members[wkm].parent_marker.unwrap());
                wkm = self.member_parent(wkm).expect("walk to root");
            }
            wkm = k1;
            while wkm != root {
                marker_stack.push(self.members[wkm].parent_marker.unwrap());
                wkm = self.member_parent(wkm).expect("walk to root");
            }
            let u0 = self.u[0].unwrap();
            let u1 = self.u[1].unwrap();
            self.squeeze_node(k0, u0);
            if k1 == root {
                self.root_squeeze_node(k1, u1, e1);
            } else {
                self.squeeze_node(k1, u1);
            }
            let mut wke = self.markers[self.members[k0].parent_marker.unwrap()].parent;
            for &wkm in r1.iter().skip(1) {
                self.squeeze_edge(wkm, wke);
                wke = self.markers[self.members[wkm].parent_marker.unwrap()].parent;
            }
            let mut wke = self.markers[self.members[k1].parent_marker.unwrap()].parent;
            for &wkm in r2.iter().skip(1) {
                self.squeeze_edge(wkm, wke);
                wke = self.markers[self.members[wkm].parent_marker.unwrap()].parent;
            }
            while let Some(mk) = marker_stack.pop() {
                self.marker_union(mk, GraphType::Prime);
            }
            self.edges[f].head = self.u[0];
            self.edges[f].tail = self.u[1];
            let f_member = self.edge_member(e1);
            self.edges[f].member = Some(f_member);
        }
        true
    }
}

impl TDecomp {
    /// Java `TDecomposition.typing(Member)`: determines the type (1-4) of a
    /// non-root member with respect to the current path, rearranging its
    /// edges and propagating the result to its parent member
    fn typing(&mut self, m: MemberId) -> bool {
        let marker = self
            .parent_marker_edge(m)
            .expect("typed member without parent marker");
        let parent_marker = self.members[m].parent_marker.unwrap();
        let pm_parent_edge = self.markers[parent_marker].parent;
        let parent = self.member_parent(m).expect("typed member without parent");

        match self.members[m].graph_type {
            GraphType::Polygon => {
                // gather the path edges into a run following the marker
                let mut wke = self.next_edge(marker);
                let mut path_edge_number = self.members[m].path_edges.len();
                while path_edge_number > 0 {
                    if self.edges[wke].path_edge_flag {
                        self.edges[wke].path_edge_flag = false;
                        wke = self.next_edge(wke);
                        path_edge_number -= 1;
                        continue;
                    }
                    let tedge = self.members[m].path_edges.remove(0);
                    if !self.edges[tedge].path_edge_flag {
                        continue;
                    }
                    self.swap_edges(wke, tedge);
                    wke = tedge;
                    self.edges[wke].path_edge_flag = false;
                    wke = self.next_edge(wke);
                    path_edge_number -= 1;
                }
                if wke == marker {
                    // type 1
                    self.add_path_edge(parent, pm_parent_edge);
                } else if self.members[m].type4.len() == 1 {
                    // type 4
                    if self.next_edge(wke) != marker {
                        return false;
                    }
                    self.members[parent].type4.push(pm_parent_edge);
                } else if self.members[m].type2or3.len() == 2 {
                    // type 4
                    let end0 = self.members[m].type2or3[0];
                    self.swap_edges(wke, end0);
                    let mk0 = self.edges[end0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let end1 = self.members[m].type2or3[1];
                    let pre = self.pre_edge(marker);
                    self.swap_edges(end1, pre);
                    let mk1 = self.edges[end1].marker.unwrap();
                    self.set_parent_orient(mk1, TAIL);
                    self.members[parent].type4.push(pm_parent_edge);
                } else if self.members[m].type2or3.len() == 1 {
                    // type 2 or 3
                    let end0 = self.members[m].type2or3[0];
                    self.swap_edges(wke, end0);
                    let mk0 = self.edges[end0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    self.set_child_orient(parent_marker, TAIL);
                    self.members[parent].type2or3.push(pm_parent_edge);
                } else {
                    // type 2 or 3
                    self.set_child_orient(parent_marker, TAIL);
                    self.members[parent].type2or3.push(pm_parent_edge);
                    let n = self.edges[wke].head.unwrap();
                    if !self.rule1(m, n) {
                        return false;
                    }
                }
                true
            }

            GraphType::Bond => {
                if self.members[m].type4.len() == 1 {
                    // type 4
                    if !self.members[m].path_edges.is_empty() {
                        return false;
                    }
                    self.members[parent].type4.push(pm_parent_edge);
                } else if self.members[m].type2or3.len() == 2 {
                    // type 4
                    if !self.members[m].path_edges.is_empty() {
                        return false;
                    }
                    let end0 = self.members[m].type2or3[0];
                    let wke = self.members[m].type2or3[1];
                    let mk0 = self.edges[end0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let h0 = self.node_find(self.edges[end0].head.unwrap());
                    let h1 = self.node_find(self.edges[wke].head.unwrap());
                    let mk1 = self.edges[wke].marker.unwrap();
                    if h0 == h1 {
                        self.set_parent_orient(mk1, TAIL);
                    } else {
                        self.set_parent_orient(mk1, HEAD);
                    }
                    self.members[parent].type4.push(pm_parent_edge);
                } else if self.members[m].type2or3.len() == 1 {
                    // type 2 or 3
                    if self.members[m].path_edges.len() > 1 {
                        return false;
                    }
                    let end0 = self.members[m].type2or3[0];
                    let mk0 = self.edges[end0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let marker_head = self.node_find(self.edges[marker].head.unwrap());
                    let end0_head = self.node_find(self.edges[end0].head.unwrap());
                    if self.members[m].path_edges.len() == 1 {
                        if marker_head == end0_head {
                            self.set_child_orient(parent_marker, TAIL);
                        } else {
                            self.set_child_orient(parent_marker, HEAD);
                        }
                    } else if marker_head == end0_head {
                        self.set_child_orient(parent_marker, HEAD);
                    } else {
                        self.set_child_orient(parent_marker, TAIL);
                    }
                    self.members[parent].type2or3.push(pm_parent_edge);
                } else {
                    // type 1
                    if self.members[m].path_edges.len() != 1 {
                        return false;
                    }
                    self.add_path_edge(parent, pm_parent_edge);
                }
                true
            }

            GraphType::Prime => {
                let mut path_nodes: Vec<NodeId> = Vec::new();
                for i in 0..self.members[m].path_edges.len() {
                    let wke = self.members[m].path_edges[i];
                    let h = self.node_find(self.edges[wke].head.unwrap());
                    if self.nodes[h].end_edges.is_empty() {
                        path_nodes.push(h);
                    }
                    self.nodes[h].end_edges.push(wke);
                    let t = self.node_find(self.edges[wke].tail.unwrap());
                    if self.nodes[t].end_edges.is_empty() {
                        path_nodes.push(t);
                    }
                    self.nodes[t].end_edges.push(wke);
                }
                let mut path_ends: Vec<NodeId> = Vec::new();
                for &wkn in &path_nodes {
                    let count = self.nodes[wkn].end_edges.len();
                    if count > 2 {
                        return false;
                    }
                    if count == 1 {
                        path_ends.push(wkn);
                    }
                }
                let marker_head = self.node_find(self.edges[marker].head.unwrap());
                let marker_tail = self.node_find(self.edges[marker].tail.unwrap());
                match path_ends.len() {
                    0 => {
                        for &wkn in &path_nodes {
                            self.nodes[wkn].end_edges.clear();
                        }
                        self.add_type234(m, marker_head, marker_tail, TAIL, 0)
                    }
                    2 => {
                        for &wkn in &path_nodes {
                            self.nodes[wkn].end_edges.clear();
                        }
                        let end0 = path_ends[0];
                        let end1 = path_ends[1];
                        if self.edge_is_end(marker, end0, end1) {
                            self.add_type234(m, marker_head, marker_tail, TAIL, 2)
                        } else if marker_head == end0 {
                            self.add_type234(m, end1, marker_tail, TAIL, 2)
                        } else if marker_tail == end0 {
                            self.add_type234(m, end1, marker_head, HEAD, 2)
                        } else if marker_head == end1 {
                            self.add_type234(m, end0, marker_tail, TAIL, 2)
                        } else if marker_tail == end1 {
                            self.add_type234(m, end0, marker_head, HEAD, 2)
                        } else {
                            false
                        }
                    }
                    4 => {
                        let mut end = [path_ends[0], path_ends[1], path_ends[2], path_ends[3]];
                        if self.is_another_end(end[0], end[2]) {
                            end.swap(1, 2);
                        }
                        if self.is_another_end(end[0], end[3]) {
                            end.swap(1, 3);
                        }
                        for &wkn in &path_nodes {
                            self.nodes[wkn].end_edges.clear();
                        }
                        if self.edge_is_end(marker, end[0], end[2]) {
                            self.add_type234(m, end[1], end[3], HEAD, 4)
                        } else if self.edge_is_end(marker, end[0], end[3]) {
                            self.add_type234(m, end[1], end[2], HEAD, 4)
                        } else if self.edge_is_end(marker, end[1], end[2]) {
                            self.add_type234(m, end[0], end[3], HEAD, 4)
                        } else if self.edge_is_end(marker, end[1], end[3]) {
                            self.add_type234(m, end[0], end[2], HEAD, 4)
                        } else {
                            false
                        }
                    }
                    _ => false,
                }
            }
        }
    }

    /// Java `TDecomposition.addType234(...)`: records a member typed 2, 3 or 4
    /// in its parent, orienting any child markers already typed 2 or 3 so that
    /// they line up with the path ends `target1`/`target2`
    fn add_type234(
        &mut self,
        m: MemberId,
        target1: NodeId,
        target2: NodeId,
        t2_type: bool,
        count: u8,
    ) -> bool {
        let parent_marker = self.members[m].parent_marker.unwrap();
        let pm_parent_edge = self.markers[parent_marker].parent;
        let parent = self.member_parent(m).expect("typed member without parent");

        if self.members[m].type4.is_empty() && self.members[m].type2or3.is_empty() {
            if count == 0 {
                return false;
            } else if count == 2 {
                let pme = self.parent_marker_edge(m).unwrap();
                if self.edge_is_end(pme, target1, target2) {
                    self.add_path_edge(parent, pm_parent_edge);
                } else {
                    self.members[parent].type2or3.push(pm_parent_edge);
                    self.set_child_orient(parent_marker, !t2_type);
                    if !self.rule1(m, target1) {
                        return false;
                    }
                }
            } else if count == 4 {
                self.members[parent].type4.push(pm_parent_edge);
                if !self.rule2(m, target1, target2) {
                    return false;
                }
            }
        } else if self.members[m].type4.len() == 1 {
            let t4 = self.members[m].type4[0];
            if self.edge_is_end(t4, target1, target2) {
                self.members[parent].type4.push(pm_parent_edge);
            } else {
                return false;
            }
        } else if self.members[m].type2or3.len() == 1 {
            let wke = self.members[m].type2or3[0];
            let mk = self.edges[wke].marker.unwrap();
            let wke_head = self.node_find(self.edges[wke].head.unwrap());
            let wke_tail = self.node_find(self.edges[wke].tail.unwrap());
            if wke_head == target1 || wke_tail == target1 {
                if wke_head == target1 {
                    self.set_parent_orient(mk, HEAD);
                } else {
                    self.set_parent_orient(mk, TAIL);
                }
                if count == 0 {
                    self.members[parent].type2or3.push(pm_parent_edge);
                    self.set_child_orient(parent_marker, !t2_type);
                } else if count == 2 {
                    self.members[parent].type2or3.push(pm_parent_edge);
                    let pme = self.parent_marker_edge(m).unwrap();
                    if self.edge_is_end(pme, target1, target2) {
                        self.set_child_orient(parent_marker, t2_type);
                    } else {
                        self.set_child_orient(parent_marker, !t2_type);
                    }
                } else if count == 4 {
                    self.members[parent].type4.push(pm_parent_edge);
                    if !self.rule3(m, target2) {
                        return false;
                    }
                }
            } else if wke_head == target2 || wke_tail == target2 {
                if wke_head == target2 {
                    self.set_parent_orient(mk, HEAD);
                }
                if wke_tail == target2 {
                    self.set_parent_orient(mk, TAIL);
                }
                if count == 0 {
                    self.members[parent].type2or3.push(pm_parent_edge);
                    self.set_child_orient(parent_marker, t2_type);
                } else if count == 2 {
                    let pme = self.parent_marker_edge(m).unwrap();
                    if self.edge_is_end(pme, target1, target2) {
                        self.members[parent].type2or3.push(pm_parent_edge);
                        self.set_child_orient(parent_marker, !t2_type);
                    } else {
                        self.members[parent].type4.push(pm_parent_edge);
                        if !self.rule3(m, target1) {
                            return false;
                        }
                    }
                } else if count == 4 {
                    self.members[parent].type4.push(pm_parent_edge);
                    if !self.rule3(m, target1) {
                        return false;
                    }
                }
            } else {
                return false;
            }
        } else if self.members[m].type2or3.len() == 2 {
            self.members[parent].type4.push(pm_parent_edge);
            let wke = [self.members[m].type2or3[0], self.members[m].type2or3[1]];
            for i in 0..2 {
                let head_i = self.node_find(self.edges[wke[i]].head.unwrap());
                let tail_i = self.node_find(self.edges[wke[i]].tail.unwrap());
                let mk_i = self.edges[wke[i]].marker.unwrap();
                let mk_other = self.edges[wke[1 - i]].marker.unwrap();
                let head_other = self.node_find(self.edges[wke[1 - i]].head.unwrap());
                let tail_other = self.node_find(self.edges[wke[1 - i]].tail.unwrap());
                if head_i == target1 || tail_i == target1 {
                    if head_i == target1 {
                        if tail_i == target2 {
                            continue;
                        }
                        self.set_parent_orient(mk_i, HEAD);
                    } else {
                        if head_i == target2 {
                            continue;
                        }
                        self.set_parent_orient(mk_i, TAIL);
                    }
                    if head_other == target2 {
                        self.set_parent_orient(mk_other, HEAD);
                    } else if tail_other == target2 {
                        self.set_parent_orient(mk_other, TAIL);
                    } else {
                        return false;
                    }
                    return true;
                } else if head_i == target2 || tail_i == target2 {
                    if head_i == target2 {
                        if tail_i == target1 {
                            continue;
                        }
                        self.set_parent_orient(mk_i, HEAD);
                    } else {
                        if head_i == target1 {
                            continue;
                        }
                        self.set_parent_orient(mk_i, TAIL);
                    }
                    if head_other == target1 {
                        self.set_parent_orient(mk_other, HEAD);
                    } else if tail_other == target1 {
                        self.set_parent_orient(mk_other, TAIL);
                    } else {
                        return false;
                    }
                    return true;
                } else {
                    return false;
                }
            }
            // both child markers span target1-target2: orient them oppositely
            let head0 = self.node_find(self.edges[wke[0]].head.unwrap());
            let head1 = self.node_find(self.edges[wke[1]].head.unwrap());
            let mk0 = self.edges[wke[0]].marker.unwrap();
            let mk1 = self.edges[wke[1]].marker.unwrap();
            if head0 == head1 {
                self.set_parent_orient(mk0, HEAD);
                self.set_parent_orient(mk1, TAIL);
            } else {
                self.set_parent_orient(mk0, HEAD);
                self.set_parent_orient(mk1, HEAD);
            }
        }
        true
    }

    /// Java `TDecomposition.isPath(Member)`: like `typing` but for the root
    /// member of the hypopath, where the path need not reach the parent marker
    fn is_path(&mut self, m: MemberId) -> bool {
        match self.members[m].graph_type {
            GraphType::Polygon => {
                if self.members[m].type4.len() == 1 {
                    return false;
                }
                let mut path_edge_number = self.members[m].path_edges.len();
                if path_edge_number == 0 {
                    if self.members[m].type2or3.len() == 2 {
                        let end_edge0 = self.members[m].type2or3[0];
                        let end_edge1 = self.members[m].type2or3[1];
                        let next = self.next_edge(end_edge0);
                        self.swap_edges(next, end_edge1);
                        let mk0 = self.edges[end_edge0].marker.unwrap();
                        self.set_parent_orient(mk0, TAIL);
                        let mk1 = self.edges[end_edge1].marker.unwrap();
                        self.set_parent_orient(mk1, HEAD);
                    } else {
                        return false;
                    }
                    return true;
                }
                let ini_edge = self.members[m].path_edges.remove(0);
                self.edges[ini_edge].path_edge_flag = false;
                path_edge_number -= 1;
                let mut wke = self.next_edge(ini_edge);
                while path_edge_number > 0 {
                    if self.edges[wke].path_edge_flag {
                        self.edges[wke].path_edge_flag = false;
                        wke = self.next_edge(wke);
                        path_edge_number -= 1;
                        continue;
                    }
                    let tedge = self.members[m].path_edges.remove(0);
                    if !self.edges[tedge].path_edge_flag {
                        continue;
                    }
                    self.swap_edges(wke, tedge);
                    wke = tedge;
                    self.edges[wke].path_edge_flag = false;
                    wke = self.next_edge(wke);
                    path_edge_number -= 1;
                }
                if wke == ini_edge {
                    false
                } else if self.members[m].type2or3.len() == 2 {
                    let end_edge0 = self.members[m].type2or3[0];
                    self.swap_edges(wke, end_edge0);
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let end_edge1 = self.members[m].type2or3[1];
                    let pre = self.pre_edge(ini_edge);
                    self.swap_edges(end_edge1, pre);
                    let mk1 = self.edges[end_edge1].marker.unwrap();
                    self.set_parent_orient(mk1, TAIL);
                    true
                } else if self.members[m].type2or3.len() == 1 {
                    let end_edge0 = self.members[m].type2or3[0];
                    self.swap_edges(wke, end_edge0);
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let n = self.edges[ini_edge].head.unwrap();
                    self.rule3(m, n)
                } else {
                    let n1 = self.edges[wke].head.unwrap();
                    let n2 = self.edges[ini_edge].head.unwrap();
                    self.rule2(m, n1, n2)
                }
            }

            GraphType::Bond => {
                if self.members[m].type2or3.len() == 2 {
                    let end_edge0 = self.members[m].type2or3[0];
                    let wke = self.members[m].type2or3[1];
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let h0 = self.node_find(self.edges[end_edge0].head.unwrap());
                    let h1 = self.node_find(self.edges[wke].head.unwrap());
                    let mk1 = self.edges[wke].marker.unwrap();
                    if !self.members[m].path_edges.is_empty() {
                        if h0 == h1 {
                            self.set_parent_orient(mk1, TAIL);
                        } else {
                            self.set_parent_orient(mk1, HEAD);
                        }
                    } else if h0 == h1 {
                        self.set_parent_orient(mk1, HEAD);
                    } else {
                        self.set_parent_orient(mk1, TAIL);
                    }
                    true
                } else if self.members[m].type2or3.len() == 1 {
                    let end_edge0 = self.members[m].type2or3[0];
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    self.set_parent_orient(mk0, HEAD);
                    let n = self.edges[end_edge0].tail.unwrap();
                    self.rule3(m, n)
                } else if self.members[m].type4.len() == 1 {
                    true
                } else {
                    let wke = self
                        .parent_marker_edge(m)
                        .expect("bond without parent marker edge");
                    let t = self.node_find(self.edges[wke].tail.unwrap());
                    let h = self.node_find(self.edges[wke].head.unwrap());
                    self.rule2(m, t, h)
                }
            }

            GraphType::Prime => {
                if self.members[m].path_edges.is_empty() {
                    if self.members[m].type2or3.len() != 2 {
                        return false;
                    }
                    let end_edge0 = self.members[m].type2or3[0];
                    let wke = self.members[m].type2or3[1];
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    let mk1 = self.edges[wke].marker.unwrap();
                    let h0 = self.node_find(self.edges[end_edge0].head.unwrap());
                    let t0 = self.node_find(self.edges[end_edge0].tail.unwrap());
                    let h1 = self.node_find(self.edges[wke].head.unwrap());
                    let t1 = self.node_find(self.edges[wke].tail.unwrap());
                    if h0 == h1 {
                        self.set_parent_orient(mk0, HEAD);
                        self.set_parent_orient(mk1, HEAD);
                    } else if h0 == t1 {
                        self.set_parent_orient(mk0, HEAD);
                        self.set_parent_orient(mk1, TAIL);
                    } else if t0 == h1 {
                        self.set_parent_orient(mk0, TAIL);
                        self.set_parent_orient(mk1, HEAD);
                    } else if t0 == t1 {
                        self.set_parent_orient(mk0, TAIL);
                        self.set_parent_orient(mk1, TAIL);
                    } else {
                        return false;
                    }
                    return true;
                }
                let mut path_nodes: Vec<NodeId> = Vec::new();
                for i in 0..self.members[m].path_edges.len() {
                    let wke = self.members[m].path_edges[i];
                    let h = self.node_find(self.edges[wke].head.unwrap());
                    if self.nodes[h].degree == 0 {
                        path_nodes.push(h);
                    }
                    self.nodes[h].degree += 1;
                    let t = self.node_find(self.edges[wke].tail.unwrap());
                    if self.nodes[t].degree == 0 {
                        path_nodes.push(t);
                    }
                    self.nodes[t].degree += 1;
                }
                let mut path_ends: Vec<NodeId> = Vec::new();
                while let Some(wkn) = path_nodes.pop() {
                    if self.nodes[wkn].degree > 2 {
                        return false;
                    }
                    if self.nodes[wkn].degree == 1 {
                        path_ends.push(wkn);
                    }
                    self.nodes[wkn].degree = 0;
                }
                if path_ends.len() != 2 {
                    return false;
                }
                let end0 = path_ends[0];
                let end1 = path_ends[1];
                if self.members[m].type4.len() == 1 {
                    let t4 = self.members[m].type4[0];
                    self.edge_is_end(t4, end0, end1)
                } else if self.members[m].type2or3.len() == 2 {
                    let end_edge = [self.members[m].type2or3[0], self.members[m].type2or3[1]];
                    let end = [end0, end1];
                    for i in 0..2 {
                        let head_i = self.node_find(self.edges[end_edge[i]].head.unwrap());
                        let tail_i = self.node_find(self.edges[end_edge[i]].tail.unwrap());
                        let mk_i = self.edges[end_edge[i]].marker.unwrap();
                        let head_other = self.node_find(self.edges[end_edge[1 - i]].head.unwrap());
                        let tail_other = self.node_find(self.edges[end_edge[1 - i]].tail.unwrap());
                        let mk_other = self.edges[end_edge[1 - i]].marker.unwrap();
                        for j in 0..2 {
                            if head_i == end[j] || tail_i == end[j] {
                                if head_i == end[j] {
                                    if tail_i == end[1 - j] {
                                        continue;
                                    }
                                    self.set_parent_orient(mk_i, HEAD);
                                } else {
                                    if head_i == end[1 - j] {
                                        continue;
                                    }
                                    self.set_parent_orient(mk_i, TAIL);
                                }
                                if head_other == end[1 - j] {
                                    self.set_parent_orient(mk_other, HEAD);
                                } else if tail_other == end[1 - j] {
                                    self.set_parent_orient(mk_other, TAIL);
                                } else {
                                    return false;
                                }
                                return true;
                            }
                        }
                    }
                    // both child markers are parallel and span the path ends
                    if !self.edges_same_end(end_edge[0], end_edge[1]) {
                        return false;
                    }
                    if !self.edge_is_end(end_edge[0], end0, end1) {
                        return false;
                    }
                    let head0 = self.node_find(self.edges[end_edge[0]].head.unwrap());
                    let head1 = self.node_find(self.edges[end_edge[1]].head.unwrap());
                    let tail1 = self.node_find(self.edges[end_edge[1]].tail.unwrap());
                    let mk0 = self.edges[end_edge[0]].marker.unwrap();
                    let mk1 = self.edges[end_edge[1]].marker.unwrap();
                    if head0 == head1 {
                        self.set_parent_orient(mk0, HEAD);
                        self.set_parent_orient(mk1, TAIL);
                    } else if head0 == tail1 {
                        self.set_parent_orient(mk0, HEAD);
                        self.set_parent_orient(mk1, HEAD);
                    } else {
                        return false;
                    }
                    true
                } else if self.members[m].type2or3.len() == 1 {
                    let end_edge0 = self.members[m].type2or3[0];
                    let mk0 = self.edges[end_edge0].marker.unwrap();
                    let head0 = self.node_find(self.edges[end_edge0].head.unwrap());
                    let tail0 = self.node_find(self.edges[end_edge0].tail.unwrap());
                    let node_for_r3;
                    if tail0 == end0 {
                        node_for_r3 = end1;
                        self.set_parent_orient(mk0, TAIL);
                    } else if head0 == end0 {
                        node_for_r3 = end1;
                        self.set_parent_orient(mk0, HEAD);
                    } else if tail0 == end1 {
                        self.set_parent_orient(mk0, TAIL);
                        node_for_r3 = end0;
                    } else if head0 == end1 {
                        self.set_parent_orient(mk0, HEAD);
                        node_for_r3 = end0;
                    } else {
                        return false;
                    }
                    self.rule3(m, node_for_r3)
                } else {
                    self.rule2(m, end0, end1)
                }
            }
        }
    }

    /// Port of Java `TDecomposition.getGraph()`: contracts every remaining
    /// marker (merging each block's t-decomposition into a single member),
    /// assigns vertex numbers, and records the endpoints of the edge realizing
    /// each matrix column. Returns the number of vertices.
    fn extract_endpoints(
        &mut self,
        pcols: &[usize],
        endpoints: &mut [Option<(usize, usize)>],
    ) -> usize {
        for i in 0..self.tree_edge.len() {
            let Some(e) = self.tree_edge[i] else { continue };
            loop {
                let m = self.edge_member(e);
                match self.members[m].parent_marker {
                    Some(mk) => self.marker_union(mk, GraphType::Prime),
                    None => break,
                }
            }
        }
        for i in 0..self.cotree_edges.len() {
            let e = self.cotree_edges[i];
            loop {
                let m = self.edge_member(e);
                match self.members[m].parent_marker {
                    Some(mk) => self.marker_union(mk, GraphType::Prime),
                    None => break,
                }
            }
        }

        let mut real_edges: Vec<EdgeId> = self.tree_edge.iter().copied().flatten().collect();
        real_edges.extend(self.cotree_edges.iter().copied());

        let mut node_count = 0;
        for &e in &real_edges {
            let h = self.node_find(self.edges[e].head.unwrap());
            if self.nodes[h].name.is_none() {
                self.nodes[h].name = Some(node_count);
                node_count += 1;
            }
            let t = self.node_find(self.edges[e].tail.unwrap());
            if self.nodes[t].name.is_none() {
                self.nodes[t].name = Some(node_count);
                node_count += 1;
            }
        }

        for i in 0..self.tree_edge.len() {
            if let Some(e) = self.tree_edge[i] {
                let h = self.node_find(self.edges[e].head.unwrap());
                let t = self.node_find(self.edges[e].tail.unwrap());
                endpoints[pcols[i]] =
                    Some((self.nodes[h].name.unwrap(), self.nodes[t].name.unwrap()));
            }
        }
        for idx in 0..self.cotree_edges.len() {
            let e = self.cotree_edges[idx];
            let h = self.node_find(self.edges[e].head.unwrap());
            let t = self.node_find(self.edges[e].tail.unwrap());
            endpoints[self.cotree_cols[idx]] =
                Some((self.nodes[h].name.unwrap(), self.nodes[t].name.unwrap()));
        }
        node_count
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{rngs::SmallRng, Rng, SeedableRng};

    fn column_weights_ok(m: &BitMatrix) -> bool {
        (0..m.cols()).all(|j| (0..m.rows()).filter(|&i| m.bit(i, j)).count() <= 2)
    }

    fn same_rowspace(a: &BitMatrix, b: &BitMatrix) -> bool {
        let rank = a.rank();
        rank == b.rank() && a.vstack(b).rank() == rank
    }

    /// incidence matrix of a random multigraph (parallel edges allowed)
    fn random_graph_incidence(rng: &mut SmallRng, vertices: usize, edges: usize) -> BitMatrix {
        let mut m = BitMatrix::zeros(vertices, edges);
        for j in 0..edges {
            let u = rng.random_range(0..vertices);
            let mut v = rng.random_range(0..vertices - 1);
            if v >= u {
                v += 1;
            }
            m.set_bit(u, j, true);
            m.set_bit(v, j, true);
        }
        m
    }

    /// incidence matrix of the complete graph on `n` vertices
    fn complete_graph_incidence(n: usize) -> BitMatrix {
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|u| (u + 1..n).map(move |v| (u, v)))
            .collect();
        BitMatrix::build(n, pairs.len(), |i, j| pairs[j].0 == i || pairs[j].1 == i)
    }

    /// a representation of the dual matroid: a basis of the nullspace
    fn dual_basis(m: &BitMatrix) -> BitMatrix {
        BitMatrix::vstack_from_iter(&m.nullspace())
    }

    fn scramble_rows(rng: &mut SmallRng, m: &BitMatrix) -> BitMatrix {
        let s = BitMatrix::random_invertible(rng, m.rows());
        &s * m
    }

    fn check_graphic(m: &BitMatrix) {
        let (n, b, skipped) = m.graphic_form_with_options(false, true);
        let n = n.expect("matrix should be graphic");
        let b = b.expect("basis change should be computed");
        assert!(skipped.is_empty());
        assert!(column_weights_ok(&n), "column weight > 2 in\n{}", n);
        assert_eq!(n.rows(), m.rank(), "output does not have full row rank");
        assert!(same_rowspace(m, &n), "rowspace changed:\n{}\n->\n{}", m, n);
        assert_eq!(&b * m, n, "basis change does not map input to output");
    }

    #[test]
    fn round_trip_random_graphs() {
        let mut rng = SmallRng::seed_from_u64(42);
        for vertices in [2, 3, 4, 5, 8, 12, 20] {
            for _ in 0..20 {
                let edges = rng.random_range(1..=3 * vertices);
                let inc = random_graph_incidence(&mut rng, vertices, edges);
                check_graphic(&inc);
                check_graphic(&scramble_rows(&mut rng, &inc));
            }
        }
    }

    #[test]
    fn round_trip_large_random_graph() {
        let mut rng = SmallRng::seed_from_u64(1234);
        for _ in 0..3 {
            let inc = random_graph_incidence(&mut rng, 40, 120);
            check_graphic(&scramble_rows(&mut rng, &inc));
        }
    }

    #[test]
    fn complete_graphs_are_graphic() {
        for n in 2..=6 {
            check_graphic(&complete_graph_incidence(n));
        }
    }

    // the example input from the reference implementation's README (10 tree
    // edges, 5 fundamental circuits), encoded as [I | B]
    #[test]
    fn readme_example_is_graphic() {
        let circuits: [&[usize]; 5] = [&[0, 1, 2, 3], &[4, 5, 6], &[7, 8], &[1, 2, 8], &[7, 0]];
        let m = BitMatrix::build(10, 15, |i, j| {
            if j < 10 {
                i == j
            } else {
                circuits[j - 10].contains(&i)
            }
        });
        check_graphic(&m);
    }

    // the Fano plane F7 is the smallest non-graphic binary matroid
    #[test]
    fn fano_is_not_graphic() {
        let f7 = BitMatrix::build(3, 7, |i, j| (j + 1) & (1 << i) != 0);
        assert!(f7.graphic_form().is_none());
        let mut rng = SmallRng::seed_from_u64(7);
        assert!(scramble_rows(&mut rng, &f7).graphic_form().is_none());
    }

    // the dual of a graphic matroid is graphic iff the graph is planar, so
    // the duals of K5 and K3,3 are not graphic while the dual of K4 is
    #[test]
    fn dual_k5_is_not_graphic() {
        let dual = dual_basis(&complete_graph_incidence(5));
        assert_eq!(dual.rank(), 6);
        assert!(dual.graphic_form().is_none());
        let mut rng = SmallRng::seed_from_u64(5);
        assert!(scramble_rows(&mut rng, &dual).graphic_form().is_none());
    }

    #[test]
    fn dual_k33_is_not_graphic() {
        let pairs: Vec<(usize, usize)> = (0..3).flat_map(|u| (3..6).map(move |v| (u, v))).collect();
        let k33 = BitMatrix::build(6, 9, |i, j| pairs[j].0 == i || pairs[j].1 == i);
        check_graphic(&k33);
        let dual = dual_basis(&k33);
        assert_eq!(dual.rank(), 4);
        assert!(dual.graphic_form().is_none());
    }

    #[test]
    fn dual_k4_is_graphic() {
        check_graphic(&dual_basis(&complete_graph_incidence(4)));
    }

    #[test]
    fn empty_and_zero_matrices() {
        let (n, b, _) = BitMatrix::zeros(0, 0).graphic_form_with_options(false, true);
        let (n, b) = (n.unwrap(), b.unwrap());
        assert_eq!((n.rows(), n.cols()), (0, 0));
        assert_eq!((b.rows(), b.cols()), (0, 0));

        let (n, b, _) = BitMatrix::zeros(3, 4).graphic_form_with_options(false, true);
        let (n, b) = (n.unwrap(), b.unwrap());
        assert_eq!((n.rows(), n.cols()), (0, 4));
        assert_eq!((b.rows(), b.cols()), (0, 3));
    }

    #[test]
    fn identity_is_graphic() {
        check_graphic(&BitMatrix::identity(5));
    }

    #[test]
    fn simple_graphic_ex() {
        let m = BitMatrix::from_int_vec(&vec![
            vec![1, 1, 0, 1, 1],
            vec![0, 1, 1, 1, 0],
            vec![1, 1, 1, 0, 0],
        ]);

        let n = m.graphic_form().unwrap();
        println!("graphic form:\n{}", n);
        check_graphic(&m);
    }

    #[test]
    fn single_entry_and_row() {
        check_graphic(&BitMatrix::build(1, 1, |_, _| true));
        // a single all-ones row realizes as a bundle of parallel edges
        check_graphic(&BitMatrix::build(1, 6, |_, _| true));
    }

    #[test]
    fn duplicate_and_zero_columns() {
        let m = BitMatrix::from_int_vec(&vec![
            vec![1, 1, 0, 0, 1, 0],
            vec![0, 1, 1, 0, 1, 0],
            vec![0, 0, 1, 0, 1, 0],
        ]);
        check_graphic(&m);
        // zero columns of the input stay zero columns of the output
        let n = m.graphic_form().unwrap();
        for i in 0..n.rows() {
            assert!(!n.bit(i, 3));
            assert!(!n.bit(i, 5));
        }
    }

    #[test]
    fn rank_deficient_input() {
        let mut rng = SmallRng::seed_from_u64(9);
        let inc = random_graph_incidence(&mut rng, 6, 10);
        check_graphic(&inc.vstack(&inc));
    }

    fn column_weight(m: &BitMatrix, j: usize) -> usize {
        (0..m.rows()).filter(|&i| m.bit(i, j)).count()
    }

    /// checks the partial-realization contract and returns the skipped columns
    fn check_partial(m: &BitMatrix) -> Vec<usize> {
        let (n, b, skipped) = m.graphic_form_with_options(true, true);
        let (n, b) = (n.unwrap(), b.unwrap());
        assert_eq!(n.rows(), m.rank(), "output does not have full row rank");
        assert!(same_rowspace(m, &n), "rowspace changed:\n{}\n->\n{}", m, n);
        assert_eq!(&b * m, n, "basis change does not map input to output");
        for j in 0..n.cols() {
            if !skipped.contains(&j) {
                assert!(
                    column_weight(&n, j) <= 2,
                    "unskipped column {} has weight > 2 in\n{}",
                    j,
                    n
                );
            }
        }
        assert_eq!(
            skipped.is_empty(),
            m.graphic_form().is_some(),
            "skipped columns must be empty exactly for graphic inputs"
        );
        skipped
    }

    #[test]
    fn partial_matches_on_graphic_inputs() {
        let mut rng = SmallRng::seed_from_u64(11);
        for _ in 0..10 {
            let inc = random_graph_incidence(&mut rng, 8, 20);
            let m = scramble_rows(&mut rng, &inc);
            let (n, b, skipped) = m.graphic_form_with_options(true, true);
            let (n, b) = (n.unwrap(), b.unwrap());
            assert!(skipped.is_empty());
            assert_eq!(Some(n.clone()), m.graphic_form());
            assert_eq!(&b * &m, n);
        }
    }

    #[test]
    fn simple_partial_graphic_ex() {
        let m = BitMatrix::from_int_vec(&vec![
            vec![1, 1, 0, 1, 1, 0, 0],
            vec![0, 1, 1, 1, 0, 1, 0],
            vec![1, 1, 1, 0, 0, 0, 1],
        ]);

        let (n, b, hyper) =
            m.graphic_form_with_options(/*partial=*/ true, /*basis_change=*/ true);
        let n = n.unwrap();
        let b = b.unwrap();
        println!("partial graphic form:\n{}", n);
        println!("basis-change matrix:\n{}", b);
        println!("hyperedge columns: {:?}", hyper);
        assert_eq!(hyper.len(), 1);
        assert_eq!(&b * &m, n);
    }

    #[test]
    fn partial_on_fano() {
        let f7 = BitMatrix::build(3, 7, |i, j| (j + 1) & (1 << i) != 0);
        assert!(!check_partial(&f7).is_empty());
    }

    #[test]
    fn partial_on_dual_k5() {
        let dual = dual_basis(&complete_graph_incidence(5));
        assert!(!check_partial(&dual).is_empty());
        let mut rng = SmallRng::seed_from_u64(55);
        assert!(!check_partial(&scramble_rows(&mut rng, &dual)).is_empty());
    }

    // dense random matrices are almost never graphic, so this exercises the
    // greedy skip-and-restart loop heavily
    #[test]
    fn partial_on_random_matrices() {
        let mut rng = SmallRng::seed_from_u64(77);
        for _ in 0..30 {
            let rows = rng.random_range(1..10);
            let cols = rng.random_range(1..30);
            let m = BitMatrix::random(&mut rng, rows, cols);
            check_partial(&m);
        }
    }

    #[test]
    fn partial_on_trivial_inputs() {
        let (n, b, skipped) = BitMatrix::zeros(0, 0).graphic_form_with_options(true, true);
        let (n, b) = (n.unwrap(), b.unwrap());
        assert_eq!((n.rows(), n.cols()), (0, 0));
        assert_eq!((b.rows(), b.cols()), (0, 0));
        assert!(skipped.is_empty());

        let (n, b, skipped) = BitMatrix::zeros(3, 4).graphic_form_with_options(true, true);
        let (n, b) = (n.unwrap(), b.unwrap());
        assert_eq!((n.rows(), n.cols()), (0, 4));
        assert_eq!((b.rows(), b.cols()), (0, 3));
        assert!(skipped.is_empty());
    }
}
