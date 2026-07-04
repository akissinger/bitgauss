// Graph Realization Problem
// (t-decomposition)
// ver 3.0.0 2/24/2003
//
//

import java.io.*;
import java.util.*;
////////////////////////////////////////////////////////////
public class GRP{
	public static void main(String args[]){
		int i, treeSize;
		boolean graphic = true;
		ArrayList block, subblock;
		MatrixDivide md;
		TDecomposition t;
		
		if(args.length != 1) {
			System.out.println("Error : >java GPR [filename]");
			System.exit(1);
		}
		md = new MatrixDivide(args[0]);
		treeSize = md.getTreeSize();
		block = md.getBlock();
		t = new TDecomposition(treeSize);
		for(i=0; i<block.size(); i++) {
			if(!t.realization((ArrayList)block.get(i))) {
				System.out.println("NonGraphic");
				graphic = false;
			}
		}
		if(graphic) {
			System.out.println("Graphic. BlockSize is " + block.size() + ".");
			t.disp();
		}
	}
}
//////////////////////////////////////////////////////////// 
class MatrixDivide {
	int treeSize, coTreeIndex = 0;
	ArrayList block = new ArrayList();
//----------------------------------------------------------
	public MatrixDivide(String fileName) {
		int i, j, k, l, num, len;
		int[] rq;
		String s;
		StringTokenizer st;
		ArrayList[] tree;
		ArrayList coTree = new ArrayList(),
		                   queue = new ArrayList(), subBlock;
		boolean[] treeFlg, coTreeFlg;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			s = br.readLine();
			treeSize = Integer.parseInt(s);
			tree = new ArrayList[treeSize+1];
			for(i=0; i<=treeSize; i++){
				tree[i] = new ArrayList();
			}

			while((s = br.readLine()) != null) {
				st = new StringTokenizer(s, " ");
				num = st.countTokens();
				rq = new int[num];
				for(i=0; i<num; i++) {
					rq[i] = Integer.valueOf(st.nextToken()).intValue();
				}			
				coTree.add(rq);
				for(i=0; i<num; i++){
					tree[rq[i]].add(new Integer(coTreeIndex));
				}
				coTreeIndex++;
			}

			treeFlg = new boolean[treeSize+1];
			for(i=0; i<=treeSize; i++) {
				treeFlg[i] = false;
			}
			coTreeFlg = new boolean[coTreeIndex];
			for(i=0; i<coTreeIndex; i++) {
				coTreeFlg[i] = false;
			}
			
			for(i=0; i<coTreeIndex; i++) {
				if(coTreeFlg[i]) continue;
				coTreeFlg[i] = true;
				queue.add((int[])coTree.get(i));
				subBlock = new ArrayList();
				subBlock.add((int[])coTree.get(i));
				while(queue.size() >0) {
					rq = (int[])queue.get(0); queue.remove(0);
					for(j=0; j<rq.length; j++) {
						if(treeFlg[rq[j]]) continue;
						treeFlg[rq[j]] = true;
						for(k=0; k<tree[rq[j]].size(); k++) {
							l = ((Integer)tree[rq[j]].get(k)).intValue();
							if(coTreeFlg[l]) continue;
							coTreeFlg[l] = true;
							queue.add((int[])coTree.get(l));
							subBlock.add((int[])coTree.get(l));
						}
					}
				}
				block.add(subBlock);
			}
			
/*
 // Test for check blocked
			for(i=0; i<block.size(); i++){
				System.out.println("block " + i);
				subBlock = (ArrayList)block.get(i);
				for(j=0; j<subBlock.size(); j++) {
					rq = (int[])subBlock.get(j);
					for(k=0; k<rq.length; k++)
						System.out.print(rq[k] + " ");
					System.out.println();
				}
			}
*/	
		}
		catch(Exception e) {
			System.out.println("Execption: " +e);
		}
	}
//----------------------------------------------------------
	public ArrayList getBlock() {
		return block;
	}
//----------------------------------------------------------
	public int getTreeSize() {
		return treeSize;
	}
//----------------------------------------------------------
}
//////////////////////////////////////////////////////////// 
class TDecomposition{
	Member[] k = new Member[2];
	Node[] u = new Node[2];
	Edge[] treeEdge;
	ArrayList coTreeEdge; // Edge
	int coTreeIndex;
//----------------------------------------------------------
	public TDecomposition(int numberOfTreeEdge) {
		treeEdge = new Edge[numberOfTreeEdge+1];
		coTreeEdge = new ArrayList();
		coTreeIndex = numberOfTreeEdge+1;
	}
//----------------------------------------------------------
	public boolean realization(ArrayList subblock) {
		int[] rq;
		int i, firstEdge;
		
		rq = (int[])subblock.get(0);
		firstEdge = rq[0];
		treeEdge[0] = new Edge(0);
		treeEdge[firstEdge] = new Edge(firstEdge);
		new Member(treeEdge[0], treeEdge[firstEdge]);
		addCircuit(rq);
		for(i=1; i<subblock.size(); i++) {
			rq = (int[])subblock.get(i);
			if(!addCircuit(rq)) return false;
		}
		return true;
	}
//----------------------------------------------------------
	boolean addCircuit(int[] s) {
		int i;
		ArrayList p = new ArrayList(); // Integer. p is in the td
		ArrayList c = new ArrayList(); // Integer. c isn't in the td
		
		for(i=0; i<2; i++) {
			k[i] = null;
			u[i] = null;
		}
		for(i=0; i<s.length; i++) {
			if(treeEdge[s[i]] == null) {
				c.add(new Integer(s[i]));
			}
			else {
				p.add(new Integer(s[i]));
				treeEdge[s[i]].pathEdgeFlag = true;
			}
		}
		c.add(new Integer(coTreeIndex++));
		return (hypopath(p) && update(c));
	}
//----------------------------------------------------------
	Stack getGraph() {
		int i, nodeCount = 0;
		Edge wke;
		Node wkn;
		Stack graphStack = new Stack();
		int[] edgeA;
		
		for(i=1; i < treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			while(wke.getMember().parentMarker != null) {
				wke.getMember().parentMarker.union(Member.prime);
			}
		}
		for(i=0; i < coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			while(wke.getMember().parentMarker != null) {
				wke.getMember().parentMarker.union(Member.prime);
			}
		}
		for(i=1; i < treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			wkn = wke.head.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
			wkn = wke.tail.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
		}
		for(i=0; i < coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			wkn = wke.head.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
			wkn = wke.tail.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
		}
		for(i=1; i<treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			edgeA = new int[3];
			edgeA[0] = i;
			edgeA[1] = wke.head.find().name;
			edgeA[2] = wke.tail.find().name;
			graphStack.push(edgeA);
		}
		for(i=0; i<coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			edgeA = new int[3];
			edgeA[0] = wke.name;
			edgeA[1] = wke.head.find().name;
			edgeA[2] = wke.tail.find().name;
			graphStack.push(edgeA);			
		}
		return graphStack;
	}
//----------------------------------------------------------
	void disp() {
		int i, nodeCount = 0;
		Edge wke;
		Node wkn;
		
		for(i=1; i < treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			while(wke.getMember().parentMarker != null) {
				wke.getMember().parentMarker.union(Member.prime);
			}
		}
		for(i=0; i < coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			while(wke.getMember().parentMarker != null) {
				wke.getMember().parentMarker.union(Member.prime);
			}
		}
		for(i=1; i < treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			wkn = wke.head.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
			wkn = wke.tail.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
		}
		for(i=0; i < coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			wkn = wke.head.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
			wkn = wke.tail.find();
			if(wkn.name == 0) wkn.name = ++nodeCount;
		}
		for(i = 1; i < treeEdge.length; i++) {
			wke = treeEdge[i];
			if(wke == null) continue;
			System.out.println(i + ": (" + wke.head.find().name + ", " + wke.tail.find().name
			 + ")");
		}
		for(i = 0; i < coTreeEdge.size(); i++) {
			wke = (Edge)coTreeEdge.get(i);
			System.out.println(wke.name + ": (" + wke.head.find().name + ", " + wke.tail.find().name
			 + ")");
		}
	}
//----------------------------------------------------------
	private boolean hypopath(ArrayList p){
		int i, j, type2or3Size, type4Size;
		Member wkm, root;
		Edge e;
		ArrayList pathMember = new ArrayList(); // Member
		Stack wks; // Member
		ArrayList depthPartition; // Stack
		
		for(i=0; i<p.size(); i++){
			e = treeEdge[((Integer)p.get(i)).intValue()];
			wkm = e.getMember();
			if(wkm.pathEdge.size() == 0) {
				pathMember.add(wkm);
			}
			wkm.addPathEdge(e);
		}
		
		root = findRoot(pathMember);
		depthPartition = makeDepthPartition(pathMember, root);
		i = depthPartition.size();
		
		while(--i > 0) {
			wks = (Stack)depthPartition.get(i);
			while(!wks.empty()) {
				wkm = (Member)wks.pop();
				type2or3Size = wkm.type2or3.size();
				type4Size = wkm.type4.size();
				if((type4Size > 1) || (type2or3Size > 2) || 
		       ((type4Size == 1)&&(type2or3Size > 0))) return false;
				if(!typing(wkm)) return false;
				for(j=0; j<wkm.pathEdge.size(); j++) {
					((Edge)wkm.pathEdge.get(j)).pathEdgeFlag = false;
				}
				wkm.pathEdge.clear();
				wkm.type2or3.clear();
				wkm.type4.clear();
			}
		}
		type2or3Size = root.type2or3.size();
		type4Size = root.type4.size();
		if((type4Size > 1) || (type2or3Size > 2) || 
		   ((type4Size == 1)&&(type2or3Size > 0))) return false;
		if(!isPath(root)) return false;
		for(j=0; j<root.pathEdge.size(); j++) {
			((Edge)root.pathEdge.get(j)).pathEdgeFlag = false;
		}
		root.pathEdge.clear();
		root.type2or3.clear();
		root.type4.clear();
		rule5();
		return true;
	}
//----------------------------------------------------------
	boolean update(ArrayList c) {
		Edge f, f1, f2, e1, e2, wke;
		Member[] wkmA = new Member[2];
		Marker m, m1, m2;
		int i;
		Member root, wkm;
		ArrayList r1 = new ArrayList(); // Member
		ArrayList r2 = new ArrayList(); // Member
		Stack markerStack = new Stack(); // Marker
		
		if(c.size() == 1){
			f = new Edge(((Integer)c.get(0)).intValue());
			coTreeEdge.add(f);
		}
		else{
			m = new Marker();
			f = m.parent;
			new Member(m.child, c, treeEdge, coTreeEdge);
		}
		if(k[0] == k[1]){
			if(k[0].graphType != Member.polygon){
				if(k[0].graphType == Member.bond) {
					k[0].edgeSize++;
					wke = k[0].parentMarkerEdge();
					u[0] = wke.head.find();
					u[1] = wke.tail.find();
				}
				f.setEndNode(u[0], u[1]);// (head, tail)
				f.setMember(k[0]);
			} else if(u[0].isJoined(u[1])){
				f1 = u[0].jointEdge(u[1]);
				m2 = new Marker();
				m2.parent.setMember(k[0]);
				f1.head.tail = m2.parent;
				f1.tail.head = m2.parent;
				m2.parent.tail = f1.tail;
				m2.parent.head = f1.head;
				new Member(m2.child, f, f1);
			} else {
				k[0].divide(u[1], u[0], f);
			}
		} else {
			root = findRoot(k[0], k[1]);
			wkm = k[0];
			while(wkm != root) {
				r1.add(wkm);
				wkm = wkm.getParent();
			}
			wkm = k[1];
			while(wkm != root) {
				r2.add(wkm);
				wkm = wkm.getParent();
			}
			wkmA[0] = (Member)r1.get(r1.size()-1);
			e1 = wkmA[0].parentMarker.parent;
			if((root != k[0]) && (root != k[1])) {
				wkmA[1] = (Member)r2.get(r2.size()-1);
				e2 = wkmA[1].parentMarker.parent;
			  if((root.graphType == Member.prime)&&(e1.isSameEnd(e2))){
					m1 = new Marker();
					f1 = m1.parent;
					f1.memberName = root;
					f1.tail = e1.tail.find();
					f1.head = e1.head.find();
					root = new Member(m1.child, e1, e2);
					if(wkmA[0].graphType == Member.bond) {
						wkmA[0].parentMarker.union(Member.bond);
						r1.remove(r1.size() -1);
						root = root.find();
						e1 = ((Member)r1.get(r1.size()-1)).parentMarker.parent;
					}
					if(wkmA[1].graphType == Member.bond) {
						wkmA[1].parentMarker.union(Member.bond);
						r2.remove(r2.size() -1);
						root = root.find();
						e2 = ((Member)r2.get(r2.size()-1)).parentMarker.parent;
					}
				}
				if((root.graphType == Member.bond)&&(root.edgeSize > 3)) {
					root.edgeSize--;
					m1 = new Marker();
					f1 = m1.parent;
					f1.memberName = root;
					f1.tail = e1.tail.find();
					f1.head = e1.head.find();
					root = new Member(m1.child, e1, e2);
				}
				root.rootSqueeze(e1, e2);
			}
			root = findRoot(k[0], k[1]);
			wkm = k[0];
			while(wkm != root) {
				markerStack.push(wkm.parentMarker);
				wkm = wkm.getParent();
			}
			wkm = k[1];
			while(wkm != root) {
				markerStack.push(wkm.parentMarker);
				wkm = wkm.getParent();
			}
			k[0].squeeze(u[0]);
			if(k[1] == root) k[1].rootSqueeze(u[1], e1);
			else k[1].squeeze(u[1]);
			wke = k[0].parentMarker.parent;
			for(i=1; i<r1.size(); i++) {
				wkm = (Member)r1.get(i);
				wkm.squeeze(wke);
				wke = wkm.parentMarker.parent;
			}
			wke = k[1].parentMarker.parent;
			for(i=1; i<r2.size(); i++) {
				wkm = (Member)r2.get(i);
				wkm.squeeze(wke);
				wke = wkm.parentMarker.parent;
			}
			while(!markerStack.isEmpty()) {
				((Marker)markerStack.pop()).union(Member.prime);
			}
			f.setEndNode(u[0], u[1]);// (head, tail)
			f.setMember(e1.getMember().find());
		}
		return true;
	}
//----------------------------------------------------------
	Member findRoot(ArrayList p) {
		int i;
		Member[] m = new Member[2];
		Member root;
		ArrayList copy = new ArrayList(p.size());
		
		for(i=0; i<p.size(); i++) {
			copy.add(p.get(i));
		}
		while(copy.size() > 1) {
			for(i=0; i<2; i++) {
				m[i] = (Member)copy.get(0);
				copy.remove(0);
			}
			copy.add(findRoot(m[0], m[1]));
		}
		return (Member)copy.get(0);
	}
//----------------------------------------------------------
	Member findRoot(Member m1, Member m2) {
		Member root = null;
		Member[] wk = new Member[2];
		Stack ini = new Stack(); // Member
		int i;

		wk[0] = m1; wk[1] = m2;
		while((wk[0] != null) && (wk[1] != null)) {
			for(i=0; i<2; i++) {
				if(wk[i].findRootFlag) {
					root = wk[i];
					while(!ini.isEmpty()) {
						((Member)ini.pop()).findRootFlag = false;
					}
					return root.find();
				}
				wk[i].findRootFlag = true;
				ini.push(wk[i]);
				wk[i] = wk[i].getParent();
			}
		}
		for(i=0; i<2; i++){
			while(wk[i] != null) {
				if(wk[i].findRootFlag) {
					root = wk[i];
					while(!ini.isEmpty()) {
						((Member)ini.pop()).findRootFlag = false;
					}
					return root.find();
				}
				wk[i].findRootFlag = true;
				ini.push(wk[i]);
				wk[i] = wk[i].getParent();
			}
		}
		return null;
	}
//----------------------------------------------------------
	ArrayList makeDepthPartition(ArrayList path, Member root) {
		int i,d;
		Member wk;
		ArrayList depthPartition = new ArrayList(); // Stack
		Stack wks = new Stack(); // Member
		Stack iniStack = new Stack(); // Member
		Stack copy = new Stack(); // Member
		Stack depthStack = new Stack(); // Member
		
		for(i=0; i<path.size(); i++) {
			copy.push(path.get(i));
		}
		
		root.depth = 0;
		iniStack.push(root);
		wks = new Stack();
		wks.push(root);
		depthPartition.add(wks);
		d = 0;
		
		while(!copy.empty()) {
			wk = (Member)copy.pop();
			while(wk.depth == -1) {
				depthStack.push(wk);
				wk = wk.getParent();
			}
			i = wk.depth;
			while(!depthStack.empty()) {
				if(++i > d) {
					wks = new Stack();
					depthPartition.add(wks);
				}
				wk = (Member)depthStack.pop();
				wk.depth = i;
				iniStack.push(wk);
				((Stack)depthPartition.get(i)).push(wk);
			}
			if(i > d) d = i;
		}
		while(!iniStack.isEmpty()) {
			((Member)iniStack.pop()).depth = -1;
		}
		return depthPartition;
	}
//----------------------------------------------------------
	boolean typing(Member m) {
		int pathEdgeNumber, i, j, count;
		Edge marker, wke, tedge, type2or3Edge;
		ArrayList pathNode = new ArrayList(); // Node
		Node wkn, marker_h, marker_t, nodeForR3;
		ArrayList pathEnd = new ArrayList(); // Node
		Node[] end = new Node[4];
		
		marker = m.parentMarkerEdge();
		switch(m.graphType) {
			case Member.polygon:
				wke = marker.nextEdge();
				pathEdgeNumber = m.pathEdge.size();
				while(pathEdgeNumber > 0) {
					if(wke.pathEdgeFlag) {
						wke.pathEdgeFlag = false;
						wke = wke.nextEdge();
						pathEdgeNumber--;
						continue;
					}
					
					tedge = (Edge)m.pathEdge.get(0);
					m.pathEdge.remove(0);
					if(!tedge.pathEdgeFlag) continue;
					wke.swap(tedge);
					wke = tedge;
					wke.pathEdgeFlag = false;
					wke = wke.nextEdge();
					pathEdgeNumber--;
				}
				if(wke == marker) { // type = 1
					m.getParent().addPathEdge(m.parentMarker.parent);
				} else if(m.type4.size() == 1) { // type = 4
					if(wke.nextEdge() != marker) return false;
					m.getParent().type4.add(m.parentMarker.parent);
				} else if(m.type2or3.size() == 2) { // type = 4
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke.swap(type2or3Edge);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					type2or3Edge = (Edge)m.type2or3.get(1);
					type2or3Edge.swap(marker.preEdge());
					type2or3Edge.markerName.setParentOrient(Marker.tail);
					m.getParent().type4.add(m.parentMarker.parent);
				} else if(m.type2or3.size() == 1) { // type = 2or3
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke.swap(type2or3Edge);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					m.parentMarker.setChildOrient(Marker.tail);
					m.getParent().type2or3.add(m.parentMarker.parent);
				} else { // type = 2or3
					m.parentMarker.setChildOrient(Marker.tail);
					m.getParent().type2or3.add(m.parentMarker.parent);
					if(!rule1(m, wke.head)) return false;
				}
				return true;
				
			case Member.bond:
				if(m.type4.size() == 1) { // type 4
					if(m.pathEdge.size() > 0) return false;
					m.getParent().type4.add(m.parentMarker.parent);
				} else if(m.type2or3.size() == 2) { // type 4
					if(m.pathEdge.size() > 0) return false;
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke = (Edge)m.type2or3.get(1);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					if(type2or3Edge.head.find() == wke.head.find()) 
						wke.markerName.setParentOrient(Marker.tail);
					else
						wke.markerName.setParentOrient(Marker.head);
					m.getParent().type4.add(m.parentMarker.parent);
				} else if(m.type2or3.size() == 1) { // type 2or3
					if(m.pathEdge.size() > 1) return false;
					type2or3Edge = (Edge)m.type2or3.get(0);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					if(m.pathEdge.size() == 1) {
						if(marker.head.find() == type2or3Edge.head.find())
							m.parentMarker.setChildOrient(Marker.tail);
						else
							m.parentMarker.setChildOrient(Marker.head);
					} else {
						if(marker.head.find() == type2or3Edge.head.find())
							m.parentMarker.setChildOrient(Marker.head);
						else
							m.parentMarker.setChildOrient(Marker.tail);
					}
					m.getParent().type2or3.add(m.parentMarker.parent);
				} else { // type 1
					if(m.pathEdge.size() != 1) return false;
					m.getParent().addPathEdge(m.parentMarker.parent);
				}
				return true;
				
			case Member.prime:
				for(i = 0; i < m.pathEdge.size(); i++) {
					wke = (Edge)m.pathEdge.get(i);
					wkn = wke.head.find();
					if(wkn.endEdgeArrayList.size() == 0) pathNode.add(wkn);
					wkn.endEdgeArrayList.add(wke);
					wkn = wke.tail.find();
					if(wkn.endEdgeArrayList.size() == 0) pathNode.add(wkn);
					wkn.endEdgeArrayList.add(wke);
				}
				for(i=0; i < pathNode.size(); i++) {
					wkn = (Node)pathNode.get(i);
					count = wkn.endEdgeArrayList.size();
					if(count > 2) return false;
					if(count == 1) pathEnd.add(wkn);
				}
				count = pathEnd.size();
				if(count == 0) {
					for(i=0; i<pathNode.size(); i++) ((Node)pathNode.get(i)).endEdgeArrayList.clear();
					return addType234(m, marker.head.find(), marker.tail.find(), Marker.tail, (byte)count);
				} else if(count == 2) {
					for(i=0; i<pathNode.size(); i++) ((Node)pathNode.get(i)).endEdgeArrayList.clear();
					end[0] = (Node)pathEnd.get(0);
					end[1] = (Node)pathEnd.get(1);
					if(marker.isEnd(end[0], end[1])) 
						return addType234(m, marker.head.find(), marker.tail.find(), Marker.tail, (byte)count);
					else if(marker.head.find() == end[0])
						return addType234(m, end[1], marker.tail.find(), Marker.tail, (byte)count);
					else if(marker.tail.find() == end[0])
						return addType234(m, end[1], marker.head.find(), Marker.head, (byte)count);
					else if(marker.head.find() == end[1]) 
						return addType234(m, end[0], marker.tail.find(), Marker.tail, (byte)count);
					else if(marker.tail.find() == end[1]) 
						return addType234(m, end[0], marker.head.find(), Marker.head, (byte)count);
					else return false;
				} else if(count == 4) {
					for(i=0; i<4; i++) {
						end[i] = (Node)pathEnd.get(i);
					}
					if(end[0].isAnotherEnd(end[2])) {
						wkn = end[1];
						end[1] = end[2];
						end[2] = wkn;
					}
					if(end[0].isAnotherEnd(end[3])) {
						wkn = end[1];
						end[1] = end[3];
						end[3] = wkn;
					}
					for(i=0; i<pathNode.size(); i++) ((Node)pathNode.get(i)).endEdgeArrayList.clear();
					if(marker.isEnd(end[0], end[2]))
						return addType234(m, end[1], end[3], Marker.head, (byte)count);
					if(marker.isEnd(end[0], end[3]))
						return addType234(m, end[1], end[2], Marker.head, (byte)count);
					if(marker.isEnd(end[1], end[2]))
						return addType234(m, end[0], end[3], Marker.head, (byte)count);
					if(marker.isEnd(end[1], end[3]))
						return addType234(m, end[0], end[2], Marker.head, (byte)count);
					return false;
				} else return false;
		}
		return false;
	}
//----------------------------------------------------------
	boolean isPath(Member m) {
		int pathEdgeNumber = m.pathEdge.size();
		int i,j;
		Edge wke, iniEdge, tedge, type2or3Edge;
		Edge[] endEdge = new Edge[2];
		Stack pathNode = new Stack(); // Node
		ArrayList pathEnd = new ArrayList(); // Node
		Node[] end = new Node[2];
		Node nodeForR3, wkn;
		
		switch(m.graphType) {
			case Member.polygon:
				if(m.type4.size() == 1) return false;
				if(m.pathEdge.size() == 0) {
					if(m.type2or3.size() == 2) {
						endEdge[0] = (Edge)m.type2or3.get(0);
						endEdge[1] = (Edge)m.type2or3.get(1);
						endEdge[0].nextEdge().swap(endEdge[1]);
						endEdge[0].markerName.setParentOrient(Marker.tail);
						endEdge[1].markerName.setParentOrient(Marker.head);
					} else return false;
					return true;
				}
				iniEdge = (Edge)m.pathEdge.get(0);
				m.pathEdge.remove(0);
				iniEdge.pathEdgeFlag = false;
				pathEdgeNumber--;
				wke = iniEdge.nextEdge();
				while(pathEdgeNumber > 0) {
					if(wke.pathEdgeFlag) {
						wke.pathEdgeFlag = false;
						wke = wke.nextEdge();
						pathEdgeNumber--;
						continue;
					}
					
					tedge = (Edge)m.pathEdge.get(0);
					m.pathEdge.remove(0);
					if(!tedge.pathEdgeFlag) continue;
					wke.swap(tedge);
					wke = tedge;
					wke.pathEdgeFlag = false;
					wke = wke.nextEdge();
					pathEdgeNumber--;
				}
				if(wke == iniEdge) return false;
				else if(m.type2or3.size() == 2) {
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke.swap(type2or3Edge);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					type2or3Edge = (Edge)m.type2or3.get(1);
					type2or3Edge.swap(iniEdge.preEdge());
					type2or3Edge.markerName.setParentOrient(Marker.tail);
				} else if(m.type2or3.size() == 1) {
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke.swap(type2or3Edge);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					if(!rule3(m, iniEdge.head)) return false;
				} else {
					if(!rule2(m, wke.head, iniEdge.head)) return false;
				}
				return true;
				
			case Member.bond:
				if(m.type2or3.size() == 2) {
					type2or3Edge = (Edge)m.type2or3.get(0);
					wke = (Edge)m.type2or3.get(1);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					if(m.pathEdge.size() > 0) {
						if(type2or3Edge.head.find() == wke.head.find())
							wke.markerName.setParentOrient(Marker.tail);
						else
							wke.markerName.setParentOrient(Marker.head);
					} else {
						if(type2or3Edge.head.find() == wke.head.find())
							wke.markerName.setParentOrient(Marker.head);
						else
							wke.markerName.setParentOrient(Marker.tail);
					}
				} else if(m.type2or3.size() == 1) {
					type2or3Edge = (Edge)m.type2or3.get(0);
					type2or3Edge.markerName.setParentOrient(Marker.head);
					if(!rule3(m, type2or3Edge.tail)) return false;
				} else if(m.type4.size() == 1) {
					;
				} else {
					wke = m.parentMarkerEdge();
					if(!rule2(m, wke.tail.find(), wke.head.find())) return false;
				}
				return true;
				
			case Member.prime:
				if(m.pathEdge.size() == 0) {
					if(m.type2or3.size() == 2) {
						type2or3Edge = (Edge)m.type2or3.get(0);
						wke = (Edge)m.type2or3.get(1);
						if(type2or3Edge.head.find() == wke.head.find() ) {
							type2or3Edge.markerName.setParentOrient(Marker.head);
							wke.markerName.setParentOrient(Marker.head);
						} else if(type2or3Edge.head.find() == wke.tail.find() ) {
							type2or3Edge.markerName.setParentOrient(Marker.head);
							wke.markerName.setParentOrient(Marker.tail);
						} else if(type2or3Edge.tail.find() == wke.head.find() ) {
							type2or3Edge.markerName.setParentOrient(Marker.tail);
							wke.markerName.setParentOrient(Marker.head);
						} else 	if(type2or3Edge.tail.find() == wke.tail.find() ) {
							type2or3Edge.markerName.setParentOrient(Marker.tail);
							wke.markerName.setParentOrient(Marker.tail);
						} else return false;
						return true;
					} else return false;
				}
				for(i=0; i<m.pathEdge.size(); i++) {
					wke = (Edge)m.pathEdge.get(i);
					wkn = wke.head.find();
					if(wkn.degree++ == 0) pathNode.push(wkn);
					wkn = wke.tail.find();
					if(wkn.degree++ == 0) pathNode.push(wkn);
				}
				while(!pathNode.isEmpty()) {
					wkn = (Node)pathNode.pop();
					if(wkn.degree > 2) return false;
					if(wkn.degree == 1) {
						wkn.degree = 0;
						pathEnd.add(wkn);
					} else
						wkn.degree = 0;
				}
				if(pathEnd.size() == 2) {
					end[0] = (Node)pathEnd.get(0);
					end[1] = (Node)pathEnd.get(1);
					if(m.type4.size() == 1) {
						if(!((Edge)m.type4.get(0)).isEnd(end[0], end[1])) return false;
						return true;
					} else if(m.type2or3.size() == 2) {
						endEdge[0] = (Edge)m.type2or3.get(0);
						endEdge[1] = (Edge)m.type2or3.get(1);
						for(i=0; i<2; i++) {
							for(j=0; j<2; j++) {
								if(endEdge[i].head.find() == end[j] || endEdge[i].tail.find() == end[j]) {
									if(endEdge[i].head.find() == end[j]) {
										if(endEdge[i].tail.find() == end[1-j]) continue;
										endEdge[i].markerName.setParentOrient(Marker.head);
									} else {
										if(endEdge[i].head.find() == end[1-j]) continue;
										endEdge[i].markerName.setParentOrient(Marker.tail);
									}
									if(endEdge[1-i].head.find() == end[1-j])
										endEdge[1-i].markerName.setParentOrient(Marker.head);
									else if(endEdge[1-i].tail.find() == end[1-j])
										endEdge[1-i].markerName.setParentOrient(Marker.tail);
									else return false;
									return true;
								}
							}
						}
						if(!endEdge[0].isSameEnd(endEdge[1])) return false;
						if(!endEdge[0].isEnd(end[0], end[1])) return false;
						if(endEdge[0].head.find() == endEdge[1].head.find()) {
							endEdge[0].markerName.setParentOrient(Marker.head);
							endEdge[1].markerName.setParentOrient(Marker.tail);
						} else if(endEdge[0].head.find() == endEdge[1].tail.find()) {
							endEdge[0].markerName.setParentOrient(Marker.head);
							endEdge[1].markerName.setParentOrient(Marker.head);
						} else return false;
						return true;
					}
					else if(m.type2or3.size() == 1) {
						type2or3Edge = (Edge)m.type2or3.get(0);
						if(type2or3Edge.tail.find() == end[0]) {
							nodeForR3 = end[1];
							type2or3Edge.markerName.setParentOrient(Marker.tail);
						} else if(type2or3Edge.head.find() == end[0]) {
							nodeForR3 = end[1];
							type2or3Edge.markerName.setParentOrient(Marker.head);
						} else if(type2or3Edge.tail.find() == end[1]) {
							type2or3Edge.markerName.setParentOrient(Marker.tail);
							nodeForR3 = end[0];
						} else if(type2or3Edge.head.find() == end[1]) {
							type2or3Edge.markerName.setParentOrient(Marker.head);
							nodeForR3 = end[0];
						} else return false;
						if(!rule3(m, nodeForR3)) return false;
						return true;
					} else {
						if(!rule2(m, end[0], end[1])) return false;
						return true;
					}
				} else return false;
		}
		return false;
	}
//----------------------------------------------------------
	boolean rule1(Member m, Node n) {
		if(u[0] == null) {
			k[0] = m;
			u[0] = n;
		} else if(u[1] == null) {
			k[1] = m;
			u[1] = n;
		} else return false;
		return true;
	}
//----------------------------------------------------------
	boolean rule2(Member m, Node n1, Node n2) {
		if(u[0] != null) return false;
		else {
			k[0] = m;
			k[1] = m;
			u[0] = n1;
			u[1] = n2;
			return true;
		}
	}
//----------------------------------------------------------
	boolean rule3(Member m, Node n) {
		Stack testMember = new Stack(); // Member
		Member wkm = k[0];
		Node wkn = n;
		if(u[1] != null) return false;
		while(wkm != m) {
			testMember.push(wkm);
			wkm = wkm.getParent();
		}
		while(!testMember.isEmpty()) {
			wkm = (Member)testMember.pop();
			if(wkm.parentMarker.parent.head.find() == wkn) {
				if(wkm.parentMarker.orient)
					wkn = wkm.parentMarkerEdge().tail.find();
				else
					wkn = wkm.parentMarkerEdge().head.find();
			}else if(wkm.parentMarker.parent.tail.find() == wkn) {
				if(wkm.parentMarker.orient)
					wkn = wkm.parentMarkerEdge().head.find();
				else
					wkn = wkm.parentMarkerEdge().tail.find();			
			}else {
				k[1] = wkm.getParent();
				u[1] = wkn;
				return true;
			}
		}
		k[1] = wkm;
		u[1] = wkn;
		return true;
	}
//----------------------------------------------------------
	void rule5() {
		Edge wke;
		if(k[0] == k[1]) {
			if(k[0].parentMarkerEdge().isEnd(u[0], u[1]) 
			   && k[0].graphType != Member.bond) {
				u[0] = k[0].parentMarker.parent.tail.find();
				u[1] = k[0].parentMarker.parent.head.find();
				k[0] = k[0].getParent();
				k[1] = k[0];
			} else if((k[0].graphType == Member.polygon)
			          && (u[0].isJoined(u[1]))
			          && (u[0].jointEdge(u[1]).name == -1)
			          && (u[0].jointEdge(u[1]).markerName.child.getMember().graphType == Member.bond)) {
			  wke = u[0].jointEdge(u[1]).markerName.child;
				u[0] = wke.head.find();
				u[1] = wke.tail.find();
				k[0] = wke.getMember();
				k[1] = k[0];
			}
		}
	}
//----------------------------------------------------------
	boolean addType234(Member m, Node target1, Node target2, boolean t2Type, byte count) {
		int i;
		Edge wke;
		Edge[] wkeA = new Edge[2];
		Node wkn;
		if(m.type4.size() + m.type2or3.size() == 0) {
			if(count == 0) return false;
			else if(count == 2) {
				if(m.parentMarkerEdge().isEnd(target1, target2)) {
					m.getParent().addPathEdge(m.parentMarker.parent);
				} else {
					m.getParent().type2or3.add(m.parentMarker.parent);
					m.parentMarker.setChildOrient(!t2Type);
					if(!rule1(m, target1)) return false;
				}
			}
			else if(count == 4) {
				m.getParent().type4.add(m.parentMarker.parent);
				if(!rule2(m, target1, target2)) return false;
			}
		}
		else if(m.type4.size() == 1) {
			if(((Edge)m.type4.get(0)).isEnd(target1, target2)) {
				m.getParent().type4.add(m.parentMarker.parent);
			} else return false;
		}
		else if(m.type2or3.size() == 1) {
			wke = (Edge)m.type2or3.get(0);
			if(wke.head.find() == target1 || wke.tail.find() == target1) {
				if(wke.head.find() == target1)
					wke.markerName.setParentOrient(Marker.head);
				else
					wke.markerName.setParentOrient(Marker.tail);
				if(count == 0) {
					m.getParent().type2or3.add(m.parentMarker.parent);
					m.parentMarker.setChildOrient(!t2Type);
				}
				else if(count == 2) {
					m.getParent().type2or3.add(m.parentMarker.parent);
					if(m.parentMarkerEdge().isEnd(target1, target2)) {
						m.parentMarker.setChildOrient(t2Type);
					} else {
						m.parentMarker.setChildOrient(!t2Type);
					}
				}
				else if(count == 4) {
					m.getParent().type4.add(m.parentMarker.parent);
					if(!rule3(m, target2)) return false;
				}
			} else if(wke.head.find() == target2 || wke.tail.find() == target2) {
				if(wke.head.find() == target2)
					wke.markerName.setParentOrient(Marker.head);
				if(wke.tail.find() == target2)
					wke.markerName.setParentOrient(Marker.tail);
				if(count == 0) {
					m.getParent().type2or3.add(m.parentMarker.parent);
					m.parentMarker.setChildOrient(t2Type);
				} else if(count == 2) {
					if(m.parentMarkerEdge().isEnd(target1, target2)) {
						m.getParent().type2or3.add(m.parentMarker.parent);
						m.parentMarker.setChildOrient(!t2Type);
					} else {
						m.getParent().type4.add(m.parentMarker.parent);
						if(!rule3(m, target1)) return false;
					}
				} else if(count == 4) {
					m.getParent().type4.add(m.parentMarker.parent);
					if(!rule3(m, target1)) return false;
				}
			} else return false;
		}
		else if(m.type2or3.size() == 2) {
			m.getParent().type4.add(m.parentMarker.parent);
			wkeA[0] = (Edge)m.type2or3.get(0);
			wkeA[1] = (Edge)m.type2or3.get(1);
			for(i=0; i<2; i++) {
				if(wkeA[i].head.find() == target1 || wkeA[i].tail.find() == target1) {
					if(wkeA[i].head.find() == target1) {
						if(wkeA[i].tail.find() == target2) continue;
						wkeA[i].markerName.setParentOrient(Marker.head);
					}
					else {
						if(wkeA[i].head.find() == target2) continue;
						wkeA[i].markerName.setParentOrient(Marker.tail);
					}
					if(wkeA[1-i].head.find() == target2)
						wkeA[1-i].markerName.setParentOrient(Marker.head);
					else if(wkeA[1-i].tail.find() == target2)
						wkeA[1-i].markerName.setParentOrient(Marker.tail);
					else return false;
					return true;
				} else if(wkeA[i].head.find() == target2 || wkeA[i].tail.find() == target2) {
					if(wkeA[i].head.find() == target2) {
						if(wkeA[i].tail.find() == target1) continue;
						wkeA[i].markerName.setParentOrient(Marker.head);
					}
					else {
						if(wkeA[i].head.find() == target1) continue;
						wkeA[i].markerName.setParentOrient(Marker.tail);
					}
					if(wkeA[1-i].head.find() == target1)
						wkeA[1-i].markerName.setParentOrient(Marker.head);
					else if(wkeA[1-i].tail.find() == target1)
						wkeA[1-i].markerName.setParentOrient(Marker.tail);
					else return false;
					return true;
				} 
				else return false;
			}
			if(wkeA[0].head.find() == wkeA[1].head.find()) {
				wkeA[0].markerName.setParentOrient(Marker.head);
				wkeA[1].markerName.setParentOrient(Marker.tail);
			} else {
				wkeA[0].markerName.setParentOrient(Marker.head);
				wkeA[1].markerName.setParentOrient(Marker.head);
			}
		}
		return true;
	}
//----------------------------------------------------------
}
////////////////////////////////////////////////////////////
class Member{
	int sizeForUF = 1;
	Member parentForUF;
	int depth = -1;
	byte graphType;
		static final byte polygon = 0;
		static final byte bond = 1;
		static final byte prime = 2;
	int edgeSize; // bond
	Marker parentMarker;
	ArrayList pathEdge = new ArrayList(); // Edge
	ArrayList type2or3 = new ArrayList(); // Edge
	ArrayList type4 = new ArrayList(); // Edge
	boolean findRootFlag = false;
	Edge rootParentMarkerEdge;
	int name; //debag
//----------------------------------------------------------
	public Member(Edge markerEdge, ArrayList p, Edge[] tree, ArrayList coTreeEdge){ // U0
		int i,j;
		Edge e;
		Node n, wn;
		
		graphType = Member.polygon;
		markerEdge.setMember(this);
		parentMarker = markerEdge.markerName;
		n = new Node();
		markerEdge.head = n;
		n.tail = markerEdge;
		for(i=0; i<p.size(); i++){
			wn = n;
			j = ((Integer)p.get(i)).intValue();
			e = new Edge(j, this);
			if(j < tree.length) tree[j] = e;
			else coTreeEdge.add(e);
			n = new Node();
			e.head = n;
			e.tail = wn;
			wn.head = e;
			n.tail = e;
		}
		markerEdge.tail = n;
		n.head = markerEdge;
		parentForUF = this;
	}
//----------------------------------------------------------	
	public Member(Edge f, Edge f1, Edge f2) {
		Node n1 = new Node();
		Node n2 = new Node();
		boolean headToTail;
		
		if(f1.head != null && f2.head != null && f1.head.find() == f2.tail.find()) 
			headToTail = true;
		else 
			headToTail = false;
		f.memberName = this;
		f1.memberName = this;
		f2.memberName = this;
		f.head = n1;
		f.tail = n2;
		f1.head = n1;
		f1.tail = n2;
		if(headToTail) {
			f2.head = n2;
			f2.tail = n1;
		} else {
			f2.head = n1;
			f2.tail = n2;	
		}
		graphType = Member.bond;
		edgeSize = 3;
		parentMarker = f.markerName;
		parentForUF = this;
	}
//----------------------------------------------------------	
	public Member(Edge f1, Edge f2) {
		Node n1 = new Node();
		Node n2 = new Node();
		
		f1.memberName = this;
		f2.memberName = this;
		f1.head = n1;
		f2.head = n1;
		f1.tail = n2;
		f2.tail = n2;
		graphType = Member.bond;
		edgeSize = 2;
		rootParentMarkerEdge = f1;
		parentForUF = this;
	}
//----------------------------------------------------------	
	public Member() {
		graphType = Member.polygon;
		parentForUF = this;
	}
//----------------------------------------------------------
	Edge parentMarkerEdge() {
		if(parentMarker == null)
		 return rootParentMarkerEdge;
		return parentMarker.child;
	}
//----------------------------------------------------------
	Member find(){
		Member wk,j, root;
		wk = this;
		while(wk.parentForUF != wk)
			wk = wk.parentForUF;
		root = wk;
		wk = this;
		while(wk != root){
			j = wk.parentForUF;
			wk.parentForUF = root;
			wk = j;
		}
		
		return root;
	}
//----------------------------------------------------------
	void union(Member m, byte gtype){
		Member m1,m2;
		m1 = this.find();
		m2 = m.find();
		if(m1.sizeForUF > m2.sizeForUF) {
			m2.parentForUF = m1;
			m1.sizeForUF += m2.sizeForUF;
			m1.parentMarker = m2.parentMarker;
			if(gtype == Member.prime) 
				m1.graphType = Member.prime;
			else if(gtype == Member.bond) 
				m1.edgeSize += m2.edgeSize -2;
		}
		else{
			m1.parentForUF = m2;
			m2.sizeForUF += m1.sizeForUF;
			if(gtype == Member.prime)
				m2.graphType = Member.prime;
			else if(gtype == Member.bond) 
				m2.edgeSize += m1.edgeSize -2;
		}
	}
//----------------------------------------------------------
	Member getParent() {
		if(parentMarker == null) return null;
		return parentMarker.parent.getMember();
	}
//----------------------------------------------------------
	void addPathEdge(Edge e) {
		pathEdge.add(e);
		e.pathEdgeFlag = true;
	}
//----------------------------------------------------------
	boolean isParentMarker(Node s, Node t) {
		return isParentMarker(s.head, t.tail);
	}
//----------------------------------------------------------
	boolean isParentMarker(Node s, Edge t) {
		return isParentMarker(s.head, t);
	}
//----------------------------------------------------------
	boolean isParentMarker(Edge s, Edge t) {
		Edge wke, pm;
		wke = s.nextEdge();
		pm = t.getMember().parentMarker.child;
		while(wke != t) {
			if(wke == pm) return true;
			wke = wke.nextEdge();
		}
		return false;
	}
//----------------------------------------------------------
	void divide(Node s, Node t, Edge f) {
		Member newMember = new Member();
		Marker m1 = new Marker(), m2 = new Marker();
		
		if(isParentMarker(s, t)) {
			makeNewPolygon(s, t, m2.parent, m1.child, newMember, true);
			newMember.parentMarker = this.parentMarker;
			this.parentMarker = m1;
			new Member(m2.child, f, m1.parent);
		} else {
			makeNewPolygon(s, t, m2.child, m1.parent, newMember, true);
			newMember.parentMarker = m2;
			new Member(m1.child, f, m2.parent);		
		}
	}
//----------------------------------------------------------
	void squeeze(Node n) {
		Edge pMarker = parentMarkerEdge();
		Marker newMarker;
		Member newMember;
		Edge endEdge;
		
		if(graphType == Member.polygon) {
			if(!((pMarker.tail == n) || (pMarker.nextEdge().tail == n))) {
				newMember = new Member();
				newMarker = new Marker();
				makeNewPolygon(pMarker.tail, n, newMarker.child, newMarker.parent, newMember, false);
				newMember.parentMarker = newMarker;
			}
			if(!((pMarker.head == n) || (pMarker.preEdge().head == n))) {
				newMember = new Member();
				newMarker = new Marker();
				endEdge = n.head;
				makeNewPolygon(pMarker.head, n, newMarker.parent, newMarker.child, newMember, true);
				newMember.parentMarker = this.parentMarker;
				this.parentMarker = newMarker;
				n = endEdge.tail;
			}
		}
	}
//----------------------------------------------------------
	void rootSqueeze(Node n, Edge e1) {
		Edge pMarker = parentMarkerEdge();
		Marker newMarker;
		Member newMember;
		boolean isParent;
		
		if(graphType == Member.polygon) {
			isParent = this.isParentMarker(n, e1);
			if(!((n.tail == e1) || (n.tail.nextEdge() == e1))) {
				newMember = new Member();
				newMarker = new Marker();
				if(isParent) {
					makeNewPolygon(n, e1.head, newMarker.parent, newMarker.child, newMember, false);
					newMember.parentMarker = this.parentMarker;
					this.parentMarker = newMarker;
				} else {
					makeNewPolygon(n, e1.head, newMarker.child, newMarker.parent, newMember, false);
					newMember.parentMarker = newMarker;
				}
			}
			if(!((n.head == e1) || (n.head.preEdge() == e1))) {;
				newMember = new Member();
				newMarker = new Marker();
				if(isParent) {
					makeNewPolygon(n, e1.tail, newMarker.parent, newMarker.child, newMember, true);
					newMember.parentMarker = this.parentMarker;
					this.parentMarker = newMarker;
				} else {
					makeNewPolygon(n, e1.tail, newMarker.child, newMarker.parent, newMember, true);
					newMember.parentMarker = newMarker;
				}
			}
		}
	}
//----------------------------------------------------------
	void squeeze(Edge cMarker) {
		Edge pMarker = parentMarkerEdge();
		Marker newMarker;
		Member newMember;
		Edge wke;int i=5;//debag
		if(graphType == Member.polygon) {
			if(!((pMarker.nextEdge() == cMarker) || (pMarker.nextEdge().nextEdge() == cMarker))) {
				newMember = new Member();
				newMarker = new Marker();
				makeNewPolygon(pMarker.tail, cMarker.head, 
												newMarker.child, newMarker.parent, newMember, false);
				newMember.parentMarker = newMarker;
			}
			if(!((pMarker.preEdge() == cMarker) || (pMarker.preEdge().preEdge() == cMarker))) {
				newMember = new Member();
				newMarker = new Marker();
				makeNewPolygon(pMarker.head, cMarker.tail, 
												newMarker.parent, newMarker.child, newMember, true);
				newMember.parentMarker = this.parentMarker;
				this.parentMarker = newMarker;
			}
		}
	}
//----------------------------------------------------------
	void rootSqueeze(Edge e1, Edge e2) {
		Edge pMarker = parentMarkerEdge(), s, t, wke1, wke2;
		Marker newMarker;
		Member newMember;
		boolean isParent;
		
		if(graphType == Member.polygon) {
			wke1 = e1; wke2 = e2;
			while((wke1 != e2) && (wke2 != e1)) {
				wke1 = wke1.nextEdge();
				wke2 = wke2.nextEdge();
			}
			if(wke1 == e2) {
				s = e1;
				t = e2;
			} else {
				s = e2;
				t = e1;
			}
			isParent = this.isParentMarker(s, t);
			if(!((s.nextEdge() == t) || (s.nextEdge().nextEdge() == t))) {
				newMember = new Member();
				newMarker = new Marker();
				if(isParent) {
					makeNewPolygon(s.tail,t.head, newMarker.parent, newMarker.child, 
													newMember, false);
					newMember.parentMarker = this.parentMarker;
					this.parentMarker = newMarker;
				} else {
					makeNewPolygon(s.tail, t.head, newMarker.child, newMarker.parent, 
													newMember, false);
					newMember.parentMarker = newMarker;
				}
			}
			if(!((s.preEdge() == t) || (s.preEdge().preEdge() == t))) {
				newMember = new Member();
				newMarker = new Marker();
				if(isParent) {
					makeNewPolygon(s.head, t.tail, newMarker.parent, newMarker.child, 
													newMember, true);
					newMember.parentMarker = this.parentMarker;
					this.parentMarker = newMarker;
				} else {
					makeNewPolygon(s.head, t.tail, newMarker.child, newMarker.parent, 
													newMember, true);
					newMember.parentMarker = newMarker;
				}
			}
		}
	}
//----------------------------------------------------------
	void makeNewPolygon(Node h, Node t, Edge edgeInNewMember,	Edge edgeInOldMember, 
											Member newMember, boolean changeEnd) {
		Edge h1, h2, t1, t2, wke;
		Node newNode1 = new Node();
		Node newNode2 = new Node();
		
		h1 = h.tail;
		t2 = h.head;
		t1 = t.head;
		h2 = t.tail;
		wke = h1;
		while(wke != h2) {
			wke.memberName = newMember;
			wke = wke.nextEdge();
		}
		if(changeEnd){
			t1.setNextEdge(edgeInNewMember);
			h1.setPreEdge(edgeInNewMember);
			t2.setNextEdge(newNode1, edgeInOldMember);
			h2.setPreEdge(newNode2, edgeInOldMember);
		} else {
			t1.setNextEdge(newNode1, edgeInNewMember);
			h1.setPreEdge(newNode2, edgeInNewMember);
			t2.setNextEdge(edgeInOldMember);
			h2.setPreEdge(edgeInOldMember);			
		}
		edgeInNewMember.memberName = newMember;
		edgeInOldMember.memberName = this;
	}
//----------------------------------------------------------	
}
////////////////////////////////////////////////////////////
class Marker{
	Edge parent = new Edge();
	Edge child = new Edge();
	boolean orient; // false:head to head, true:head to tail
	boolean childOrient;
		static final boolean head = false;
		static final boolean tail = true;
//----------------------------------------------------------
	public Marker(){
		parent.markerName = this;
		child.markerName = this;
	}
//----------------------------------------------------------
	void setChildOrient(boolean b) {
		childOrient = b;
	}
//----------------------------------------------------------
	void setParentOrient(boolean parentOrient) {
		if(childOrient ^ parentOrient) orient = true;
		else orient = false;
	}
//----------------------------------------------------------
	void union(byte gtype) {
		child.memberName.union(parent.memberName, gtype);
		if(orient) {
			child.head.union(parent.tail);
			child.tail.union(parent.head);
		} else {
			child.head.union(parent.head);
			child.tail.union(parent.tail);
		}
	}
//----------------------------------------------------------
}
//////////////////////////////////////////////////////////// 
class Edge{
	Member memberName;
	Node head;
	Node tail;
	int name; // -1: marker
	Marker markerName;
	boolean pathEdgeFlag = false; // false:No, true:tree or type 1
	int edgeName;
//----------------------------------------------------------
	public Edge(){
		name = -1;
	}
//----------------------------------------------------------
	public Edge(int n){
		name = n;
	}
//----------------------------------------------------------
	public Edge(Member m){
		memberName = m;
	}
//----------------------------------------------------------
	public Edge(int n, Member m){
		name = n;
		memberName = m;
	}
//----------------------------------------------------------
	void setEndNode(Node h, Node t){
		head = h;
		tail = t;
	}
//----------------------------------------------------------
	void setMember(Member m){
			memberName = m;
	}
//----------------------------------------------------------
	Member getMember() {
		return memberName.find();
	}
//----------------------------------------------------------
	Edge nextEdge() {
		return this.tail.tail;
	}
//----------------------------------------------------------
	Edge preEdge() {
		return this.head.head;
	}
//----------------------------------------------------------
	void swap(Edge m) { // for polygon
		Node this_h, this_t, m_h, m_t;
		
		this_h = this.head;
		this_t = this.tail;
		m_h = m.head;
		m_t = m.tail;
		
		this_h.tail = m;
		this_t.head = m;
		m_h.tail = this;
		m_t.head = this;
		m.head = this_h;
		m.tail = this_t;
		this.head = m_h;
		this.tail = m_t;
	}
//----------------------------------------------------------
	boolean isEnd(Node n1, Node n2) {
		Node wkn1 = tail.find();
		Node wkn2 = head.find();
		return (((wkn1 == n1) && (wkn2 == n2)) ||
		        ((wkn1 == n2) && (wkn2 == n1)));
	}
//----------------------------------------------------------
	boolean isSameEnd(Edge e) {
		Node n1t = this.tail.find();
		Node n1h = this.head.find();
		Node n2t = e.tail.find();
		Node n2h = e.head.find();
		if(((n1t == n2t) && (n1h == n2h))
		   ||((n1t == n2h) && (n1h == n2t)))
		  return true;
		return false;
	}
//----------------------------------------------------------
	void setNextEdge(Edge targetEdge) {
		
		this.tail.tail = targetEdge;
		targetEdge.head = this.tail;
	}
//----------------------------------------------------------
	void setNextEdge(Node targetNode, Edge targetEdge) {
		targetNode.head = this;
		targetNode.tail = targetEdge;
		targetEdge.head = targetNode;
		this.tail = targetNode;
	}
//----------------------------------------------------------
	void setPreEdge(Edge targetEdge) {
		this.head.head = targetEdge;
		targetEdge.tail = this.head;
	}
//----------------------------------------------------------
	void setPreEdge(Node targetNode, Edge targetEdge) {
		targetNode.tail = this;
		targetNode.head = targetEdge;
		targetEdge.tail = targetNode;
		this.head = targetNode;
	}
//----------------------------------------------------------
}
////////////////////////////////////////////////////////////
class Node{
	int sizeForUF, name = 0;
	int degree = 0;
	Node parentForUF;
	Edge head, tail;
	Stack edgeStack = new Stack(); // Edge
	ArrayList endEdgeArrayList = new ArrayList(); // Edge
//----------------------------------------------------------
	public Node(){
		sizeForUF = 1;
		parentForUF = this;
	}
//----------------------------------------------------------
	Node find(){
		Node wk, j, root;
		wk = this;
		while(wk.parentForUF != wk) {
			wk = wk.parentForUF;
		}
		root = wk;
		wk = this;
		while(wk != root){
			j = wk.parentForUF;
			wk.parentForUF = root;
			wk = j;
		}
		
		return root;
	}
//----------------------------------------------------------
	void union(Node m){
		Node m1,m2;
		m1 = this.find();
		m2 = m.find();
		if(m1.sizeForUF > m2.sizeForUF) {
			m2.parentForUF = m1;
			m1.sizeForUF += m2.sizeForUF;
		}
		else{
			m1.parentForUF = m2;
			m2.sizeForUF += m1.sizeForUF;
		}
	}
//----------------------------------------------------------
	boolean isJoined(Node n) { // Test if two vertices are adjacent.
		if(((this.head.head) == n) || (this.tail.tail) == n)
			return true;
		return false;
	}
//----------------------------------------------------------
	Edge jointEdge(Node n) { // Return an edge between this and n.
		if(this.tail.tail == n) 
			return this.tail;
		if(this.head.head == n)
			return this.head;
		return null;
	}
//----------------------------------------------------------
	boolean isAnotherEnd(Node testNode) {
		Node wkn, pren;
		Edge wke;
		
		wkn = this;
		wke = (Edge)wkn.endEdgeArrayList.get(0);
		if(wke.head.find() == wkn) wkn = wke.tail.find();
		else wkn = wke.head.find();
		pren = wkn;
		while(wkn.endEdgeArrayList.size() == 2) {
			if((Edge)wkn.endEdgeArrayList.get(0) == wke)
				wke = (Edge)wkn.endEdgeArrayList.get(1);
			else wke = (Edge)wkn.endEdgeArrayList.get(0);
			if(wke.head.find() == wkn) wkn = wke.tail.find();
			else wkn = wke.head.find();
			pren = wkn;
		}
		if(wkn == testNode) return true;
		else return false;
	}
//----------------------------------------------------------
}
////////////////////////////////////////////////////////////
