This is the reference implementation of the graph realisation
algorithm used in `graphic.rs`. It was retrieved from the
author's website in 2026, and included in this repo for
posterity. See:

https://www.math.keio.ac.jp/~kakimura/GRP/

****************************************************************


<How to Use Graph Realization Program>  


This is a JAVA implementation of the Bixby-Wagner algorithm for 
the graph relaization problem. See the following paper for the 
algorithm and its verification. 

R. Bixby and D. Wagner: 
An almost linear-time algorithm for graph realization, 
Mathematics of Operations Research, Vol.13 (1988), pp.99-123. 


1. How to Compile the Program 

In the directory where GRP.java exists, execute  

>javac GRP.java

to produce the following seven files: 

GRP.class, TDecomposition.class, MatrixDivide.class,
Member.class, Node.class, Edge.class, Marker.class.

* if you use jdk 1.5 or later, compile this program as follows:
> javac -target 1.4 -source 1.4 GRP.java


2. How to Make an Input File.

Write the maximum number of tree edges in the first line.
The subsequent lines correspond to fundamental circuits. 
In each line, enumerate all the tree arcs in the fundamental 
circuit. Their should be a blank space between tree arcs. 
You can put the fundamental circuits in an arbitrary order. 

Example: 
------------------------------------------
10
1 2 3 4
5 6 7
8 9
2 3 9
8 1
------------------------------------------


3. How to Run the Program  

Put the input file (e.g., input.txt) in the same directory 
as GRP.class, and execute 

>java GRP input.txt

to obtain an output in console. If you want to 
save the result, execute 
 
>java GRP input.txt > output.txt

to obtain an output file (i.e., output.txt). 



4. How to Read the Output

The output contains the number of sequentially connected components
and the adjacency lists. 

Example 1 (graphic case):
------------------------------------------
Graphic. BlockSize is 2.
1: (1, 2)
2: (3, 4)
3: (4, 1)
4: (5, 3)
5: (6, 7)
6: (8, 6)
7: (9, 8)
8: (1, 10)
9: (1, 11)
11: (2, 5)
12: (10, 2)
13: (11, 3)
14: (11, 10)
15: (7, 9)
------------------------------------------

The first line: The number of sequentially connected components. 

Each subsequent line: 
[name of arc]: ([name of end-vertex], [name of end-vertex])
Remark: The order of the non-tree arcs can be different from 
the order of in the input file. 



Example 2 (nongraphic case): 
------------------------------------------
NonGraphic. 
------------------------------------------


5. LICENSE

Copyright is reserved by Takahiro Ohto.
You may redistribute with or without modification under (modified)
BSD license.  See 'LICENSE' file for more detailed information.


=======================================================================
The code was written by Takahiro Ohto in 2001-2003, when he was 
a graduate course student at University of Tokyo. For the bug 
information and other questions, please take contact to: 

Satoru Iwata 
Research Institute for Mathematical Sciences 
Kyoto University 
Tel: +81 75 753 7236 Fax: +81 75 753 7272 
E-mail: iwata@kurims.kyoto-u.ac.jp 
=======================================================================












