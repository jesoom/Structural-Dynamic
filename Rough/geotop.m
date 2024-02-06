%
% Function geotop: defines geometry and topology data
%
function [nInc,nElements,dXY,nNodes]=geotop

% Matrix of the nodal coordinates dXY:
% the n-th row in dXY collects the nodal coordinates of node i-th 
% dXY(n,:)=[x coordinate of n-th node,  y coordinate of n-th node]

% Force -> N    , Length= m


 
    
   dXY=[
       0.000	0.000	;
5.721	0.000	;
11.410	0.000	;
17.115	0.000	;
22.820	0.000	;
28.525	0.000	;
0.854	3.417	;
5.845	3.417	;
10.829	3.417	;
15.817	3.417	;
20.804	3.417	;
25.791	3.417	;
1.709	6.834	;
5.981	6.834	;
10.248	6.834	;
14.518	6.834	;
18.788	6.834	;
23.058	6.834	;
2.563	10.250	;
6.117	10.250	;
9.668	10.250	;
13.220	10.250	;
16.773	10.250	;
20.325	10.250	;
3.417	13.667	;
6.253	13.667	;
9.087	13.667	;
11.922	13.667	;
14.757	13.667	;
17.591	13.667	;
4.271	17.084	;
6.389	17.084	;
8.506	17.084	;
10.623	17.084	;
12.740	17.084	;
14.858	17.084	;
5.125	20.500	;
6.525	20.500	;
7.925	20.500	;
9.325	20.500	;
10.725	20.500	;
12.125	20.500	;
5.125	23.500	;
6.358	23.500	;
7.592	23.500	;
8.825	23.500	;
10.058	23.500	;
11.292	23.500	;
5.125	26.500	;
6.191	26.500	;
7.258	26.500	;
8.325	26.500	;
9.392	26.500	;
10.458	26.500	;
5.125	29.500	;
6.026	29.500	;
6.925	29.500	;
7.825	29.500	;
8.725	29.500	;
9.625	29.500	;
5.125	32.000	;
6.026	32.000	;
6.925	32.000	;
7.825	32.000	;
8.725	32.000	;
9.625	32.000	;
5.125	34.500	;
6.026	34.500	;
6.925	34.500	;
7.825	34.500	;
8.725	34.500	;
9.625	34.500	;

        ];

      
    
        
 %Total number of nodes
  [nNodes,~]=size(dXY);  

 % Connection matrix nInc
 % the n-th row of nInc contains the node numbers for the n-th finite element (INPUT DATA)
 % and the corresponding dofs (evaluated from the former INPUT)
 % nInc(ne,:)=[n1, n2, n3, n4, n1u, n1v, n2u, n2v, n3u, n3v, n4u, n4v] 
 
 
nInc=zeros(55,4);


% nInc(1:68,:)=[ nout(1:end-1)' , ninn2(1:end-1)' , ninn2(2:end)',  nout(2:end)'];
% nInc(69:136,:)=[ninn2(1:end-1)' , ninn1(1:end-1)' , ninn1(2:end)', ninn2(2:end)'];
% nInc(137:204,:)=[ninn1(1:end-1)' ,  ninn(1:end-1)' ,  ninn(2:end)', ninn1(2:end)'];


  nInc=[
      1	2	8	7	;
2	3	9	8	;
3	4	10	9	;
4	5	11	10	;
5	6	12	11	;
7	8	14	13	;
8	9	15	14	;
9	10	16	15	;
10	11	17	16	;
11	12	18	17	;
13	14	20	19	;
14	15	21	20	;
15	16	22	21	;
16	17	23	22	;
17	18	24	23	;
19	20	26	25	;
20	21	27	26	;
21	22	28	27	;
22	23	29	28	;
23	24	30	29	;
25	26	32	31	;
26	27	33	32	;
27	28	34	33	;
28	29	35	34	;
29	30	36	35	;
31	32	38	37	;
32	33	39	38	;
33	34	40	39	;
34	35	41	40	;
35	36	42	41	;
37	38	44	43	;
38	39	45	44	;
39	40	46	45	;
40	41	47	46	;
41	42	48	47	;
43	44	50	49	;
44	45	51	50	;
45	46	52	51	;
46	47	53	52	;
47	48	54	53	;
49	50	56	55	;
50	51	57	56	;
51	52	58	57	;
52	53	59	58	;
53	54	60	59	;
55	56	62	61	;
56	57	63	62	;
57	58	64	63	;
58	59	65	64	;
59	60	66	65	;
61	62	68	67	;
62	63	69	68	;
63	64	70	69	;
64	65	71	70	;
65	66	72	71	;

];

   nInc=[
      1	2	8	7	;
2	3	9	8	;
3	4	10	9	;
4	5	11	10	;
5	6	12	11	;
7	8	14	13	;
8	9	15	14	;
9	10	16	15	;
10	11	17	16	;
11	12	18	17	;
13	14	20	19	;
14	15	21	20	;
15	16	22	21	;
16	17	23	22	;
17	18	24	23	;
19	20	26	25	;
20	21	27	26	;
21	22	28	27	;
22	23	29	28	;
23	24	30	29	;
25	26	32	31	;
26	27	33	32	;
27	28	34	33	;
28	29	35	34	;
29	30	36	35	;
31	32	38	37	;
32	33	39	38	;
33	34	40	39	;
34	35	41	40	;
35	36	42	41	;
37	38	44	43	;
38	39	45	44	;
39	40	46	45	;
40	41	47	46	;
41	42	48	47	;
43	44	50	49	;
44	45	51	50	;
45	46	52	51	;
46	47	53	52	;
47	48	54	53	;
49	50	56	55	;
50	51	57	56	;
51	52	58	57	;
52	53	59	58	;
53	54	60	59	;
55	56	62	61	;
56	57	63	62	;
57	58	64	63	;
58	59	65	64	;
59	60	66	65	;
61	62	68	67	;
62	63	69	68	;
63	64	70	69	;
64	65	71	70	;
65	66	72	71	;

];
   nInc=[nInc,nInc(:,1)*2-1,nInc(:,1)*2,nInc(:,2)*2-1,nInc(:,2)*2,nInc(:,3)*2-1,nInc(:,3)*2,nInc(:,4)*2-1,nInc(:,4)*2];
 %Total number of plane elements in the structure
  [nElements,~]=size(nInc);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%