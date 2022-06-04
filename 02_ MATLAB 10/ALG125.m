% Finite Element Algorithm 12.5
%
% To approximate the solution to an elliptic partial-differential
% equation subject to Dirichlet, mixed, or Neumann boundary
% conditions:
%
% Input:   see STEP 0
%
% Output:  description of triangles, nodes, line integrals, basis
%          functions, linear system to be solved, and the
%          coefficients of the basis functions
%
%
% Step 0
% General outline
%
%    1. Triangles numbered: 1 to K for triangles with no edges on
%       Script-S-1 or Script-S-2, K+1 to N for triangles with
%       edges on script-S-2, N+1 to M for remaining triangles.
%       Note: K=0 implies that no triangle is interior to D.
%       Note: M=N implies that all triangles have edges on
%       script-S-2.
%
%    2. Nodes numbered: 1 to LN for interior nodes and nodes on
%       script-S-2, LN+1 to LM for nodes on script-S-1.
%       Note: LM and LN represent lower case m and n resp.
%       Note: LN=LM implies that script-S-1 contains no nodes.
%       Note: If a node is on both script-S-1 and script-S-2, then
%       it should be treated as if it were only on script-S-1.
%
%    3. NL=number of line segments on script-S-2
%       line(I,J) is an NL by 2 array where
%       line(I,1)= first node on line I and
%       line(I,2)= second node on line I taken
%       in positive direction along line I
%
%    4. For the node labelled KK,KK=1,...,LM we have:
%       A) coordinates XX(KK),YY(KK)
%       B) number of triangles in which KK is a vertex= LL(KK)
%       C) II(KK,J) labels the triangles KK is in and
%       NV(KK,J) labels which vertex node KK is for
%       each J=1,...,LL(KK)
%
%    5. NTR is an M by 3 array where
%       NTR(I,1)=node number of vertex 1 in triangle I
%       NTR(I,2)=node number of vertex 2 in triangle I
%       NTR(I,3)=node number of vertex 3 in triangle I
%
%    6. Function subprograms:
%       A) P,Q,R,F,G,G1,G2 are the functions specified by
%          the particular differential equation
%       B) RR is the integrand
%          R*N(J)*N(K) on triangle I in step 4
%       C) FFF is the integrand F*N(J) on triangle I in step 4
%       D) GG1=G1*N(J)*N(K)
%          GG2=G2*N(J)
%          GG3=G2*N(K)
%          GG4=G1*N(J)*N(J)
%          GG5=G1*N(K)*N(K)
%          integrands in step 5
%       E) QQ(FF) computes the double integral of any
%          integrand FF over a triangle with vertices given by
%          nodes J1,J2,J3 - the method is an O(H**2) approximation
%          for triangles
%       F) SQ(PP) computes the line integral of any integrand PP
%          along the line from (XX(J1),YY(J1)) to (XX(J2),YY(J2))
%          by using a parameterization given by:
%          X=XX(J1)+(XX(J2)-XX(J1))*T
%          Y=YY(J1)+(YY(J2)-YY(J1))*T
%          for 0 <= t <= 1
%          and applying Simpson's composite method with H=.01
%
%    7. Arrays:
%       A) A,B,C are M by 3 arrays where the basis function N
%          for the Ith triangle, Jth vertex is
%          N(X,Y)=A(I,J)+B(I,J)*X+C(I,J)*Y
%          for J=1,2,3 and I=1,2,...,M
%       B) XX,YY are LM by 1 arrays to hold coordinates of nodes
%       C) line,LL,II,NV,NTR have been explained above
%       D) Gamma and Alpha are clear
%
%    8. Note that A,B,C,XX,YY,I,I1,I2,J1,J2,J3,Delta are reserved
%       storage so that in any subprogram we know that
%       triangle I has vertices (XX(J1),YY(J1)),(XX(J2),YY(J2)),
%       (XX(J3),YY(J3)). That is, vertex 1 is node J1, vertex 2 is
%       node J2, vertex 3 is node J3 unless a line integral is
%       involved in which case the line integral goes from node J1
%       to node J2 in triangle I or unless vertex I1 is node J1
%       and vertex I2 is node J2 - the basis functions involve
%       A(I,I1)+B(I,I1)*X+C(I,I1)*Y for vertex I1 in triangle I
%       and A(I,I2)+B(I,I2)*X+C(I,I2)*Y for vertex I2 in triangle I
%       delta is 1/2 the area of triangle I
%
%  To change problems:
%    1) change functions P,Q,R,F,G,G1,G2
%    2) change data input for K,N,M,LN,LM,NL.
%    3) change data input for XX,YY,LLL,II,NV.
%    4) change data input for line.
 syms('OK', 'AA', 'NAME', 'INP', 'K', 'N', 'M', 'LN1', 'LM', 'NL');
 syms('KK', 'XX', 'YY', 'LLL', 'J', 'II', 'NV', 'LL', 'N1', 'N2');
 syms('NTR', 'I', 'LINE', 'FLAG', 'OUP', 'K1', 'L', 'GAMMA');
 syms('ALPHA', 'J1', 'J2', 'J3', 'DELTA', 'A', 'B', 'C', 'I1');
 syms('JJ1', 'I2', 'JJ2', 'ZZ', 'HH', 'JI', 'XJ', 'XJ1', 'XJ2');
 syms('XI1', 'XI2', 'INN', 'CCC', 'x', 's', 'y');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Finite Element Method.\n');
 OK = FALSE;
 fprintf(1,'The input will be from a text file in the following form:\n');
 fprintf(1,'K  N  M  n  m  nl\n\n');
 fprintf(1,'Each of the above is an integer -');
 fprintf(1,'separate with at least one blank.\n\n');
 fprintf(1,'Follow with the input for each node in the form:\n');
 fprintf(1,'x-coord., y-coord., number of triangles in which the');
 fprintf(1,' node is a vertex.\n\n');
 fprintf(1,'Continue with the triangle number and vertex number for\n');
 fprintf(1,'each triangle in which the node is a vertex.\n');
 fprintf(1,'Separate each entry with at least one blank.\n\n');
 fprintf(1,'After all nodes have been entered follow with information\n');
 fprintf(1,'on the lines over which line integrals must be computed.\n');
 fprintf(1,'The format of this data will be the node number of the\n');
 fprintf(1,'starting node, followed by the node number of the ending\n');
 fprintf(1,'node for each line, taken in the positive direction.\n');
 fprintf(1,'There should be 2 * nl such entries, each an integer\n');
 fprintf(1,'separated by a blank.\n\n');
 fprintf(1,'Functions can be input or coded as procedures.\n');
 fprintf(1,'The example has all functions as procedures.\n');
 fprintf(1,'However, the functions could be input by uncommenting\n');
 fprintf(1,'the code as indicated.\n');
 fprintf(1,'Have the functions P,Q,R,F,G,G1,G2 been created and\n');
 fprintf(1,'has the input file been created?  Answer Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y'
% Remove "%" from the following statements if functions are to be input.
% fprintf(1,'Input function P(X,Y) in terms of x and y\n');
% s = input(' ');
%  P = inline(s,'x','y');
% fprintf(1,'Input function Q(X,Y) in terms of x and y\n');
% s = input(' ');
%  Q = inline(s,'x','y');
% fprintf(1,'Input function R(X,Y) in terms of x and y\n');
% s = input(' ');
%  R = inline(s,'x','y');
% fprintf(1,'Input function F(X,Y) in terms of x and y\n');
% s = input(' ');
%  F = inline(s,'x','y');
% fprintf(1,'Input function G(X,Y) in terms of x and y\n');
% s = input(' ');
% G = inline(s,'x','y');
% fprintf(1,'Input function G1(X,Y) in terms of x and y\n');
% s = input(' ');
% G1 = inline(s,'x','y');
% fprintf(1,'Input function G2(X,Y) in terms of x and y\n');
% s = input(' ');
% G2 = inline(s,'x','y');
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:  A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = TRUE;
 K = fscanf(INP, '%d',1);
 N = fscanf(INP, '%d',1);
 M = fscanf(INP, '%d',1);
 LN1 = fscanf(INP, '%d',1);
 LM = fscanf(INP, '%d',1);
 NL = fscanf(INP, '%d',1);
 XX = zeros(1,LM);
 YY = zeros(1,LM);
 LL = zeros(1,LM);
 II = zeros(LM,M);
 NV = zeros(LM,M);
 LINE = zeros(NL,2);
 NTR = zeros(M,3);
 A = zeros(M,3);
 B = zeros(M,3);
 C = zeros(M,3);
 GAMMA = zeros(1,LM);
 ALPHA = zeros(LN1,LN1+1);
 for KK = 1 : LM
 XX(KK) = fscanf(INP, '%f',1);
 YY(KK) = fscanf(INP, '%f',1);
 LLL = fscanf(INP, '%d',1);
 for J = 1 : LLL 
 II(KK,J) = fscanf(INP, '%d',1);
 NV(KK,J) = fscanf(INP, '%d',1);
 end;
 LL(KK) = LLL;
 for J = 1 : LLL 
 N1 = II(KK,J);
 N2 = NV(KK,J);
 NTR(N1,N2) = KK;
 end;
 end;
 if NL > 0 
 for I = 1 : NL 
 for J = 1 : 2 
 LINE(I,J) = fscanf(INP, '%d',1);
 end;
 end;
 end;
 fclose(INP);
 else
 fprintf(1,'The program will end so that the preparations\n');
 fprintf(1,'can be completed.\n');
 end;
 if OK == TRUE 
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:  A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'FINITE ELEMENT METHOD\n\n\n');
 K1 = K+1;
 N1 = LN1+1;
 fprintf(OUP, 'Vertices and Nodes of Triangles\n');
 fprintf(OUP, 'Triangle-node number for vertex 1 to 3\n');
 for I = 1 : M 
 fprintf(OUP, '%5d', I);
 for J = 1 : 3 
 fprintf(OUP, '%4d', NTR(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 fprintf(OUP, 'x and y coordinates of nodes\n');
 for KK = 1 : LM 
 fprintf(OUP, '%5d%11.8f%11.8f\n', KK, XX(KK), YY(KK));
 end;
 fprintf(OUP, 'Lines of the Domain\n');
 for KK = 1 : NL 
 fprintf(OUP, '%5d%4d%4d\n', KK, LINE(KK,1), LINE(KK,2));
 end;
% STEP 1
 if LM ~= LN1 
 for L = N1 : LM 
 GAMMA(L) = G(XX(L),YY(L));
 end;
 end;
% STEP 2 - initialization of ALPHA - note that ALPHA(I,LN1+1) = BETA(I)
 for I = 1 : LN1 
 for J = 1 : N1 
 ALPHA(I,J) = 0;
 end;
 end;
% STEPS 3, 4, and 6-12 are within the next loop
% for each triangle I let node J1 be vertex 1, node J2 be vertex 2
% and node J3 be vertex 3
% STEP 3
 for I = 1 : M 
 J1 = NTR(I,1);
 J2 = NTR(I,2);
 J3 = NTR(I,3);
 DELTA =  XX(J2)*YY(J3)-XX(J3)*YY(J2)-XX(J1)*(YY(J3)-YY(J2));
 DELTA =  DELTA + YY(J1)*(XX(J3)-XX(J2));
 A(I,1) = (XX(J2)*YY(J3)-YY(J2)*XX(J3))/DELTA;
 B(I,1) = (YY(J2)-YY(J3))/DELTA;
 C(I,1) = (XX(J3)-XX(J2))/DELTA;
 A(I,2) = (XX(J3)*YY(J1)-YY(J3)*XX(J1))/DELTA;
 B(I,2) = (YY(J3)-YY(J1))/DELTA;
 C(I,2) = (XX(J1)-XX(J3))/DELTA;
 A(I,3) = (XX(J1)*YY(J2)-YY(J1)*XX(J2))/DELTA;
 B(I,3) = (YY(J1)-YY(J2))/DELTA;
 C(I,3) = (XX(J2)-XX(J1))/DELTA;
% STEP 4
% I1 = J for  STEP 4 and I1 = K for  STEP 7
 for I1 = 1 : 3 
% STEP 8
 JJ1 = NTR(I,I1);
% I2 = K for  STEP 4 and I2 = J for  STEP 9
 for I2 = 1 : I1 
% STEP 4 and  STEP 10
 JJ2 = NTR(I,I2);
 ZZ = B(I,I1)*B(I,I2)*QQ(1,XX,YY,DELTA,J1,J2,J3,I1,I2,A,B,C);
 ZZ = ZZ + C(I,I1)*C(I,I2)*QQ(2,XX,YY,DELTA,J1,J2,J3,I1,I2,A,B,C);
 ZZ = ZZ - QQ(3,XX,YY,DELTA,J1,J2,J3,I1,I2,A,B,C);
% STEPS 11 and 12
 if JJ1 <= LN1 
 if JJ2 <= LN1 
 ALPHA(JJ1,JJ2) = ALPHA(JJ1,JJ2)+ZZ;
 if JJ1 ~= JJ2 
 ALPHA(JJ2,JJ1) = ALPHA(JJ2,JJ1)+ZZ;
 end;
 else
 ALPHA(JJ1,N1) = ALPHA(JJ1,N1)-ZZ*GAMMA(JJ2);
 end;
 else
 if JJ2 <= LN1 
 ALPHA(JJ2,N1) = ALPHA(JJ2,N1)-ZZ*GAMMA(JJ1);
 end;
 end;
 end;
 HH = -QQ(4,XX,YY,DELTA,J1,J2,J3,I1,I2,A,B,C);
 if JJ1 <= LN1 
 ALPHA(JJ1,N1) = ALPHA(JJ1,N1)+HH;
 end;
 end;
 end;
% output the basis functions
 fprintf(OUP, 'Basis Functions\n');
 fprintf(OUP, 'Triangle - Vertex - Node - Function\n');
 for I = 1 : M 
 for J = 1 : 3 
 fprintf(OUP,'%4d%3d%3d%13.8f%13.8f',I,J,NTR(I,J),A(I,J),B(I,J));
 fprintf(OUP,'%13.8f\n',C(I,J));
 end;
 end;
% STEP 5
% for each line segment JI = 1 ,..., NL and for each triangle
% I, I = K1 ,..., N which may contain line JI, search all 3 vertices
% for possible correspondences.
% STEP 5 and  STEPS 13-19
 if NL ~= 0 & N ~= K 
 for JI = 1 : NL 
 for I = K1 : N 
 for I1 = 1 : 3 
% I1 = J in  STEP 5 and I1 = K in  STEP 14
% STEP 15
 J1 = NTR(I,I1);
 if LINE(JI,1) == J1 
 for I2 = 1 : 3 
% I2 = K in  STEP 5 and I2 = J in  STEP 16
% STEP 17
 J2 = NTR(I,I2);
 if LINE(JI,2) == J2 
% There is a correspondence of vertex I1 in triangle I with
% node J1 as the start of line JI and vertex I2 with 
% node J2 as the end of line JI
% STEP 5
 XJ = SQ(1,XX,YY,J1,J2,I,I1,I2,A,B,C);
 XJ1 = SQ(4,XX,YY,J1,J2,I,I1,I2,A,B,C);
 XJ2 = SQ(5,XX,YY,J1,J2,I,I1,I2,A,B,C);
 XI1 = SQ(2,XX,YY,J1,J2,I,I1,I2,A,B,C);
 XI2 = SQ(3,XX,YY,J1,J2,I,I1,I2,A,B,C);
% STEPS 8 and 19
 if J1 <= LN1 
 if J2 <= LN1 
 ALPHA(J1,J1) = ALPHA(J1,J1)+XJ1;
 ALPHA(J1,J2) = ALPHA(J1,J2)+XJ;
 ALPHA(J2,J2) = ALPHA(J2,J2)+XJ2;
 ALPHA(J2,J1) = ALPHA(J2,J1)+XJ;
 ALPHA(J1,N1) = ALPHA(J1,N1)+XI1;
 ALPHA(J2,N1) = ALPHA(J2,N1)+XI2;
 else
 ALPHA(J1,N1) = ALPHA(J1,N1)-XJ*GAMMA(J2)+XI1;
 ALPHA(J1,J1) = ALPHA(J1,J1)+XJ1;
 end;
 else
 if J2 <= LN1 
 ALPHA(J2,N1) = ALPHA(J2,N1)-XJ*GAMMA(J1)+XI2;
 ALPHA(J2,J2) = ALPHA(J2,J2)+XJ2;
 end;
 end;
 end;
 end;
 end;
 end;
 end;
 end;
 end;
% output ALPHA
 fprintf(OUP, 'Matrix ALPHA follows:\n');
 for I = 1 : LN1 
 fprintf(OUP, 'Row%4d\n', I);
 for J = 1 : N1 
 fprintf(OUP, '%13.10e\n ', ALPHA(I,J));
 end;
 end;
 fprintf(OUP, '\n');
% STEP 20
 if LN1 > 1 
 INN = LN1-1;
 for I = 1 : INN 
 I1 = I+1;
 for J = I1 : LN1 
 CCC = ALPHA(J,I)/ALPHA(I,I);
 for J1 = I1 : N1 
 ALPHA(J,J1) = ALPHA(J,J1)-CCC*ALPHA(I,J1);
 end;
 ALPHA(J,I) = 0;
 end;
 end;
 end;
 GAMMA(LN1) = ALPHA(LN1,N1)/ALPHA(LN1,LN1);
 if LN1 > 1 
 for I = 1 : INN 
 J = LN1-I;
 CCC = ALPHA(J,N1);
 J1 = J+1;
 for KK = J1 : LN1 
 CCC = CCC-ALPHA(J,KK)*GAMMA(KK);
 end;
 GAMMA(J) = CCC/ALPHA(J,J);
 end;
 end;
% STEP 21
% output gamma
 fprintf(OUP, 'Coefficients of Basis Functions:\n');
 for I = 1 : LM 
 LLL = LL(I);
 fprintf(OUP, '%3d %14.8f %d', I, GAMMA(I), I);
 for J = 1 : LLL 
 fprintf(OUP, '%4d', II(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
% STEP 23
 end;
