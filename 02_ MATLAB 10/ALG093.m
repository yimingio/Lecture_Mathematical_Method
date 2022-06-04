% INVERSE POWER METHOD ALGORITHM 9.3
%
% To approximate an eigenvalue and an associated eigenvector of the
% n by n matrix A given a nonzero vector x:
%
% INPUT:   Dimension n; matrix A; vector x; tolerance TOL;
%          maximum number of iterations N.
%
% OUTPUT:  Approximate eigenvalue MU; approximate eigenvector x
%          or a message that the maximum number of iterations was
%          exceeded.
 syms('OK', 'AA', 'NAME', 'INP', 'N', 'I');
 syms('J', 'A', 'X', 'TOL', 'NN', 'FLAG', 'OUP', 'Q');
 syms('S', 'K', 'LP', 'AMAX', 'B', 'YMU', 'ERR', 'T');
 syms('M','IMAX','IP','L1','L2','JJ','I1','J1');
 syms('N1','L','N2','KK');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Inverse Power Method.\n');
 OK = FALSE;
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,n), \n');
 fprintf(1,'A(2,1), A(2,2), ..., A(2,n),\n');
 fprintf(1,'..., A(n,1), A(n,2), ..., A(n,n) \n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n');
 fprintf(1,'The initial approximation should follow in same format.\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y'
 fprintf(1,'Input the file name in the form - drive:name.ext \n');
 fprintf(1,'for example:   A:DATA.DTA \n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the dimension n.\n');
 N = input(' ');
 if N > 0
 A = zeros(N,N);
 X = zeros(1,N);
 Y = zeros(1,N);
 B = zeros(1,N);
 NROW = zeros(1,N);
 for I = 1 : N
 for J = 1 : N
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 for I = 1 : N
 X(I) = fscanf(INP, '%f',1);
 end;
 fclose(INP);
 while OK == FALSE
 fprintf(1,'Input the tolerance.\n');
 TOL = input(' ');
 if TOL > 0
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be positive number.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input maximum number of iterations ');
 fprintf(1,'- integer.\n');
 NN = input(' ');
% use NN for n
 if NN > 0
 OK = TRUE;
 else
 fprintf(1,'Number must be positive integer.\n');
 end;
 end;
 else
 fprintf(1,'The dimension must be a positive integer.\n');
 end;
 end;
 else
 fprintf(1,'The program will end so the input file can be created.\n');
 end;
 if OK == TRUE
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2
 fprintf(1,'Input the file name in the form - drive:name.ext\n');
 fprintf(1,'for example   A:OUTPUT.DTA \n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'INVERSE POWER METHOD\n');
% STEP 1
% Q could be input instead of computed.
 Q = 0;
 S = 0;
 for I = 1 : N
 S = S + X(I) * X(I);
 for J = 1 : N
 Q = Q + A(I,J) * X(I) * X(J);
 end;
 end;
 Q = Q / S;
 fprintf(1,'Q is %.8e\n', Q);
 fprintf(1,'Input new Q? Enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' |  AA == 'y'
 fprintf(1,'input new Q\n');
 Q = input(' ');
 end;
 fprintf(OUP, 'Iteration  Eigenvalue  Eigenvector\n');
% STEP 2
 K3 = 1;
 for I = 1 : N
 A(I,I) = A(I,I) - Q;
 end;
% Determine the row ordering and multipliers
% for the matrix A-Q*I for use in Gaussian elimination
% with Partial Pivoting.
% NROW holds the ordering of the rows for interchanges.
 for I = 1 : N
 NROW(I) = I;
 end;
 OK = TRUE;
 I = 1;
 M = N - 1;
 while I <= M & OK == TRUE
 IMAX = I;
 J = I+1;
 for IP = J : N
 L1 = NROW(IMAX);
 L2 = NROW(IP);
 if abs(A(L2,I)) > abs(A(L1,I))
 IMAX = IP;
 end;
 end;
 if abs(A(NROW(IMAX),I)) <= 1.0e-20
 OK = FALSE;
 fprintf(1,'A - Q * I is singular, Q = %.8e is an eigenvalue\n', Q);
 else
 JJ = NROW(I);
 NROW(I) = NROW(IMAX);
 NROW(IMAX) = JJ;
 I1 = NROW(I);
 for JJ = J : N
 J1 = NROW(JJ);
 A(J1,I) = A(J1,I) / A(I1,I);
 for K = J : N
 A(J1,K) = A(J1,K) - A(J1,I) * A(I1,K);
 end;
 end;
 end;
 I = I+1;
 end;
 if abs(A(NROW(N),N)) <= 1.0e-20
 OK = FALSE;
 fprintf(1,'A - Q * I is singular, Q = %.8e is an eigenvalue\n', Q);
 end;
 if OK == TRUE
% STEP 3
 LP = 1;
 for I = 2 : N
 if abs(X(I)) > abs(X(LP))
 LP = I;
 end;
 end;
% STEP 4
 AMAX = X(LP);
 for I = 1 : N
 X(I) = X(I) / (AMAX);
 end;
% STEP 5
 while K3 <= NN & OK == TRUE
% STEPS 6 AND 7
 for I = 1 : N
 B(I) = X(I);
 end;
% Solve the linear system (A-Q*I) Y = B given a new
% vector X and the row ordering and multipliers from procedure MULTIP.
 M = N - 1;
 for I = 1 : M
 J = I+1;
 I1 = NROW(I);
 for JJ = J : N
 J1 = NROW(JJ);
 B(J1) = B(J1) - A(J1,I) * B(I1);
 end;
 end;
 N1 = NROW(N);
 Y(N) = B(N1) / A(N1,N);
 L = N - 1;
 for K1 = 1 : L
 J = L - K1 + 1;
 JJ = J + 1;
 N2 = NROW(J);
 Y(J) = B(N2);
 for KK = JJ : N
 Y(J) = Y(J) - A(N2,KK) * Y(KK);
 end;
 Y(J) = Y(J) / A(N2,J);
 end;
% STEP 8
 YMU = Y(LP);
% STEPS 9 AND 10
 LP = 1;
 for I = 2 : N
 if abs(Y(I)) > abs(Y(LP))
 LP = I;
 end;
 end;
 AMAX = Y(LP);
 ERR = 0;
 for I = 1 : N
 T = Y(I) / AMAX;
 if abs(X(I) - T) > ERR
 ERR = abs(X(I) - T);
 end;
 X(I) = T;
 end;
 YMU = 1 / YMU + Q;
% STEP 11
 fprintf(OUP, '%3d %12.8f\n', K3, YMU);
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
 if ERR < TOL
 OK = FALSE;
 fprintf(OUP, 'Eigenvalue = %12.8f', YMU);
 fprintf(OUP, ' to tolerance = %.10e\n', TOL);
 fprintf(OUP, 'obtained on iteration number %d\n\n', K3);
 fprintf(OUP, 'Unit eigenvector is :\n');
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
 else
% STEP 12
 K3 = K3+1;
 end;
 end;
 if K3 > NN
 fprintf(OUP,  'No convergence in %d iterations\n',NN);
 end;
 end;
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
