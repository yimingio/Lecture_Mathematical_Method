% GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING ALGORITHM 6.2
%
% To solve the n by n linear system
%
% E1:  A(1,1) X(1) + A(1,2) X(2) +...+ A(1,n) X(n) = A(1,n+1)
% E2:  A(2,1) X(1) + A(2,2) X(2) +...+ A(2,n) X(n) = A(2,n+1)
% :
% .
% EN:  A(n,1) X(1) + A(n,2) X(2) +...+ A(n,n) X(n) = A(n,n+1)
%
% INPUT:   number of unknowns and equations n; augmented
%          matrix A = (A(I,J)) where 1<=I<=n and 1<=J<=n+1.
%
% OUTPUT:  solution x(1), x(2),...,x(n) or a message that the
%         linear system has no unique solution.
 syms('AA', 'NAME', 'INP', 'OK', 'N', 'I', 'J', 'A');
 syms('M', 'NROW', 'NN', 'ICHG', 'IMAX', 'AMAX', 'JJ');
 syms('IP', 'JP', 'NCOPY', 'I1', 'J1', 'XM', 'K', 'N1');
 syms('X', 'N2', 'SUM', 'KK', 'FLAG', 'OUP');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Gaussian Elimination with Partial Pivoting.\n');
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,N+1) \n');
 fprintf(1,'A(2,1), A(2,2), ..., A(2,N+1),\n');
 fprintf(1,'..., A(N,1), A(N,2), ..., A(N,N+1)\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the number of equations - an integer.\n');
 N = input(' ');
 if N > 0 
 A = zeros(N,N+1);
 X = zeros(1,N);
 NROW = zeros(1,N);
 for I = 1:N
 for J = 1:N+1 
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 OK = TRUE;
 fclose(INP);
 else fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 else 
 fprintf(1,'The program will end so the input file can be created.\n');
 end;
 if OK == TRUE 
 M = N+1;
% STEP 1
 for I = 1:N
 NROW(I) = I;
 end;
% initialize row pointer
 NN = N-1;
 ICHG = 0;
 I = 1;
% STEP 2
 while OK == TRUE & I <= NN 
% STEP 3
 IMAX = NROW(I);
 AMAX = abs(A(IMAX,I));
 IMAX = I;
 JJ = I+1;
 for IP = JJ:N
 JP = NROW(IP);
 if abs(A(JP,I)) > AMAX 
 AMAX = abs(A(JP,I));
 IMAX = IP;
 end;
 end;
% STEP 4
 if AMAX <= 1.0e-20 
 OK = FALSE;
 else
% STEP 5
% simulate row interchange
 if NROW(I) ~= NROW(IMAX) 
 ICHG = ICHG+1;
 NCOPY = NROW(I);
 NROW(I) = NROW(IMAX);
 NROW(IMAX) = NCOPY;
 end;
 I1 = NROW(I);
% STEP 6
 for J = JJ:N
 J1 = NROW(J);
% STEP 7
 XM = A(J1,I)/A(I1,I);
% STEP 8
 for K = JJ:M
 A(J1,K) = A(J1,K)-XM*A(I1,K);
 end;
% Multiplier XM could be saved in A(J1,I)
 A(J1,I) = 0;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
% STEP 9
 N1 = NROW(N);
 if abs(A(N1,N)) <= 1.0e-20 
 OK = FALSE;
% system has no unique solution
 else
% STEP 10
% start backward substitution
 X(N) = A(N1,M) / A(N1,N);
% STEP 11
 for K = 1:NN
 I = NN - K + 1;
 JJ = I + 1;
 N2 = NROW(I);
 SUM = 0;
 for KK = JJ:N
 SUM = SUM-A(N2,KK)*X(KK);
 end;
 X(I) = (A(N2,M) + SUM) / A(N2,I);
 end;
% STEP 12
% procedure completed successfully
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
 fprintf(OUP, 'GAUSSIAN ELIMINATION - PARTIAL PIVOTING\n\n');
 fprintf(OUP, 'The reduced system - output by rows:\n');
 for I = 1:N
 for J = 1:M
 fprintf(OUP, ' %11.8f', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 fprintf(OUP, '\n\nHas solution vector:\n');
 for I = 1:N
 fprintf(OUP, '  %12.8f', X(I));
 end;
 fprintf (OUP, '\n\nwith %d row interchange(s)\n', ICHG);
 fprintf(OUP, '\nThe rows have been logically re-ordered to:\n');
 for I = 1:N
 fprintf(OUP, ' %2d', NROW(I)); 
 end;
 fprintf(OUP,'\n');
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
 end;
 if OK == FALSE 
 fprintf(1,'System has no unique solution\n');
 end;
 end;
