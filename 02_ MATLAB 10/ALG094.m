% WIELANDT'S DEFLATION ALGORITHM 9.4
%
% To approximate the second most dominant eigenvalue and an
% associated eigenvector of the n by n matrix A given an
% approximation LAMBDA to the dominant eigenvalue, an
% approximation V to a corresponding eigenvector and a vector X
% belonging to R^(n-1), tolerance TOL, maximum number of
% iterations N.
%
% INPUT:   Dimension n; matrix A; approximate eigenvalue LAMBDA;
%          approximate eigenvector V belonging to R^n; vector X
%          belonging to R^(n-1).
%
% OUTPUT:  Approximate eigenvalue MU; approximate eigenvector U or
%          a message that the method fails.
 syms('OK', 'AA', 'NAME', 'INP', 'N', 'TOL', 'NN', 'I', 'J');
 syms('A', 'V', 'XMU', 'YMU', 'M', 'X', 'FLAG', 'OUP');
 syms('AMAX', 'K', 'B', 'W', 'S', 'VV', 'L1', 'L2', 'Y');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Wielandt Deflation.\n');
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,n), \n');
 fprintf(1,'A(2,1), A(2,2), ..., A(2,n),\n');
 fprintf(1,'..., A(n,1), A(n,2), ..., A(n,n)\n\n');
 fprintf(1,'Next place the approximate eigenvector V(1), ..., ');
 fprintf(1,'V(n) and follow it\n');
 fprintf(1,'by the approximate eigenvalue. Finally, an ');
 fprintf(1,'initial approximate\n');
 fprintf(1,'eigenvector of dimension n-1: X(1), ..., X(n-1) ');
 fprintf(1,'should follow.\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y'
 fprintf(1,'Input the file name in the form - drive:name.ext\n');
 fprintf(1,'for example:   A:DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the dimension n.\n');
 N = input(' ');
 if N > 1
 A = zeros(N,N);
 B = zeros(N,N);
 X = zeros(1,N);
 W = zeros(1,N);
 V = zeros(1,N);
 Y = zeros(1,N);
 VV = zeros(1,N);
 OK = TRUE;
 else
 fprintf(1,'Dimension must be greater than 1.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input a positive tolerance for the power method.\n');
 TOL = input(' ');
 if TOL > 0
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be a positive number.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the maximum number of iterations for the ');
 fprintf(1,'power method.\n');
 NN = input(' ');
 if NN > 0
 OK = TRUE;
 else
 fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 for I = 1 : N
 for J = 1 : N
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 OK = FALSE;
 for I = 1 : N
 V(I) = fscanf(INP, '%f',1);
 if abs(V(I)) > 0
 OK = TRUE;
 end;
 end;
 XMU = fscanf(INP, '%f',1);
 M = N-1;
 if OK == TRUE
 OK = FALSE;
 for I = 1 : M
 X(I) = fscanf(INP, '%f',1);
 if abs(X(I)) > 0
 OK = TRUE;
 end;
 end;
 end;
 if OK == FALSE
 fprintf(1,'Input Error - All vectors must be nonzero.\n');
 end;
 fclose(INP);
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
 fprintf(1,'for example   A:OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'WIELANDT DEFLATION\n\n');
% STEP 1
 I = 1;
 AMAX = abs(V(1));
 for J = 2 : N
 if abs(V(J)) > AMAX
 I = J;
 AMAX = abs(V(J));
 end;
 end;
% STEP 2
 if I ~= 1
 for K = 1 : I-1
 for J = 1 : I-1
 B(K,J) = A(K,J)-V(K)*A(I,J)/V(I);
 end;
 end;
 end;
% STEP 3
 if I ~= 1 & I ~= N
 for K = I : N-1
 for J = 1 : I-1
 B(K,J) = A(K+1,J)-V(K+1)*A(I,J)/V(I);
 B(J,K) = A(J,K+1)-V(J)*A(I,K+1)/V(I);
 end;
 end;
 end;
% STEP 4
 if I ~= N
 for K = I : N-1
 for J = I : N-1
 B(K,J) = A(K+1,J+1)-V(K+1)*A(I,J+1)/V(I);
 end;
 end;
 end;
 I3 = I;
 K = 1;
 LP = 1;
 AMAX = abs(X(1));
 for I = 2 : M
 if abs(X(I)) > AMAX
 AMAX = abs(X(I));
 LP = I;
 end;
 end;
 DONE = FALSE;
 for I = 1 : M
 X(I) = X(I) / AMAX;
 end;
 while K <= NN & OK == TRUE & DONE == FALSE
 for I = 1 : M
 Y(I) = 0;
 for J = 1 : M
 Y(I) = Y(I) + B(I,J) * X(J);
 end;
 end;
 YMU = Y(LP);
 LP = 1;
 AMAX = abs(Y(1));
 for I = 2 : M
 if abs(Y(I)) > AMAX
 AMAX = abs(Y(I));
 LP = I;
 end;
 end;
 if AMAX <= 1.0e-20
 fprintf(1,'Zero eigenvalue - B is singular\n');
 OK = FALSE;
 else
 ERR = 0;
 for I = 1 : M
 T = Y(I)/Y(LP);
 if abs(X(I)-T) > ERR
 ERR = abs(X(I)-T);
 end;
 X(I) = T;
 end;
 if ERR < TOL
 for I = 1 : M
 Y(I) = X(I);
 end;
 DONE = TRUE;
 else
 K = K+1;
 end;
 end;
 end;
 if K > NN & OK == TRUE
 fprintf(1,'Power Method did not converge in %d iterations.\n',NN);
 OK = FALSE;
 else
 fprintf(OUP, 'Number Iterations for Power Method = %d\n', K);
 end;
 I = I3;
 if OK == TRUE
% STEP 6
 if I ~= 1
 for K = 1 : I-1
 W(K) = Y(K);
 end;
 end;
% STEP 7
 W(I) = 0;
% STEP 8
 if I ~=  N
 for K = I+1 : N
 W(K) = Y(K - 1);
 end;
 end;
% STEP 9
 S = 0;
 for J = 1 : N
 S = S + A(I,J) * W(J);
 end;
 S = S/V(I);
 for K = 1 : N
% Compute eigenvector
% VV is used in place of u here
 VV(K) = (YMU-XMU)*W(K)+S*V(K);
 end;
 fprintf(OUP, 'The reduced matrix B:\n');
 for L1 = 1 : M
 for L2 = 1 : M
 fprintf(OUP, '%.10e  ', B(L1,L2));
 end;
 fprintf(OUP, '\n');
 end;
 fprintf(OUP, '\nThe Eigenvalue = %12.8f', YMU);
 fprintf(OUP, ' to Tolerance = %.10e\n\n', TOL);
 fprintf(OUP, 'Eigenvector is:\n');
 for I = 1 : N
 fprintf(OUP,' %11.8f', VV(I));
 end;
 fprintf(OUP, '\n');
 end;
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
