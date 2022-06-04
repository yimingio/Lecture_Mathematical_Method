% HOUSEHOLDER'S ALGORITHM 9.5
%
% To obtain a symmetric tridiagonal matrix A(n-1) similar
% to the symmetric matrix A = A(1), construct the following
% matrices A(2),A(3),...,A(n-1) where A(K) = A(I,J)**K, for
% each K = 1,2,...,n-1:
%
% INPUT:   Dimension n; matrix A.
%
% OUTPUT:  A(n-1) (At each step, A can be overwritten.)
 syms('OK', 'AA', 'NAME', 'INP', 'N', 'I', 'J', 'A', 'K');
 syms('Q', 'KK', 'S', 'RSQ', 'V', 'U', 'PROD', 'Z');
 syms('L', 'FLAG', 'OUP');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Householder Method.\n');
 OK = FALSE;
 fprintf(1,'The symmetric array A will be input from a text file\n');
 fprintf(1,'in the order:\n');
 fprintf(1,'              A(1,1), A(1,2), A(1,3), ..., A(1,n),\n');
 fprintf(1,'                      A(2,2), A(2,3), ..., A(2,n),\n');
 fprintf(1,'                              A(3,3), ..., A(3,n),\n');
 fprintf(1,'                                      ..., A(n,n)\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - drive:name.ext\n');
 fprintf(1,'for example:    A:DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the dimension n.\n');
 N = input(' ');
 if N > 1
 A = zeros(N,N);
 U = zeros(1,N);
 V = zeros(1,N);
 Z = zeros(1,N);
 for I = 1 : N 
 for J = I : N 
 A(I,J) = fscanf(INP, '%f',1);
 A(J,I) = A(I,J);
 end;
 end;
 fclose(INP);
 OK = TRUE;
 else
 fprintf(1,'Dimension must be greater than 1.\n');
 end;
 end;
 else
 fprintf(1,'The program will end so the input file can be created.\n');
 end;
 if OK == TRUE 
% STEP 1
 for K = 1 : N-2 
 Q = 0;
 KK = K+1;
% STEP 2
 for I = KK : N 
 Q = Q+A(I,K)*A(I,K);
 end;
% STEP 3
% S is used in place of alpha.
 if abs(A(K+1,K)) <= 1.0e-20 
 S = sqrt(Q);
 else
 S = A(K+1,K)/abs(A(K+1,K))*sqrt(Q);
 end;
% STEP 4
 RSQ = (S+A(K+1,K))*S;
% STEP 5
 V(K) = 0;
 V(K+1) = A(K+1,K)+S;
 for J = K+2 : N 
 V(J) = A(J,K);
 end;
% STEP 6
 for J = K : N 
 U(J) = 0;
 for I = KK : N 
 U(J) = U(J)+A(J,I)*V(I);
 end;
 U(J) = U(J)/RSQ;
 end;
% STEP 7
 PROD = 0;
 for I = K+1 : N 
 PROD = PROD + V(I)*U(I);
 end;
% STEP 8
 for J = K : N 
 Z(J) = U(J) - 0.5*PROD*V(J)/RSQ;
 end;
% STEP 9
 for L = K+1 : N-1 
% STEP 10
 for J = L+1 : N 
 A(J,L) = A(J,L)-V(L)*Z(J)-V(J)*Z(L);
 A(L,J) = A(J,L);
 end;
% STEP 11
 A(L,L) = A(L,L) - 2*V(L)*Z(L);
 end;
% STEP 12
 A(N,N) = A(N,N)-2*V(N)*Z(N);
% STEP 13
 for J = K+2 : N 
 A(K,J) = 0;
 A(J,K) = 0;
 end;
% STEP 14
 A(K+1,K) = A(K+1,K)-V(K+1)*Z(K);
 A(K,K+1) = A(K+1,K);
 end;
% STEP 15
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:name.ext\n');
 fprintf(1,'for example    A:OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'HOUSEHOLDER METHOD\n\n');
 fprintf(OUP, 'The similar tridiagonal matrix follows - output by rows\n\n');
 for I = 1 : N 
 for J = 1 : N 
 fprintf(OUP, ' %11.8f', A(I,J));
 end;
 fprintf(OUP, '\n\n');
 end;
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
