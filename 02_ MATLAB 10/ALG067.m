% CROUT FACTORIZATION FOR TRIDIAGONAL LINEAR SYSTEMS ALGORITHM 6.7
%
% To solve the n x n linear system
%
% E1:  A(1,1) X(1) + A(1,2) X(2)                  = A(1,n+1)
% E2:  A(2,1) X(1) + A(2,2) X(2) + A(2,3) X(3)    = A(2,n+1)
% :
% .
% E(n):          A(n,n-1) X(n-1) + A(n,n) X(n)    = A(n,n+1)
%
% INPUT:   the dimension n; the entries of A.
%
% OUTPUT:  the solution X(1), ..., X(N).
 syms('AA', 'OK', 'NAME', 'INP', 'N', 'I', 'A', 'B', 'NN');
 syms('C', 'BB', 'Z', 'X', 'II', 'FLAG', 'OUP', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Crout Method for tridiagonal linear systems.\n');
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'all diagonal entries, all lower sub-diagonal entries, all ');
 fprintf(1,'upper sub-diagonal\n');
 fprintf(1,'entries, inhomogeneous term.\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 OK = FALSE;
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
 A = zeros(1,N);
 B = zeros(1,N);
 C = zeros(1,N);
 BB = zeros(1,N);
 X = zeros(1,N);
 Z = zeros(1,N);
% A(I,I) is stored in A(I), 1 <= I <= n */
 for I = 1 : N 
 A(I) = fscanf(INP, '%f',1);
 end;
% the lower sub-diagonal A(I,I-1) is stored
%in B(I), 2 <= I <= n */
 for I = 2 : N 
 B(I) = fscanf(INP, '%f',1);
 end;
% the upper sub-diagonal A(I,I+1) is stored
%in C(I), 1 <= I <= n-1 */
 NN = N-1;
 for I = 1 : NN 
 C(I) = fscanf(INP, '%f',1);
 end;
% A(I,N+1) is stored in BB(I), 1 <= I <= n */
 for I = 1 : N 
 BB(I) = fscanf(INP, '%f',1);
 end;
 OK = TRUE;
 fclose(INP);
 else
 fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 else
 fprintf(1,'The program will end so the input file can be created.\n')
 end;
 if OK == TRUE 
% Steps 1-3 set up and solve LZ = B
% STEP 1
% the entries of U overwrite C and the entries of  L overwrite A
 C(1) = C(1)/A(1);
 Z(1) = BB(1)/A(1);
% STEP 2
 for I = 2 : NN 
 A(I) = A(I)-B(I)*C(I-1);
 C(I) = C(I)/A(I);
 Z(I) = (BB(I)-B(I)*Z(I-1))/A(I);
 end;
% STEP 3
 A(N) = A(N)-B(N)*C(N-1);
 Z(N) = (BB(N)-B(N)*Z(N-1))/A(N);
% STEP 4
% STEPS 4, 5 solve UX = Z
 X(N) = Z(N);
% STEP 5
 for II = 1 : NN 
 I = NN-II+1;
 X(I) = Z(I)-C(I)*X(I+1);
 end;
% STEP 6
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'CROUT METHOD FOR TRIDIAGONAL LINEAR SYSTEMS\n\n');
 fprintf(OUP, 'The solution is\n');
 for I = 1 : N 
 fprintf(OUP, '  %12.8f', X(I));
 end;
 fprintf(OUP, '\n');
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
