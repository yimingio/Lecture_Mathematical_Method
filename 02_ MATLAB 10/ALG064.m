% DIRECT FACTORIZATION ALGORITHM 6.4
%
% To factor the n by n matrix A = (A(I,J)) into the product of the
% lower triangular matrix L = (L(I,J)) and the upper triangular
% matrix U = (U(I,J)), that is A = LU, where the main diagonal of 
% either L or U consists of all ones:
%
% INPUT:   dimension n; the entries A(I,J), 1<=I, J<=n, of A;
%          the diagonal L(1,1), ..., L(N,N) of L or the diagonal
%          U(1,1), ..., U(N,N) of U.
%
% OUTPUT:  the entries L(I,J), 1<=J<=I, 1<=I<=n of L and the entries
%          U(I,J), I<=J<=n, 1<=I<=n of U.
 syms('AA', 'NAME', 'INP', 'OK', 'N', 'I', 'J', 'A');
 syms('FLAG', 'ISW', 'XL', 'M', 'KK', 'S', 'K', 'JJ');
 syms('SS', 'OUP', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the general LU factorization method.\n');
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,N), \n') 
 fprintf(1,'A(2,1), A(2,2), ..., A(2,N),\n');
 fprintf(1,'..., A(N,1), A(N,2), ..., A(N,N)\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate\n');
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
 fprintf(1,'Input the dimension n - an integer.\n');
 N = input(' ');
 if N > 0 
 A = zeros(N,N);
 XL = zeros(1,N);
 for I = 1 : N 
 for J = 1 : N 
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 OK = TRUE;
 fclose(INP);
 else fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 fprintf(1,'Choice of diagonals:\n');
 fprintf(1,'1. Diagonal of L consists of ones\n'); 
 fprintf(1,'2. Diagonal of U consists of ones\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 1 
 ISW = 0;
 else
 ISW = 1;
 end
 else 
 fprintf(1,'The program will end so the input file can be created.\n');
 OK = FALSE;
 end;
 if OK == TRUE 
 for I = 1 : N 
 XL(I) = 1;
 end;
% STEP 1
 if abs(A(1,1)) <= 1.0e-20 
 OK = FALSE;
 else
% the entries of L below the main diagonal will be placed 
% in the corresponding entries of A; the entries of U 
% above the main diagonal will be placed in the 
% corresponding entries of A; the main diagonal which 
% was not input will become the main diagonal of A; 
% the input main diagonal of L or U is, 
% of course, placed in XL
 A(1,1) = A(1,1)/XL(1);
% STEP 2
 for J = 2 : N 
 if ISW == 0 
% first row of U
 A(1,J) = A(1,J)/XL(1);
% first column of L
 A(J,1) = A(J,1)/A(1,1);
 else
% first row of U
 A(1,J) = A(1,J)/A(1,1);
% first column of L
 A(J,1) = A(J,1)/XL(1);
 end;
 end;
% STEP 3
 M = N-1;
 I = 2;
 while I <= M & OK == TRUE 
% STEP 4
 KK = I-1;
 S = 0;
 for K = 1 : KK 
 S = S-A(I,K)*A(K,I);
 end;
 A(I,I) = (A(I,I)+S)/XL(I);
 if abs(A(I,I)) <= 1.0e-20 
 OK = FALSE;
 else
% STEP 5
 JJ = I+1;
 for J = JJ : N 
 SS = 0;
 S = 0;
 for K = 1 : KK 
 SS = SS-A(I,K)*A(K,J);
 S = S-A(J,K)*A(K,I);
 end;
 if ISW == 0 
% Ith row of U
 A(I,J) = (A(I,J)+SS)/XL(I);
% Ith column of L
 A(J,I) = (A(J,I)+S)/A(I,I);
 else
% Ith row of U
 A(I,J) = (A(I,J)+SS)/A(I,I);
% Ith column of L
 A(J,I) = (A(J,I)+S)/XL(I);
 end;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
% STEP 6
 S = 0;
 for K = 1 : M 
 S = S-A(N,K)*A(K,N);
 end;
 A(N,N) = (A(N,N)+S)/XL(N);
% If A(N,N) = 0 then A = LU but the matrix is singular.
% Process is complete, all entries of A have been determined.
% STEP 7
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'For example   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'GENERAL LU FACTORIZATION\n\n');
 if ISW == 0  
 fprintf(OUP, 'The diagonal of L consists of all entries = 1.0\n');
 else
 fprintf(OUP, 'The diagonal of U consists of all entries = 1.0\n');
 end;
 fprintf(OUP, '\nEntries of L below/on diagonal and entries of U above');
 fprintf(OUP, '/on diagonal\n');
 fprintf(OUP, '- output by rows in overwrite format:\n');
 for I = 1 : N 
 for J = 1 : N 
 fprintf(OUP, ' %11.8f', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
 end;
 if OK == FALSE 
 fprintf(1,'The matrix does not have an LU factorization.\n');
 end;
 end;
