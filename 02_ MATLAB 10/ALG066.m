% CHOLESKI'S ALGORITHM 6.6
%
% To factor the positive definite n by n matrix A into LL**T,
% where L is lower triangular.
%
% INPUT:   the dimension n; entries A(I,J), 1<=I, J<=n of A.
%
% OUTPUT:  the entries L(I,J), 1<=J<=I, 1<=I<=n of L.
%
 syms('AA', 'NAME', 'INP', 'OK', 'N', 'I', 'J', 'A', 'NN');
 syms('KK', 'S', 'K', 'JJ', 'FLAG', 'OUP', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Choleski Factorization Method.\n');
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,N), \n');
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
 else 
 fprintf(1,'The program will end so the input file can be created.\n');
 OK = FALSE;
 end;
 if OK == TRUE 
% STEP 1
 A(1,1) = sqrt(A(1,1));
% STEP 2
 for J = 2 : N 
 A(J,1) = A(J,1)/A(1,1);
 end;
% STEP 3
 NN = N-1;
 for I = 2 : NN 
% STEP 4
 KK = I-1;
 S = 0;
 for K = 1 : KK 
 S = S-A(I,K)*A(I,K);
 end;
 A(I,I) = sqrt(A(I,I)+S);
% STEP 5
 JJ = I+1;
 for J = JJ : N 
 S = 0;
 KK = I-1;
 for K = 1 : KK 
 S = S - A(J,K)*A(I,K);
 end;
 A(J,I) = (A(J,I)+S)/A(I,I);
 end;
 end;
% STEP 6
 S = 0;
 for K = 1 : NN 
 S = S-A(N,K)*A(N,K);
 end;
 A(N,N) = sqrt(A(N,N)+S);
% STEP 7
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
 fprintf(OUP, 'CHOLESKI FACTORIZATION\n\n');
 fprintf(OUP, 'The matrix L output by rows:\n');
 for I = 1 : N 
 for J = 1 : I 
 fprintf(OUP, '  %12.8f', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
