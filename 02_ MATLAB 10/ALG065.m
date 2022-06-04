% LDL^t ALGORITHM 6.5
%
% To factor the positive definite n by n matrix A into LDL**T,
% where L is a lower triangular matrix with ones along the diagonal
% and D is a diagonal matrix with positive entries on the
% diagonal.
%
% INPUT:   the dimension n; entries A(I,J), 1<=I, J<=n of A.
%
% OUTPUT:  the entries L(I,J), 1<=J<I, 1<=I<=N of L and D(I),
%          1<=I<=n of D.
 syms('AA', 'NAME', 'INP', 'OK', 'N', 'I', 'J', 'A');
 syms('V', 'D', 'K', 'FLAG', 'OUP', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the LDL^t Method for Positive Definite Matrices.\n');
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
 V = zeros(1,N);
 D = zeros(1,N);
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
 for I = 1 : N 
% STEP 2
 for J = 1 : I - 1 
 V(J) = A(I,J)*D(J);
 end;
% STEP 3
 D(I) = A(I,I);
 for J = 1 : I-1 
 D(I) = D(I)-A(I,J)*V(J);
 end;
% STEP 4
 for J = I+1 : N 
 for K = 1 : I-1 
 A(J,I) = A(J,I)-A(J,K)*V(K);
 end;
 A(J,I) = A(J,I)/D(I);
 end;
 end;
% STEP 5
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
 fprintf(OUP, 'LDL^t FACTORIZATION\n\n');
 fprintf(OUP, 'The matrix L output by rows:\n');
 for I = 1 : N 
 for J = 1 : I - 1
 fprintf(OUP, '  %12.8f', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 fprintf(OUP, 'The diagonal of D:\n');
 for I = 1 : N 
 fprintf(OUP, '  %12.8f', D(I));
 end;
 fprintf(OUP, '\n');
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
