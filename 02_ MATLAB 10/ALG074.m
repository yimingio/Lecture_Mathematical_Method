% ITERATIVE REFINEMENT ALGORITHM 7.4
%
% To approximate the solution to the linear system Ax=b when A is
% suspected to be ill-conditioned:
%
% INPUT:  The number of equations and unknowns n; the entries
%         A(i,j), 1<=i, j<=n, of the matrix A; the entries b(i),
%         1<=i<=n, of the inhomogeneous term b; the maximum number
%         of iterations N.
%
% OUTPUT: The approximation XX(1),...,XX(n) or a message that the
%         number of iterations was exceeded. 
 syms('A1', 'OK', 'NAME', 'INP', 'N', 'I', 'J', 'A', 'NN');
 syms('RND', 'D', 'TOL', 'FLAG', 'OUP', 'M', 'NROW', 'B', 'KK');
 syms('IS', 'C', 'L', 'X', 'S', 'K', 'XX', 'LL', 'R', 'I1', 'J1');
 syms('XXMAX', 'YMAX', 'ERR1', 'TEMP', 'COND');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Iterative Refinement Method.\n');
 fprintf(1,'This program used the file CHIP.m.\n');
 fprintf(1,'The array will be input from a text file in the order\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,n+1) \n');
 fprintf(1,' A(2,1), A(2,2), ..., A(2,n+1) \n');
 fprintf(1,',..., A(n,1), A(n,2), ..., A(n,n+1)\n');
 fprintf(1,'Place as many entries as desired on each line, but separate\n');
 fprintf(1,'entries with ');
 fprintf(1,'at least one blank.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 A1 = input(' ','s');
 OK = FALSE;
 if A1 == 'Y' | A1 == 'y' 
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
 B = zeros(N,N+1);
 NROW = zeros(1,N);
 X = zeros(1,N);
 XX = zeros(1,N);
 R = zeros(1,N);
 for I = 1 : N 
 for J = 1 : N+1 
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 OK = TRUE;
 fclose(INP);
 else
 fprintf(1,'The number must be a positive integer\n');
 end;
 end;
% NN is used for the maximun number of iterations
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input maximum number of iterations.\n');
 NN = input(' ');
 if NN > 0 
 OK = TRUE;
 else
 fprintf(1,'Number must be a positive integer.\n');
 end;
 end;
 OK = FALSE;
 fprintf(1,'Choice of rounding or chopping:\n');
 fprintf(1,'1. Rounding\n');
 fprintf(1,'2. Chopping\n');
 fprintf(1,'Enter 1 or 2.\n');
 RND = input(' ');
 while OK == FALSE 
 fprintf(1,'Input number of digits D <= 8 of rounding\n');
 D = input(' ');
 if D > 0 & D < 9
 OK = TRUE;
 else
 fprintf(1,'D must be a positive integer < 9\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input tolerance, which is usually 10^(-D).\n');
 TOL = input(' ');
 if TOL > 0 
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be a positive.\n');
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
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'ITERATIVE REFINEMENT METHOD\n\n');
 M = N+1;
 fprintf(OUP, 'Original system\n');
 for I = 1 : N 
 for J = 1 : M 
 fprintf(OUP,' %.10e',A(I,J));
 end;
 fprintf(OUP,'\n');
 end;
 if RND == 1 
 fprintf(OUP,'Rounding to %d Digits.\n',D);
 else fprintf(OUP,'Chopping to %d Digits.\n',D);
 end;
 fprintf(OUP,'\n Modified System \n');
 for I = 1 : N 
 NROW(I) = I;
 for J = 1 : M
 A(I,J) = CHIP(RND,D,A(I,J));
 B(I,J) = A(I,J);
 fprintf(OUP,'  %.10e', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
% NROW and B have been initialized, Gauss elimination will begin
% STEP 0
 I = 1;
 while I <= N-1 & OK == TRUE 
 KK = I;
 while abs(A(KK,I)) < 1.0e-20 & KK <= N 
 KK = KK+1;
 end;
 if KK > N 
 OK = false;
 fprintf(OUP, 'System does not have a unique solution.\n');
 else 
 if KK ~= I 
% Row interchange necessary
 IS = NROW(I);
 NROW(I) = NROW(KK);
 NROW(KK) = IS;
 for J = 1 : M 
 C = A(I,J);
 A(I,J) = A(KK,J);
 A(KK,J) = C;
 end;
 end;
 for J = I+1 : N 
 A(J,I) = CHIP(RND,D,A(J,I)/A(I,I));
 for L = I+1 : M 
 A(J,L) = CHIP(RND,D,A(J,L)-CHIP(RND,D,A(J,I)*A(I,L)));
 end;
 end;
 end;
 I = I+1;
 end;
 if abs(A(N,N)) < 1.0e-20 & OK == TRUE 
 OK = FALSE;
 fprintf(OUP, 'System has singular matrix\n');
 end;
 if OK == TRUE 
 fprintf(OUP, 'Reduced system\n');
 for I = 1 : N 
 for J = 1 : M 
 fprintf(OUP, '  %.10e', A(I,J));
 end;
 fprintf(OUP, '\n');
 end;
 X(N) = CHIP(RND,D,A(N,M)/A(N,N));
 for I = 1 : N-1 
 J = N-I;
 S = 0.0;
 for L = J+1 : N 
 S = CHIP(RND,D,S-CHIP(RND,D,A(J,L)*X(L)));
 end;
 S = CHIP(RND,D,A(J,M)+S);
 X(J) = CHIP(RND,D,S/A(J,J));
 end;
 end;
 fprintf(OUP, 'Initial solution\n');
 for I = 1 : N 
 fprintf(OUP,'  %.10e', X(I));
 end;
 fprintf(OUP, '\n');
% Refinement begins
% STEP 1
 if OK == TRUE 
 K = 1;
 for I = 1 : N 
 XX(I) = X(I);
 end;
 end;
% STEP 2
 while OK == TRUE & K <= NN 
% LL is set to 1 if the desired accuracy in any component is not 
% achieved. Thus LL is initially 0 for each iteration.
 LL = 0;
% STEP 3
 for I = 1 : N 
 R(I) = 0;
 for J = 1 : N 
 R(I) = CHIP(RND,2*D,R(I)-CHIP(RND,2*D,B(I,J)*XX(J)));
 end;
 R(I) = CHIP(RND,2*D,B(I,M)+R(I));
 end;
 fprintf(OUP, 'Residual number %d\n', K);
 for I = 1 : N 
 R(I) = CHIP(RND,D,R(I));
 fprintf(OUP, '%18.10e ', R(I));
 end;
 fprintf(OUP, '\n');
% STEP 4
% Solve the linear system in the same order as in step 0.  
% The solution will be placed in X instead of Y.
 for I = 1 : N-1 
 I1 = NROW(I);
 for J = I+1 : N 
 J1 = NROW(J);
 R(J1) = CHIP(RND,D,R(J1)-CHIP(RND,D,A(J,I)*R(I1)));
 end;
 end;
 X(N) = CHIP(RND,D,R(NROW(N))/A(N,N));
 for I = 1 : N-1 
 J = N-I;
 S = 0;
 for L = J+1 : N 
 S = CHIP(RND,D,S-CHIP(RND,D,A(J,L)*X(L)));
 end;
 S = CHIP(RND,D,S+R(NROW(J)));
 X(J) = CHIP(RND,D,S/A(J,J));
 end;
 fprintf(OUP, 'Vector Y\n');
 for I = 1 : N 
 fprintf(OUP,'%18.10e ', X(I));
 end;
 fprintf(OUP, '\n');
% Steps 5 and 6
 XXMAX = 0;
 YMAX = 0;
 ERR1 = 0;
 for I = 1 : N 
% If not accurate set LL to 1
 if abs(X(I)) > TOL 
 LL = 1;
 end;
 if K == 1 
 if abs(X(I)) > YMAX 
 YMAX = abs(X(I));
 end;
 if abs(XX(I)) > XXMAX 
 XXMAX = abs(XX(I));
 end;
 end;
 TEMP = XX(I);
 XX(I) = CHIP(RND,D,XX(I)+X(I));
 TEMP = abs(TEMP-XX(I));
 if TEMP > ERR1 
 ERR1 = TEMP;
 end;
 end;
 if ERR1 <= TOL 
 LL = 2;
 end;
 if K == 1 
 COND = YMAX/XXMAX*10^D;
 end;
 fprintf(OUP, 'New approximation\n');
 for I = 1 : N 
 fprintf(OUP, '%18.10e ', XX(I));
 end;
 fprintf(OUP, '\n');
% STEP 7
 if LL == 0 
 fprintf(OUP, 'The above vector is the solution.\n');
 OK = FALSE;
 else
 if LL == 2 
 fprintf(OUP,'The above vector is the best possible\n');
 fprintf(OUP,'with TOL = %18.10e \n',TOL);
 OK = FALSE;
 else
 K = K+1;
 end
 end;
% STEP 8 is not used in this implementation
 end;
 if K > NN 
 fprintf(OUP, 'Maximum Number of Iterations Exceeded.\n');
 end;
 fprintf(OUP, 'Condition number is %.10e\n', COND);
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
