% CONJUGATE GRADIENT ALGORITHM 7.5
%
% To solve Ax = b given the preconditioning matrix C inverse
% and an initial approximation
% x(0):
%
% INPUT:   the number of equations and unknowns n; the entries
%          A(I,J), 1<=I, J<=n, of the matrix A; the entries
%          B(I), 1<=I<=n, of the inhomogeneous term b; the
%          entries C({I,J), 1<=I, J<=n, of the preconditioning
%          matrix C inverse, entries XO(I), 1<=I<=n, of x(0);
%          tolerance TOL; maximum number of iterations N.
%
% OUTPUT:  the approximate solution X(1),...,X(n) and its
%          residual vector R(1),...,R(N) or a message
%          that the number of iterations was exceeded.
 syms('OK','AA','NAME','INP','N','I','J','A','X1','TOL','NN','W');
 syms('K','ERR','S','FLAG','OUP','R','T','ALPHA','BETA','U','V');
 syms('CI','QERR','ERR1','CT','SS','Z');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Conjugate Gradient Method for Linear Systems.\n');
 OK = FALSE;
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,n+1), A(2,1), A(2,2),\n');
 fprintf(1,'..., A(2,n+1),\n');
 fprintf(1,'..., A(n,1), A(n,2), ..., A(n,n+1)\n\n');
 fprintf(1,'Place as many entries as desired on each line,\n');
 fprintf(1,'but separate ');
 fprintf(1,'entries with at least one blank.\n');
 fprintf(1,'Do the same for the input of the inverse of C.\n');
 fprintf(1,'The initial approximation should follow \n');
 fprintf(1,'in same format.\n\n\n');
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
    CI = zeros(N,N);
    Q = zeros(N,N);
    X1 = zeros(1,N);
    R = zeros(1,N);
    W = zeros(1,N);
    V = zeros(1,N);
    U = zeros(1,N);
    Z = zeros(1,N);
 for I = 1 : N 
 for J = 1 : N+1
 A(I,J) = fscanf(INP, '%f',1);
 end;
 end;
 for I = 1 : N 
 for J = 1 : N 
 CI(I,J) = fscanf(INP, '%f',1);
 CT(J,I) = CI(I,J);
 end;
 end;
 for I = 1 : N 
 X1(I) = fscanf(INP, '%f',1);
 end;
% use X1 for X0
 OK = TRUE;
 fclose(INP);
 else
 fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the tolerance.\n');
 TOL = input(' ');
 if TOL > 0 
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be a positive number.\n');
 end;
 end;
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
 OUP = fopen(NAME, 'wt');
 else
 OUP = 1;
 end;
% STEP 1
 for I = 1 : N
 R(I) = A(I,N+1);
 for J = 1 : N
 R(I) = R(I)-A(I,J)*X1(J);
 end;
 end;
 for I = 1 : N
 W(I) = 0;
 for J = 1 : N
 W(I) = W(I)+CI(I,J)*R(J);
 end;
 end;
 for I = 1 : N
 V(I) = 0;
 for J = 1 : N
 V(I) = V(I)+CT(I,J)*W(J);
 end;
 end;
 ALPHA = 0.0;
 for I = 1 : N
 ALPHA = ALPHA + W(I)*W(I);
 end;
% Step 2
 K = 1;
 OK = FALSE;
% STEP 3
 while (OK == FALSE) & (K <= NN)
% ERR is used to test accuracy - it measures the 2 norm
 ERR = 0;
 for I = 1 : N
 ERR = ERR + V(I)*V(I);
 end;
% STEP 4
 if sqrt(ERR) < TOL
 K = K -1;
 OK = TRUE;
 else
% Step 5
 for I = 1 : N
 U(I) = 0.0;
 for J = 1 : N
 U(I) = U(I)+A(I,J)*V(J);
 end;
 end;
 S = 0.0;
 for I = 1 : N
 S = S + V(I)*U(I);
 end;
 T = ALPHA/S;
 for I = 1 : N
 X1(I) = X1(I)+T*V(I);
 R(I) = R(I) - T*U(I);
 end;
 fprintf(OUP, 'The approximation :\n');
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X1(I));
 end;
 fprintf(OUP, '\nafter %d iterations with\n', K);
 fprintf(OUP, 'The residual vector is :\n');
 for I = 1 : N
 fprintf(OUP, ' %11.8f', R(I));
 end;
 fprintf(OUP,' \n');
 for I = 1 : N
 W(I) = 0.0;
 for J = 1 : N
 W(I) = W(I)+CI(I,J)*R(J);
 end;
 end;
 BETA = 0.0;
 for I = 1 : N
 BETA = BETA + W(I)*W(I);
 end;
 ERR1 = sqrt(BETA);
% Step 6
 if ERR1 <= TOL
 ERR = 0.0;
 for I = 1 : N
 ERR = ERR + R(I)*R(I)
 end;
 if sqrt(ERR) < TOL
 OK = TRUE;
 end;
 end;
 if OK == FALSE
% Step 7
 K = K + 1;
 S = BETA/ALPHA;
 for I = 1 : N
 Z(I) = 0;
 for J = 1 : N
 Z(I) = Z(I) + CT(I,J)*W(J);
 end;
 end;
 for I = 1 : N
 V(I) = Z(I)+S*V(I);
 end;
 ALPHA = BETA;
 end;
 end;
 end;
 if OK == FALSE
% Step 8
 fprintf(1,'Maximum Number of Iterations Exceeded.\n');
% procedure completed unsuccessfully
 else
 fprintf(OUP, 'CONJUGATE GRADIENT METHOD FOR LINEAR SYSTEMS\n\n');
 fprintf(OUP, 'The solution vector is :\n');
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X1(I));
 end;
 fprintf(OUP, '\nusing %d iterations with\n', K);
 fprintf(OUP, 'Tolerance  %.10e in infinity-norm\n', TOL);
 fprintf(OUP, 'The residual vector is :\n');
 for I = 1 : N
 fprintf(OUP, ' %11.8f', R(I));
 end;
 fprintf(OUP,' \n');
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
 end;
