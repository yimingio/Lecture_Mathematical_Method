% POWER METHOD ALGORITHM 9.1
%
% To approximate the dominant eigenvalue and an associated
% eigenvector of the n by n matrix A given a nonzero vector x:
%
% INPUT:   Dimension n; matrix A; vector x; tolerance TOL; maximum
%          number of iterations N.
%
% OUTPUT:  Approximate eigenvalue MU; approximate eigenvector x
%          or a message that the maximum number of iterations was
%          exceeded.
 syms('OK', 'AA', 'NAME', 'INP', 'N', 'I', 'J', 'A', 'X', 'TOL');
 syms('NN', 'FLAG', 'OUP', 'K', 'LP', 'AMAX', 'Y', 'YMU');
 syms('ERR', 'T');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Power Method.\n');
 OK = FALSE;
 fprintf(1,'The array will be input from a text file in the order:\n');
 fprintf(1,'A(1,1), A(1,2), ..., A(1,n), \n');
 fprintf(1,'A(2,1), A(2,2), ..., A(2,n),\n');
 fprintf(1,'..., A(n,1), A(n,2), ..., A(n,n)\n\n');
 fprintf(1,'Place as many entries as desired on each line, but ');
 fprintf(1,'separate entries with\n');
 fprintf(1,'at least one blank.\n');
 fprintf(1,'The initial approximation should follow in same format.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input (' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - drive:name.ext\n');
 fprintf(1,'for example: A:DATA.DTA \n');
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
% use NN for N
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
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'POWER METHOD\n\n');
 fprintf(OUP, 'iter  approx       approx eigenvector\n');
 fprintf(OUP, '      eigenvalue\n');
% STEP 1
 K = 1;
% STEP 2
 LP = 1;
 AMAX = abs(X(1));
 for I = 2 : N 
 if abs(X(I)) > AMAX 
 AMAX = abs(X(I));
 LP = I;
 end;
 end;
% STEP 3
 for I = 1 : N 
 X(I) = X(I)/AMAX;
 end;
% STEP 4
 while K <= NN & OK == TRUE 
% STEP 5
 for I = 1 : N 
 Y(I) = 0;
 for J = 1 : N 
 Y(I) = Y(I) + A(I,J) * X(J);
 end;
 end;
% STEP 6
 YMU = Y(LP);
% STEP 7
 LP = 1;
 AMAX = abs(Y(1));
 for I = 2 : N 
 if abs(Y(I)) > AMAX 
 AMAX = abs(Y(I));
 LP = I;
 end;
 end;
% STEP 8
 if AMAX <= 0 
 fprintf(1,'0 eigenvalue - select another ');
 fprintf(1,'initial vector and begin again\n');
 OK = FALSE;
 else
% STEP 9
 ERR = 0;
 for I = 1 : N 
 T = Y(I)/Y(LP);
 if abs(X(I)-T) > ERR 
 ERR = abs(X(I)-T);
 end;
 X(I) = T;
 end;
 fprintf(OUP, '%d    %12.8f', K, YMU);
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
% STEP 10
 if ERR <= TOL 
 fprintf(OUP, '\n\nThe eigenvalue = %12.8f',YMU);
 fprintf(OUP, ' to tolerance = %.10e\n', TOL);
 fprintf(OUP, 'obtained on iteration number = %d\n\n', K);
 fprintf(OUP, 'Unit eigenvector is :\n\n');
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
 OK = FALSE;
 end;
% STEP 11
 K = K+1;
 end;
 end;
% STEP 12
 if K > NN 
 fprintf(1,'Method did not converge within %d iterations\n', NN);
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
