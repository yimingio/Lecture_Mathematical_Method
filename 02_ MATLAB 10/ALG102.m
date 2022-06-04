% BROYDEN ALGORITHM 10.2
%
% To approximate the solution of the nonlinear system F(X) = 0
% given an initial approximation X.
%
% INPUT:   Number n of equations and unknowns; initial
%          approximation X = (X(1),...,X(n)); tolerance TOL;
%          maximum number of iterations N.
%
% OUTPUT:  Approximate solution X = (X(1),...,X(n)) or a message
%          that the number of iterations was exceeded.
 syms('OK', 'N', 'I', 'ZZ', 'J', 's1', 's2', 'TOL', 'NN', 'X');
 syms('FLAG', 'NAME', 'OUP', 'A', 'V', 'B', 'I1', 'I2','kk');
 syms('C', 'K', 'SN', 'S', 'VV', 'Y', 'ZN', 'Z', 'P', 'U', 'KK');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Broyden Method for Nonlinear Systems.\n');
 fprintf(1,'The functions could be input or defined in code.\n');
 fprintf(1,'This code assumes input of functions - see \n');
 fprintf(1,'comments in code for alternate version.\n');
 fprintf(1,'This program also uses M-files JAC.M and FN.M \n');
 fprintf(1,'If the number of equations exceeds 7 then JAC.M\n');
 fprintf(1,'and FN.M must be changed.\n');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the number n of equations.\n');
 N = input(' ');
 if N >= 2 & N < 8
 OK = TRUE;
 else
 fprintf(1,'N must be an integer, 1 < N < 8.\n');
 end;
 end;
 s1 = cell(N,1);
 s2 = cell(N*N,1);
 for I = 1 : N
 fprintf(1,'Input the function F_(%d) in terms of y1...y%d \n',I,N);
 kk = input(' ');
 s1{I} = kk;
 end;
% Define components of F as follows:
% s1{1} = '3*y1-cos(y2*y3)-0.5';
% s1{2} = 'y1^2-81*(y2+0.1)^2+sin(y3)+1.06';
% s1{3} = 'exp(-y1*y2)+20*y3+(10*pi-3)/3';
 for I = 1 : N
 for J = 1 : N
 fprintf(1,'Input the partial of F_(%d) with respect to x_%d \n',I,J);
 fprintf(1,'in terms of y1, ..., y%d \n',N);
 kk = input(' ');
 s2{(I-1)*N+J} = kk;
 end;
 end;
% Define the entries of the Jacobian in row major ordering.
% s2{1} = '3';
% s2{2} = 'y3*sin(y2*y3)';
% s2{3} = 'y2*sin(y2*y3)';
% s2{4} = '2*y1';
% s2{5} = '-162*(y2+0.1)';
% s2{6} = 'cos(y3)';
% s2{7} = '-y2*exp(-y1*y2)';
% s2{8} = '-y1*exp(-y1*y2)';
% s2(9} = '20';
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input tolerance\n');
 TOL = input(' ');
 if TOL > 0
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be positive.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the maximum number of iterations.\n');
 NN = input(' ');
 if NN > 0
 OK = TRUE;
 else
 fprintf(1,'Must be a positive integer.\n');
 end;
 end;
 X = zeros(1,N);
 A = zeros(N,N);
 B = zeros(N,N);
 V = zeros(1,N);
 S = zeros(1,N);
 Y = zeros(1,N);
 U = zeros(1,N);
 Z = zeros(1,N);
 for I = 1 : N
 fprintf(1,'Input initial approximation X(%d).\n', I);
 X(I) = input(' ');
 end;
 if OK == TRUE
 fprintf(1,'Select output destination\n');
 fprintf(1,'1. Screen\n');
 fprintf(1,'2. Text file\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example  A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(1,'Select amount of output\n');
 fprintf(1,'1. Answer only\n');
 fprintf(1,'2. All intermediate approximations\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 fprintf(OUP, 'BROYDENS METHOD FOR NONLINEAR SYSTEMS\n\n');
 if FLAG == 2
 fprintf(OUP, 'Iteration, Approximation, Error\n');
 end;
% STEP 1
% A will hold the Jacobian for the initial approximation.
 for I = 1 : N
 for J = 1 : N
 ZZ = JAC(I,J,N,X,s2);
 A(I,J) = ZZ;
 end;
% Compute V = F(x(0))
 V(I) = FN(I,N,X,s1);
 end;
% STEP 2
% Invert the Jacobian.
 for I = 1 : N
 for J = 1 : N
 B(I,J) = 0;
 end;
 B(I,I) = 1;
 end;
 I = 1;
 while I <= N & OK == TRUE
 I1 = I+1;
 I2 = I;
 if I ~= N
 C = abs(A(I,I));
 for J = I1 : N
 if abs(A(J,I)) > C
 I2 = J;
 C = abs(A(J,I));
 end;
 end;
 if C <= 1.0e-20 
 OK = FALSE;
 else
 if I2 ~= I 
 for J = 1 : N
 C = A(I,J);
 A(I,J) = A(I2,J);
 A(I2,J) = C;
 C = B(I,J);
 B(I,J) = B(I2,J);
 B(I2,J) = C;
 end;
 end;
 end;
 else
 if abs(A(N,N)) <= 1.0e-20 
 OK = FALSE;
 end;
 end;
 if OK == TRUE 
 for J = 1 : N 
 if J ~= I 
 C = A(J,I)/A(I,I);
 for K = 1 : N 
 A(J,K) = A(J,K)-C*A(I,K);
 B(J,K) = B(J,K)-C*B(I,K);
 end;
 end;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
 for I = 1 : N 
 C = A(I,I);
 for J = 1 : N 
 A(I,J) = B(I,J)/C;
 end;
 end;
 else
 fprintf(1,'Jacobian has no inverse\n');
 end;
 if OK == TRUE 
% STEP 3
 K = 2;
% Note: S = S(1)
% Compute the product S = -Av and the L2 norm SN of S
 SN = 0;
 for I = 1 : N
 S(I) = 0;
 for J = 1 : N
 S(I) = S(I)-A(I,J)*V(J);
 end;
 SN = SN+S(I)^2;
 end;
 SN = sqrt(SN);
 for I = 1 : N 
 X(I) = X(I)+S(I);
 end;
 if FLAG == 2 
 fprintf(OUP,' %d',K-1);
 for I = 1 : N 
 fprintf(OUP,' %11.8f',X(I));
 end;
 fprintf(OUP,'\n %12.6e\n',SN);
 end;
% STEP 4
 while K <= NN & OK == TRUE 
% STEP 5
 for I = 1 : N
 VV = FN(I,N,X,s1);
 Y(I) = VV-V(I);
 V(I) = VV;
 end;
% Note: V = F(X(K)) and Y = Y(K)
% STEP 6
% Form Z = -Ay and norm ZN of Z.
 ZN = 0;
 for I = 1 : N
 Z(I) = 0;
 for J = 1 : N
 Z(I) = Z(I)-A(I,J)*Y(J);
 end;
 ZN = ZN+Z(I)*Z(I);
 end;
 ZN = sqrt(ZN);
% Note = Z = -A(K-1)^(-1)*Y(K)
% STEP 7
 P = 0;
% P will be S(K)^T*A(K)^(-1)*Y(K)
 for I = 1 : N 
 P = P-S(I)*Z(I);
 end;
% STEP 8
 for I = 1 : N 
 U(I) = 0;
 for J = 1 : N 
 U(I) = U(I)+S(J)*A(J,I);
 end;
 end;
% STEP 9
 for I = 1 : N 
 for J = 1 : N 
 A(I,J) = A(I,J)+(S(I)+Z(I))*U(J)/P;
 end;
 end;
% STEP 10
% Form S = -Av and norm SN of S.
 SN = 0;
 for I = 1 : N
 S(I) = 0;
 for J = 1 : N 
 S(I) = S(I)-A(I,J)*V(J);
 end;
 SN = SN+S(I)^2;
 end;
 SN = sqrt(SN);
% Note = A = A(K)^(-1) and S = -A(K)^(-1)*F(X(K))
% STEP 11
 for I = 1 : N 
 X(I) = X(I)+S(I);;
 end;
% Note: X = X(K+1)
 KK = K+1;
 if FLAG == 2 
 fprintf(OUP, ' %2d', K);
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n%12.6e\n', SN);
 end;
 if SN <= TOL 
% procedure completed successfully
 OK = FALSE;
 fprintf(OUP, 'Iteration number %d', K);
 fprintf(OUP, ' gives solution:\n\n');
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n\nto within tolerance %.10e\n\n', TOL);
 fprintf(OUP, 'Process is complete\n');
 else
% STEP 13
 K = KK;
 end;
 end;
 if K >= NN 
% STEP 14
 fprintf(OUP, 'Procedure does not converge in %d iterations\n', NN);
 end;
 end;
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
