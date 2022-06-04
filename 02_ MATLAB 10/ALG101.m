% NEWTON'S METHOD FOR SYSTEMS ALGORITHM 10.1
%
% To approximate the solution of the nonlinear system F(X)=0 given
% an initial approximation X:
%
% INPUT:   Number n of equations and unknowns; initial approximation
%          X=(X(1),...,X(n)); tolerance TOL; maximum number of
%          iterations N.
%
% OUTPUT:  Approximate solution X=(X(1),...,X(n)) or a message
%          that the number of iterations was exceeded.
 syms('OK', 'N', 'I', 'J', 'P', 'TOL', 'NN', 'X','ZZ','kk');
 syms('FLAG', 'NAME', 'OUP', 'K', 'A', 'R', 's1', 's2');
 syms('K1', 'I1', 'Z1', 'IR1', 'IA1', 'J1', 'C1', 'L1', 'JA1');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Newton Method for Nonlinear Systems.\n');
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
 fprintf(1,'N must be an integer, 1 < N < 8 .\n');
 end;
 end;
 s1 = cell(N,1);
 s2 = cell(N*N,1);
 A = zeros(N,N+1);
 X = zeros(1,N);
 Y = zeros(1,N);
% Define components of F as follows:
% s1{1} = '3*y1-cos(y2*y3)-0.5';
% s1{2} = 'y1^2-81*(y2+0.1)^2+sin(y3)+1.06';
% s1{3} = 'exp(-y1*y2)+20*y3+(10*pi-3)/3';
 for I = 1 : N
 fprintf(1,'Input the function F_(%d) in terms of y1 ... y%d \n' ,I ,N);
 kk = input(' ');
 s1{I} = kk;
 end;
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
% s2{9}) = '20';
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the Tolerance.\n');
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
 for I = 1 : N
 fprintf(1,'Input initial approximation X(%d).\n', I);
 X(I) = input(' ');
 end;
 fprintf(1,'Select output destination\n');
 fprintf(1,'1. Screen\n');
 fprintf(1,'2. Text file\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2
 fprintf(1,'Input the file name in the form - drive\\:name.ext\n');
 fprintf(1,'for example   A:\\OUTPUT.DTA\n');
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
 fprintf(OUP, 'NEWTONS METHOD FOR NONLINEAR SYSTEMS\n\n');
 if FLAG == 2
 fprintf(OUP, 'Iteration, Approximation, Error\n');
 end;
% STEP 1
 K = 1;
% STEP 2
 while OK == TRUE & K <= NN
% STEP 3
 for I = 1 : N
 for J = 1 : N
 ZZ = JAC(I,J,N,X,s2);
 A(I,J) = ZZ;
 end;
 ZZ = -FN(I,N,X,s1);
 A(I,N+1) = ZZ;
 end;
% STEP 4
 K1 = N-1;
 OK = TRUE;
 I1 = 1;
 while OK == TRUE & I1 <= K1
 Z1 = abs(A(I1,I1));
 IR1 = I1;
 IA1 = I1+1;
 for J1 = IA1 : N
 if abs(A(J1,I1)) > Z1
 IR1 = J1;
 Z1 = abs(A(J1,I1));
 end;
 end;
 if Z1 <= 1.0e-20
 OK = FALSE;
 else
 if IR1 ~= I1
 for J1 = I1 : N+1
 C1 = A(I1,J1);
 A(I1,J1) = A(IR1,J1);
 A(IR1,J1) = C1;
 end;
 end;
 for J1 = IA1 : N
 C1 = A(J1,I1)/A(I1,I1);
 if abs(C1) <= 1.0e-20
 C1 = 0;
 end;
 for L1 = I1 : N+1
 A(J1,L1) = A(J1,L1)-C1*A(I1,L1);
 end;
 end;
 end;
 I1 = I1+1;
 end;
 if OK == TRUE
 if abs(A(N,N)) <= 1.0e-20
 OK = FALSE;
 else
 Y(N) = A(N,N+1)/A(N,N);
 for I1 = 1 : K1
 J1 = N-I1;
 JA1 = J1+1;
 C1 = A(J1,N+1);
 for L1 = JA1 : N
 C1 = C1-A(J1,L1)*Y(L1);
 end;
 Y(J1) = C1/A(J1,J1);
 end;
 end;
 end;
 if OK == FALSE
 fprintf(1,'Linear system is singular\n');
 end;
 if OK == TRUE
% STEP 5
 R = 0;
 for I = 1 : N
 if abs(Y(I)) > R
 R = abs(Y(I));
 end;
 X(I) = X(I)+Y(I);
 end;
 if FLAG == 2
 fprintf(OUP, ' %2d', K);
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n%12.6e\n', R);
 end;
% STEP 6
 if R < TOL
 OK = FALSE;
 fprintf(OUP, 'Iteration %d gives solution:\n\n', K);
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n\nto within tolerance %.10e\n', TOL);
% STEP 7
 else
 K = K+1;
 end;
 end;
 end;
 if K > NN
% STEP 8
 fprintf(OUP, 'Procedure does not converge in %d iterations\n', NN);
 end;
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
