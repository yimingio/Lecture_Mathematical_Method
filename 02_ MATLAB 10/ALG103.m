% STEEPEST DESCENT ALGORITHM 10.3
%
% To approximate a solution P to the minimization problem
%                G(P) = MIN( G(X) : X in R(n) )
% given an initial approximation X:
%
% INPUT:   Number n of variables; initial approximation X;
%          tolerance TOL; maximum number of iterations N.
%
% OUTPUT:  Approximate solution X or a message of failure.
 syms('OK', 'N', 'I', 'P', 'J', 'TOL', 'NN', 'X', 'FLAG1');
 syms('NAME', 'OUP', 'K', 'G', 'Z0', 'Z', 'A', 'X0', 'C', 'AA');
 syms('G0', 'FLAG', 'H1', 'H2', 'H3', 'A0','s1','s2','ZZ','kk');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Steepest Descent Method.\n');
 fprintf(1,'The functions could be input or defined in code.\n');
 fprintf(1,'This code assumes functions are defined in code - see \n');
 fprintf(1,'comments in code for alternate version.\n');
 fprintf(1,'This program also uses M-files JAC.M, FF.M and CF.M\n');
 fprintf(1,'If the number of equations exceeds 7 then JAC.M,\n');
 fprintf(1,'FF.M and CF.M must be changed.\n');
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
% The next segment of code allows for input of functions.
% fprintf(1,'The function CF_(I) is the Ith component of F\n');
% for I = 1 : N
% fprintf(1,'Input the function CF_(%d) in terms of y1 ... y%d\n',I,N);
% kk = input(' ');
% S1{I} = kk;
% end;
% for I = 1 : N
% for J = 1 : N
% fprintf(1,'Input the partial of CF_(%d) with respect to x_%d \n',I,J);
% fprintf(1,'in terms of y1 ... y%d \n',N);
% kk = input(' ');
% s2{(I-1)*N+J} = kk;
% end;
% end;
% Define the components of F as follows:
 s1{1} = '3*y1-cos(y2*y3)-0.5';
 s1{2} = 'y1^2-81*(y2+0.1)^2+sin(y3)+1.06';
 s1{3} = 'exp(-y1*y2)+20*y3+(10*pi-3)/3';
% Define the entries of the Jacobian in row major ordering.
 s2{1} = '3';
 s2{2} = 'y3*sin(y2*y3)';
 s2{3} = 'y2*sin(y2*y3)';
 s2{4} = '2*y1';
 s2{5} = '-162*(y2+0.1)';
 s2{6} = 'cos(y3)';
 s2{7} = '-y2*exp(-y1*y2)';
 s2{8} = '-y1*exp(-y1*y2)';
 s2{9} = '20';
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
 for I = 1 : N
 fprintf(1,'Input initial approximation X(%d).\n', I);
 X(I) = input(' ');
 end;
 if OK == TRUE
 fprintf(1,'Select output destination\n');
 fprintf(1,'1. Screen\n');
 fprintf(1,'2. Text file\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG1 = input(' ');
 if FLAG1 == 2
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
 FLAG1 = input(' ');
 fprintf(OUP, 'STEEPEST DESCENT METHOD FOR NONLINEAR SYSTEMS\n\n');
 if FLAG1 == 2
 fprintf(OUP, 'Iteration, Approximation\n');
 end;
% STEP 1
 K = 1;
 G = zeros(1,3);
 Z = zeros(1,N);
 A = zeros(1,3);
 C = zeros(1,N);
 AA = zeros(N,N);
% STEP 2
 while OK == TRUE & K <= NN
% STEP 3
 G(1) = FF(N,X,s1);
% AA is the  Jacobian
 for I = 1 : N
 for J = 1 : N
 ZZ = JAC(I,J,N,X,s2);
 AA(I,J) = ZZ;
 end;
 end;
 Z0 = 0;
 for I = 1 : N
 ZZ = 0;
 for J = 1 : N
 ZZ = ZZ + 2*CF(J,N,X,s1)*AA(J,I);
 end;
 Z(I) = ZZ;
 Z0 = Z0+(Z(I))*(Z(I));
 end;
 Z0 = sqrt(Z0);
% STEP 4
 if Z0 <= 1.0e-20
 OK = FALSE;
 fprintf(OUP, '0 qradient - may have a minimum\n');
 else
% STEP 5
 for I = 1 : N
 Z(I) = Z(I) / Z0;
 end;
 A(1) = 0;
 X0 = 1;
 for I = 1 : N
 C(I) = X(I)-X0*Z(I);
 end;
 G0 = FF(N,C,s1);
% STEP 6
 FLAG = TRUE;
 if G0 < G(1)
 FLAG = FALSE;
 end;
 while FLAG == TRUE & OK == TRUE 
% STEPS 7 and 8
 X0 = 0.5*X0;
 if X0 <= 1.0e-20 
 OK = FALSE;
 fprintf(OUP, 'No likely improvement - may\n');
 fprintf(OUP, 'have a minimum\n');
 else
 for I = 1 : N 
 C(I) = X(I)-X0*Z(I);
 end;
 G0 = FF(N,C,s1);
 end;
 if G0 < G(1) 
 FLAG = FALSE;
 end;
 end;
 if OK == TRUE 
 A(3) = X0;
 G(3) = G0;
% STEP 9
 X0 = 0.5*X0;
 for I = 1 : N 
 C(I) = X(I)-X0*Z(I);
 end;
 A(2) = X0;
 G(2) = FF(N,C,s1);
% STEP 10
 H1 = (G(2)-G(1))/(A(2)-A(1));
 H2 = (G(3)-G(2))/(A(3)-A(2));
 H3 = (H2-H1)/(A(3)-A(1));
% STEP 11
 X0 = 0.5*(A(1)+A(2)-H1/H3);
 for I = 1 : N 
 C(I) = X(I)-X0*Z(I);
 end;
 G0 = FF(N,C,s1);
% STEP 12
 A0 = X0;
 for I = 1 : N 
 if abs(G(I)) < abs(G0) 
 A0 = A(I);
 G0 = G(I);
 end;
 end;
 if abs(A0) <= 1.0e-20 
 OK = FALSE;
 fprintf(OUP, 'No change likely\n');
 fprintf(OUP, '- probably rounding error problems\n');
 else
% STEP 13
 for I = 1 : N 
 X(I) = X(I)-A0*Z(I);
 end;
% STEP 14
 if FLAG1 == 2 
 fprintf(OUP, ' %2d', K);
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
 end;
 if abs(G0) < TOL | abs(G0-G(1)) < TOL 
 OK = FALSE;
 fprintf(OUP, 'Iteration number %d\n', K);
 fprintf(OUP, 'gives solution\n\n');
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n\nto within %.10e\n\n', TOL);
 fprintf(OUP, 'Process is complete\n');
 else
% STEP 15
 K = K+1;
 end;
 end;
 end;
 end;
 end;
 if K > NN 
% STEP 16
 fprintf(OUP, 'Process does not converge in %d\n', NN);
 fprintf(OUP, ' iterations\n');
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
