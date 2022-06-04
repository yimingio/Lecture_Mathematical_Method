% CONTINUATION METHOD FOR SYSTEMS ALGORITHM 10.4
%
% To approximate the solution of the nonlinear system F(X)=0 given
% an initial approximation X:
%
% INPUT:   Number n of equations and unknowns; initial approximation
%          X=(X(1),...,X(n)); number of RK4 steps N.
%
% OUTPUT:  Approximate solution X=(X(1),...,X(n)).
%         
 syms('OK', 'N', 'I', 'J', 'P', 'NN', 'X','ZZ');
 syms('FLAG', 'NAME', 'OUP', 'K', 'A', 's', 'ss','mm','s1','s2');
 syms('KK1', 'I1', 'Z1', 'IR1', 'IA1', 'J1', 'C1', 'L1', 'JA1');
 syms('K1','K2','K3','K4','X1','KK','H','b');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Continuation Method for Nonlinear Systems.\n');
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
 fprintf(1,'N must be an integer greater, 1 < N < 8.\n');
 end;
 end;
 s1 = cell(N,1);
 s2 = cell(N*N,1);
 A = zeros(N,N+1);
 K1 = zeros(4,N);
 X = zeros(1,N);
 Y = zeros(1,N);
 X1 = zeros(1,N);
 b = zeros(1,N);
 mm = zeros(1,4);
 mm(1) = 0.5;
 mm(2) = 0.5;
 mm(3) = 1.0;
 mm(4) = 0;
% Define components of F as follows:
 s1{1} = '3*y1-cos(y2*y3)-0.5';
 s1{2} = 'y1^2-81*(y2+0.1)^2+sin(y3)+1.06';
 s1{3} = 'exp(-y1*y2)+20*y3+(10*pi-3)/3';
% for I = 1 : N
% fprintf(1,'Input the function F_(%d) in terms of y1 ... y%d \n' ,I ,N);
% kkk = input(' ');
% s{I} = kkk;
% end;
% for I = 1 : N
% for J = 1 : N
% fprintf(1,'Input the partial of F_(%d) with respect to x_%d \n',I,J);
% fprintf(1,'in terms of y1, ..., y%d \n',N);
% kkk = input(' ');
% s2{(I-1)*N+J} = kkk;
% end;
% end;
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
 fprintf(1,'Input the number of RK4 steps.\n');
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
 fprintf(OUP, 'Continuation Method for Nonlinear Systems \n');
 fprintf(OUP, 'Iteration  Approximation\n');
% STEP 1
 H = 1/NN;
 for I = 1 : N
 ZZ = -H*FN(I,N,X,s1);
 b(I) = ZZ;
 end;
% STEP 2
 K = 1;
 while OK == TRUE & K <= NN
% STEPS 3 - 6
 for I = 1 : N
 X1(I) = X(I);
 end;
 KK = 1;
 while OK == TRUE & KK <= 4
 for I = 1 : N
 for J = 1 : N
 ZZ = JAC(I,J,N,X1,s2);
 A(I,J) = ZZ;
 end;
 A(I,N+1) = b(I);
 end;
 KK1 = N-1;
 OK = TRUE;
 I1 = 1;
 while OK == TRUE & I1 <= KK1
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
 for I1 = 1 : KK1
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
 for I = 1 : N
 K1(KK,I) = Y(I);
 X1(I) =  X(I) + mm(KK)*K1(KK,I);
 end;
 KK = KK + 1;
 end;
 end;
% STEP 7
 if OK == TRUE
 for I = 1 : N
 X(I) = X(I) + (K1(1,I)+2*K1(2,I)+2*K1(3,I)+K1(4,I))/6;
 end;
 fprintf(OUP, ' %2d',  K);
 for I = 1 : N
 fprintf(OUP, ' %11.8f', X(I));
 end;
 fprintf(OUP, '\n');
 end;
 K = K + 1;
 end;
% STEP 8
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
