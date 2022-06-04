% RUNGE-KUTTA FOR SYSTEMS OF DIFFERENTIAL EQUATIONS ALGORITHM 5.7
%
% TO APPROXIMATE THE SOLUTION OF THE MTH-ORDER SYSTEM OF FIRST-
% ORDER INITIAL-VALUE PROBLEMS
%            UJ' = FJ( T, U1, U2, ..., UM ), J = 1, 2, ..., M
%            A <= T <= B, UJ(A) = ALPHAJ, J = 1, 2, ..., M
% AT (N+1) EQUALLY SPACED NUMBERS IN THE INTERVAL (A,B).
%
% INPUT:   ENDPOINTS A,B; NUMBER OF EQUATIONS M; INITIAL
%          CONDITIONS ALPHA1, ..., ALPHAM; INTEGER N.
%
% OUTPUT:  APPROXIMATION WJ TO UJ(T) AT THE (N+1) VALUES OF T.
 syms('OK', 'M', 'I', 'A', 'B', 'ALPHA', 'N', 'FLAG');
 syms('NAME', 'OUP', 'H', 'T', 'J', 'W', 'L', 'K','ss');
 syms('K1','K2','K3','K4','Z','kk');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Runge-Kutta Method for Systems of m equations\n');
 fprintf(1,'This program uses the file F.m.  If the number of equations\n');
 fprintf(1,'exceeds 7, then F.m must be changed.\n');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the number of equations\n');
 M = input(' ');
 if M <= 0 | M > 7
 fprintf(1,'Number must be a positive integer < 8\n');
 else
 OK = TRUE;
 end;
 end;
 ss = cell(M,1);
 for I = 1:M
 fprintf(1,'Input the function F_(%d) in terms of t and y1 ... y%d\n', I,M);
 fprintf(1,'For example: y1-t^2+1 \n');
 kk = input(' ');
 ss{I} = kk;
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 A = input(' ');
 B = input(' ');
 if A >= B
 fprintf(1,'Left endpoint must be less than right endpoint\n');
 else
 OK = TRUE;
 end;
 end;
 ALPHA = zeros(1,M);
 for I = 1:M
 fprintf(1,'Input the initial condition alpha(%d)\n', I);
 ALPHA(I) = input(' ');
 end;
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input a positive integer for the number of subintervals\n');
 N = input(' ');
 if N <= 0
 fprintf(1,'Number must be a positive integer\n');
 else
 OK = TRUE;
 end;
 end;
 if OK == TRUE
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
 fprintf(OUP,'RUNGE-KUTTA METHOD FOR SYSTEMS OF DIFFERENTIAL EQUATIONS\n\n');
 fprintf(OUP, '    T');
 for I = 1:M
 fprintf(OUP, '          W%d', I);
 end;
% STEP 1
 W = zeros(1,M);
 V = zeros(1,M+1);
 K1 = zeros(1,M);
 K2 = zeros(1,M);
 K3 = zeros(1,M);
 K4 = zeros(1,M);
 H = (B-A)/N;
 T = A;
% STEP 2
 for J = 1:M
 W(J) = ALPHA(J);
 end;
% STEP 3
 fprintf(OUP, '\n%5.3f', T);
 for I = 1:M
 fprintf(OUP, ' %11.8f', W(I));
 end;
 fprintf(OUP, '\n');
% STEP 4
 for L = 1:N
% STEP 5
 V(1) = T;
 for J = 2:M+1
 V(J) = W(J-1);
 end;
 for J = 1:M
 Z = H*F(J,M,V,ss);
 K1(J) = Z;
 end;
% STEP 6
 V(1) = T+H/2;
 for J = 2:M+1
 V(J) = W(J-1)+K1(J-1)/2;
 end;
 for J = 1:M
 Z = H*F(J,M,V,ss);
 K2(J) = Z;
 end;
% STEP 7
 for J = 2:M+1
 V(J) = W(J-1)+K2(J-1)/2;
 end;
 for J = 1:M
 Z = H*F(J,M,V,ss);
 K3(J) = Z;
 end;
% STEP 8
 V(1) = T + H;
 for J = 2:M+1
 V(J) = W(J-1)+K3(J-1);
 end;
 for J = 1:M
 Z = H*F(J,M,V,ss);
 K4(J) = Z;
 end;
% STEP 9
 for J = 1:M
 W(J) = W(J)+(K1(J)+2.0*K2(J)+2.0*K3(J)+K4(J))/6.0;
 end;
% STEP 10
 T = A+L*H;
% STEP 11
 fprintf(OUP, '%5.3f', T);
 for I = 1:M
 fprintf(OUP, ' %11.8f', W(I));
 end;
 fprintf(OUP, '\n');
 end;
% STEP 12
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
