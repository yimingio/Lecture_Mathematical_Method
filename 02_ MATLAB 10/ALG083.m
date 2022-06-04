% FAST FOURIER TRANSFORM ALGORITHM 8.3
%
% To compute the coefficients in the discrete approximation
% for the data (x(J),y(J)), 0<=J<=2m-1 where m=2^p and
% x(J)=-pi+J*pi/m for 0<=J<=2m-1.
%
% INPUT:  m; y(0),y(1),...y(2m-1).
%
% OUTPUT: complex numbers c(0),...,c(2m-1); real numbers
%         a(0),...,a(m); b(1),...,b(m-1).
 syms('OK', 'FLAG', 'M', 'N', 'JJ', 'J', 'Y', 'A');
 syms('NAME', 'INP', 'F', 'Z', 'OUP', 'TW', 'N2', 'C');
 syms('NG', 'NU1', 'YY', 'WW', 'X', 'W', 'K', 'L', 'M1');
 syms('NP', 'T1', 'T3','x','s');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Fast Fourier Transform.\n\n');
 fprintf(1,'The file IBR.m is used by this program.\n');
 fprintf(1,'The user must make provisions if the\n');
 fprintf(1,'interval is not (-pi,pi).\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Choice of input method:\n');
 fprintf(1,'1. Input entry by entry from keyboard\n');
 fprintf(1,'2. Input data from a text file\n');
 fprintf(1,'3. Generate data using a function F\n');
 fprintf(1,'Choose 1, 2, or 3 please\n');
 FLAG = input(' ');
 if FLAG == 1 | FLAG == 2 | FLAG == 3 
 OK = TRUE;
 end;
 end;
 if FLAG == 1 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input m\n');
 M = input(' ');
 if M > 0 
 OK = TRUE;
 N = 2*M;
 Y = zeros(1,N);
 for JJ = 1 : N 
 J = JJ-1;
 fprintf(1,'Input y(%d).\n', J);
 Y(JJ) = input(' ');
 end;
 else
 fprintf(1,'Number must be a positive integer.\n');
 end;
 end;
 end;
 if FLAG == 2 
 fprintf(1,'Has a text file been created with the ');
 fprintf(1,'entries y(0),...,y(2m-1)\n');
 fprintf(1,'separated by a blank?\n');
 fprintf(1,'Enter Y or N\n');
 A = input(' ','s');
 if A == 'Y' | A == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input number m.\n');
 M = input(' ');
 N = 2*M;
 Y = zeros(1,N);
 if N > 0 
 for JJ = 1 : N 
 Y(JJ) = fscanf(INP, '%f',1);
 end;
 fclose(INP);
 OK = TRUE;
 else fprintf(1,'Number must be a positive integer.\n');
 end;
 end;
 else
 fprintf(1,'The program will end so the input file can ');
 fprintf(1,'be created.\n');
 OK = FALSE;
 end;
 end;
 if FLAG == 3 
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'for example:   cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input the number m.\n');
 M = input(' ');
 N = 2*M;
 Y = zeros(1,N);
 if N > 0 
 for JJ = 1 : N 
 Z = -pi+(JJ-1)*pi/M;
 Y(JJ) = F(Z);
 end;
 OK = TRUE;
 else fprintf(1,'Number must be a postive integer.\n');
 end;
 end; 
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
 fprintf(OUP, 'FAST FOURIER TRANSFORM\n\n');
 TW = log(2);
% STEP 1
% use N2 for m, NG for p, NU1 for q, WW for zeta
 N2 = floor(N/2);
% STEP 2
 C = zeros(1,N);
 for JJ = 1 : N 
 C(JJ) = Y(JJ);
 end;
 Z = N;
 NG = round(log(Z)/TW);
 NU1 = NG-1;
 YY = 2*pi*sqrt(-1)/N;
 WW = exp(YY);
% STEP 3
 W = zeros(1,N);
 for JJ = 1 : N2 
 X = 1;
 YY = 1;
 for J = 1 : JJ 
 YY = X*WW;
 X = YY;
 end;
 W(JJ) = X;
 W(N2+JJ) = -X;
 end;
% STEP 4
 K = 0;
% STEP 5
 for L = 1 : NG 
% STEP 6
 while K < N-1 
% STEP 7
 for JJ = 1 : N2 
% STEP 8
 Z = exp(NU1*TW);
 M1 = round(Z);
 M1 = floor(K/M1);
% IBR does the bit reversal
 NP = IBR(M1,NG);
% T1 is eta
 T1 = C(K+N2+1);
% STEP 9
 if NP ~= 0 
 X = T1;
 T1 = X*W(NP);
 end;
 C(K+N2+1) = C(K+1)-T1;
 C(K+1) = C(K+1)+T1;
% STEP 10
 K = K+1;
 end;
% STEP 11
 K = K+N2;
 end;
% STEP 12
 K = 0;
 N2 = floor(N2/2);
 NU1 = NU1-1;
 end;
% STEP 13
 while K < N-1 
% STEP 14
 JJ = IBR(K,NG);
% STEP 15
 if JJ > K 
 T3 = C(K+1);
 C(K+1) = C(JJ+1);
 C(JJ+1) = T3;
 end;
% STEP 16
 K = K+1;
 end;
% STEPS 17 and 18
 fprintf(OUP, 'Coefficients c(0), ... , c(2m-1)\n\n');
 for JJ = 1 : N 
 YY = -(JJ-1)*pi*sqrt(-1);
 X = exp(YY);
 YY = X*C(JJ);
 C(JJ) = YY/(0.5*N);
 K = JJ-1;
 fprintf(OUP, '%3d %.8f %.8f \n', K, real(C(JJ)),imag(C(JJ)));
 end;
 fprintf(OUP, '\nCoefficients a(0), ..., a(m)\n\n');
 for JJ = 1 : M+1 
 fprintf(OUP, '%.8f\n', real(C(JJ)));
 end;
 fprintf(OUP, '\nCoefficients b(1), ..., b(m-1)\n\n');
 for JJ = 2 : M  
 fprintf(OUP, '%.8f\n', imag(C(JJ))); 
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
