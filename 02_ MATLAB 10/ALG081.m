% PADE RATIONAL APPROXIMATION ALGORITHM 8.1
%
% To obtain the rational approximation
%
%     r(x) = p(x) / q(x)
%          = (p0 + p1*x + ... + Pn*x^n) / (q0 + q1*x + ... + qm*x^m)
%
% for a given function f(x):
%
% INPUT  nonnegative integers m and n.
%
% OUTPUT  coefficients q0, q1, ... , qm, p0, p1, ... , pn.
%
% The coefficients of the Maclaurin polynomial a0, a1,  ... could
% be calculated instead of input as is assumed in this program.
 syms('OK', 'LM', 'LN', 'BN', 'FLAG', 'I', 'AA', 'AAA');
 syms('NAME', 'INP', 'N', 'M', 'NROW', 'NN', 'Q', 'P', 'J');
 syms('A', 'IMAX', 'AMAX', 'JJ', 'IP', 'JP', 'NCOPY', 'I1');
 syms('J1', 'XM', 'K', 'N1', 'PP', 'N2', 'SUM', 'KK', 'LL', 'OUP');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Pade Approximation.\n\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input m and n on separate lines.\n');
 LM = input(' ');
 LN = input(' ');
 BN = LM+LN;
 if LM >= 0 & LN >= 0 
 OK = TRUE;
 else
 fprintf(1,'m and n must both be nonnegative.\n');
 end;
 if LM == 0 & LN == 0 
 OK = FALSE;
 fprintf(1,'Not both m and n can be zero\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'The MacLaurin coefficients a(0), a(1), ... , a(N)\n');
 fprintf(1,'are to be input.\n');
 fprintf(1,'Choice of input method:\n');
 fprintf(1,'1. Input entry by entry from keyboard\n');
 fprintf(1,'2. Input data from a text file\n');
 fprintf(1,'Choose 1 or 2 please\n');
 FLAG = input(' ');
 if FLAG == 1 | FLAG == 2 
 OK = TRUE;
 end;
 end;
 AA = zeros(1,BN+1);
 if FLAG == 1 
 fprintf(1,'Input in order a(0) to a(N)\n');
 for I = 0 : BN 
 fprintf(1,'Input A( %d ) \n',I);
 AA(I+1) = input(' ');
 end;
 end;
 if FLAG == 2 
 fprintf(1,'As many entries as desired can be placed\n');
 fprintf(1,'on each line of the file each separated by blank.\n');
 fprintf(1,'Has such a text file been created?\n');
 fprintf(1,'Enter Y or N\n');
 AAA = input(' ','s');
 if AAA == 'Y' | AAA == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 for I = 0 : BN 
 AA(I+1) = fscanf(INP, '%f',1);
 end;
 fclose(INP);
 else
 fprintf(1,'Please create the input file.\n');
 fprintf(1,'The program will end so the input file can ');
 fprintf(1,'be created.\n');
 OK = FALSE;
 end;
 end;
 if OK == TRUE 
% STEP 1
 N = BN;
 M = N+1;
% STEP 2 - performed in input
 NROW = zeros(1,N);
 for I = 1 : N 
 NROW(I) = I;
 end;
% initialize row pointer for linear system
 NN = N-1;
% STEP 3
 Q = zeros(1, LM + 1);
 P = zeros(1, LN + 1);
 A = zeros(N,N+1);
 Q(1) = 1;
 P(1) = AA(1);
% STEP 4
% Set up a linear system, but use A(i,j) instead of B(i,j).
 for I = 1 : N 
% STEP 5
 for J = 1 : I-1 
 if J <= LN 
 A(I,J) = 0;
 end;
 end;
% STEP 6
 if I <= LN 
 A(I,I) = 1;
 end;
% STEP 7
 for J = I+1 : LN 
 A(I,J) = 0;
 end;
% STEP 8
 for J = 1 : I 
 if  J <= LM 
 A(I,LN+J) = -AA(I-J+1);
 end;
 end;
% STEP 9
 for J = LN+I+1 : N 
 A(I,J) = 0;
 end;
% STEP 10
 A(I,N+1) = AA(I+1);
 end;
% Solve the linear system using partial pivoting.
 I = LN+1;
% STEP 11
 while OK == TRUE & I <= NN 
% STEP 12
 IMAX = NROW(I);
 AMAX = abs(A(IMAX,I));
 IMAX = I;
 JJ = I+1;
 for IP = JJ : N 
 JP = NROW(IP);
 if abs(A(JP,I)) > AMAX 
 AMAX = abs(A(JP,I));
 IMAX = IP;
 end;
 end;
% STEP 13
 if AMAX <= 1.0e-20 
 OK = false;
 else
% STEP 14
% simulate row interchange
 if NROW(I) ~= NROW(IMAX) 
 NCOPY = NROW(I);
 NROW(I) = NROW(IMAX);
 NROW(IMAX) = NCOPY;
 end;
 I1 = NROW(I);
% STEP 15
% Perform elimination.
 for J = JJ : N 
 J1 = NROW(J);
% STEP 16
 XM = A(J1,I)/A(I1,I);
% STEP 17
 for K = JJ : M 
 A(J1,K) = A(J1,K)-XM * A(I1,K);
 end;
% STEP 18
 A(J1,I) = 0;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
% STEP 19
 N1 = NROW(N);
 if abs(A(N1,N)) <= 1.0e-20 
 OK = FALSE;
% system has no unique solution
 else
% STEP 20
% Start backward substitution.
 if LM > 0 
 Q(LM+1) = A(N1,M)/A(N1,N);
 A(N1,M) = Q(LM+1);
 end;
 PP = 1;
% STEP 21
 for K = LN+1 : NN 
 I = NN-K+LN+1;
 JJ = I+1;
 N2 = NROW(I);
 SUM = A(N2,N+1);
 for KK = JJ : N 
 LL = NROW(KK);
 SUM = SUM-A(N2,KK)*A(LL,M);
 end;
 A(N2,M) = SUM/A(N2,I);
 Q(LM-PP+1) = A(N2,M);
 PP = PP+1;
 end;
% STEP 22
 for K = 1 : LN 
 I = LN-K+1;
 N2 = NROW(I);
 SUM = A(N2,N+1);
 for KK = LN+1 : N 
 LL = NROW(KK);
 SUM = SUM-A(N2,KK)*A(LL,M);
 end;
 A(N2,M) = SUM;
 P(LN-K+2) = A(N2,M);
 end;
% STEP 23
% procedure completed successfully
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'PADE RATIONAL APPROXIMATION\n\n');
 fprintf(OUP, 'Denominator Coefficients Q(0), ..., Q(M) \n');
 for I = 0 : LM 
 fprintf(OUP, ' %11.8f', Q(I+1));
 end;
 fprintf(OUP, '\n');
 fprintf(OUP, 'Numerator Coefficients P(0), ..., P(N)\n');
 for I = 0 : LN 
 fprintf(OUP, ' %11.8f', P(I+1));
 end;
 fprintf(OUP, '\n');
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
 end;
 if OK == FALSE 
 fprintf(1,'System has no unique solution\n');
 end;
 end;
