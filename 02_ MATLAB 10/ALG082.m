% CHEBYSHEV RATIONAL APPROXIMATION ALGORITHM 8.2
%
% To obtain the rational approximation
%
% rT(x) = (p0*T0 + p1*T1 +...+ pn*Tn) / (q0*T0 + q1*T1 +...+ qm*Tm)
%
% for a given function f(x):
%
% INPUT  nonnegative integers m and n.
%
% OUTPUT  coefficients q0, q1, ... , qm, p0, p1, ... , pn.
%
% The coefficients of the Chebyshev expansion a0, a1, ..., aN could
% be calculated instead of input as is assumed in this program.
 syms('OK', 'LM', 'LN', 'BN', 'FLAG', 'I', 'AA', 'AAA', 'NAME');
 syms('INP', 'N', 'M', 'NROW', 'NN', 'Q', 'J', 'A', 'PP', 'IMAX');
 syms('AMAX', 'JJ', 'IP', 'JP', 'NCOPY', 'I1', 'J1', 'XM', 'K');
 syms('N1', 'N2', 'SUM', 'KK', 'LL', 'P', 'OUP');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Chebyshev Rational Approximation.\n\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input m and n on separate lines.\n');
 LM = input(' ');
 LN = input(' ');
 BN = LM+LN;
 if LM >= 0 & LN >= 0 
 OK = TRUE;
 else fprintf(1,'m and n must both be nonnegative.\n');
 end;
 if LM == 0 & LN == 0 
 OK = FALSE;
 fprintf(1,'Not both m and n can be zero\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'The Chebyshev coefficients a(0), a(1), ... , a(N+m)\n');
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
 AA = zeros(1,BN+LM+1);
 NROW =  zeros(1,BN+1);
 P = zeros(1,LN+1);
 Q = zeros(1,LM+1);
 A = zeros(BN+1,BN+2);
 if FLAG == 1 
 fprintf(1,'Input in order a(0) to a(N+m)\n');
 for I = 0 : BN+LM 
 fprintf(1,'Input A(%d)\n', I);
 AA(I+1) = input(' ');
 end;
 end;
 if FLAG == 2 
 fprintf(1,'The text file may contain as many entries\n');
 fprintf(1,'per line as desired each separated by blank.\n');
 fprintf(1,'Has such a text file been created?\n');
 fprintf(1,'Enter Y or N\n');
 AAA = input(' ','s');
 if AAA == 'Y' | AAA == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'for example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 for I = 0 : BN+LM 
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
% STEP 2 - performed on input
 for I = 1 : M 
 NROW(I) = I;
 end;
% initialize row pointer
 NN = N-1;
% STEP 3
 Q(1) = 1.0;
% STEP 4
% set up a linear system with matrix A instead of B
 for I = 0 : N 
% STEP 5
 for J = 0 : I 
 if J <= LN 
 A(I+1,J+1) = 0;
 end;
 end;
% STEP 6
 if I <= LN 
 A(I+1,I+1) = 1.0;
 end;
% STEP 7
 for J = I+1 : LN 
 A(I+1,J+1) = 0;
 end;
% STEP 8
 for J = LN+1 : N 
 if I ~= 0 
 PP = I-J+LN;
 if PP < 0 
 PP = -PP;
 end;
 A(I+1,J+1) = -(AA(I+J-LN+1)+AA(PP+1))/2.0;
 else
 A(I+1,J+1) = -AA(J-LN+1)/2.0;
 end;
 end;
 A(I+1,N+2) = AA(I+1);
 end;
% STEP 9
 A(1,N+2) = A(1,N+2)/2.0;
% STEPS 10 -21 solve the linear system using partial pivoting
 I = LN+2;
% STEP 10
 while OK == TRUE & I <= N 
% STEP 11
 IMAX = NROW(I);
 AMAX = abs(A(IMAX,I));
 IMAX = I;
 JJ = I+1;
 for IP = JJ : N + 1 
 JP = NROW(IP);
 if abs(A(JP,I)) > AMAX 
 AMAX = abs(A(JP,I));
 IMAX = IP;
 end;
 end;
% STEP 12
 if AMAX <= 1.0e-20 
 OK = false;
 else
% STEP 13
% simulate row interchange
 if NROW(I) ~= NROW(IMAX) 
 NCOPY = NROW(I);
 NROW(I) = NROW(IMAX);
 NROW(IMAX) = NCOPY;
 end;
 I1 = NROW(I);
% STEP 14
% perform elimination
 for J = JJ : M  
 J1 = NROW(J);
% STEP 15
 XM = A(J1,I)/A(I1,I);
% STEP 16
 for K = JJ : M + 1
 A(J1,K) = A(J1,K)-XM*A(I1,K);
 end;
% STEP 17
 A(J1,I) = 0;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE
% STEP 18
 N1 = NROW(N+1);
 if abs(A(N1,N+1)) <= 1.0e-20 
 OK = false;
% system has no unique solution
 else
% STEP 19
% start backward substitution
 if LM > 0 
 Q(LM+1) = A(N1,M+1)/A(N1,N+1);
 A(N1,M+1) = Q(LM+1);
 end;
 PP = 1;
% STEP 20
 for K = LN+2 : N 
 I = N-K+LN+2;
 JJ = I+1;
 N2 = NROW(I);
 SUM = A(N2,M+1);
 for KK = JJ : N + 1
 LL = NROW(KK);
 SUM = SUM - A(N2,KK) * A(LL,M+1);
 end;
 A(N2,M+1) = SUM / A(N2,I);
 Q(LM-PP+1) = A(N2,M+1);
 PP = PP+1;
 end;
% STEP 21
 for K = 1 : LN + 1 
 I = LN+1-K+1;
 N2 = NROW(I);
 SUM = A(N2,M+1);
 for KK = LN+2 : N + 1 
 LL = NROW(KK);
 SUM = SUM-A(N2,KK)*A(LL,M+1);
 end;
 A(N2,M+1) = SUM ;
 P(LN-K+2) = A(N2,M+1);
 end;
% STEP 22
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
 fprintf(OUP, 'CHEBYSHEV RATIONAL APPROXIMATION\n\n');
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
