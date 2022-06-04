% CRANK-NICOLSON ALGORITHM 12.3
%
% To approximate the solution of the parabolic partial-differential
% equation subject to the boundary conditions
%            u(0,t) = u(l,t) = 0, 0 < t < T = max t
% and the initial conditions
%             u(x,0) = F(x), 0 <= x <= l:
%
% INPUT:   endpoint l; maximum time T; constant ALPHA; integers m, N:
%
% OUTPUT:  approximations W(I,J) to u(x(I),t(J)) for each
%          I = 1,..., m-1 and J = 1,..., N.
 syms('OK', 'FX', 'FT', 'ALPHA', 'M', 'N', 'M1', 'M2', 'H', 'K');
 syms('VV', 'V', 'I', 'L', 'U', 'J', 'T', 'Z', 'I1');
 syms('FLAG', 'NAME', 'OUP', 'X', 's', 'x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Crank-Nicolson Method.\n');
 fprintf(1,'Input the function F(X) in terms of x.\n');
 fprintf(1,'for example, sin(pi*x) \n');
 s = input(' ');
 F = inline(s,'x');
 fprintf(1,'The lefthand endpoint on the X-axis is 0.\n');
 OK =FALSE;
 while OK == FALSE 
 fprintf(1,'Input the righthand endpoint on the X-axis.\n');
 FX = input(' ');
 if FX <= 0 
 fprintf(1,'Must be positive number.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the maximum value of the time variable T.\n');
 FT = input(' ');
 if FT <= 0 
 fprintf(1,'Must be positive number.\n');
 else
 OK = TRUE;
 end;
 end;
 fprintf(1,'Input the constant alpha.\n');
 ALPHA = input(' ');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input integer m = number of intervals on X-axis\n');
 fprintf(1,'and N = number of time intervals - on separate lines.\n');
 fprintf(1,'Note that m must be 3 or larger.\n');
 M = input(' ');
 N = input(' ');
 if M <= 2 | N <= 0 
 fprintf(1,'Numbers are not within correct range.\n');
 else
 OK = TRUE;
 end;
 end;
 if OK == TRUE
 V = zeros(1,M);
 L = zeros(1,M);
 U = zeros(1,M);
 Z = zeros(1,M);
 M1 = M-1;
 M2 = M-2;
% STEP1
 H = FX/M;
 K = FT/N;
% VV is used for lambda
 VV = ALPHA^2*K/(H^2);
% set V(M) = 0
 V(M) = 0;
% STEP 2
 for I = 1 : M1 
 V(I) = F(I*H);
 end;
% STEP 3
% STEPS 3 through 11 solve a tridiagonal linear system
% using Crout reduction
 L(1) = 1+VV;
 U(1) = -VV/(2*L(1));
% STEP 4
 for I = 2 : M2 
 L(I) = 1+VV+VV*U(I-1)/2;
 U(I) = -VV/(2*L(I));
 end;
% STEP 5
 L(M1) = 1+VV+0.5*VV*U(M2);
% STEP 6
 for J = 1 : N 
% STEP 7
% current t(j)
 T = J*K;
 Z(1) = ((1-VV)*V(1)+VV*V(2)/2)/L(1);
% STEP 8
 for I = 2 : M1 
 Z(I) = ((1-VV)*V(I)+0.5*VV*(V(I+1)+V(I-1)+Z(I-1)))/L(I);
 end;
% STEP 9
 V(M1) = Z(M1);
% STEP 10
 for I1 = 1 : M2 
 I = M2-I1+1;
 V(I) = Z(I)-U(I)*V(I+1);
 end;
 end;
% STEP 11
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:  A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'CRANK-NICOLSON METHOD\n\n');
 fprintf(OUP, '  I     X(I)         W(X(I),%12.6e)\n', FT);
 for I = 1 : M1 
 X = I*H;
 fprintf(OUP, '%3d %11.8f %13.8f\n', I, X, V(I));
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
% STEP 12
