% WAVE EQUATION FINITE-DIFFERENCE ALGORITHM 12.4
%
% To approximate the solution to the wave equation:
% subject to the boundary conditions
%              u(0,t) = u(l,t) = 0, 0 < t < T = max t
% and the initial conditions
%              u(x,0) = F(x) and Du(x,0)/Dt = G(x), 0 <= x <= l:
%
% INPUT:   endpoint l; maximum time T; constant ALPHA; integers m, N.
%
% OUTPUT:  approximations W(I,J) to u(x(I),t(J)) for each I = 0, ..., m
%          and J=0,...,N.
 syms('OK', 'FX', 'FT', 'ALPHA', 'M', 'N', 'M1', 'M2', 'N1', 'N2');
 syms('H', 'K', 'V', 'J', 'W', 'I', 'FLAG', 'NAME', 'OUP', 'X', 's', 'x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Finite-Difference Method for the Wave Equation.\n');
 fprintf(1,'Input the functions F(X) and G(X) in terms of x, \n');
 fprintf(1,'on separate lines. \n');
 fprintf(1,'For example:     sin(pi*x)\n');
 fprintf(1,'                     0    \n');
 s = input(' ');
 F = inline(s,'x');
 s = input(' ');
 G = inline(s,'x');
 fprintf(1,'The lefthand endpoint on the X-axis is 0.\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the righthand endpoint on the X-axis.\n');
 FX = input(' ');
 if FX <= 0 
 fprintf(1,'Must be a positive number.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the maximum value of the time variable T.\n');
 FT = input(' ');
 if FT <= 0 
 fprintf(1,'Must be a positive number.\n');
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
 W = zeros(M+1,N+1);
 M1 = M+1;
 M2 = M-1;
 N1 = N+1;
 N2 = N-1;
% STEP 1
% V is used for lambda
 H = FX/M;
 K = FT/N;
 V = ALPHA*K/H;
% STEP 2
 for J = 2 : N1 
 W(1,J) = 0;
 W(M1,J) = 0;
 end;
% STEP 3
 W(1,1) = F(0);
 W(M1,1) = F(FX);
% STEP 4
 for I = 2 : M 
 W(I,1) = F(H*(I-1));
 W(I,2) = (1-V^2)*F(H*(I-1))+V^2*(F(I*H)+F(H*(I-2)))/2+K*G(H*(I-1));
 end;
% STEP 5
 for J = 2 : N 
 for I = 2 : M 
 W(I,J+1) = 2*(1-V^2)*W(I,J)+V^2*(W(I+1,J)+W(I-1,J))-W(I,J-1);
 end;
 end;
% STEP 6
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
 fprintf(OUP, 'FINITE DIFFERENCE METHOD FOR THE WAVE EQUATION\n\n');
 fprintf(OUP, '  I    X(I)     W(X(I),%12.6e)\n', FT);
 for I = 1 : M1 
 X = (I-1)*H;
 fprintf(OUP, '%3d %11.8f %13.8f\n', I, X, W(I,N1));
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
% STEP 7
