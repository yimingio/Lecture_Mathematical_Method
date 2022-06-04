% NONLINEAR SHOOTING ALGORITHM 11.2
%
% To approximate the solution of the nonlinear boundary-value problem
%
%          Y'' = F(X,Y,Y'), A<=X<=B, Y(A) = ALPHA, Y(B) = BETA:
%
%
% INPUT:   Endpoints A,B; boundary conditions ALPHA, BETA; number of
%          subintervals N; tolerance TOL; maximum number of iterations M.
%
% OUTPUT:  Approximations W(1,I) TO Y(X(I)); W(2,I) TO Y'(X(I))
%          for each I=0,1,...,N or a message that the maximum
%          number of iterations was exceeded.
 syms('OK', 'A', 'B', 'ALPHA', 'BETA', 'TK', 'AA', 'N');
 syms('TOL', 'NN', 'FLAG', 'NAME', 'OUP', 'H', 'K', 'W1');
 syms('W2', 'U1', 'U2', 'I', 'X', 'T', 'K11', 'K12', 'K21');
 syms('K22', 'K31', 'K32', 'K41', 'K42', 'J', 's', 'x', 'y', 'z');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Nonlinear Shooting Method.\n');
 fprintf(1,'Input the function F(X,Y,Z) in terms of x, y, z.\n');
 fprintf(1,'followed by the partial of F with respect to y on the \n');
 fprintf(1,'next line followed by the partial of F with respect to \n');
 fprintf(1,'z or y-prime on the next line. \n');
 fprintf(1,'For example:   (32+2*x^3-y*z)/8 \n');
 fprintf(1,'               -z/8 \n');
 fprintf(1,'               -y/8 \n');
 s = input(' ');
 F = inline(s,'x','y','z');
 s = input(' ');
 FY = inline(s,'x','y','z');
 s = input(' ');
 FYP = inline(s,'x','y','z');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 A = input(' ');
 B = input(' ');
 if A >= B 
 fprintf(1,'Left endpoint must be less than right endpoint.\n');
 else OK = TRUE;
 end;
 end;
 fprintf(1,'Input Y(%.10e).\n', A);
 ALPHA = input(' ');
 fprintf(1,'Input Y(%.10e).\n', B);
 BETA = input(' ');
 TK = (BETA-ALPHA)/(B-A);
 fprintf(1,'TK = %.8e\n', TK);
 fprintf(1,'Input new TK? Enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'input new TK\n');
 TK = input(' ');
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input an integer > 1 for the number of subintervals.\n');
 N = input(' ');
 if N <= 1 
 fprintf(1,'Number must exceed 1.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input Tolerance.\n');
 TOL = input(' ');
 if TOL <= 0 
 fprintf(1,'Tolerance must be positive.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input maximum number of iterations.\n');
 NN = input(' ');
 if NN <= 0 
 fprintf(1,'Must be positive integer.\n');
 else
 OK = TRUE;
 end;
 end;
 if OK == TRUE 
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text File\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example  A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'NONLINEAR SHOOTING METHOD\n\n');
 fprintf(OUP, '  I    X(I)         W1(I)        W2(I)\n');
% STEP 1
 W1 = zeros(1,N+1);
 W2 = zeros(1,N+1);
 H = (B-A)/N;
 K = 1;
% TK already computed
 OK = FALSE;
% STEP 2
 while K <= NN & OK == FALSE 
% STEP 3
 W1(1) = ALPHA;
 W2(1) = TK;
 U1 = 0 ;
 U2 = 1;
% STEP 4
% Rung-Kutta method for systems is used in STEPS 5 and 6
 for I = 1 : N 
%  STEP 5
 X = A+(I-1)*H;
 T = X+0.5*H;
% STEP 6
 K11 = H*W2(I);
 K12 = H*F(X,W1(I),W2(I));
 K21 = H*(W2(I)+0.5*K12);
 K22 = H*F(T,W1(I)+0.5*K11,W2(I)+0.5*K12);
 K31 = H*(W2(I)+0.5*K22);
 K32 = H*F(T,W1(I)+0.5*K21,W2(I)+0.5*K22);
 K41 = H*(W2(I)+K32);
 K42 = H*F(X+H,W1(I)+K31,W2(I)+K32);
 W1(I+1) = W1(I)+(K11+2*(K21+K31)+K41)/6;
 W2(I+1) = W2(I)+(K12+2*(K22+K32)+K42)/6;
 K11 = H*U2;
 K12 = H*(FY(X,W1(I),W2(I))*U1+FYP(X,W1(I),W2(I))*U2);
 K21 = H*(U2+0.5*K12);
 K22 = H*(FY(T,W1(I),W2(I))*(U1+0.5*K11)+FYP(T,W1(I),W2(I))*(U2+0.5*K21));
 K31 = H*(U2+0.5*K22);
 K32 = H*(FY(T,W1(I),W2(I))*(U1+0.5*K21)+FYP(T,W1(I),W2(I))*(U2+0.5*K22));
 K41 = H*(U2+K32);
 K42 = H*(FY(X+H,W1(I),W2(I))*(U1+K31)+FYP(X+H,W1(I),W2(I))*(U2+K32));
 U1 = U1+(K11+2*(K21+K31)+K41)/6;
 U2 = U2+(K12+2*(K22+K32)+K42)/6;
 end;
% STEP 7
% test for accuracy
 if abs(W1(N+1)-BETA) < TOL 
% STEP 8
 I = 0;
 fprintf(OUP, '%3d %13.8f %13.8f %13.8f\n', I, A, ALPHA, TK);
 for I = 1 : N 
 J = I+1;
 X = A+I*H;
 fprintf(OUP, '%3d %13.8f %13.8f %13.8f\n', I, X, W1(J), W2(J));
 end;
 fprintf(OUP, 'Convergence in %d iterations\n', K);
 fprintf(OUP, ' t = %14.7e\n', TK);
% STEP 9
 OK = TRUE;
 else
% STEP 10
% Newton's method applied to improve TK
 TK = TK-(W1(N+1)-BETA)/U1;
 K = K+1;
 end;
 end;
% STEP 11
% method failed
 if OK == FALSE 
 fprintf(OUP, 'Method failed after %d iterations\n', NN);
 end;
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
