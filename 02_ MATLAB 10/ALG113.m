% LINEAR FINITE-DIFFERENCE ALGORITHM 11.3
%
% To approximate the solution of the boundary-value problem
%
%    Y'' = P(X)Y' + Q(X)Y + R(X), A<=X<=B, Y(A) = ALPHA, Y(B) = BETA:
%
% INPUT:   Endpoints A, B; boundary conditions ALPHA, BETA;
%          integer N.
%
% OUTPUT:  Approximations W(I) to Y(X(I)) for each I=0,1,...,N+1.
 syms('OK', 'AA', 'BB', 'ALPHA', 'BETA', 'N', 'FLAG', 'NAME');
 syms('OUP', 'H', 'X', 'A', 'B', 'D', 'M', 'I', 'C', 'L', 'U');
 syms('Z', 'W', 'J', 's', 'x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Linear Finite-Difference Method.\n');
 fprintf(1,'Input the functions P(X), Q(X) and R(X) in terms of x, \n');
 fprintf(1,'on separate lines.\n');
 fprintf(1,'For example: -2/x \n');
 fprintf(1,'              2/(x^2) \n');
 fprintf(1,'              sin(log(x))/(x^2)\n');
 s = input(' ');
 P = inline(s,'x');
 s = input(' ');
 Q = inline(s,'x');
 s = input(' ');
 R = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 AA = input(' ');
 BB = input(' ');
 if AA >= BB 
 fprintf(1,'Left endpoint must be less than right endpoint.\n');
 else
 OK = TRUE;
 end;
 end;
 fprintf(1,'Input Y(  %.10e).\n', AA);
 ALPHA = input(' ');
 fprintf(1,'Input Y(  %.10e).\n', BB);
 BETA = input(' ');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input an integer > 1 for the number of\n');
 fprintf(1,'subintervals.  Note that h = (b-a)/(n+1)\n');
 N = input(' ');
 if N <= 1 
 fprintf(1,'Number must exceed 1.\n');
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
 fprintf(OUP, 'LINEAR FINITE DIFFERENCE METHOD\n\n');
 fprintf(OUP, '  I      X(I)           W(I)\n');
% STEP 1 */
 H = (BB-AA)/(N+1);
 A = zeros(1,N+1);
 B = zeros(1,N+1);
 C = zeros(1,N+1);
 D = zeros(1,N+1);
 L = zeros(1,N+1);
 U = zeros(1,N+1);
 Z = zeros(1,N+1);
 W = zeros(1,N+1);
 X = AA+H;
 A(1) = 2+H^2*Q(X);
 B(1) = -1+0.5*H*P(X);
 D(1) = -H^2*R(X)+(1+0.5*H*P(X))*ALPHA;
 M = N-1;
% STEP 2 */
 for I = 2 : M 
 X = AA+I*H;
 A(I) = 2+H^2*Q(X);
 B(I) = -1+0.5*H*P(X);
 C(I) = -1-0.5*H*P(X);
 D(I) = -H^2*R(X);
 end;
% STEP 3 */
 X = BB-H;
 A(N) = 2+H^2*Q(X);
 C(N) = -1-0.5*H*P(X);
 D(N) = -H^2*R(X)+(1-0.5*H*P(X))*BETA;
% STEP 4 */
% STEPS 4 through 8 solve a tridiagonal linear system using
% Crout factorization */
 L(1) = A(1);
 U(1) = B(1)/A(1);
 Z(1) = D(1)/L(1);
% STEP 5 */
 for I = 2 : M 
 L(I) = A(I)-C(I)*U(I-1);
 U(I) = B(I)/L(I);
 Z(I) = (D(I)-C(I)*Z(I-1))/L(I);
 end;
% STEP 6 */
 L(N) = A(N)-C(N)*U(N-1);
 Z(N) = (D(N)-C(N)*Z(N-1))/L(N);
% STEP 7 */
 W(N) = Z(N);
% STEP 8 */
 for J = 1 : M 
 I = N-J;
 W(I) = Z(I)-U(I)*W(I+1);
 end;
 I = 0;
% STEP 9 */
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, AA, ALPHA);
 for I = 1 : N 
 X = AA+I*H;
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, X, W(I));
 end;
 I = N+1;
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, BB, BETA);
% STEP 12 */
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
