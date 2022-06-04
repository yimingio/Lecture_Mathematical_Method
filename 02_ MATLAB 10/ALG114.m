% NONLINEAR FINITE-DIFFERENCE ALGORITHM 11.4
%
% To approximate the solution to the nonlinear boundary-value problem
%
%    Y'' = F(X,Y,Y'), A<=X<=B, Y(A) = ALPHA, Y(B) = BETA:
%
% INPUT:   Endpoints A,B; boundary conditions ALPHA, BETA;
%          integer N; tolerance TOL; maximum number of iterations M.
%
% OUTPUT:  Approximations W(I) TO Y(X(I)) for each I=0,1,...,N+1
%          or a message that the maximum number of iterations was
%          exceeded.
 syms('OK', 'AA', 'BB', 'ALPHA', 'BETA', 'N', 'TOL', 'NN');
 syms('FLAG', 'NAME', 'OUP', 'N1', 'H', 'I', 'W', 'K', 'X');
 syms('T', 'A', 'B', 'D', 'C', 'L', 'U', 'Z', 'V');
 syms('VMAX', 'J', 'x', 'y', 'z', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Nonlinear Finite-Difference Method.\n');
 fprintf(1,'Input the function F(X,Y,Z) in terms of x, y, z\n');
 fprintf(1,'followed by the partial of F with respect to y on \n');
 fprintf(1,'the next line and the partial of F with respect \n');
 fprintf(1,'to z = y-prime on the third line. \n');
 fprintf(1,'For example:   (32+2*x^3-y*z)/8 \n');
 fprintf(1,'                        -z/8    \n');
 fprintf(1,'                        -y/8    \n');
 s = input(' ');
 F = inline(s,'x','y','z');
 s = input(' ');
 FY = inline(s,'x','y','z');
 s = input(' ');
 FYP = inline(s,'x','y','z');
 OK = FALSE;
 while  OK == FALSE 
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
 OK = FALSE;
 while  OK == FALSE 
 fprintf(1,'Input Tolerance.\n');
 TOL = input(' ');
 if TOL <= 0 
 fprintf(1,'Tolerance must be positive.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while  OK == FALSE 
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
 fprintf(1,'for example   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'NONLINEAR FINITE-DIFFERENCE METHOD\n\n');
 fprintf(OUP, '  I    X(I)         W(I)\n');
% STEP 1
 A = zeros(1,N);
 B = zeros(1,N);
 C = zeros(1,N);
 D = zeros(1,N);
 W = zeros(1,N);
 V = zeros(1,N);
 Z = zeros(1,N);
 U = zeros(1,N);
 L = zeros(1,N);
 N1 = N-1;
 H = (BB-AA)/(N+1);
% STEP 2
 for I = 1 : N 
 W(I) = ALPHA+I*H*(BETA-ALPHA)/(BB-AA);
 end;
% STEP 3
 K = 1;
% STEP 4
 while K <= NN & OK == TRUE 
% STEP 5
 X = AA+H;
 T = (W(2)-ALPHA)/(2*H);
 A(1) = 2+H*H*FY(X,W(1),T);
 B(1) = -1+H*FYP(X,W(1),T)/2;
 D(1) = -(2*W(1)-W(2)-ALPHA+H*H*F(X,W(1),T));
% STEP 6
 for I = 2 : N1 
 X = AA+I*H;
 T = (W(I+1)-W(I-1))/(2*H);
 A(I) = 2+H*H*FY(X,W(I),T);
 B(I) = -1+H*FYP(X,W(I),T)/2;
 C(I) = -1-H*FYP(X,W(I),T)/2;
 D(I) = -(2*W(I)-W(I+1)-W(I-1)+H*H*F(X,W(I),T));
 end;
% STEP 7
 X = BB - H;
 T = (BETA-W(N-1))/(2*H);
 A(N) = 2+H*H*FY(X,W(N),T);
 C(N) = -1-H*FYP(X,W(N),T)/2;
 D(N) = -(2*W(N)-W(N-1)-BETA+H*H*F(X,W(N),T));
% STEP 8
% STEPS 8 through 12 solve a tridiagonal linear system using 
% Crout reduction
 L(1) = A(1);
 U(1) = B(1)/A(1);
 Z(1) = D(1)/L(1);
% STEP 9
 for I = 2 : N1 
 L(I) = A(I)-C(I)*U(I-1);
 U(I) = B(I)/L(I);
 Z(I) = (D(I)-C(I)*Z(I-1))/L(I);
 end;
% STEP 10
 L(N) = A(N)-C(N)*U(N-1);
 Z(N) = (D(N)-C(N)*Z(N-1))/L(N);
% STEP 11
 V(N) = Z(N);
 VMAX = abs(V(N));
 W(N) = W(N)+V(N);
% STEP 12
 for J = 1 : N1 
 I = N-J;
 V(I) = Z(I)-U(I)*V(I+1);
 W(I) = W(I)+V(I);
 if abs(V(I)) > VMAX 
 VMAX = abs(V(I));
 end;
 end;
% STEP 13
% test for accuracy
 if VMAX <= TOL 
 I = 0;
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, AA, ALPHA);
 for I = 1 : N 
 X = AA+I*H;
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, X, W(I));
 end;
 I = N+1;
 fprintf(OUP, '%3d %13.8f %13.8f\n', I, BB, BETA);
 OK = FALSE;
 else
% STEP 18
 K = K+1;
 end;
 end;
% STEP 19
 if K > NN 
 fprintf(OUP, 'No convergence in %d iterations\n', NN);
 end;
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
