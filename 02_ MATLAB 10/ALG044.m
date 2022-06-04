 % DOUBLE INTEGAL ALGORITHM 4.4
 %
 % To approximate I = double integral ( ( f(x,y) dy dx ) ) with limits
 % of integration from a to b for x and from c(x) to d(x) for y:
 %
 % INPUT:    endpoints a, b; even positive integers m, n.
 %
 % OUTPUT:   approximation J to I. 
 syms('OK', 'A', 'B', 'N', 'M', 'NN', 'MM', 'H', 'AN', 'AE', 'AO');
 syms('I', 'X', 'YA', 'YB', 'HX', 'BN', 'BE', 'BO', 'J', 'Y', 'Z');
 syms('A1', 'AC','x','s','y');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Simpsons Method for double integrals.\n\n');
 fprintf(1,'Input the functions F(X,Y), C(X), and D(X) in terms of x\n');
 fprintf(1,'and y on separate lines.\n');
 fprintf(1,'For example: cos(x+y)\n');
 fprintf(1,'             x^3     \n');
 fprintf(1,'             x       \n');
 s = input(' ');
 F = inline(s,'x','y');
 s = input(' ');
 C = inline(s,'x');
 s = input(' ');
 D = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input lower limit of integration and ');
 fprintf(1,'upper limit of integration\n');
 fprintf(1,'on separate lines\n');
 A = input(' ');
 B = input(' ');
 if A > B 
 fprintf(1,'Lower limit must be less than upper limit\n');
 else
 OK = TRUE;
 end
 end 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input two even positive integer N and M.\n');
 fprintf(1,'N corresponds to the outer integral and M \n');
 fprintf(1,'corresponds to the inner integral.  Place \n');
 fprintf(1,'on separate lines.\n');
 N = input(' ');
 M = input(' ');
 if N > 0 & rem(N,2) == 0 & M > 0 & rem(M,2) == 0 
 OK = TRUE;
 else
 fprintf(1,'N and M must both be even and positive\n');
 end
 end
 if OK == TRUE 
% STEP 1
 H = (B-A)/N;
% use AN, AE, AO, for J(1), J(2), J(3) resp.
% end terms
 AN = 0;
% even terms
 AE = 0;
% odd terms
 AO = 0;
% STEP 2
 for I = 0:N
% STEP 3
% Composite Simpson's Method for X
 X = A+I*H;
 YA = C(X);
 YB = D(X);
 HX = (YB-YA)/(M);
% use BN, BE, BO for K(1), K(2), K(3) resp.
% end terms
 BN = F(X,YA)+F(X,YB);
% even terms
 BE = 0;
% odd terms
 BO = 0;
% STEP 4
 for J = 1:M-1
% STEP 5
 Y = YA+J*HX;
 Z = F(X, Y);
% STEP 6
 if rem(J,2) == 0 
 BE = BE+Z;
 else
 BO = BO+Z;
 end
 end
% STEP 7
% use A1 for L, which is the integral of F(X(I), Y) from
% C(X(I)) to D(X(I)) by Composite Simpson's Method
 A1 = (BN+2*BE+4*BO)*HX/3;
% STEP 8
 if I == 0 | I == N 
 AN = AN+A1;
 else
 if rem(I,2) == 0 
 AE = AE + A1;
 else
 AO = AO + A1;
 end
 end
 end
% STEP 9
% Use AC for J
 AC = (AN + 2 * AE + 4 * AO) * H /3;
% STEP 10
 fprintf(1,'\nThe double integral of F from %12.8f to %12.8f is\n', A, B);
 fprintf(1,'%12.8f', AC);
 fprintf(1,' obtained with N = %3d and M = %3d\n', N, M);
 end
