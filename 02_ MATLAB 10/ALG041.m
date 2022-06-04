 % SIMPSON'S COMPOSITE ALGORITHM 4.1
 %
 % To approximate I = integral ( ( f(x) dx ) ) from a to b:
 %
 % INPUT:   endpoints a, b; even positive integer n.
 %
 % OUTPUT:  approximation XI to I.
 syms('OK','A','B','N','H','XI0','XI1','XI2','NN','I','X','XI','s','x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Simpsons Method.\n\n');
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
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
 fprintf(1,'Input an even positive integer N.\n');
 N = input(' ');
 if N > 0 & rem(N,2) == 0 
 OK = TRUE;
 else
 fprintf(1,'Input must be even and positive\n');
 end
 end
 if OK == TRUE 
% STEP 1
 H = (B-A)/N;
% STEP 2
 XI0 = F(A) + F(B);
% summation of f(x(2*I-1))
 XI1 = 0.0;
% summation of f(x(2*I))
 XI2 = 0.0;
% STEP 3
 NN = N - 1;
 for I = 1:NN
% STEP 4
 X = A + I * H;
% STEP 5
 if rem(I,2) == 0  
 XI2 = XI2 + F(X);
 else
 XI1 = XI1 + F(X);      
 end
 end
% STEP 6
 XI = (XI0 + 2.0 * XI2 + 4.0 * XI1) * H / 3.0;
% STEP 7
 fprintf(1,'\nThe integral of F from %12.8f to %12.8f is\n', A, B);
 fprintf(1,'%12.8f\n', XI);
 end
