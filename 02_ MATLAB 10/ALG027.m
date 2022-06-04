 %  HORNER'S ALGORITHM 2.7
 %
 %  To evaluate the polynomial
 %  p(x) = a(n) * x^n + a(n-1) * x^(n-1) + ... + a(1) * x + a(0)
 %  and its derivative p'(x) at x = x0;
 %
 %  INPUT:   degree n; coefficients aa(0),aa(1),...,aa(n);
 %           value of x0.
 %
 %  OUTPUT:  y = p(x0), z = p'(x0).
 syms('OK', 'N', 'I', 'AA', 'X0', 'Y', 'Z', 'MM', 'J');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Horners Method\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input degree n of polynomial - no decimal point\n');
 N = input(' ');
 if N < 0 
 fprintf(1,'Incorrect input - degree must be positive.\n');
 else
 OK = TRUE;
 end
 end
 AA = zeros(1,N+1);
 fprintf(1,'Input coefficients of P(X) in ascending order\n');
 for I = 0:N 
 fprintf(1,'Input coefficient of X^%d\n', I);
 AA(I+1) = input(' ');
 end
 fprintf(1,'Input argument X0 at which to evaluate P(X)\n');
 X0 = input(' ');
% STEP 1
% compute b(n) for p(x)
 Y = AA(N+1);
% compute b(n-1) for q(x) = p'(x)
 if N == 0 
 Z = 0;
 else
 Z = AA(N+1);
 end
 MM = N-1;
% STEP 2
 for I = 1:MM
 J = N-I;
% compute b(j) for p(x)
 Y = Y*X0+AA(J+1);
% compute b(j-1) for q(x)
 Z = Z*X0+Y;
 end
% STEP 3
% compute b(0) for p(x)
 if N ~= 0 
 Y = Y*X0+AA(1);
 end
 fprintf(1,'Coefficients of polynominal P :\n');
% STEP 4
 for I = 0:N
 fprintf(1,'Exponent = %3d Coefficient = %12.8f\n', I, AA(I+1));
 end
 fprintf(1,'\n P ( %.10e ) = %12.8f\n', X0, Y);
 fprintf(1,' P` ( %.10e ) = %12.8f\n', X0, Z);
