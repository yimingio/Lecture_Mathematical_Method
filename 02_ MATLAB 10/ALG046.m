 % GAUSSIAN TRIPLE INTEGRAL ALGORITHM 4.6
 % 
 % To approximate I = triple integral ( ( f(x,y,z) dz dy dx ) ) with
 % limits of integration from a to b for x, from c(x) to d(x) for y, and
 % from alpha(x,y) to beta(x,y) for z.
 %
 % INPUT:   endpoints a, b; positive integers m, n, p. (Assume that the
 %          roots r(i,j) and coefficients c(i,j) are available for i
 %          equals m, n, and p and for 1 <= j <= i.
 syms('OK', 'A', 'B', 'M', 'N', 'P', 'r', 'co', 'H1', 'H2', 'AJ', 'I');
 syms('X', 'JX', 'C1', 'D1', 'K1', 'K2', 'J', 'Y', 'JY', 'Z1', 'Z2');
 syms('L1', 'L2', 'K', 'Z', 'Q','x','s','y','z');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Gaussian Quadrature for triple integrals.\n');  
 fprintf(1,'Input the function F(x,y,z) in terms of x, y, and z\n');
 fprintf(1,'For example: sqrt(x^2+y^2)*z\n');
 s = input(' ');
 F = inline(s,'x','y','z');
 fprintf(1,'Input the functions C(x), and D(x) in terms of x \n');
 fprintf(1,'on separate lines\n');
 fprintf(1,'For example: cos(x) \n');
 fprintf(1,'             sin(x) \n');
 s = input(' ');
 C = inline(s,'x');
 s = input(' ');
 D = inline(s,'x');
 fprintf(1,'Input the functions alpha(x,y) and beta(x,y) in \n');
 fprintf(1,'terms of x and y on separate lines\n'); 
 s = input(' ');
 alpha = inline(s,'x','y');
 s = input(' ');
 beta = inline(s,'x','y');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input lower limit of integration and upper limit of \n'); 
 fprintf(1,'integration on separate lines\n');
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
 fprintf(1,'Input three integers M > 1, N > 1, P > 1. This\n');
 fprintf(1,'implementation of Gaussian quadrature requires\n');
 fprintf(1,'that they are less than or equal to 5.\n');
 fprintf(1,'They will be used in ');  
 fprintf(1,'first, second, and third dimensions, resp.\n');
 fprintf(1,' on separate lines.\n');
 M = input(' ');
 N = input(' ');
 P = input(' ');
 if M <= 1 | N <= 1 | P <= 1  
 fprintf(1,'Integers must be greater than 1.\n');
 else
 if M > 5 | N > 5 | P > 5 
 fprintf(1,'Integers must be less than or equal to 5.\n');
 else
 OK = TRUE;
 end
 end
 end
 r = zeros(4,5);
 co = zeros(4,5);
 if OK == TRUE
    r(1,1) = 0.5773502692;
    r(1,2) = -r(1,1);
    co(1,1) = 1.0;
    co(1,2) = 1.0; 
    r(2,1) = 0.7745966692;
    r(2,2) = 0.0;
    r(2,3) = -r(2,1);
    co(2,1) = 0.5555555556;
    co(2,2) = 0.8888888889;
    co(2,3) = co(2,1);
    r(3,1) = 0.8611363116;
    r(3,2) = 0.3399810436;
    r(3,3) = -r(3,2);
    r(3,4) = -r(3,1);
    co(3,1) = 0.3478548451;
    co(3,2) = 0.6521451549;
    co(3,3) = co(3,2);
    co(3,4) = co(3,1);
    r(4,1) = 0.9061798459;
    r(4,2) = 0.5384693101; 
    r(4,3) = 0.0;
    r(4,4) = -r(4,2);
    r(4,5) = -r(4,1);
    co(4,1) = 0.2369268850;
    co(4,2) = 0.4786286705;
    co(4,3) = 0.5688888889;
    co(4,4) = co(4,2);
    co(4,5) = co(4,1);
% STEP 1
 H1 = (B-A)/2;
 H2 = (B+A)/2;
% use AJ instead of J
 AJ = 0;
% STEP 2
 for I = 1:M
% STEP 3
 X = H1*r(M-1,I)+H2;
 JX = 0;
 C1 = C(X);
 D1 = D(X);
 K1 = (D1-C1)/2;
 K2 = (D1+C1)/2;
% STEP 4
 for J = 1:N
% STEP 5
 Y = K1*r(N-1,J)+K2;
 JY = 0;
% use Z1 for Beta and Z2 for Alpha
 Z1 = beta(X,Y);
 Z2 = alpha(X,Y);
 L1 = (Z1-Z2)/2;
 L2 = (Z1+Z2)/2;
% STEP 6
 for K = 1:P
 Z = L1*r(P-1,K)+L2;
 Q = F(X,Y,Z);
 JY = JY+co(P-1,K)*Q;
 end
% STEP 7
 JX = JX+co(N-1,J)*L1*JY;
 end
% STEP 8
 AJ = AJ+co(M-1,I)*K1*JX;
 end
% STEP 9
 AJ = AJ*H1;
% STEP 10
 fprintf(1,'\nThe triple integral of F from %12.8f to %12.8f is   ', A, B);
 fprintf(1,'%.10e\n', AJ);
 fprintf(1,' obtained with M = %3d , N = %3d and P = %3d\n', M, N, P);
 end
