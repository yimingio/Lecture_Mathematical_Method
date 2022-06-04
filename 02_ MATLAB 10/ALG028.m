 % MULLER'S ALGORITHM 2.8
 %
 % To find a solution to f(x) = 0 given three approximations x0, x1
 % and x2:
 %
 % INPUT:  x0,x1,x2; tolerance TOL; maximum number of iterations NO.
 %
 % OUTPUT: approximate solution p or message of failure.
 %
 % This implementation allows for a switch to complex arithmetic.
 % The coefficients are stored in the vector A, so the dimension
 % of A may have to be changed.
 syms('P', 'OK', 'TOL', 'M', 'X', 'FLAG', 'NAME', 'OUP', 'F', 'H');
 syms('r','DEL1', 'DEL', 'I', 'B', 'D', 'E', 'J','x','s','N');
 TRUE = 1;
 FALSE = 0;
 F = zeros(1,4);
 X = zeros(1,4);
 H = zeros(1,3);
 DEL1 = zeros(1,2);
 fprintf(1,'This is Mullers Method.\n');
 fprintf(1,'Input the Polynomial P(x)\n');
 fprintf(1,'For example: to input x^3-2*x+4 enter \n');
 fprintf(1,' [ 1 0 -2 4 ] \n'); 
 P = input(' ');
 OK = TRUE;
 N = length(P);
 if N == 2 
 r = -P(N)/P(N-1);
 fprintf(1,'Polynomial is linear: root is %11.8f\n', r);
 OK = FALSE;
 end
 if OK == TRUE 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input tolerance\n');
 TOL = input(' ');
 if TOL <= 0 
 fprintf(1,'Tolerance must be positive\n');
 else 
 OK = TRUE;
 end
 end
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input maximum number of iterations - no decimal point\n');
 M = input(' ');
 if M <= 0 
 fprintf(1,'Must be positive integer\n');
 else 
 OK = TRUE;
 end
 end
 fprintf(1,'Input the first of three starting values\n');
 X(1) = input(' ');
 fprintf(1,'Input the second of three starting values\n');
 X(2) = input(' ');
 fprintf(1,'Input the third starting value\n');
 X(3) = input(' ');
 end
 if OK == TRUE 
 fprintf(1,'Select output destination\n');
 fprintf(1,'1. Screen\n');
 fprintf(1,'2. Text file\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'For example:   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end
 fprintf(OUP, 'MULLERS METHOD\n');
 fprintf(OUP, 'The output is i, approximation x(i), f(x(i))\n\n');
 fprintf(OUP,'The real and imaginary parts of x(i) are\n');
 fprintf(OUP,'followed by real and imaginary parts of f(x(i)).\n\n');
 F(1) = polyval(P,X(1));
 F(2) = polyval(P,X(2));
 F(3) = polyval(P,X(3));
% STEP 1
 H(1) = X(2)-X(1);
 H(2) = X(3)-X(2);
 DEL1(1) = (F(2)-F(1))/H(1);
 DEL1(2) = (F(3)-F(2))/H(2);
 DEL = (DEL1(2)-DEL1(1))/(H(2)+H(1));
 I = 3;
% STEP 2
 while I <= M & OK == TRUE 
% STEP 3
 B = DEL1(2)+H(2)*DEL;
 D = B*B-4*F(3)*DEL;
% test to see if straight line
 if abs(DEL) <= 1.0e-20 
% straight line - test if horizontal line
 if abs(DEL1(2)) <= 1.0e-20 
 fprintf(1,'Horizontal Line\n');
 OK = FALSE;
 else
% straight line but not horizontal
 X(4) = (F(3)-DEL1(2)*X(3))/DEL1(2);
 H(3) = X(4)-X(3);
 end
 else
% not a straight line
 D = sqrt(D);
% STEP 4
 E = B+D;
 if abs(B-D) > abs(E) 
 E = B-D;
 end
% STEP 5
 H(3) = -2*F(3)/E;
 X(4) = X(3)+H(3);
 end
 if OK == TRUE 
 F(4) = polyval(P,X(4));
 fprintf(OUP, '%d %f %f %f %f\n',I,X(4),imag(X(4)),F(4),imag(F(4)));
 end
% STEP 6
 if abs(H(3)) < TOL 
% Procedure completed successfully. 
 fprintf(OUP, '\nMethod Succeeds\n');
 fprintf(OUP, 'Approximation is within %.10e\n', TOL);
 fprintf(OUP, 'in %d iterations\n', I);
 OK = FALSE;
 else
% STEP 7
 for J = 1:2
 H(J) = H(J+1);
 X(J) = X(J+1);
 F(J) = F(J+1);
 end
 X(3) = X(4);
 F(3) = F(4);
 DEL1(1) = DEL1(2);
 DEL1(2) = (F(3)-F(2))/H(2);
 DEL = (DEL1(2)-DEL1(1))/(H(2)+H(1));
 end
 I = I+1;
 end
% STEP 8
 if I > M & OK == TRUE 
% Procedure completed unsuccessfully.
 fprintf (OUP, 'Method Failed\n');
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created sucessfully\n',NAME);
 end
 end
