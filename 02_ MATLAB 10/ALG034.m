 % NATURAL CUBIC SPLINE ALGORITHM 3.4                                     %
 % To construct the cubic spline interpolant S for the function f,
 % defined at the numbers x(0) < x(1) < ... < x(n), satisfying
 % S''(x(0)) = S''(x(n)) = 0:
 %
 % INPUT:   n; x(0), x(1), ..., x(n); either generate A(I) = f(x(I))
 %          for I = 0, 1, ..., n or input A(I) for I = 0, 1, ..., n.
 %
 % OUTPUT:  A(J), B(J), C(J), D(J) for J = 0, 1, ..., n - 1.
 %
 % NOTE:    S(x) = A(J) + B(J)*( x - x(J) ) + C(J)*( x - x(J) )**2 +
 %          D(J) * ( x - x(J) )**3 for x(J) <= x < x(J + 1) 
 syms('OK','FLAG','N','I','X','A','AA','NAME','INP','F','M','H','XA')
 syms('XL','XU','XZ','C','J','B','D','OUP','x','s');
 TRUE=1;
 FALSE=0;
 fprintf(1,'This is the natural cubic spline interpolation.\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Choice of input method:\n');
 fprintf(1,'1. Input entry by entry from keyboard\n');
 fprintf(1,'2. Input data from a text file\n');
 fprintf(1,'3. Generate data using a function F with nodes entered ');
 fprintf(1,'from keyboard\n');
 fprintf(1,'4. Generate data using a function F with nodes from ');
 fprintf(1,'a text file\n');
 fprintf(1,'Choose 1, 2, 3, or 4 please\n');
 FLAG = input(' ');
 if FLAG >= 1 & FLAG <= 4 
 OK = TRUE;
 end
 end
 if FLAG == 1 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0 
 OK = TRUE;
 X = zeros(1,N+1);
 A = zeros(1,N+1);
 for I = 0:N 
 fprintf(1,'Input X(%d) and F(X(%d)) ', I, I);
 fprintf(1,'on separate lines.\n');
 X(I+1) = input(' ');
 A(I+1) = input(' ');
 end
 else fprintf(1,'Number must be a positive integer\n');
 end
 end
 end
 if FLAG == 2 
 fprintf(1,'Has a text file been created with the data in two ');
 fprintf(1,'columns ?\n');
 fprintf(1,'Enter Y or N\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'For example: A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0
 X = zeros(1,N+1);
 A = zeros(1,N+1);
 for I = 0:N
 X(I+1) = fscanf(INP, '%f',1);
 A(I+1) = fscanf(INP, '%f',1);
 end
 fclose(INP);
 OK = TRUE;
 else 
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 else
 fprintf(1,'Please create the input file in two column ');
 fprintf(1,'form with the\n');
 fprintf(1,'X values and F(X) values in the corresponding columns.\n');
 fprintf(1,'The program will end so the input file can be created.\n');
 OK = FALSE;
 end
 end
 if FLAG == 3
 fprintf(1,'Input the function F(x) in terms of x.\n');
 fprintf(1,'For example: cos(x) \n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0
 X = zeros(1,N+1);
 A = zeros(1,N+1);
 for I = 0:N 
 fprintf(1,'Input X(%d)\n', I);
 X(I+1) = input(' ');
 A(I+1) = F(X(I+1));
 end
 OK = TRUE;
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 end
 if FLAG == 4 
 fprintf(1,'Has the text file with X-values been created?\n');
 fprintf(1,'Enter Y or N\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'For example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 fprintf(1,'Input the function F(x) in terms of x.\n');
 fprintf(1,'For example: cos(x) \n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0
 OK = TRUE;
 X = zeros(1,N+1);
 A = zeros(1,N+1);
 for I = 0:N
 X(I+1) = fscanf(INP, '%f',1);
 A(I+1) = F(X(I+1));
 end
 fclose(INP);
 else fprintf(1,'Number must be a positive integer\n');
 end
 end
 else 
 fprintf(1,'The program will end so the input file can be created\n');
 OK = FALSE;
 end
 end
 if OK == TRUE 
 M = N - 1;
% STEP 1
 H = zeros(1,M+1);
 for I = 0:M
 H(I+1) = X(I+2) - X(I+1);
 end
% STEP 2
% Use XA in place of ALPHA
 XA = zeros(1,M+1);
 for I = 1:M
 XA(I+1) = 3.0*(A(I+2)*H(I)-A(I+1)*(X(I+2)-X(I))+A(I)*H(I+1))/(H(I+1)*H(I));
 end
% STEP 3
% STEPs 3, 4, 5 and part of 6 solve the tridiagonal system using 
% Crout reduction.
% use XL, XU, XZ in place of L, MU, Z resp.
 XL = zeros(1,N+1);
 XU = zeros(1,N+1);
 XZ = zeros(1,N+1);
 XL(1) = 1;
 XU(1) = 0;
 XZ(1) = 0;
% STEP 4
 for I = 1:M
 XL(I+1) = 2*(X(I+2)-X(I))-H(I)*XU(I);
 XU(I+1) = H(I+1)/XL(I+1);
 XZ(I+1) = (XA(I+1)-H(I)*XZ(I))/XL(I+1);
 end
% STEP 5
 XL(N+1) = 1;
 XZ(N+1) = 0;
 B = zeros(1,N+1);
 C = zeros(1,N+1);
 D = zeros(1,N+1);
 C(N+1) = XZ(N+1);
% STEP 6
 for I = 0:M
 J = M-I;
 C(J+1) = XZ(J+1)-XU(J+1)*C(J+2);
 B(J+1) = (A(J+2)-A(J+1))/H(J+1) - H(J+1) * (C(J+2) + 2.0 * C(J+1)) / 3.0;
 D(J+1) = (C(J+2) - C(J+1)) / (3.0 * H(J+1));
 end
% STEP 7
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
 fprintf(OUP, 'NATURAL CUBIC SPLINE INTERPOLATION\n\n');
 fprintf(OUP, 'The numbers X(0), ..., X(N) are:\n');
 for I = 0:N
 fprintf(OUP, ' %11.8f', X(I+1));
 end
 fprintf(OUP, '\n\nThe coefficients of the spline on the subintervals '); 
 fprintf(OUP, 'are:\n');
 fprintf(OUP, 'for I = 0, ..., N-1\n');
 fprintf(OUP, '     A(I)          B(I)           C(I)         D(I)\n');
 for I = 0:M
 fprintf(OUP,'%13.8f %13.8f %13.8f %13.8f \n',A(I+1),B(I+1),C(I+1),D(I+1));
 end
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end
 end
   
