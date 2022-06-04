% TRAPEZOIDAL WITH NEWTON ITERATION ALGORITHM 5.8
%
% TO APPROXIMATE THE SOLUTION OF THE INITIAL VALUE PROBLEM:
%            Y' = F(T,Y), A <= T <= B, Y(A) = ALPHA,
% AT (N+1) EQUALLY SPACED NUMBERS IN THE INTERVAL [A,B].
%
% INPUT:   ENDPOINTS A,B; INITIAL CONDITION ALPHA; INTEGER N;
%          TOLERANCE TOL; MAXIMUM NUMBER OF ITERATIONS M AT ANY ONE STEP.
%
% OUTPUT:  APPROXIMATION W TO Y AT THE (N+1) VALUES OF T
%          OR A MESSAGE OF FAILURE.
 syms('F', 'FYP', 'OK', 'A', 'B', 'ALPHA', 'N', 'TOL', 'M');
 syms('FLAG', 'NAME', 'OUP', 'W', 'T', 'H', 'I', 'XK1', 'W0');
 syms('J', 'IFLAG','y','t');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Implicit Trapezoidal Method.\n');
 fprintf(1,'Input the function F(t,y) in terms of t and y\n');
 fprintf(1,'For example: y^2-y*t^2+1 \n');
 s = input(' ');
 F = inline(s,'t','y');
 fprintf(1,'Input the partial derivative of F(t,y) with respect to y \n');
 fprintf(1,'in terms of t and y.\n');
 fprintf(1,'for example:    2*y-t^2 \n');
 s = input(' ');
 FYP = inline(s,'t','y');
 OK = FALSE;
 while  OK == FALSE 
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 A = input(' ');
 B = input(' ');
 if A >= B 
 fprintf(1,'Left endpoint must be less than right endpoint.\n');
 else
 OK = TRUE;
 end;
 end;
 fprintf(1,'Input the initial condition.\n');
 ALPHA = input(' ');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input a positive integer for the number of subintervals.\n');
 N = input(' ');
 if N <= 0 
 fprintf(1,'Number must be a postiive integer.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input tolerance.\n');
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
 M = input(' ');
 if M > 0 
 OK = TRUE;
 else
 fprintf(1,'Number of iterations must be positive.\n');
 end;
 end;
 if OK == TRUE 
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
 fprintf(1,'Please enter 1 or 2\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'For example   A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'IMPLICIT TRAPEZOIDAL METHOD USING NEWTONS METHOD\n\n');
 fprintf(OUP, '    t           w #iter\n');
% STEP 1
 W = ALPHA;
 T = A;
 H = (B-A)/N;
 fprintf(OUP, '%5.3f %11.8f   0\n', T, W);
 I = 1;
 OK = TRUE;
% STEP 2
 while I <= N & OK == TRUE 
% STEP 3
 XK1 = W+0.5*H*F(T, W);
 W0 = XK1;
 J = 1;
 IFLAG = 0;
% STEP 4
 while IFLAG == 0 & OK == TRUE 
% STEP 5
 W = W0-(W0-XK1-0.5*H*F(T+H, W0))/(1-0.5*H*FYP(T+H, W0));
% STEP 6
 if abs(W-W0) < TOL 
 IFLAG = 1;
% STEP 7
 T = A+I*H;
 fprintf(OUP,'%5.3f %11.8f %3d\n', T, W, J);
 I = I+1;
 else
 J = J+1;
 W0 = W;
 if J > M 
 OK = FALSE;
 end;
 end;
 end;
 end;
 if OK == FALSE 
 fprintf(OUP, 'Maximum Number of Iterations Exceeded\n');
 end;
% STEP 8
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
