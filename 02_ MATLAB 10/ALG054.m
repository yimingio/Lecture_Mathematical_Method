 % ADAMS-FOURTH ORDER PREDICTOR-CORRECTOR ALGORITHM 5.4
 %
 % To approximate the solution of the initial value problem
 %        y' = f(t,y), a <= t <= b, y(a) = alpha,
 % at N+1 equally spaced points in the interval [a,b].
 %
 % INPUT:   endpoints a,b; initial condition alpha; integer N.
 %
 % OUTPUT:  approximation w to y at the (N+1) values of t.
 syms('F', 'OK', 'A', 'B', 'ALPHA', 'N', 'FLAG', 'NAME', 'OUP');
 syms('H', 'T', 'W', 'I', 'K1', 'K2', 'K3', 'K4', 'T0', 'W0', 'J');
 syms('t','y', 's','Part1','Part2');
 TRUE = 1;
 FALSE = 0;
 T = zeros(1,4);
 W = zeros(1,4);
 fprintf(1,'This is Adams-Bashforth Predictor Corrector Method\n');
 fprintf(1,'Input the function F(t,y) in terms of t and y\n');
 fprintf(1,'For example: y-t^2+1 \n');
 s = input(' ');
 F = inline(s,'t','y');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 A = input(' ');
 B = input(' ');
 if A >= B  
 fprintf(1,'Left endpoint must be less than right endpoint\n');
 else
 OK = TRUE;
 end;
 end;
 fprintf(1,'Input the initial condition\n');
 ALPHA = input(' ');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input an integer > 3 for the number of subintervals\n');
 N = input(' ');
 if N <= 3 
 fprintf(1,'Number must be at least 4.\n');
 else
 OK = TRUE;
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
 fprintf(OUP, 'ADAMS-BASHFORTH FOURTH ORDER PREDICTOR CORRECTOR METHOD\n\n');
 fprintf(OUP, '    t           w\n');
% STEP 1
 H = (B-A)/N;
 T(1) = A;
 W(1) = ALPHA;
 fprintf(OUP, '%5.3f %11.7f\n', T(1), W(1));
% STEP 2
 for I = 1:3 
% STEP 3 AND 4
% compute starting values using Runge-Kutta method
 T(I+1) = T(I)+H;
 K1 = H*F(T(I), W(I));
 K2 = H*F(T(I)+0.5*H, W(I)+0.5*K1);
 K3 = H*F(T(I)+0.5*H, W(I)+0.5*K2);
 K4 = H*F(T(I+1), W(I)+K3);
 W(I+1) = W(I)+(K1+2.0*(K2+K3)+K4)/6.0;
% STEP 5
 fprintf(OUP, '%5.3f %11.7f\n', T(I+1), W(I+1));
 end;
% STEP 6
 for I = 4:N 
% STEP 7
% T0, W0 will be used in place of t, w resp.
 T0 = A+I*H;
% predict W(I)
 Part1 = 55.0*F(T(4),W(4))-59.0*F(T(3),W(3))+37.0*F(T(2),W(2));
 Part2 = -9.0*F(T(1),W(1));
 W0 = W(4)+H*(Part1+Part2)/24.0;
% correct W(I)
 Part1 = 9.0*F(T0,W0)+19.0*F(T(4),W(4))-5.0*F(T(3),W(3))+F(T(2),W(2));
 W0 = W(4)+H*(Part1)/24.0;
% STEP 8
 fprintf(OUP, '%5.3f %11.7f\n', T0, W0);
% STEP 9
% prepare for next iteration
 for J = 1:3 
 T(J) = T(J+1);
 W(J) = W(J+1);
 end;
% STEP 10
 T(4) = T0;
 W(4) = W0;
 end;
 end;
% STEP 11
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
