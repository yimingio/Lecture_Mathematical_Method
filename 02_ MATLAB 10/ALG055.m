% ADAMS VARIABLE STEP-SIZE PREDICTOR-CORRECTOR ALGORITHM 5.5
%
% To approximate the solution of the initial value problem
%        y' = f( t, y ), a <= t <= b, y(a) = ALPHA,
%
% with local truncation error within a given tolerance:
%
% INPUT:   endpoints a, b; initial condition ALPHA; tolerance TOL;
%          maximum step size HMAX; minimum step size HMIN.
%
% OUTPUT:  I, T(I), W(I), H where at the Ith step W(I) approximates
%          y(T(I)) and step size H was used or a message that the
%          minimum step size was exceeded.
 syms('OK', 'A', 'B', 'ALPHA', 'TOL', 'HMIN', 'HMAX');
 syms('FLAG', 'NAME', 'OUP',  'DONE', 'KK', 'NFLAG');
 syms('I', 'WP', 'WC', 'SIG', 'K', 'J', 'Q');
 syms('W','T','t','y','s','TT','WW','K1','K2','K3','K4');
 syms('P1','P2');
 TRUE = 1;
 FALSE = 0;
% STEP 1  Runge-Kutta Order 4 Method is implemented within the
%         following code.
 fprintf(1,'This is the Adams Variable Step-size Predictor-');
 fprintf(1,'Corrector Method\n');
 fprintf(1,'Input the function F(t,y) in terms of t and y\n');
 fprintf(1,'For example: y-t^2+1 \n');
 s = input(' ');
 F = inline(s,'t','y');
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
 OK = FALSE;
 fprintf(1,'Input the initial condition.\n');
 ALPHA = input(' ');
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
 fprintf(1,'Input minimum and maximum mesh spacing  ');
 fprintf(1,'on separate lines.\n');
 HMIN = input(' ');
 HMAX = input(' ');
 if HMIN < HMAX & HMIN > 0 
 OK = TRUE;
 else
 fprintf(1,'Minimum mesh spacing must be a  ');
 fprintf(1,'positive real number and less than\n');
 fprintf(1,'the maximum mesh spacing.\n');
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
 fprintf(OUP, 'ADAMS VARIABLE STEP-SIZE PREDICTOR CORRECTOR METHOD\n\n');
 fprintf(OUP, '           t           w           h       sigma\n');
% STEP 2
 T = zeros(1,100);
 W = zeros(1,100);
 T(1) = A;
 W(1) = ALPHA;
 H = HMAX;
% OK is used in place of FLAG to exit the loop in Step 4.
 OK = TRUE;
% DONE is used in place of last to indicate when last value 
% is calculated
 DONE = FALSE;
% STEP 3
 for KK = 1:3
 X = T(KK);
 Y = W(KK);
 K1 = H*F(X,Y);
 K2 = H*F(X+0.5*H,Y+0.5*K1);
 K3 = H*F(X+0.5*H,Y+0.5*K2);
 K4 = H*F(X+H,Y+K3);
 WW = Y+(K1+2.0*(K2+K3)+K4)/6.0;
 TT = X+H;
 T(KK+1) = TT;
 W(KK+1) = WW;
 end;
% NFLAG indicates the computation from RK4
 NFLAG = 1;
 I = 5;
% use TT in place of t
 TT = T(4) + H;
% STEP 4
 while DONE == FALSE
% STEP 5
% predict W(I)
 P1 = 55.0*F(T(I-1),W(I-1))-59.0*F(T(I-2),W(I-2));
 P2 = 37.0*F(T(I-3),W(I-3))-9.0*F(T(I-4),W(I-4));
 WP = W(I-1)+H*(P1+P2)/24.0;
% correct W(I)
 P1 = 9.0*F(TT,WP)+19.0*F(T(I-1),W(I-1))-5.0*F(T(I-2),W(I-2));
 P2 = F(T(I-3),W(I-3));
 WC = W(I-1)+H*(P1+P2)/24.0;
 SIG = 19.0*abs(WC-WP)/(270.0*H);
% STEP 6
 if SIG <= TOL 
% STEP 7
% result accepted
 W(I) = WC;
 T(I) = TT;
% STEP 8
 if NFLAG == 1 
 K = I-3;
 KK = I-1;
% Previous results are also accepted.
 for J = K:KK
 fprintf(OUP, '%12.8f %11.8f %11.8f %11.8f\n', T(J), W(J), H, SIG);
 end;
 fprintf(OUP, '%12.8f %11.8f %11.8f %11.8f\n', T(I), W(I), H, SIG);
 else
% Previous results were already accepted.
 fprintf(OUP, '%12.8f %11.8f %11.8f %11.8f\n', T(I), W(I), H, SIG);
 end;
% STEP 9
 if OK == FALSE 
% Next step is 20.
 DONE = TRUE;
 else
% STEP 10
 I = I+1;
 NFLAG = 0;
% STEP 11
 if SIG <= 0.1*TOL | T(I-1)+H > B 
% Increase H if more accuracy than required has been obtained, 
% or decrease H to include b as a mesh point.
% STEP 12
% to avoid underflow
 if SIG <= 1.0e-20 
 Q = 4.0;
 else
 Q = (0.5*TOL/SIG)^(1/4);
 end;
% STEP 13
 if Q > 4.0 
 H = 4.0*H;
 else
 H = Q * H;
 end;
% STEP 14
 if H > HMAX 
 H = HMAX;
 end;
% STEP 15
 if T(I-1)+4.0*H > B 
 H = 0.25*(B-T(I-1));
 if H < TOL 
 DONE = TRUE;
 end;
 OK = FALSE;
 end;
% STEP 16
 for KK = I-1:I+2
 X = T(KK);
 Y = W(KK);
 K1 = H*F(X,Y);
 K2 = H*F(X+0.5*H,Y+0.5*K1);
 K3 = H*F(X+0.5*H,Y+0.5*K2);
 K4 = H*F(X+H,Y+K3);
 WW = Y+(K1+2.0*(K2+K3)+K4)/6.0;
 TT = X+H;
 T(KK+1) = TT;
 W(KK+1) = WW;
 end;
 NFLAG = 1;
 I = I+3;
 end;
 end;
 else
% FALSE branch for Step 6 - result rejected.
% STEP 17
 Q = (0.5*TOL/SIG)^(1/4);
% STEP 18
 if Q < 0.1 
 H = 0.1 * H;
 else
 H = Q * H;
 end;
% STEP 19
 if H < HMIN 
 fprintf(OUP, 'HMIN exceeded\n');
 DONE = TRUE;
 else
 if T(I-1)+4.0*H > B 
 H = 0.25*(B-T(I-1));
 end;
 if NFLAG == 1 
% Previous results also rejected.
 I = I-3;
 end;
 for KK = I-1:I+2
 X = T(KK);
 Y = W(KK);
 K1 = H*F(X,Y);
 K2 = H*F(X+0.5*H,Y+0.5*K1);
 K3 = H*F(X+0.5*H,Y+0.5*K2);
 K4 = H*F(X+H,Y+K3);
 WW = Y+(K1+2.0*(K2+K3)+K4)/6.0;
 TT = X+H;
 T(KK+1) = TT;
 W(KK+1) = WW;
 end;
 I = I+3;
 NFLAG = 1;
 end;
 end;
% STEP 20
 TT = T(I-1) + H;
 end;
% STEP 21
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
