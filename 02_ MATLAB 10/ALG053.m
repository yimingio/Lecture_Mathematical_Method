 % RUNGE-KUTTA-FEHLBERG ALGORITHM 5.3
 %
 % TO APPROXIMATE THE SOLUTION OF THE INITIAL VALUE PROBLEM:
 %            Y' = F(T,Y), A<=T<=B, Y(A) = ALPHA,
 % WITH LOCAL TRUNCATION ERROR WITHIN A GIVEN TOLERANCE.
 %
 % INPUT:   ENDPOINTS A,B; INITIAL CONDITION ALPHA; TOLERANCE TOL;
 %          MAXIMUM STEPSIZE HMAX; MINIMUM STEPSIZE HMIN.
 %
 % OUTPUT:  T, W, H WHERE W APPROXIMATES Y(T) AND STEPSIZE H WAS 
 % USED OR A MESSAGE THAT THE MINIMUM STEPSIZE WAS EXCEEDED. 
 syms('F', 'OK', 'A', 'B', 'ALPHA', 'TOL', 'HMIN', 'HMAX', 'FLAG');
 syms('NAME', 'OUP', 'H', 'T', 'W', 'K1', 'K2', 'K3', 'K4', 'K5', 'K6');
 syms('R', 'DELTA', 't', 's');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Runge-Kutta-Fehlberg Method.\n');
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
 fprintf(1,'Input tolerance\n');
 TOL = input(' ');
 if TOL <= 0 
 fprintf(1,'Tolerance must be a positive.\n');
 else
 OK = TRUE;
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input minimum and maximum mesh spacing on separate lines.\n');
 HMIN = input(' ');
 HMAX = input(' ');
 if HMIN < HMAX & HMIN > 0 
 OK = TRUE;
 else
 fprintf(1,'Minimum mesh spacing must be a positive real\n');
 fprintf(1,'number and less than the maximum mesh spacing\n');
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
 fprintf(OUP, 'RUNGE-KUTTA-FEHLBERG METHOD\n\n');
 fprintf(OUP, '    T(I)           W(I)          H           R\n\n');
% STEP 1
 H = HMAX;
 T = A;
 W = ALPHA;
 fprintf(OUP, '%12.7f %11.7f    0           0\n', T, W);
 OK = TRUE;
% STEP 2
 while T < B & OK == TRUE 
% STEP 3
 K1 = H*F(T,W);
 K2 = H*F(T+H/4,W+K1/4);K3 = H*F(T+3*H/8,W+(3*K1+9*K2)/32);
 K4 = H*F(T+12*H/13,W+(1932*K1-7200*K2+7296*K3)/2197);
 K5 = H*F(T+H,W+439*K1/216-8*K2+3680*K3/513-845*K4/4104);
 K6 = H*F(T+H/2,W-8*K1/27+2*K2-3544*K3/2565+1859*K4/4104-11*K5/40);
% STEP 4
 R = abs(K1/360-128*K3/4275-2197*K4/75240.0+K5/50+2*K6/55)/H;
% STEP 5
 if R <= TOL 
% STEP 6
% Approximation accepted
 T = T+H;
 W = W+25*K1/216+1408*K3/2565+2197*K4/4104-K5/5;
% STEP 7
 fprintf(OUP, '%12.7f %11.7f %11.7f %11.7f\n', T, W, H, R);
 end;
% STEP 8
% To avoid underflow
 if R > 1.0E-20 
 DELTA = 0.84 * exp(0.25 * log(TOL / R));
 else
 DELTA = 10.0;
 end;
% STEP 9
% Calculate new H
 if DELTA <= 0.1 
 H = 0.1*H;
 else
 if DELTA >= 4 
 H = 4.0 * H;
 else
 H = DELTA * H;
 end;
 end;
% STEP 10
 if H > HMAX 
 H = HMAX;
 end;
% STEP 11
 if H < HMIN 
 OK = FALSE;
 else
 if T+H > B 
 if abs(B-T) < TOL 
 T = B;
 else
 H = B-T;
 end;
 end;
 end;
 end;
 if OK == FALSE 
 fprintf(OUP, 'Minimal H exceeded\n');
 end;
% STEP 12
% Process is complete
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
