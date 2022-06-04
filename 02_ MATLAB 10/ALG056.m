% EXTRAPOLATION ALGORITHM 5.6
%
% To approximate the solution of the initial value problem:
%               y' = f(t,y), a <= t <= b, y(a) = ALPHA,
% with local truncation error within a given tolerance:
%
% INPUT:   endpoints a,b; initial condition ALPHA; tolerance TOL;
%          maximum stepsize HMAX; minimum stepsize HMIN.
%
% OUTPUT:  T, W, H where W approximates y(T) and stepsize H was
%          used or a message that minimum stepsize was exceeded.
 syms('F', 'OK', 'A', 'B', 'ALPHA', 'TOL', 'HMIN', 'HMAX');
 syms('FLAG', 'NAME', 'OUP', 'NK', 'J', 'I', 'T0', 'W0');
 syms('H', 'DONE', 'Q', 'K', 'NFLAG', 'HK', 'T', 'W2');
 syms('W3', 'M', 'W1', 'Y', 'V','t','y');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Gragg Extrapolation\n');
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
 fprintf(OUP, 'GRAGG EXTRAPOLATION\n\n');
 fprintf(OUP, '           T           W           H      K\n');
% STEP 1 */
 NK = zeros(1,8);
 NK(1) = 2;
 NK(2) = 4;
 for J = 1:3
 I = 2*J;
 NK(I+1) = 3*NK(I)/2;
 NK(I+2) = 2*NK(I);
 end;
% STEP 2 */
 T0 = A;
 W0 = ALPHA;
 H = HMAX;
%  DONE is used in place of FLAG to exit the loop in Step 4  */
 DONE = FALSE;
% STEP 3 */
 Q = zeros(8,8);
 for I = 1:7
 for J = 1:I
 Q(I,J) = (NK(I+1)*1/NK(J))*(NK(I+1)*1/NK(J));
 end;
 end;
% STEP 4 */
 while DONE == FALSE 
% STEP 5 */
 K = 1;
% when desired accuracy achieved, NFLAG is set to 1 */
 NFLAG = 0;
 Y = zeros(1,10);
% STEP 6 */
 while K <= 8 & NFLAG == 0 
% STEP 7 */
 HK = H/NK(K);
 T = T0;
 W2 = W0;
% Euler first step *
 W3 = W2+HK*F(T, W2);
 T = T0+HK;
% STEP 8 */
 M = NK(K)-1;
 for J = 1:M
 W1 = W2;
 W2 = W3;
% midpoint method */
 W3 = W1+2*HK*F(T, W2);
 T = T0+(J+1)*HK;
 end;
% STEP 9 */
% endpoint correction to compute Y(K,1) */
 Y(K) = (W3+W2+HK*F(T, W3))/2;
% STEP 10 */
% NOTE: Y(K-1)=Y(K-1,1),Y(K-2)=Y(K-2,2),..., */
% Y(1)=Y(K-1,K-1) since only previous row of table */
% is saved  */
 if K >= 2 
% STEP 11 */
 J = K;
%  save Y(K-1,K-1)  */
 V = Y(1);
% STEP 12 */
 while J >= 2 
% extrapolation to compute */
% Y(J-1) = Y(K,K-J+2) */
 Y(J-1) = Y(J)+(Y(J)-Y(J-1))/(Q(K-1,J-1)-1);
 J = J-1;
 end;
% STEP 13 */
 if abs(Y(1) - V) <= TOL 
 NFLAG = 1;
 end;
% Y(1) accepted as new w */
 end;
% STEP 14 */
 K = K+1;
 end;
% STEP 15 */
 K = K-1;
% STEP 16 */
 if NFLAG == 0 
% STEP 17 */
% new value for w rejected, decrease H */
 H = H/2;
% STEP 18 */
 if H < HMIN 
 fprintf(OUP, 'HMIN exceeded\n');
 DONE = TRUE;
 end;
 else
% STEP 19 */
% new value for w accepted */
 W0 = Y(1);
 T0 = T0 + H;
 fprintf(OUP, '%12.8f %11.8f %11.8f %6d\n', T0, W0, H, K);
% STEP 20 */
% increase H if possible */
 if T0 >= B 
 DONE = TRUE;
 else
 if T0 + H > B 
 H = B - T0;
 else
 if K <= 3 
 H = 2*H;
 end;
 end;
 end;
 if H > HMAX 
 H = H/2;
 end;
 end;
 end;
% STEP 21 */
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
