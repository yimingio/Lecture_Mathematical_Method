% LINEAR SHOOTING ALGORITHM 11.1
%
% To approximate the solution of the boundary-value problem
%
% -Y'' + P(X)Y' + Q(X)Y + R(X) = 0, A<=X<=B, Y(A)=ALPHA, Y(B)=BETA:
%
%
% INPUT: Endpoints A,B; boundary conditions ALPHA, BETA; number of
%        subintervals N.
%
% OUTPUT: Approximations W(1,I) to Y(X(I)); W(2,I) to Y'(X(I))
%         for each I=0,1,...,N.
 syms('OK', 'A', 'B', 'ALPHA', 'BETA', 'N', 'FLAG', 'NAME');
 syms('OUP', 'H', 'U1', 'U2', 'V1', 'V2', 'I', 'X', 'T');
 syms('K11', 'K12', 'K21', 'K22', 'K31', 'K32', 'K41', 'K42');
 syms('U', 'V', 'W1', 'Z', 'W2', 's', 'x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Linear Shooting Method.\n');
 fprintf(1,'Input the functions P(X), Q(X) and R(X) in terms of x,\n'); 
 fprintf(1,'on separate lines.\n');
 fprintf(1,'For example: -2/x \n');
 fprintf(1,'2/(x^2) \n');
 fprintf(1,'sin(log(x))/(x^2)\n');
 s = input(' ');
 P = inline(s,'x');
 s = input(' ');
 Q = inline(s,'x');
 s = input(' ');
 R = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input left and right endpoints on separate lines.\n');
 A = input(' ');
 B = input(' ');
 if A >= B 
 fprintf(1,'Left endpoint must be less than right endpoint.\n');
 else
 OK = TRUE;
 end;
 end;
 fprintf(1,'Input Y(  %.10e).\n', A);
 ALPHA = input(' ');
 fprintf(1,'Input Y(  %.10e).\n', B);
 BETA = input(' ');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input a positive integer for the number of subintervals.\n');
 N = input(' ');
 if N <= 0 
 fprintf(1,'Number must be a positive integer.\n');
 else
 OK = TRUE;
 end;
 end;
 if OK == TRUE 
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text File\n');
 fprintf(1,'Please enter 1 or 2.\n');
 FLAG = input(' ');
 if FLAG == 2 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example  A:\\OUTPUT.DTA\n');
 NAME = input(' ','s');
 OUP = fopen(NAME,'wt');
 else
 OUP = 1;
 end;
 fprintf(OUP, 'LINEAR SHOOTING METHOD\n\n');
 fprintf(OUP, '  I  X(I)        W(1,I)      W(2,I)\n');
% STEP 1
 H = (B-A)/N;
 U1 = ALPHA;
 U2 = 0;
 V1 = 0;
 V2 = 1;
 U = zeros(2,N);
 V = zeros(2,N);
% STEP 2
 for I = 1 : N 
% STEP 3
 X = A+(I-1)*H;
 T = X+0.5*H;
% STEP 4
 K11 = H*U2;
 K12 = H*(P(X)*U2+Q(X)*U1+R(X));
 K21 = H*(U2+0.5*K12);
 K22 = H*(P(T)*(U2+0.5*K12)+Q(T)*(U1+0.5*K11)+R(T));
 K31 = H*(U2+0.5*K22);
 K32 = H*(P(T)*(U2+0.5*K22)+Q(T)*(U1+0.5*K21)+R(T));
 T = X+H;
 K41 = H*(U2+K32);
 K42 = H*(P(T)*(U2+K32)+Q(T)*(U1+K31)+R(T));
 U1 = U1+(K11+2*(K21+K31)+K41)/6;
 U2 = U2+(K12+2*(K22+K32)+K42)/6;
 K11 = H*V2;
 K12 = H*(P(X)*V2+Q(X)*V1);
 T = X+0.5*H;
 K21 = H*(V2+0.5*K12);
 K22 = H*(P(T)*(V2+0.5*K12)+Q(T)*(V1+0.5*K11));
 K31 = H*(V2+0.5*K22);
 K32 = H*(P(T)*(V2+0.5*K22)+Q(T)*(V1+0.5*K21));
 T = X+H;
 K41 = H*(V2+K32);
 K42 = H*(P(T)*(V2+K32)+Q(T)*(V1+K31));
 V1 = V1+(K11+2*(K21+K31)+K41)/6;
 V2 = V2+(K12+2*(K22+K32)+K42)/6;
 U(1,I) = U1;
 U(2,I) = U2;
 V(1,I) = V1;
 V(2,I) = V2;
 end;
% STEP 5
 W1 = ALPHA;
 Z = (BETA-U(1,N))/V(1,N);
 X = A;
 I = 0;
 fprintf(OUP, '%3d %11.8f %11.8f %11.8f\n', I, X, W1, Z);
 for I = 1 : N 
 X = A+I*H;
 W1 = U(1,I)+Z*V(1,I);
 W2 = U(2,I)+Z*V(2,I);
 fprintf(OUP, '%3d %11.8f %11.8f %11.8f\n', I, X, W1, W2);
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
% STEP 7
