% PIECEWISE LINEAR RAYLEIGH-RITZ ALGORITHM 11.5
%
% To approximate the solution of the boundary-value problem
%
%      -D(P(X)Y')/DX + Q(X)Y = F(X), 0 <= X <= 1,
%             Y(0) = Y(1) = 0,
%
% with a piecewise linear function:
%
% INPUT:   integer N; mesh points X(0) = 0 < X(1) < ...
%          < X(N) < X(N+1) = 1
%
% OUTPUT:  coefficients C(1),...,C(N) of the basis functions
 syms('AA', 'OK', 'N', 'X', 'FLAG', 'HC', 'J', 'H', 'NAME');
 syms('INP', 'N1', 'Q', 'ALPHA', 'BETA', 'B', 'A', 'ZETA');
 syms('Z', 'C', 'J1', 'OUP', 'x', 'sqq', 'sp', 'sf', 'I');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Piecewise Linear Rayleigh-Ritz Method.\n');
 fprintf(1,'The M-file SIMPSON.M is used by this program.\n');
 fprintf(1,'Input F(X), Q(X), and P(X) in terms of x \n');
 fprintf(1,' on separate lines. \n');
 fprintf(1,'For example:   2*pi^2*sin(pi*x) \n');
 fprintf(1,'                        pi^2 \n');
 fprintf(1,'                           1 \n');
 sf = input(' ');
 F = inline(sf,'x');
 sqq = input(' ');
 QQ = inline(sqq,'x');
 sp = input(' ');
 P = inline(sp,'x');
 fprintf(1,'X(0), ..., X(N+1) are to be supplied.\n');
 fprintf(1,'Are the preparations complete? Answer Y or N.\n');
 AA = input(' ','s');
 OK = FALSE;
 if AA == 'Y' | AA == 'y'
 OK = FALSE;
 while OK == FALSE
 fprintf(1,'Input integer N where X(0) = 0, X(N+1) = 1.\n');
 N = input(' ');
 if N <= 1
 fprintf(1,'N must be greater than one.\n');
 else
 OK = TRUE;
 end;
 end;
 X = zeros(1,N+2);
 H = zeros(1,N+1);
 Q = zeros(6,N+1);
 A = zeros(1,N+1);
 B = zeros(1,N+1);
 C = zeros(1,N+1);
 ALPHA = zeros(1,N+1);
 BETA = zeros(1,N+1);
 ZETA = zeros(1,N+1);
 Z = zeros(1,N+1);
 X(1) = 0;
 X(N+2) = 1;
 fprintf(1,'Choice of method to input X(1), ..., X(N):\n');
 fprintf(1,'1.  Input from keyboard at the prompt\n');
 fprintf(1,'2.  Equally spaced nodes to be calculated\n');
 fprintf(1,'3.  Input from text file\n');
 fprintf(1,'Please enter 1, 2, or 3.\n');
 FLAG = input(' ');
 if FLAG == 2
 HC = 1/(N+1);
 for J = 1 : N
 X(J+1) = J*HC;
 H(J) = HC;
 end;
 H(N+1) = HC;
 else
 if FLAG == 3
 fprintf(1,'Has the input file been created?  ');
 fprintf(1,'Enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y'
 fprintf(1,'Enter the input file name using the format\n');
 fprintf(1,' - drive:\\name.ext,\n');
 fprintf(1,'for example:  A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 for J = 2 : N + 1
 X(J) = fscanf(INP, '%f',1);
 end;
 for J = 1 : N + 1
 H(J) = X(J+1)-X(J);
 end;
 fclose(INP);
 else
 fprintf(1,'The program will end so that the input ');
 fprintf(1,'file can be created.\n');
 OK = FALSE;
 end;
 else
 for J = 2 : N + 1
 I = J-1;      
 fprintf(1,'Input X(%d).\n', I);
 X(J) = input(' ');
 H(J-1) = X(J)-X(J-1);
 end;
 H(N+1) = X(N+2)-X(N+1);
 end;
 end;
 else
 fprintf(1,'The program will end so that the preparations\n');
 fprintf(1,'can be completed.\n');
 OK = FALSE;
 end;
% STEP 1 is done within the input procedure
 if OK == TRUE
 N1 = N-1;
% STEP 3
 for J = 2 : N
 Q(1,J-1) = SIMPSON(1,X(J),X(J+1),sqq,sf,sp)/((H(J))*(H(J)));
 Q(2,J-1) = SIMPSON(2,X(J-1),X(J),sqq,sf,sp)/((H(J-1))*(H(J-1)));
 Q(3,J-1) = SIMPSON(3,X(J),X(J+1),sqq,sf,sp)/((H(J))*(H(J)));
 Q(4,J-1) = SIMPSON(4,X(J-1),X(J),sqq,sf,sp)/((H(J-1))*(H(J-1)));
 Q(5,J-1) = SIMPSON(5,X(J-1),X(J),sqq,sf,sp)/H(J-1) ;
 Q(6,J-1) = SIMPSON(6,X(J),X(J+1),sqq,sf,sp)/H(J);
 end;
 Q(2,N) = SIMPSON(2,X(N),X(N+1),sqq,sf,sp)/((H(N))*(H(N)));
 Q(3,N) = SIMPSON(3,X(N+1),X(N+2),sqq,sf,sp)/((H(N+1))*(H(N+1)));
 Q(4,N) = SIMPSON(4,X(N),X(N+1),sqq,sf,sp)/((H(N))*(H(N)));
 Q(4,N+1) = SIMPSON(4,X(N+1),X(N+2),sqq,sf,sp)/((H(N+1))*(H(N+1)));
 Q(5,N) = SIMPSON(5,X(N),X(N+1),sqq,sf,sp)/H(N);
 Q(6,N) = SIMPSON(6,X(N+1),X(N+2),sqq,sf,sp)/H(N+1);
% STEP 4
 for J = 2 : N1 + 1
 ALPHA(J-1) = Q(4,J-1)+Q(4,J)+Q(2,J-1)+Q(3,J-1);
 BETA(J-1) = Q(1,J-1)-Q(4,J);
 B(J-1) = Q(5,J-1)+Q(6,J-1);
 end;
% STEP 5
 ALPHA(N) = Q(4,N)+Q(4,N+1)+Q(2,N)+Q(3,N);
 B(N) = Q(5,N)+Q(6,N);
% STEPS 6-10 solve a symmetric tridiagonal linear system using Algorithm 6.7
% STEP 6
 A(1) = ALPHA(1);
 ZETA(1) = BETA(1)/ALPHA(1);
 Z(1) = B(1)/A(1);
% STEP 7
 for J = 2 : N1
 A(J) = ALPHA(J)-BETA(J-1)*ZETA(J-1);
 ZETA(J) = BETA(J)/A(J);
 Z(J) = (B(J)-BETA(J-1)*Z(J-1))/A(J);
 end;
% STEP 8
 A(N) = ALPHA(N) - BETA(N-1) * ZETA(N-1);
 Z(N) = (B(N)-BETA(N-1)*Z(N-1))/A(N);
% STEP 9
 C(N) = Z(N);
 for J = 1 : N1
 J1 = N - J;
 C(J1) = Z(J1)-ZETA(J1)*C(J1+1);
 end;
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text file\n');
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
 fprintf(OUP, 'PIECEWISE LINEAR RAYLEIGH-RITZ METHOD\n\n');
 fprintf(OUP, ' I      X(I-1)       X(I)       X(I+1)        C(I)\n\n');
 for J = 1 : N
 fprintf(OUP,'%3d %11.8f %11.8f %11.8f %13.8f\n',J,X(J),X(J+1),X(J+2),C(J));
 end;
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
