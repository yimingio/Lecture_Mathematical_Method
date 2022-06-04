
%  BISECTION ALGORITHM 2.1
% 
%  To find a solution to f(x) = 0 given the continuous function
%  f on the interval [a,b], where f(a) and f(b) have
%  opposite signs:
% 
%  INPUT:   endpoints a,b; tolerance TOL;
%           maximum number of iterations NO.
% 
%  OUTPUT:  approximate solution p or
%           a message that the algorithm fails.
 syms('OK','A','B','X','FA','FB','TOL','NO','FLAG','NAME','OUP','I')
 syms('C','P','FP','x','s')
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Bisection Method.\n');
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n ');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input endpoints A < B on separate lines\n');
 A = input(' ');
 B = input(' ');
 if A > B 
 X = A;
 A = B;
 B = X;
 end
 if A == B 
 fprintf(1,'A cannot equal B\n');
 else
 FA = F(A);
 FB = F(B);
 if FA*FB > 0
 fprintf(1,'F(A) and F(B) have same sign\n');
 else
 OK = TRUE;
 end
 end
 end
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
 NO = input(' ');
 if NO <= 0 
 fprintf(1,'Must be positive integer\n');
 else 
 OK = TRUE;
 end
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
 fprintf(1,'Select amount of output\n');
 fprintf(1,'1. Answer only\n');
 fprintf(1,'2. All intermediate approximations\n');
 fprintf(1,'Enter 1 or 2\n');
 FLAG = input(' ');
 fprintf(OUP,'Bisection Method\n');
 if FLAG == 2 
 fprintf(OUP, '  I    P                  F(P)\n');
 end
% STEP 1
 I = 1;
% STEP 2
 OK = TRUE;
 while I <= NO & OK == TRUE 
% STEP 3
% Compute P(I)
 C = (B - A) / 2.0;
 P = A + C;
% STEP 4
 FP = F(P);
 if FLAG == 2 
 fprintf(OUP,'%3d   %15.8e   %15.7e \n',I,P,FP);
 end
 if abs(FP) < 1.0e-20 | C < TOL 
% procedure completed successfully
 fprintf(OUP,'\nApproximate solution P = %11.8f \n',P);
 fprintf(OUP,'with F(P) = %12.8f\n',FP);
 fprintf(OUP,'Number of iterations = %3d',I);
 fprintf(OUP,' Tolerance = %15.8e\n',TOL);
 OK = FALSE;
 else
% STEP 5
 I = I+1;
% STEP 6
% compute A(I) and B(I)
 if FA*FP > 0
 A = P;
 FA = FP;
 else
 B = P;
 FB = FP;
 end
 end
 end
 if OK == TRUE 
% STEP 7
% procedure completed unsuccessfully
 fprintf(OUP,'\nIteration number %3d',NO);
 fprintf(OUP,' gave approximation %12.8f\n',P);
 fprintf(OUP,'F(P) = %12.8f not within tolerance : %15.8e\n',FP,TOL);
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end
 end
 
 
