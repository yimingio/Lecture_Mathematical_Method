 % METHOD OF FALSE POSITION ALGORITHM 2.5
 %
 % To find a solution to f(x) = 0 given the continuous function
 % f on the interval [p0,p1], where f(p0) and f(p1) have
 % opposite signs:
 %
 % INPUT:   endpoints p0, p1; tolerance TOL;
 %          maximum number of iterations N0.
 %
 % OUTPUT:  approximate solution p or
 %          a message that the algorithm fails.
 syms('OK', 'P0', 'P1', 'X', 'Q0', 'Q1', 'TOL', 'NO', 'FLAG');
 syms('NAME', 'OUP', 'I', 'P', 'Q','x','s');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Method of False Position\n');
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input endpoints P0 < P1 on separate lines.\n');
 P0 = input(' '); 
 P1 = input(' ');
 if P0 > P1 
 X = P0;
 P0 = P1;
 P1 = X;
 end
 if P0 == P1 
 fprintf(1,'P0 cannot equal P1\n');
 else
 Q0 = F(P0);
 Q1 = F(P1);
 if Q0*Q1 > 0 
 fprintf(1,'F(P0) and F(P1) have the same sign.\n');
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
 fprintf(OUP, 'METHOD OF FALSE POSITION OR REGULA FALSII\n\n');
 if FLAG == 2 
 fprintf(OUP, '  I    P                 F(P)\n');
 end
% STEP 1
 I = 2;
 OK = TRUE;
 Q0 = F(P0);
 Q1 = F(P1);
% STEP 2
 while I <= NO & OK == TRUE
% STEP 3
% compute P(I)
 P = P1-Q1*(P1-P0)/(Q1-Q0);
 Q = F(P);
 if FLAG == 2 
 fprintf(OUP,'%3d %15.8e %15.8e\n',I,P,Q);
 end
% STEP 4
 if abs(P-P1) < TOL 
% procedure completed successfully
 fprintf(OUP,'\nApproximate solution P = %12.8f\n',P);
 fprintf(OUP,'with F(P) = %12.8f\n',Q);
 fprintf(OUP,'Number of iterations = %3d',I);
 fprintf(OUP,' Tolerance = %15.8e\n',TOL);
 OK = FALSE;
 else
% STEP 5
 I = I+1;
% STEP 6
% compute P0(I) and P1(I)
 if Q*Q1 < 0 
 P0 = P1;
 Q0 = Q1;
 end
% STEP 7
 P1 = P;
 Q1 = Q;
 end
 end
 if OK == TRUE 
% procedure completed unsuccessfully
 fprintf(OUP,'\nIteration number %3d',NO);
 fprintf(OUP,' gave approximation %12.8f\n',P);
 fprintf(OUP,'F(P) = %12.8f not within tolerance: %15.8e\n',Q,TOL);
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully\n',NAME);
 end
 end

