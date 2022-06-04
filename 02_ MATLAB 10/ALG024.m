 % SECANT ALGORITHM 2.4
 %
 % To find a solution to the equation f(x) = 0
 % given initial approximations p0 and p1:
 %
 % INPUT:   initial approximation p0, p1; tolerance TOL;
 %          maximum number of iterations N0.
 %
 % OUTPUT:  approximate solution p or
 %          a message that the algorithm fails.
 syms('OK', 'P0', 'P1', 'TOL', 'NO', 'FLAG', 'NAME', 'OUP', 'F0');
 syms('I', 'F1', 'P', 'FP','s','x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Secant Method\n');
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input initial approximations P0 and P1 on separate lines.\n');
 P0 = input(' '); 
 P1 = input(' ');
 if P0 == P1 
 fprintf(1,'P0 cannot equal P1\n');
 else
 OK = TRUE;
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
 fprintf(OUP, 'Secant Method\n');
 if FLAG == 2 
 fprintf(OUP, '  I    P                 F(P)\n');
 end
% STEP 1
 I = 2;
 F0 = F(P0);
 F1 = F(P1);
 OK = TRUE;        
% STEP 2
 while I <= NO & OK == TRUE   
% STEP 3
% compute P(I)
 P = P1-F1*(P1-P0)/(F1-F0);
% STEP 4
 FP = F(P);
 if FLAG == 2 
 fprintf(OUP,'%3d   %15.8e   %15.8e\n',I,P,FP);
 end
% STEP 4
 if abs(P-P1) < TOL 
% procedure completed successfully
 fprintf(OUP,'\nApproximate solution P = %12.8f\n',P);
 fprintf(OUP,'with F(P) = %12.8f\n',FP);
 fprintf(OUP,'Number of iterations = %d\n',I);
 fprintf(OUP,'Tolerance = %14.8e\n',TOL);
 OK = FALSE;
% STEP 5
 else
 I = I+1;
% STEP 6
% update P0, F0, P1, F1
 P0 = P1;
 F0 = F1;
 P1 = P;
 F1 = FP;
 end
 end
 if OK == TRUE 
% STEP 7
% procedure completed unsuccessfully
 fprintf(OUP,'\nIteration number %d',NO);
 fprintf(OUP,' gave approximation %12.8f\n',P);
 fprintf(OUP,'with F(P) = %12.8f not within tolerance  %15.8e\n',FP,TOL);
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully\n',NAME);
 end
 end
