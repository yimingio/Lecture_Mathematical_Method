 % FIXED-POINT ALGORITHM 2.2
 %
 % To find a solution to p = g(p) given an
 % initial approximation p0
 %
 % INPUT:  initial approximation p0; tolerance TOL;
 %         maximum number of iterations NO.
 %
 % OUTPUT: approximate solution p or 
 %         a message that the method fails.
 syms('OK', 'P0', 'TOL', 'NO', 'FLAG', 'NAME', 'OUP', 'I', 'P');
 syms('x','s');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Fixed-Point Method.\n');
 fprintf(1,'Input the function G(x) in terms of x\n');
 fprintf(1,'For example: cos(x) \n');
 s = input(' ');
 G = inline(s,'x');
 OK = FALSE;
 fprintf(1,'Input initial approximation\n');
 P0 = input(' '); 
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
 fprintf(OUP, 'FIXED-POINT METHOD\n');
 if FLAG == 2 
 fprintf(OUP, '  I    P\n');
 end
% STEP 1
 I = 1;
 OK = TRUE; 
% STEP 2
 while I <= NO & OK == TRUE 
% STEP 3
% compute P(I)
 P = G(P0);
 if FLAG == 2 
 fprintf(OUP, '%3d   %15.8e\n', I, P);
 end
% STEP 4
 if abs(P-P0) < TOL 
% procedure completed successfully
 fprintf(OUP, '\nApproximate solution P = %12.8f\n', P);
 fprintf(OUP, 'Number of iterations = %3d', I);
 fprintf(OUP, '    Tolerance = %14.8e\n',TOL);
 OK = FALSE;
 else
% STEP 5
 I = I+1;
% STEP 6
% update P0
 P0 = P;
 end
 end
 if OK == TRUE 
% STEP 7
% procedure completed unsuccessfully
 fprintf(OUP, '\nIteration number %3d', NO);
 fprintf(OUP, ' gave approximation %12.8f\n', P);
 fprintf(OUP, 'not within tolerance %14.8e\n',TOL);
 end
 if OUP ~= 1
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end
 end
