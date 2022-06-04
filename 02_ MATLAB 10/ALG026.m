 % STEFFENSEN'S ALGORITHM 2.6
 %
 % To find a solution to g(x) = x
 % given an initial approximation p0:
 %
 % INPUT:   initial approximation p0; tolerance TOL;
 %          maximum number of iterations N0.
 %
 % OUTPUT:  approximate solution p or
 %          a message that the method fails.
 syms('OK', 'P0', 'TOL', 'NO', 'FLAG', 'NAME', 'OUP', 'I', 'P1', 'P2');
 syms('D', 'P','x','s');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Steffensens Method.\n');
 fprintf(1,'Input the function G(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
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
 fprintf(1,'A:\\OUTPUT.DTA\n');
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
 fprintf(OUP, 'STEFFENSENS METHOD\n');
 if FLAG == 2 
 fprintf(OUP, '  I    P\n');
 end
% STEP 1
 I = 1;
 OK = TRUE;
% STEP 2
 while I <= NO & OK == TRUE 
% STEP 3
% compute P(1) with superscript (I-1)
 P1 = G(P0);
% compute P(2) with superscript (I-1)
 P2 = G(P1);
 if abs(P2-2*P1+P0) < 1.0e-20 
 FLAG = 1;
 D = 10;
 fprintf(OUP,'Denominator = 0, method fails\n');
 fprintf(OUP,'best possible is P2(%2d) = %15.8f\n',I,P2);
 OK = FALSE;
 else
 D = (P1-P0)*(P1-P0)/(P2-2*P1+P0);
 end
% compute P(0) with superscript (I-1)
 P = P0-D;
 if FLAG == 2 
 fprintf(OUP, '%3d %15.8e \n', I, P);
 end
% STEP 4
 if abs(D) < TOL 
% procedure completed successfully
 fprintf(OUP, '\nApproximate solution = %12.8f\n', P);
 fprintf(OUP, 'Number of iterations = %3d', I);
 fprintf(OUP, ' Tolerance = %15.8e\n',TOL);
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
 fprintf(OUP, 'not within tolerance %14e\n',TOL);
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully\n',NAME);
 end
 end
