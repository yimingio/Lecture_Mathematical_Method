 % NEVILLE'S ITERATED INTERPOLATION ALGORITHM 3.1
 %
 % To evaluate the interpolating polynomial P on the
 % (n+1) distinct numbers x(0), ..., x(n) at the number x
 % for the function f:
 %
 % INPUT:   numbers x(0),..., x(n) as XX(0),...,XX(N);
 %          number x; values of f as the first column of Q
 %          or may be computed if function f is supplied.
 %
 % OUTPUT:  the table Q with P(x) = Q(N+1,N+1).
 syms('TRUE', 'FALSE', 'OK', 'FLAG', 'N', 'I', 'XX', 'Q', 'A', 'NAME');
 syms('INP', 'X', 'D', 'J', 'OUP','x','s');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Nevilles Method.\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Choice of input method:\n');
 fprintf(1,'1. Input entry by entry from keyboard\n');
 fprintf(1,'2. Input data from a text file\n');
 fprintf(1,'3. Generate data using a function F\n');
 fprintf(1,'Choose 1, 2, or 3 please\n');
 FLAG = input(' ');
 if FLAG == 1 | FLAG == 2 | FLAG == 3 
 OK = TRUE;
 end
 end
 if FLAG == 1 
 OK = FALSE;
 while OK ~= TRUE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0 
 OK = TRUE;
 XX = zeros(N+1);
 Q = zeros(N+1,N+1);
 for I = 0:N
 fprintf(1,'Input X(%d) and F(X(%d)) ', I, I);
 fprintf(1,'on separate lines.\n');
 XX(I+1) = input(' ');
 Q(I+1,1) = input(' ');
 end
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 end
 if FLAG == 2 
 fprintf(1,'Has a text file been created with the data in two columns?\n');
 fprintf(1,'Enter Y or N\n');
 A = input(' ','s');
 if A == 'Y' | A == 'y' 
 fprintf(1,'Input the file name in the form - ');
 fprintf(1,'drive:\\name.ext\n');
 fprintf(1,'For example:   A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input N\n');
 N = input(' ');
 if N > 0 
 XX = zeros(N+1);
 Q = zeros(N+1,N+1);
 for I = 0:N
 XX(I+1) = fscanf(INP, '%f',1);
 Q(I+1,1) = fscanf(INP, '%f',1);
 end
 fclose(INP);
 OK = TRUE;
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 else
 fprintf(1,'Please create the input file in two column ');
 fprintf(1,'form with the X values and\n');
 fprintf(1,'F(X) values in the corresponding columns.\n');
 fprintf(1,'The program will end so the input file can ');
 fprintf(1,'be created.\n');
 OK = FALSE;
 end
 end
 if FLAG == 3 
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0 
 XX = zeros(N+1);
 Q = zeros(N+1,N+1);
 for I = 0:N
 fprintf(1,'Input X(%d)\n', I);
 XX(I+1) = input(' ');
 Q(I+1,1) = F(XX(I+1));
 end
 OK = TRUE;
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 end
 if OK == TRUE 
 fprintf(1,'Input point at which the polynomial is to be evaluated\n');
 X = input(' ');
 end
 if OK == TRUE 
% STEP 1
 D = zeros(N+1);
 D(1) = X-XX(1);
 for I = 1:N
 D(I+1) = X-XX(I+1);
 for J = 1:I
 Q(I+1,J+1) = (D(I+1)*Q(I,J)-D(I-J+1)*Q(I+1,J))/(D(I+1)-D(I-J+1));
 end
 end
% STEP 2
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
 fprintf(OUP, 'NEVILLES METHOD\n');
 fprintf(OUP, 'Table for P evaluated at X = %12.8f , follows: \n', X);
 fprintf(OUP, 'Entries are XX(I), Q(I,0), ..., Q(I,I) ');
 fprintf(OUP, 'for each I = 0, ..., N where N = %3d\n\n', N); 
 for I = 0:N
 fprintf(OUP, '%11.8f ', XX(I+1));
 for J = 0:I
 fprintf(OUP, '%11.8f ', Q(I+1,J+1));
 end
 fprintf(OUP, '\n');
 end
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully\n',NAME);
 end
 end
