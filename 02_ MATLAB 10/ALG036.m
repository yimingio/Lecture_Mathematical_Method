 % BEZIER CURVE ALGORITHM 3.6
 %
 % To construct the cubic Bezier curves C0, ..., Cn-1 in
 % parameter form, where Ci is represented by
 %
 % (xi(t),yi(t)) = ( a0(i) + a1(i)*t + a2(i)*t^2 + a3(i)*t^3,
 %                   b0(i) + b1(i)*t + b2(i)*t^2 + b3(i)*t^3)
 %
 % for 0 <= t <= 1 as determined by the left endpoint (x(i),y(i)),
 % left guidepoint (x+(i),y+(i)), right endpoint (x(i+1),y(i+1)) and
 % right guidepoint (x-(i+1),y-(i+1)) for each i = 0, 1, ... , n-1;
 %
 % INPUT  n, ( (x(i),y(i)), i = 0,...,n ),
 %           ( (x+(i),y+(i)), i = 0,...,n-1 ),
 %           ( (x-(i),y-(i)), i = 1,...,n ).
 %
 % OUTPUT coefficients ( a0(i), a1(i), a2(i), a3(i),
 %                       b0(i), b1(i), b2(i), b3(i), i = 0, ... , n-1 ).
 syms('OK', 'FLAG', 'N', 'X', 'Y', 'XPL', 'YPL', 'I', 'XMI', 'YMI');
 syms('A','NAME','INP','OUP','A0','B0','A1','B1','A2','B2','A3','B3');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Bezier Curve Algorithm.\n');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Choice of input method:\n');
 fprintf(1,'1. Input entry by entry from keyboard\n');
 fprintf(1,'2. Input data from a text file\n');
 fprintf(1,'Choose 1 or 2 please\n');
 FLAG = input(' ');
 if FLAG == 1 | FLAG == 2 
 OK = TRUE;
 end
 end
 if FLAG == 1 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0 
 OK = TRUE;
 X = zeros(1,N+1);
 Y = zeros(1,N+1);
 XPL = zeros(1,N+1);
 YPL = zeros(1,N+1);
 XMI = zeros(1,N+1);
 YMI = zeros(1,N+1);
 fprintf(1,'Input X(0),Y(0),X+(0),Y+(0)\n');
 fprintf(1,'on separate lines\n');
 X(1) = input(' ');
 Y(1) = input(' ');
 XPL(1) = input(' ');
 YPL(1) = input(' ');
 for I = 1:N-1
 fprintf(1,'Input X(%d),Y(%d)\n', I, I);
 fprintf(1,'on separate lines\n');
 X(I+1) = input(' ');
 Y(I+1) = input(' ');
 fprintf(1,'Input X-(%d),Y-(%d)\n', I, I);
 fprintf(1,'on separate lines\n');
 XMI(I) = input(' ');
 YMI(I) = input(' ');
 fprintf(1,'Input X+(%d),Y+(%d)\n', I, I);
 fprintf(1,'on separate lines\n');
 XPL(I+1) = input(' ');
 YPL(I+1) = input(' ');
 end
 fprintf(1,'Input X(n),Y(n),X-(n),Y-(n)\n');
 fprintf(1,'on separate lines\n');
 X(N+1) = input(' ');
 Y(N+1) = input(' ');
 XMI(N) = input(' ');
 YMI(N) = input(' ');
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 end
 if FLAG == 2 
 fprintf(1,'Has a text file been created with the data as follows ?\n\n');
 fprintf(1,'X(0)    Y(0)    X+(0)    Y+(0)\n');
 fprintf(1,'X(1)    Y(1)    X-(1)    Y-(1)    X+(1)    Y+(1)\n');
 fprintf(1,'...\n');
 fprintf(1,'X(n-1)  Y(n-1)  X-(n-1)  Y-(n-1)  X+(n-1)  Y+(n-1)\n');
 fprintf(1,'X(n)    Y(n)    X-(n)    Y-(n)\n\n');
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
 fprintf(1,'Input n\n');
 N = input(' ');
 if N > 0 
 OK = TRUE;
 X = zeros(1,N+1);
 Y = zeros(1,N+1);
 XPL = zeros(1,N+1);
 YPL = zeros(1,N+1);
 XMI = zeros(1,N+1);
 YMI = zeros(1,N+1);
 X(1) = fscanf(INP, '%f',1);
 Y(1) = fscanf(INP, '%f',1);
 XPL(1) = fscanf(INP, '%f',1);
 YPL(1) = fscanf(INP, '%f',1);
 for I = 1:N-1
 X(I+1) = fscanf(INP, '%f',1);
 Y(I+1) = fscanf(INP, '%f',1);
 XMI(I) = fscanf(INP, '%f',1);
 YMI(I) = fscanf(INP, '%f',1);
 XPL(I+1) = fscanf(INP, '%f',1);
 YPL(I+1) = fscanf(INP, '%f',1);
 end
 X(N+1) = fscanf(INP, '%f',1);
 Y(N+1) = fscanf(INP, '%f',1);
 XMI(N) = fscanf(INP, '%f',1);
 YMI(N) = fscanf(INP, '%f',1);
 fclose(INP);
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
 else
 fprintf(1,'Please create the input file as indicated.\n');
 fprintf(1,'The program will end so the input file can ');
 fprintf(1,'be created.\n');
 OK = FALSE;
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
 fprintf(OUP, 'BEZIER CURVE ALGORITHM\n\n');
 fprintf(OUP, '          A0          A1          A2          A3');
 fprintf(OUP,'  on the first line\n');
 fprintf(OUP, '          B0          B1          B2          B3');
 fprintf(OUP,'  on the second line\n');
% STEP1
 A0 = zeros(1,N+1);
 B0 = zeros(1,N+1);
 A1 = zeros(1,N+1);
 B1 = zeros(1,N+1);
 A2 = zeros(1,N+1);
 B2 = zeros(1,N+1);
 A3 = zeros(1,N+1);
 B3 = zeros(1,N+1);
 for I = 0:N-1
% STEP 2
 A0(I+1) = X(I+1);
 B0(I+1) = Y(I+1);
 A1(I+1) = 3*(XPL(I+1) - X(I+1));
 B1(I+1) = 3*(YPL(I+1) - Y(I+1));
 A2(I+1) = 3*(X(I+1)+XMI(I+1)-2*XPL(I+1));
 B2(I+1) = 3*(Y(I+1)+YMI(I+1)-2*YPL(I+1));
 A3(I+1) = X(I+2)-X(I+1)+3*XPL(I+1)-3*XMI(I+1);
 B3(I+1) = Y(I+2)-Y(I+1)+3*YPL(I+1)-3*YMI(I+1);
% STEP 3
 fprintf(OUP,'%11.6f %11.6f %11.6f %11.6f\n',A0(I+1),A1(I+1),A2(I+1),A3(I+1));
 fprintf(OUP,'%11.6f %11.6f %11.6f %11.6f\n',B0(I+1),B1(I+1),B2(I+1),B3(I+1));
 fprintf(OUP, '\n');
 end
% STEP 4
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully\n',NAME);
 end
 end
