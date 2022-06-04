 % ROMBERG ALGORITHM 4.2
 %     
 % To approximate I = integral ( f(x) dx ) from a to b:
 %
 % INPUT:   endpoints a, b; integer n.
 %
 % OUTPUT:  an array R. ( R(2,n) is the approximation to I. )
 %
 % R is computed by rows; only 2 rows saved in storage   
 syms('OK','A','B','N','H','R','I','SUM','M','K','J','L','s','x');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is Romberg Integration.\n\n');
 fprintf(1,'Input the function F(x) in terms of x\n');
 fprintf(1,'For example: cos(x)\n');
 s = input(' ');
 F = inline(s,'x');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input lower limit of integration and ');
 fprintf(1,'upper limit of integration\n');
 fprintf(1,'on separate lines\n');
 A = input(' ');
 B = input(' ');
 if A > B 
 fprintf(1,'Lower limit must be less than upper limit\n');
 else
 OK = TRUE;
 end
 end 
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input number of rows - no decimal point.\n');
 N = input(' ');
 if N > 0 
 OK = TRUE;
 else
 fprintf(1,'Number must be a positive integer\n');
 end
 end
% STEP 1
 if OK == TRUE 
 H = B-A;
 R = zeros(2,N+1);
 R(1,1) = (F(A)+F(B))/2*H;
% STEP 2
 fprintf(1,'Initial Data:\n');
 fprintf(1,'Limits of integration = (%12.8f, %12.8f)\n', A, B);
 fprintf(1,'Number of rows = %3d\n', N);
 fprintf(1,'\nRomberg Integration Table:\n');
 fprintf(1,'\n%12.8f\n\n', R(1,1));
% STEP 3
 for I = 2:N
% STEP 4
% approximation from Trapezoidal method
 SUM = 0;
 M = 2^(I-2);
 for K = 1:M
 SUM = SUM+F(A+(K-0.5)*H);
 end
 R(2,1) = (R(1,1)+H*SUM)/2;
% STEP 5
% extrapolation
 for J = 2:I
 L = 2^(2*(J-1));
 R(2,J) = R(2,J-1)+(R(2,J-1)-R(1,J-1))/(L-1);
 end
% STEP 6
 for K = 1:I
 fprintf(1,' %11.8f',R(2,K));
 end
 fprintf(1,'\n\n');
% STEP 7
 H = H/2;
% since only two rows are kept in storage, this step is 
% to prepare for the next row
% update row 1 of R
% STEP 8
 for J = 1:I
 R(1,J) = R(2,J);
 end
 end
 end
