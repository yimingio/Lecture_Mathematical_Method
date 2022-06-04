% QR ALGORITHM 9.6
%
% To obtain the eigenvalues of a symmetric, tridiagonal n by n matrix
%
%          a(1)   b(2)
%          b(2)   a(2)   b(3)
%             .      .      .
%               .      .      .
%                 .      .      .
%                 b(n-1)  a(n-1)  b(n)
%                            b(n)   a(n)
%
% INPUT:   n; A(1),...,A(n) (diagonal of A); B(2),...,B(n)
%          (off-diagonal of A); maximum number of iterations M, 
%          tolerance TOL.
%
% OUTPUT:  Eigenvalues of A or recommended splitting of A, or a 
%          message that the maximum number of iterations was 
%          exceeded.
 syms('OK', 'AA', 'NAME', 'INP', 'N', 'I', 'A', 'B');
 syms('TOL', 'L', 'FLAG', 'OUP', 'SHIFT', 'K', 'J');
 syms('M', 'B1', 'C1', 'D1', 'X1', 'X2', 'D', 'X', 'Y');
 syms('Z', 'C', 'S', 'Q', 'R', 'MM');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the QR Method.\n');
 OK = FALSE;
 fprintf(1,'The tridiagonal symmetric array A will be input from ');
 fprintf(1,'a text file in the order:\n');
 fprintf(1,' (diagonal): A(1), A(2), ..., A(n),\n');
 fprintf(1,' (subdiagonal): B(2), B(3), ..., B(n).\n\n');
 fprintf(1,'Place as many entries as desired on each line, but separate ');
 fprintf(1,'entries with\n');
 fprintf(1,'at least one blank.\n\n\n');
 fprintf(1,'Has the input file been created? - enter Y or N.\n');
 AA = input(' ','s');
 if AA == 'Y' | AA == 'y' 
 fprintf(1,'Input the file name in the form - drive:\\name.ext\n');
 fprintf(1,'for example:  A:\\DATA.DTA\n');
 NAME = input(' ','s');
 INP = fopen(NAME,'rt');
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the dimension n.\n');
 N = input(' ');
 if N > 1 
 A = zeros(1,N);
 B = zeros(1,N);
 C = zeros(1,N);
 D = zeros(1,N);
 Q = zeros(1,N);
 R = zeros(1,N);
 S = zeros(1,N);
 X = zeros(1,N);
 Y = zeros(1,N);
 Z = zeros(1,N);
 for I = 1 : N 
 A(I) = fscanf(INP, '%f',1);
 end;
 for I = 2 : N 
 B(I) = fscanf(INP, '%f',1);
 end;
 OK = TRUE;
 else
 fprintf(1,'Dimension must be greater then 1.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the tolerance.\n');
 TOL = input(' ');
 if TOL > 0 
 OK = TRUE;
 else
 fprintf(1,'Tolerance must be a positive real number.\n');
 end;
 end;
 OK = FALSE;
 while OK == FALSE 
 fprintf(1,'Input the maximum number of iterations.\n');
 L = input(' ');
 if L > 0 
 OK = TRUE;
 else
 fprintf(1,'The number must be a positive integer.\n');
 end;
 end;
 fclose(INP);
 else
 fprintf(1,'The program will end so the input file can be created.\n');
 end;
 if OK == TRUE 
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
 fprintf(OUP, 'QR METHOD\n\n');
% STEP 1 */
 SHIFT = 0;
 K = 1;
% STEP 2 */
 while K <= L & OK == TRUE 
 fprintf(OUP, 'Iteration number %d N = %d\n', K, N);
 fprintf(OUP, 'The array A is now as follows:\n');
 fprintf(OUP, 'Diagonal:\n');
 for I = 1 : N 
 fprintf(OUP, ' %11.8f', A(I));
 end;
 fprintf(OUP, '\nOff diagonal:\n');
 for I = 2 : N 
 fprintf(OUP, ' %11.8f', B(I));
 end;
 fprintf(OUP, '\n');
% Steps 3-7 test for success  */
% STEP 3  */
 if abs(B(N)) <= TOL 
 A(N) = A(N) + SHIFT;
 fprintf(OUP, 'Eigenvalue = %12.8f\n', A(N));
 N = N-1;
 end;
% STEP 4*/
 if abs(B(2)) <= TOL 
 A(1) = A(1)+SHIFT;
 fprintf(OUP, 'Eigenvalue = %12.8f\n', A(1));
 N = N-1;
 A(1) = A(2);
 for J = 2 : N 
 A(J) = A(J+1);
 B(J) = B(J+1);
 end;
 end;
% STEP 5  */
 if N == 0 
 OK = FALSE;
 end;
% STEP 6*/
 if N == 1 
 A(1) = A(1) + SHIFT;
 fprintf(OUP,'Eigenvalue = %12.8f\n', A(1));
 OK = FALSE;
 end;
% STEP 7  */
 if OK == TRUE 
 M = N-1;
 if M >= 2 
 for I = 2 : M 
 if abs(B(I)) <= TOL 
 OK = FALSE;
 J = I;
 end;
 end;
 if OK == FALSE 
 fprintf(OUP, 'Split the matrix into\n');
 for I = 1 : J-1 
 fprintf(OUP,'%11.8f',A(I));
 end;
 fprintf(OUP,'\n');
 for I = 2 : J-1 
 fprintf(OUP,'%11.8f',B(I));
 end;
 fprintf(OUP,'\n and \n');
 for I = J : N 
 fprintf(OUP,'%11.8f',A(I));
 end;
 fprintf(OUP,'\n');
 for I = J+1 : N 
 fprintf(OUP,'%11.8f',B(I));
 end;
 fprintf(OUP,'\n');
 end;
 end;
 end;
 if OK == TRUE 
% STEP 8 */
% compute shift */
 B1 = -(A(N)+A(N-1));
 C1 = A(N)*A(N-1)-B(N)*B(N);
 D1 = B1*B1-4*C1;
 if D1 < 0 
 fprintf(OUP, 'Problem with complex roots, D1 = %.8e\n', D1);
 OK = FALSE;
 else
 D1 = sqrt(D1);
%  STEP 9*/
 if B1 > 0 
 X1 = -2*C1/(B1+D1);
 X2 = -(B1+D1)/2;
 else
 X1 = (D1-B1)/2;
 X2 = 2*C1/(D1-B1);
 end;
% if N = 2 then the 2 eigenvalues have been computed */
%  STEP 10  */
 if N == 2 
 X1 = X1+SHIFT;
 X2 = X2+SHIFT;
 fprintf(OUP, 'The last two eigenvalues are: %12.8f%11.8f\n',X1, X2);
 OK = FALSE;
 else
% STEP 11 */
 if abs(A(N)-X1) > abs(A(N)-X2) 
 X1 = X2;
 end;
% STEP 12 */
% accumulate shift */
 SHIFT = SHIFT+X1;
% STEP 13*/
% perform shift */
 for I = 1 : N 
 D(I) = A(I)-X1;
 end;
% STEP 14 and 15 compute R(K) */
% STEP 14 */
 X(1) = D(1);
 Y(1) = B(2);
% STEP 15 */
 for J = 2 : N 
 Z(J-1) = sqrt((X(J-1)*X(J-1))+(B(J)*B(J)));
 C(J) = X(J-1)/Z(J-1);
 S(J) = B(J)/Z(J-1);
 Q(J-1) = C(J)*Y(J-1)+S(J)*D(J);
 X(J) = C(J)*D(J)-S(J)*Y(J-1);
 if J ~= N 
 R(J-1) = S(J)*B(J+1);
 Y(J) = C(J)*B(J+1);
 end;
 end;
 M = N-1;
 MM = N-2;
% Steps 16-18 compute A(K+1) */
% STEP 16 */
 Z(N) = X(N);
 A(1) = C(2)*Z(1)+S(2)*Q(1);
 B(2) = S(2)*Z(2);
% STEP 17 */
 if N > 2 
 for J = 2 : M 
 A(J) = C(J+1)*C(J)*Z(J)+S(J+1)*Q(J);
 B(J+1) = S(J+1)*Z(J+1);
 end;
 end;
% STEP 18 */         
 A(N) = C(N)*Z(N);
 end;
 end;
 end;
% STEP 19 */
 K = K+1;
 end;
% STEP 20 */
 if OK == TRUE & K > L 
 fprintf(OUP, 'Maximum Number of Iterations Exceeded.\n');
 end;
% the process is complete */
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
 end;
