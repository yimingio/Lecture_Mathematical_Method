% CUBIC SPLINE RAYLEIGH-RITZ ALGORITHM 11.6
%
% To approximate the solution to the boundary-value problem
%
%    -D(P(X)Y')/DX + Q(X)Y = F(X), 0 <= X <= 1, Y(0)=Y(1)=0
%
% With a sum of cubic splines:
%
% INPUT:  Integer n
%
% OUTPUT: Coefficients C(0),...,C(n+1) of the basis functions
%
%   GENERAL OUTLINE
%
%       1. Nodes labelled X(I)=(I-1)*H, 1 <= I <= N+2, where
%          H=1/(N+1) so that zero subscripts are avoided
%       2. The functions PHI(I) and PHI'(I) are shifted so that
%          PHI(1) and PHI'(1) are centered at X(1), PHI(2) and PHI'(2)
%          are centered at X(2), . . . , PHI(N+2) and
%          PHI'(N+2) are centered at (X(N+2)---for example,
%                   PHI(3) = S((X-X(3))/H)
%                          = S(X/H + 2)
%       3. The functions PHI(I) are represented in terms of their
%          coefficients in the following way:
%          (PHI(I))(X) = CO(I,K,1) + CO(I,K,2)*(X-X(J)) +
%              CO(I,K,3)*(X-X(J))**2 + CO(I,K,4)*(X-X(J))**3
%          for X(J) <= X <= X(J+1) where
%          K=1 IF J=I-2, K=2 IF J=I-1, K=3 IF J=I, K=4 IF J=I+1
%          since PHI(I) is nonzero only between X(I-2) and X(I+2)
%          unless I = 1, 2, N+1 or N+2
%          (see subroutine PHICO)
%       4. The derivative of PHI(I) denoted PHI'(I) is represented
%          as in 3. By its coefficients DCO(I,K,L), L = 1, 2, 3
%          (See subroutine DPHICO).
%       5. The functions P,Q and F are represented by their cubic
%          spline interpolants using clamped boundary conditions
%          (see Algorithm 3.5).  Thus, for X(I) <= X <= X(I+1) we
%          use AF(I)+BF(I)*(X-X[I])+CF(I)*(X-X[I])^2+DF(I)*(X-X[I])^3
%          to represent F(X). Similarly, AP,BP,CP,DP are used for P
%          and AQ,BQ,CQ,DQ are used for Q.  (see subroutine COEF).
%       6. The integrands in STEPS 6 and 9 are replaced by products
%          of cubic polynomial approximations on each subinterval of
%          length H and the integrals of the resulting polynomials
%          are computed exactly.  (see subroutine XINT).
%
%
 syms('s','S','SS','FPL','FPR','PPL','PPR','QPL','QPR','OK');
 syms('N','FLAG','NAME','OUP','H','N1','N2','N3','A','C','X');
 syms('CO','DCO','AF','BF','CF','DF','AP','BP','CP','DP');
 syms('AQ','BQ','CQ','DQ','AA','BB','CC','DD','AA1','BB1');
 syms('CC1','DD1','XU1','XL1','XA','XL','XZ','XU','I','J','K');
 syms('JJ','J0','J1','J2','JJ1','JJ2','KK','E','A1','B1','C1','D1');
 syms('A2','B2','C2','D2','A3','B3','C3','D3','A4','B4','C4','D4');
 syms('ZZ1','ZZ2','K2','K3');
 TRUE = 1;
 FALSE = 0;
 fprintf(1,'This is the Cubic Spline Rayleigh-Ritz Method.\n');
 OK = FALSE;
 fprintf(1,'Input functions P(X), Q(X), and F(X) on separate lines.\n');
 fprintf(1,'For example: 1 \n');
 fprintf(1,'             pi*pi \n');
 fprintf(1,'             2*pi*pi*sin(pi*x) \n');
 s = input(' ');
 P = inline(s,'x');
 s = input(' ');
 Q = inline(s,'x');
 s = input(' ');
 F = inline(s,'x');
 fprintf(1,'Input derivative of F evaluated at 0 \n');
 FPL = input(' ');
 fprintf(1,'Input derivative of F evaluated at 1 \n');
 FPR = input(' ');
 fprintf(1,'Input derivative of Q evaluated at 0 \n');
 QPL = input(' ');
 fprintf(1,'Input derivative of Q evaluated at 1 \n');
 QPR = input(' ');
 fprintf(1,'Input derivative of P evaluated at 0 \n');
 PPL = input(' ');
 fprintf(1,'Input derivative of P evaluated at 1 \n');
 PPR = input(' ');
 while OK == FALSE
 fprintf(1,'Input positive integer n, where x(0) = 0, ');
 fprintf(1,'..., x(n+1) = 1.\n');
 N = input(' ');
 if N <= 0
 fprintf(1,'Number must be a positive integer.\n');
 else
 OK = TRUE;
 end;
 end;
 if OK == TRUE
 fprintf(1,'Choice of output method:\n');
 fprintf(1,'1. Output to screen\n');
 fprintf(1,'2. Output to text File\n');
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
 fprintf(OUP, 'CUBIC SPLINE RAYLEIGH-RITZ METHOD\n\n');
% STEP 1
 H = 1/(N+1);
 N1 = N+1;
 N2 = N+2;
 N3 = N+3;
% Initialize matrix A at zero, note that A(I, N+3) =B(I)
 A = zeros(N2,N3);
 C = zeros(1,N+2);
 X = zeros(1,N+2);
 CO = zeros(N+2,4,4);
 DCO = zeros(N+2,4,3);
 AF = zeros(1,N+1);
 BF = zeros(1,N+1);
 CF = zeros(1,N+1);
 DF = zeros(1,N+1);
 AP = zeros(1,N+1);
 BP = zeros(1,N+1);
 CP = zeros(1,N+1);
 DP = zeros(1,N+1);
 AQ = zeros(1,N+1);
 BQ = zeros(1,N+1);
 CQ = zeros(1,N+1);
 DQ = zeros(1,N+1);
 AA1 = zeros(1,N+2);
 BB1 = zeros(1,N+2);
 CC1 = zeros(1,N+2);
 DD1 = zeros(1,N+2);
 XA = zeros(1,N+2);
 XL1 = zeros(1,N+2);
 XU1 = zeros(1,N+2);
 XZ = zeros(1,N+2);
% STEP 2
% X(1)=0,...,X(I) = (I-1)*H,...,X(N+1) = 1-H, X(N+2) = 1
 for I = 1 : N2
 X(I) = (I-1)*H;
 end;
% STEPS 3 and 4 are implemented in what follows. Initialize coefficients 
% CO(I,J,K), DCO(I,J,K) */
 for I = 1 : N2 
 for J = 1 : 4 
% JJ corresponds the coefficients of phi and phi' to the proper interval 
% involving J */
 JJ = I+J-3;
 CO(I,J,1) = 0;
 CO(I,J,2) = 0;
 CO(I,J,3) = 0;
 CO(I,J,4) = 0;
 E = I-1;
 OK = TRUE;
 if JJ < I-2 | JJ >= I+2 
 OK = FALSE;
 end;
 if I == 1 & JJ < I 
 OK = FALSE;
 end;
 if I == 2 & JJ < I-1 
 OK = FALSE;
 end;
 if I == N+1 & JJ > N+1 
 OK = FALSE;
 end;
 if I == N+2 & JJ >= N+2 
 OK = FALSE;
 end;
 if OK == TRUE 
 if JJ <= I-2 
 CO(I,J,1) = (((-E+6)*E-12)*E+8)/24;
 CO(I,J,2) = ((E-4)*E+4)/(8*H);
 CO(I,J,3) = (-E+2)/(8*H^2);
 CO(I,J,4) = 1/(24*H^3);
 OK = FALSE;
 else
 if JJ > I 
 CO(I,J,1) = (((E+6)*E+12)*E+8)/24;
 CO(I,J,2) = ((-E-4)*E-4)/(8*H);
 CO(I,J,3) = (E+2)/(8*H^2);
 CO(I,J,4) = -1/(24*H^3);
 OK = FALSE;
 else
 if JJ > I-1
 CO(I,J,1) = ((-3*E-6)*E*E+4)/24;
 CO(I,J,2) = (3*E+4)*E/(8*H);
 CO(I,J,3) = (-3*E-2)/(8*H^2);
 CO(I,J,4) = 1/(8*H^3);
 if I ~= 1 & I ~= N+1 
 OK = FALSE;
 end;
 else
 CO(I,J,1) = ((3*E-6)*E*E+4)/24;
 CO(I,J,2) = (-3*E+4)*E/(8*H);
 CO(I,J,3) = (3*E-2)/(8*H^2);
 CO(I,J,4) = -1/(8*H^3);
 if I ~= 2 & I ~= N+2 
 OK = FALSE;
 end;
 end;
 end;
 end;
 end;
 if OK == TRUE 
 if I <= 2 
 AA = 1/24;
 BB = -1/(8*H);
 CC = 1/(8*H^2);
 DD = -1/(24*H^3);
 if I == 2 
 CO(I,J,1) = CO(I,J,1)-AA;
 CO(I,J,2) = CO(I,J,2)-BB;
 CO(I,J,3) = CO(I,J,3)-CC;
 CO(I,J,4) = CO(I,J,4)-DD;
 else
 CO(I,J,1) = CO(I,J,1)-4*AA;
 CO(I,J,2) = CO(I,J,2)-4*BB;
 CO(I,J,3) = CO(I,J,3)-4*CC;
 CO(I,J,4) = CO(I,J,4)-4*DD;
 end;
 else
 EE = N+2;
 AA = (((-EE+6)*EE-12)*EE+8)/24;
 BB = ((EE-4)*EE+4)/(8*H);
 CC = (-EE+2)/(8*H^2);
 DD = 1/(24*H^3);
 if I == N+1 
 CO(I,J,1) = CO(I,J,1)-AA;
 CO(I,J,2) = CO(I,J,2)-BB;
 CO(I,J,3) = CO(I,J,3)-CC;
 CO(I,J,4) = CO(I,J,4)-DD;
 else
 CO(I,J,1) = CO(I,J,1)-4*AA;
 CO(I,J,2) = CO(I,J,2)-4*BB;
 CO(I,J,3) = CO(I,J,3)-4*CC;
 CO(I,J,4) = CO(I,J,4)-4*DD;
 end;
 end;
 end;

 DCO(I,J,1) = 0;
 DCO(I,J,2) = 0;
 DCO(I,J,3) = 0;
 E = I-1;
 OK = TRUE;
 if JJ < I-2 | JJ >= I+2 
 OK = FALSE;
 end;
 if I == 1 & JJ < I 
 OK = FALSE;
 end;
 if I == 2 & JJ < I-1 
 OK = FALSE;
 end;
 if I == N+1 & JJ > N+1 
 OK = FALSE;
 end;
 if I == N+2 & JJ >= N+2 
 OK = FALSE;
 end;
 if OK == TRUE 
 if JJ <= I-2 
 DCO(I,J,1) = ((E-4)*E+4)/(8*H);
 DCO(I,J,2) = (-E+2)/(4*H^2);
 DCO(I,J,3) = 1/(8*H^3);
 OK = FALSE;
 else
 if JJ > I 
 DCO(I,J,1) = ((-E-4)*E-4)/(8*H);
 DCO(I,J,2) = (E+2)/(4*H^2);
 DCO(I,J,3) = -1/(8*H^3);
 OK = FALSE;
 else
 if JJ > I-1 
 DCO(I,J,1) = (3*E+4)*E/(8*H);
 DCO(I,J,2) = (-3.0*E-2.0)/(4.0*H^2);
 DCO(I,J,3) = 3/(8*H^3);
 if I ~= 1 & I ~= N+1 
 OK = FALSE;
 end;
 else
 DCO(I,J,1) = (-3*E+4)*E/(8*H);
 DCO(I,J,2) = (3*E-2)/(4*H^2);
 DCO(I,J,3) = -3/(8*H^3);
 if I ~= 2 & I ~= N+2 
 OK = FALSE;
 end;
 end;
 end;
 end;
 end;
 if OK == TRUE 
 if I <= 2 
 AA = -1/(8*H);
 BB = 1/(4*H^2);
 CC = -1/(8*H^3);
 if I == 2 
 DCO(I,J,1) = DCO(I,J,1)-AA;
 DCO(I,J,2) = DCO(I,J,2)-BB;
 DCO(I,J,3) = DCO(I,J,3)-CC;
 else
 DCO(I,J,1) = DCO(I,J,1)-4*AA;
 DCO(I,J,2) = DCO(I,J,2)-4*BB;
 DCO(I,J,3) = DCO(I,J,3)-4*CC;
 end;
 else
 EE = N+2;
 AA = ((EE-4)*EE+4)/(8*H);
 BB = (-EE+2)/(4*H^2);
 CC = 1/(8*H^3);
 if I == N+1 
 DCO(I,J,1) = DCO(I,J,1)-AA;
 DCO(I,J,2) = DCO(I,J,2)-BB;
 DCO(I,J,3) = DCO(I,J,3)-CC;
 else
 DCO(I,J,1) = DCO(I,J,1)-4*AA;
 DCO(I,J,2) = DCO(I,J,2)-4*BB;
 DCO(I,J,3) = DCO(I,J,3)-4*CC;
 end;
 end;
 end;
 end;
 end;
% Output the basis functions. */
 fprintf(OUP, 'Basis Function: A + B*X + C*X**2 + D*X**3\n\n');
 fprintf(OUP, '  A  B  C  D\n\n');
 for I = 1 : N2 
 fprintf(OUP, 'phi( %d )\n\n', I);
 for J = 1 : 4 
 if I ~= 1 | (J ~= 1 & J ~= 2) 
 if I ~= 2 | J ~= 2 
 if I ~= N1 | J ~= 4 
 if I ~= N2 | (J ~= 3 & J ~= 4) 
 JJ1 = I+J-3;
 JJ2 = I+J-2;
 fprintf(OUP, 'On (X( %d ), X( %d )) ', JJ1, JJ2);
 for K = 1 : 4
 fprintf(OUP, ' %12.8f ', CO(I,J,K));
 end;
 fprintf(OUP, '\n');
 end;
 end;
 end;
 end;
 end;
 end;
% Obtain coefficients for F, P, Q
 for I = 1 : N2
 AA1(I) = F(X(I));
 end;
 XA(1) = 3.0*(AA1(2)-AA1(1))/H-3.0*FPL;
 XA(N2) = 3.0*FPR-3.0*(AA1(N2)-AA1(N2-1))/H;
 XL1(1) = 2.0*H;
 XU1(1) = 0.5;
 XZ(1) = XA(1)/XL1(1);
 for I = 2 : N1 
 XA(I) = 3.0*(AA1(I+1)-2.0*AA1(I)+AA1(I-1))/H;
 XL1(I) = H*(4.0-XU1(I-1));
 XU1(I) = H/XL1(I);
 XZ(I) = (XA(I)-H*XZ(I-1))/XL1(I);
 end;
 XL1(N2) = H*(2.0-XU1(N2-1));
 XZ(N2) = (XA(N2)-H*XZ(N2-1))/XL1(N2);
 CC1(N2) = XZ(N2);
 for I = 1 : N1 
 J = N2-I;
 CC1(J) = XZ(J)-XU1(J)*CC1(J+1);
 BB1(J) = (AA1(J+1)-AA1(J))/H-H*(CC1(J+1)+2.0*CC1(J))/3.0;
 DD1(J) = (CC1(J+1)-CC1(J))/(3.0*H);
 end;
 for I = 1 : N1
 AF(I) = ((-DD1(I)*X(I)+CC1(I))*X(I)-BB1(I))*X(I)+AA1(I);
 BF(I) = (3.0*DD1(I)*X(I)-2.0*CC1(I))*X(I)+BB1(I);
 CF(I) = CC1(I)-3.0*DD1(I)*X(I);
 DF(I) = DD1(I);
 end;
 for I = 1 : N2 
 AA1(I) = P(X(I));
 end;
 XA(1) = 3.0*(AA1(2)-AA1(1))/H-3.0*PPL;
 XA(N2) = 3.0*PPR-3.0*(AA1(N2)-AA1(N2-1))/H;
 XL1(1) = 2.0*H;
 XU1(1) = 0.5;
 XZ(1) = XA(1)/XL1(1);
 for I = 2 : N1 
 XA(I) = 3.0*(AA1(I+1)-2.0*AA1(I)+AA1(I-1))/H;
 XL1(I) = H*(4.0-XU1(I-1));
 XU1(I) = H/XL1(I);
 XZ(I) = (XA(I)-H*XZ(I-1))/XL1(I);
 end;
 XL1(N2) = H*(2.0-XU1(N2-1));
 XZ(N2) = (XA(N2)-H*XZ(N2-1))/XL1(N2);
 CC1(N2) = XZ(N2);
 for I = 1 : N1
 J = N2-I;
 CC1(J) = XZ(J)-XU1(J)*CC1(J+1);
 BB1(J) = (AA1(J+1)-AA1(J))/H -H*(CC1(J+1)+2.0*CC1(J))/3.0;
 DD1(J) = (CC1(J+1)-CC1(J))/(3.0*H);
 end;
 for I = 1 : N1 
 AP(I) = ((-DD1(I)*X(I)+CC1(I))*X(I)-BB1(I))*X(I)+AA1(I);
 BP(I) = (3.0*DD1(I)*X(I)-2.0*CC1(I))*X(I)+BB1(I);
 CP(I) = CC1(I)-3.0*DD1(I)*X(I);
 DP(I) = DD1(I);
 end;
 for I = 1 : N2 
 AA1(I) = Q(X(I));
 end;
 XA(1) = 3.0*(AA1(2)-AA1(1))/H-3.0*QPL;
 XA(N2) = 3.0*QPR-3.0*(AA1(N2)-AA1(N2-1))/H;
 XL1(1) = 2.0*H;
 XU1(1) = 0.5;
 XZ(1) = XA(1)/XL1(1);
 for I = 2 : N1 
 XA(I) = 3.0*(AA1(I+1)-2.0*AA1(I)+AA1(I-1))/H;
 XL1(I) = H*(4.0-XU1(I-1));
 XU1(I) = H/XL1(I);
 XZ(I) = (XA(I)-H*XZ(I-1))/XL1(I);
 end;
 XL1(N2) = H*(2.0-XU1(N2-1));
 XZ(N2) = (XA(N2)-H*XZ(N2-1))/XL1(N2);
 CC1(N2) = XZ(N2);
 for I = 1 : N1 
 J = N2-I;
 CC1(J) = XZ(J)-XU1(J)*CC1(J+1);
 BB1(J) = (AA1(J+1)-AA1(J))/H -H*(CC1(J+1)+2.0*CC1(J))/3.0;
 DD1(J) = (CC1(J+1)-CC1(J))/(3.0*H);
 end;
 for I = 1 : N1 
 AQ(I) = ((-DD1(I)*X(I)+CC1(I))*X(I)-BB1(I))*X(I)+AA1(I);
 BQ(I) = (3.0*DD1(I)*X(I)-2.0*CC1(I))*X(I)+BB1(I);
 CQ(I) = CC1(I)-3.0*DD1(I)*X(I);
 DQ(I) = DD1(I);
 end;
% STEPS 5-9 are implemented in what follows
 for I = 1 : N2 
% indices for limits of integration for A(I,J) and B(I)
 J1 = min(I+2,N+2);
 J0 = max(I-2,1);
 J2 = J1-1;
% integrate over each subinterval where phi(I) is nonzero
 for JJ = J0 : J2 
% limits of integration for each call
 XU = X(JJ+1);
 XL = X(JJ);
% coefficients of bases
 K = INTE(I,JJ);
 A1 = DCO(I,K,1);
 B1 = DCO(I,K,2);
 C1 = DCO(I,K,3);
 D1 = 0;
 A2 = CO(I,K,1);
 B2 = CO(I,K,2);
 C2 = CO(I,K,3);
 D2 = CO(I,K,4);
% call subprogram for integrations
 ZZ1 = XINT(XU,XL,AP(JJ),BP(JJ),CP(JJ),DP(JJ),A1,B1,C1,D1,A1,B1,C1,D1);
 ZZ2 = XINT(XU,XL,AQ(JJ),BQ(JJ),CQ(JJ),DQ(JJ),A2,B2,C2,D2,A2,B2,C2,D2);
 A(I,I) = A(I,I) + ZZ1 + ZZ2;
 ZZ1 = XINT(XU,XL,AF(JJ),BF(JJ),CF(JJ),DF(JJ),A2,B2,C2,D2,1,0,0,0);
 A(I,N+3) = A(I,N+3) + ZZ1;
 end;
% compute A(I,J) for J = I+1,...,Min(I+3,N+2)
 K3 = I+1;
 if K3 <= N2 
 K2 = min(I+3,N+2);
 for J = K3 : K2 
 J0 = max(J-2,1);
 for JJ = J0 : J2 
 XU = X(JJ+1);
 XL = X(JJ);
 K = INTE(I,JJ);
 A1 = DCO(I,K,1);
 B1 = DCO(I,K,2);
 C1 = DCO(I,K,3);
 D1 = 0;
 A2 = CO(I,K,1);
 B2 = CO(I,K,2);
 C2 = CO(I,K,3);
 D2 = CO(I,K,4);
 K = INTE(J,JJ);
 A3 = DCO(J,K,1);
 B3 = DCO(J,K,2);
 C3 = DCO(J,K,3);
 D3 = 0;
 A4 = CO(J,K,1);
 B4 = CO(J,K,2);
 C4 = CO(J,K,3);
 D4 = CO(J,K,4);
 ZZ1 = XINT(XU,XL,AP(JJ),BP(JJ),CP(JJ),DP(JJ),A1,B1,C1,D1,A3,B3,C3,D3);
 ZZ2 = XINT(XU,XL,AQ(JJ),BQ(JJ),CQ(JJ),DQ(JJ),A2,B2,C2,D2,A4,B4,C4,D4);
 A(I,J) = A(I,J) + ZZ1 + ZZ2;
 end;
 A(J,I) = A(I,J);
 end;
 end;
 end;
% STEP 10
 for I = 1 : N1 
 for J = I+1 : N2 
 CC = A(J,I)/A(I,I);
 for K = I+1 : N3 
 A(J,K) = A(J,K)-CC*A(I,K);
 end;
 A(J,I) = 0;
 end;
 end;
 C(N2) = A(N2,N3)/A(N2,N2);
 for I = 1 : N1 
 J = N1-I+1;
 C(J) = A(J,N3);
 for KK = J+1 : N2 
 C(J) = C(J)-A(J,KK)*C(KK);
 end;
 C(J) = C(J)/A(J,J);
 end;
% STEP 11
% output coefficients
 fprintf(OUP, '\nCoefficients:  c(1), c(2), ... , c(n+1)\n\n');
 for I = 1 : N1
 fprintf(OUP, '  %12.6e \n', C(I));
 end;
 fprintf(OUP, '\n');
% compute and output value of the approximation at the nodes
 fprintf(OUP, 'The approximation evaluated at the nodes:\n\n');
 fprintf(OUP, ' Node          Value\n\n');
 for I = 1 : N2
 S = 0;
 for J = 1 : N2 
 J0 = max(J-2,1);
 J1 = min(J+2,N+2);
 SS = 0;
 if I < J0 | I >= J1 
 S = S + C(J)*SS;
 else
 K = INTE(J,I);
 SS = ((CO(J,K,4)*X(I)+CO(J,K,3))*X(I)+CO(J,K,2))*X(I)+CO(J,K,1);
 S = S + C(J)*SS;
 end;
 end;
 fprintf(OUP, '%12.8f %12.8f\n', X(I), S);
 end;
 end;
 if OUP ~= 1 
 fclose(OUP);
 fprintf(1,'Output file %s created successfully \n',NAME);
 end;
