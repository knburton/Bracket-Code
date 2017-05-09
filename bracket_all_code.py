#Kathryn smells
#Danny does too

N=Rational(7)

T=str()
for i in range(1,N+1):
    for j in range(1,N+1):
        T=T+str('x')+str(i)+str('_')+str(j)+str(',')
T=T[:-1]

Q = PolynomialRing(QQ,T)
Q.inject_variables();
x=Q.gens()
R=matrix(QQ,N**2)
mu=Rational(2)
nu=Rational(1)
lam=1/N
X=matrix(Q,N)
XX=matrix(Q,N**2)

def FillMatEl(i,j,k,m,f,W):
    ii = N*(i-1)+k;
    jj = N*(j-1)+m;
    W[ii-1,jj-1]=f;
    
def R_stand():
    for i in range(1,N):
        for j in range(i+1,N+1):
            FillMatEl(i,j,j,i,nu,R)
            FillMatEl(j,i,i,j,-nu,R);
    
def R_CG():
    for i in range(1,N):
        for j in range(i+1,N+1):
            for m in range(1,j-i):
                FillMatEl(i,j-m,j,i+m,mu,R)
                FillMatEl(j,i+m,i,j-m,-mu,R);
                
                
def R0():
    for i in range(1,N):
        for j in range(i+1,N+1):
            FillMatEl(i,i,j,j,(N+2*(i-j))*lam,R)
            FillMatEl(j,j,i,i,-(N+2*(i-j))*lam,R);
def Xmat():
    for i in range(N):
        for j in range(N):
            X[i,j]= Q.gens()[N*i+j];
            
def XD():
    for i in range(1,N+1):
        for j in range(1,N+1):
            for k in range(1,N+1):
                for l in range(1,N+1):
                    FillMatEl(i,j,k,l,X[i-1,j-1]*X[k-1,l-1],XX);

R_stand();
R_CG();
R0();
Xmat();
XD();
ANS=R*XX-XX*R;


def linear_bracket(m,n):
    for i in range(len(Q.gens())):
        for j in range(len(Q.gens())):
            if m == Q.gens()[i]:
                if n == Q.gens()[j]:
                    return ANS[N*floor(i*lam)+floor(j*lam),N*(i-N*floor(i*lam))+(j-N*floor(j*lam))]

def monomial_bracket(m,n):
    return sum([diff(m,x[j])*sum([diff(n,x[i])*linear_bracket(x[j],x[i]) for i in range(len(x))]) for j in range(len(x))])

def bracket(f,g):
        return sum([sum([f.coefficients()[j]*g.coefficients()[i]*monomial_bracket(f.monomials()[j],g.monomials()[i]) for i in range(len(g.monomials()))]) for j in range(len(f.monomials()))])
    


U = matrix(Q,24,25)

for i in range(5):
    for j in range(7):
        U[i,j] = X[i,j]
for i in range(5,9):
    for j in range(7):
        U[i,j+1] = X[i-5,j]
    for j in range(8,15):
        U[i,j] = X.matrix_from_rows([1,2,3,4])[i-5,j-8]
for i in range(9,14):
    for j in range(9,16):
        U[i,j] = X[i-9,j-9]
for i in range(14,19):
    for j in range(10,17):
        U[i,j] = X[i-14,j-10]
    for j in range(17,24):
        U[i,j] = X[i-13,j-17]
for i in range(19,24):
    for j in range(18,25):
        U[i,j] = X[i-19,j-18]    
    
#numX will have to be changed for different N
#numX = matrix(QQ,[[1,-1,-4,-3,4,-4,2],[3,-3,1,-5,0,0,2],[2,-2,-5,-1,-3,1,-3],[-4,-1,3,-3,2,-1,-4],[-4,-3,1,-3,0,-5,3],[-2,-2,0,3,2,2,4],[-2,1,-4,-2,4,4,-3]]) 
numX = matrix(QQ,[[5,23,-1,-1,-8,43,-3],[-1,-1,2,-1,1,-1,4],[1,2,9,2,3,1,7],[-6,2,1,-3,-1,-8,3],[3,1,-1,1,1,1,81],[-25,1,-16,-1,5,1,-4],[1,-2,23,-1,1,4,1]])
        
numU = matrix(Q,24,25)

for i in range(5):
    for j in range(7):
        numU[i,j] = numX[i,j]
for i in range(5,9):
    for j in range(7):
        numU[i,j+1] = numX[i-5,j]
    for j in range(8,15):
        numU[i,j] = numX.matrix_from_rows([1,2,3,4])[i-5,j-8]
for i in range(9,14):
    for j in range(9,16):
        numU[i,j] = numX[i-9,j-9]
for i in range(14,19):
    for j in range(10,17):
        numU[i,j] = numX[i-14,j-10]
    for j in range(17,24):
        numU[i,j] = numX[i-13,j-17]
for i in range(19,24):
    for j in range(18,25):
        numU[i,j] = numX[i-19,j-18]






ma=[0]
for i in range(1,9):
    ma.append(U.matrix_from_rows_and_columns(range(9-i,9),range(8-i,8)));

mb=[0]
for i in range(1,10):
    mb.append(U.matrix_from_rows_and_columns(range(9-i,9),range(9-i,9)));

mc=[0]
for i in range(1,10):
    mc.append(U.matrix_from_rows_and_columns(range(9-i,24),range(10-i,25)));

md=[0]
for i in range(1,10):
    md.append(U.matrix_from_rows_and_columns(range(14-i,14),range(16-i,16)));

cX = matrix(Q,6,8);

for i in range(6):
    for j in range(7):
        cX[i,j] = X[i+1,j]

cY = matrix(Q,6,8);

for i in range(6):
    for j in range(1,8):
        cY[i,j] = X[i,j-1]



Chi = matrix(Q,N-1,N+1)
for i in range(N-1):
    for j in range(N):
        Chi[i,j] = X[i+1,j]

Upsilon = matrix(Q,N-1,N+1)
for i  in range(N-1):
    for j in range(N):
        Upsilon[i,j+1] = X[i,j]

pU = matrix(Q,floor((N+1)/2)*(N-1),(floor((N+1)/2)+1)*(N+1))
for k in range(floor((N+1)/2)):
    for i in range(N-1):
        for j in range(N+1):
            pU[i+k*(N-1),j+k*(N+1)] = Upsilon[i,j]
            pU[i+k*(N-1),j+N+1+k*(N+1)] = Chi[i,j]

ncX = matrix(Q,6,8);

for i in range(6):
    for j in range(7):
        ncX[i,j] = numX[i+1,j]

ncY = matrix(Q,6,8);

for i in range(6):
    for j in range(1,8):
        ncY[i,j] = numX[i,j-1]
            
npU = matrix(Q,24,40);

for k in range(4):
    for i in range(6):
        for j in range(8):
            npU[(k)*(N-1)+i,(k)*(N+1)+j]=ncY[i,j]
            npU[(k)*(N-1)+i,(k+1)*(N+1)+j]=ncX[i,j]
            
mtheta = [0];

for i in range(1,N):
    mtheta.append(X.matrix_from_rows_and_columns(range(N-i,N),range(N-i,N)));
            
            
mphi = [0];

for i in range(1,floor((N+1)/2)*(N-1)+1):
    mphi.append(pU.matrix_from_rows_and_columns(range(floor((N+1)/2)*(N-1)-i,floor((N+1)/2)*(N-1)),range(floor((N+1)/2)*(N+1)-i,floor((N+1)/2)*(N+1))));
    
mpsi = [0];

if N%2 == 0:
    M = floor((N+1)/2)*(N-1)
else:
    M = floor((N+1)/2)*(N-1)-N+1
for i in range(1,M+1):
    mpsi.append(pU.matrix_from_rows_and_columns(range(floor((N+1)/2)*(N-1)-i,floor((N+1)/2)*(N-1)),range(floor((N+1)/2)*(N+1)-i+1,floor((N+1)/2)*(N+1)+1)));

xgens=matrix(Q,1,len(Q.gens()))

def xg():
    for i in range(len(Q.gens())):
        xgens[0,i]=Q.gens()[i];
    
xg();
rXX=matrix(QQ,N**2)


def ranXX(A):
    for i in range(1,N+1):
        for j in range(1,N+1):
            for k in range(1,N+1):
                for l in range(1,N+1):
                    FillMatEl(i,j,k,l,A[i-1,j-1]*A[k-1,l-1],rXX);
                    ANS2 = R*rXX-rXX*R;
    return ANS2

NXX=ranXX(numX);

def Num_linear_bracket(m,n): # m and n are generators of Q, [removed:A random matrix(probably just use numX)]
    for i in range(len(Q.gens())):
        for j in range(len(Q.gens())):
            if m == Q.gens()[i]:
                if n == Q.gens()[j]:
                    return NXX[N*floor(i*lam)+floor(j*lam),N*(i-N*floor(i*lam))+(j-N*floor(j*lam))]
                 
    
                   

        
def Num_Pder(A,B,x): # A is in Mat(Q,n), B is in Mat(ZZ,n), x is a generator of Q
    s=0
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if A[i,j]==x:
                s += (-1)**(i+j)*det(B.matrix_from_rows_and_columns(range(i)+range(i+1,A.nrows()),range(j)+range(j+1,A.ncols())))
    return s;



subs_dict = {}
for i in range(len(Q.gens())):
    subs_dict[Q.gens()[i]] = numX[(i-i%N)/N,i%N]

def Num_bracket(A,B): # A is in Mat(Q,n), B is in Mat(Q,m),[Removed: C is in Mat(ZZ,n), D is in Mat(ZZ,m), X is the random matrix from which A,B,C,D should be created, It should be noted that this will run with A,B,C,D not consistant with X or one another but will not give the correct answer.]
    p=0
    C = A.subs(subs_dict)
    D = B.subs(subs_dict)
    for i in range(len(Q.gens())):
        factor1 = Num_Pder(A,C,Q.gens()[i])
        for j in range(len(Q.gens())):
            bracket = Num_linear_bracket(Q.gens()[i],Q.gens()[j])
            if bracket != 0:
                p += factor1*Num_Pder(B,D,Q.gens()[j])*bracket;
    return p
        
def Num_log_check(A,B): # A is in Mat(Q,n), B is in Mat(Q,m),[Removed: C is in Mat(ZZ,n), D is in Mat(ZZ,m), X is the random matrix from which A,B,C,D should be created, It should be noted that this will run with A,B,C,D not consistant with X or one another but will not give the correct answer.]
    p=0
    C = A.subs(subs_dict)
    D = B.subs(subs_dict)
    for i in range(len(Q.gens())):
        factor1 = Num_Pder(A,C,Q.gens()[i])
        if factor1 != 0:
            for j in range(len(Q.gens())):
                bracket = Num_linear_bracket(Q.gens()[i],Q.gens()[j])
                if bracket != 0:
                    p += factor1*Num_Pder(B,D,Q.gens()[j])*bracket;
    return p/(det(C)*det(D))
                   
mishaC=[]

for i in range(9):
    mishaC.append(mc[9].matrix_from_rows_and_columns(range(8-i,14),range(8-i,8)+range(9,15)));
          
