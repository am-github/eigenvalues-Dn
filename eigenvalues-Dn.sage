d7=load("d7.sobj")
d5=load("d5.sobj")
d3=load("d3.sobj")
d1=load("d1.sobj")
d83=load("d83.sobj")

def rep(n,x):
    if n==7:
        return d7[x]
    if n==5:
    	return d5[x]
    if n==3:
    	return d3[x]
    if n==1:
    	return d1[x]

def repminodd(k,D,r):
    if not (-2*D-k/4).is_integer():
        raise ValueError('Incorrect value of D!')  
    sm2 = ''
    bd = '-2*x-k/4'
    sm = ''
    pr = ''
    fnc = 's=0; x='+str(D)+'; k ='+str(k)+'; r='+str(r)+';'
    for idx in range(k):
        fnc += 'for(i%s = 0, sqrt('%idx + bd + '),'
        bd += ' - i%s^2 -i%s'%(idx,idx)
        fnc += 'if('+bd+'<0, break);'
        sm2 += '+ i%s^2 + i%s'%(idx,idx)
        sm += '+ i%s'%idx
        pr += '* (i%s+1/2)'%idx
    fnc += 'if('+sm2[1:]+'== -2*x-k/4,p=0;'
    fnc += 's += 2^k*(-1)^(r+k/2'+sm+')'+pr+';)'+')'*k
    gp(fnc)
    rk = gp('s')
    return ZZ(rk)

def CE4(n,D): #the Fourier coefficients C(D,r) of the Eisenstein series of weight 4 for D_n for any r in Z
	return rep(8-n,-2*D)

def CTE4(n,l,D): #the Fourier coefficients C(D,r) of T(l)E_{4,D_n,0} for any r in Z and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum+=kronecker((-1)**((n-1)/2)*D*8,i)*CE4(n,D*l**2/i**2)*i**(4-1-(n+1)/2)
    return sum

def lE4(n,l,D): #the eigenvalue of E_{4,D_n,0} for T(l)
    return CTE4(n,l,D)/CE4(n,D)

def CE6(n,D): #the Fourier coefficients C(D,r) of the Eisenstein series of weight 6 for D_n for any r in Z
    sum =-(D+(8-n)/24)*CE4(n,D)
    for i in xrange(1,-D+1):
        sum=sum+(8-n)*sigma(i,1)*CE4(n,D+i)
    return sum

def CTE6(n,l,D): #the Fourier coefficients C(D,r) of T(l)E_{6,D_n,0} for any r in Z and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*CE6(n,D*l**2/i**2)*i**(6-1-(n+1)/2)
    return sum

def lE6(n,l,D): #the eigenvalue of E_{6,D_n,0} for T(l)
    return CTE6(n,l,D)/CE6(n,D)

def CE8(n,D): #the Fourier coefficients C(D,r) of E_{8,D_n,0} for any r in Z and D coprime to l
    sum =-(D+(12-n)/24)*CE6(n,D)
    for i in xrange(1,-D+1):
        sum=sum+(12-n)*sigma(i,1)*CE6(n,D+i)
    return sum

def C8(n,D):  #the Fourier coefficients C(D,r) of the cusp form of weight 8 for the lattice D_n for any r in Z^n if n!=3 and r=0 otherwise
    if n==3:
        return d83[-D]
    sum = CE4(n,D)-CE8(n,D)*24**2/((8-n)*(12-n))
    for i in range (1,-D+1):
            sum=sum+240*sigma(i,3)*CE4(n,D+i)
    return sum

def CT8(n,l,D): #the Fourier coefficients C(D,r) of T(l)phi_8 for any D coprime to l if n!=3 and C(D,0) otherwise
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*C8(n,D*l**2/i**2)*i**(8-1-(n+1)/2)
    return sum

def l8(n,l,D): #the eigenvalue of phi_8 for T(l)
    return CT8(n,l,D)/C8(n,D)

def C10(n,D):  #the Fourier coefficients C(D,r) of the cusp form of weight 10 for the lattice D_n for any r in Z^n
    sum = CE4(n,D)*(8-n)/24+CE6(n,D)
    for i in range (1,-D+1):
            sum=sum-504*sigma(i,5)*CE4(n,D+i)*(8-n)/24+240*sigma(i,3)*CE6(n,D+i)
    return sum

def CT10(n,l,D): #the Fourier coefficients C(D,r) of T(l)phi_{10} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*C10(n,D*l**2/i**2)*i**(10-1-(n+1)/2)
    return sum

def l10(n,l,D): #the eigenvalue of phi_{10} for T(l)
    return CT10(n,l,D)/C10(n,D)

def Cbeta1(D): #the Fourier coefficients C(D,r) of the cusp form beta_1 of weight 12 for D_3 for any r in Z^3
    sum=-D*CE4(3,D)
    for i in xrange(1,-D+1):
        sum=sum+(5*sigma(i,1)+100*sigma(i,7)+504*sigma(i,5)*(D+i+5/24))*CE4(3,D+i)
        for j in xrange(1,-D-i+1):
            sum=sum-2520*sigma(i,5)*sigma(j,1)*CE4(3,D+i+j)
    return sum

def CTbeta1(l,D): #the Fourier coefficients C(D,r) of T(l)beta_1 for any r in Z^3 and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*(i**9)*Cbeta1(D*l**2/i**2)
    return sum

def Cbeta2(D): #the Fourier coefficients C(D,0) of the cusp form beta_2 of weight 12 for D_3
    sum=C8(3,D)
    for i in xrange(1,-D):
        sum=sum+240*sigma(i,3)*C8(3,D+i)
    return sum

def CTbeta2(l,D): #the Fourier coefficients C(D,0) of T(l)beta_2 for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*Cbeta2(D*l**2/i**2)*i**9
    return sum

def T3(l): #the matrix of T(l) on J_{12,D_3}
    G = matrix([[Cbeta1(-1),Cbeta2(-1)],[Cbeta1(-2),Cbeta2(-2)]])
    Vbeta1 = matrix([[CTbeta1(l,-1)],[CTbeta1(l,-2)]])
    Gbeta1 = G.inverse()*Vbeta1
    Vbeta2 = matrix([[CTbeta2(l,-1)],[CTbeta2(l,-2)]])
    Gbeta2 = G.inverse()*Vbeta2
    return matrix([[Gbeta1[0,0],Gbeta1[1,0]],[Gbeta2[0,0],Gbeta2[1,0]]])

def CPsi12(D): # Fourier coefficients C(D,0) of the first Hecke eigenform in J_{12,D_3}, computed by diagonalizing T(3)
    return 2/9*Cbeta1(D)+35/9*Cbeta2(D)

def CTPsi12(l,D): #the Fourier coefficients C(D,0) of T(l)Psi_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*CPsi12(D*l**2/i**2)*i**9
    return sum

def lPsi12(l,D): #the eigenvalue of Psi_{12} for T(l)
    return CTPsi12(l,D)/CPsi12(D)

def Cvp12(D): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{12,D_3}, computed by diagonalizing T(3)
    return 7/9*Cbeta1(D)-35/9*Cbeta2(D)

def CTvp12(l,D): #the Fourier coefficients C(D,0) of T(l)varphi_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*Cvp12(D*l**2/i**2)*i**9
    return sum

def lvp12(l,D): #the eigenvalue of varphi_{12} for T(l)
    return CTvp12(l,D)/Cvp12(D)

def Cminodd(n,D,r): #the Fourier coefficients C(D,r) of the cusp form of weight 12-n for the lattice D_n for any r in (1/2+Z)^n
    return repminodd(8-n,D,r)

def CTminodd(n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{12-n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and (8*D)%l==0:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd(n,D*l**2/i**2,r*l/i)*i**(12-n-1-(n+1)/2)
    return sum

def lminodd(n,l,D,r): #the eigenvalue of phi_{12-n} for T(l)
    return CTminodd(n,l,D,r)/Cminodd(n,D,r)

def Cminodd4(n,D,r): #the Fourier coefficients C(D,r) of the cusp form of weight 16-n for the lattice D_n for any r in (1/2+Z)^n
    sum=Cminodd(n,D,r)
    for i in xrange(1,-D+1):
        sum = sum+240*Cminodd(n,D+i,r)*sigma(i,3)
    return sum

def CTminodd4(n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{16-n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd4(n,D*l**2/i**2,r*l/i)*i**(16-n-1-(n+1)/2)
    return sum

def lminodd4(n,l,D,r): #the eigenvalue of phi_{16-n} for T(l)
    return CTminodd4(n,l,D,r)/Cminodd4(n,D,r)

def Cminodd6(n,D,r): #the Fourier coefficients C(D,r) of the cusp form of weight 18-n for the lattice D_n for any r in (1/2+Z)^n
    sum=Cminodd(n,D,r)
    for i in xrange(1,-D+1):
        sum = sum-504*Cminodd(n,D+i,r)*sigma(i,5)
    return sum

def CTminodd6(n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{18-n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd6(n,D*l**2/i**2,r*l/i)*i**(18-n-1-(n+1)/2)
    return sum

def lminodd6(n,l,D,r): #the eigenvalue of phi_{18-n} for T(l)
    return CTminodd6(n,l,D,r)/Cminodd6(n,D,r)

def Cminodd8(n,D,r): #the Fourier coefficients C(D,r) of the cusp form of weight 20-n for the lattice D_n for any r in (1/2+Z)^n
    sum=Cminodd(n,D,r)
    for i in xrange(1,-D+1):
        sum = sum+480*Cminodd(n,D+i,r)*sigma(i,7)
    return sum

def CTminodd8(n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{20-n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd8(n,D*l**2/i**2,r*l/i)*i**(20-n-1-(n+1)/2)
    return sum

def lminodd8(n,l,D,r): #the eigenvalue of phi_{20-n} for T(l)
    return CTminodd8(n,l,D,r)/Cminodd8(n,D,r)

def Cminodd10(n,D,r): #the Fourier coefficients C(D,r) of the cusp form of weight 22-n for the lattice D_n for any r in (1/2+Z)^n
    sum=Cminodd(n,D,r)
    for i in xrange(1,-D+1):
        sum = sum-264*Cminodd(n,D+i,r)*sigma(i,9)
    return sum

def CTminodd10(n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{22-n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd10(n,D*l**2/i**2,r*l/i)*i**(22-n-1-(n+1)/2)
    return sum

def lminodd10(n,l,D,r): #the eigenvalue of phi_{22-n} for T(l)
    return CTminodd10(n,l,D,r)/Cminodd10(n,D,r)
