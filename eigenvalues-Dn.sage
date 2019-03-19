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

def CE4(n,D): #the Fourier coefficients C(D,r) of the Eisenstein series of weight 4 for D_n for any r in Z^n
	return rep(8-n,-2*D)

def CTE4(n,l,D): #the Fourier coefficients C(D,r) of T(l)E_{4,D_n,0} for any r in Z^n and D coprime to l
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

def CE6(n,D): #the Fourier coefficients C(D,r) of the Eisenstein series of weight 6 for D_n for any r in Z^n
    sum =-(D+(8-n)/24)*CE4(n,D)
    for i in xrange(1,-D+1):
        sum=sum+(8-n)*sigma(i,1)*CE4(n,D+i)
    return sum

def CTE6(n,l,D): #the Fourier coefficients C(D,r) of T(l)E_{6,D_n,0} for any r in Z^n and D coprime to l
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

def CEis8(n,D): #the Fourier coefficients C(D,r) of E_{8,D_n,0} for any r in Z^n and D coprime to l
    sum =-(D+(12-n)/24)*CE6(n,D)
    for i in xrange(1,-D+1):
        sum=sum+(12-n)*sigma(i,1)*CE6(n,D+i)
    return sum

def C8(n,D):  #the Fourier coefficients C(D,r) of the cusp form of weight 8 for the lattice D_n for any r in Z^n if n!=3 and r=0 otherwise
    if n==3:
        return d83[-D]
    sum = CE4(n,D)-CEis8(n,D)*24**2/((8-n)*(12-n))
    for i in range (1,-D+1):
            sum=sum+240*sigma(i,3)*CE4(n,D+i)
    return sum

def CT8(n,l,D): #the Fourier coefficients C(D,r) of T(l)phi_8 for any r in Z^n and any D coprime to l if n!=3 and C(D,0) otherwise
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

def CB121(n,D): #the Fourier coefficients C(D,r) of the cusp form B_{12,n}^1 of weight 12 for D_n for any r in Z
    sum=CE4(n,D)-CEis8(n,D)*576/((12-n)*(8-n))
    for i in xrange(1,-D+1):
        sum=sum+480*CE4(n,D+i)*sigma(i,7)-240*CEis8(n,D+i)*sigma(i,3)*576/((12-n)*(8-n))
    return sum

def CB122(n,D): #the Fourier coefficients C(D,r) of the cusp form B_{12,n}^2 of weight 12 for D_1 for any r in Z
    sum=CE6(n,D)+CEis8(n,D)*24/(12-n)
    for i in xrange(1,-D+1):
        sum=sum-504*CE6(n,D+i)*sigma(i,5)+240*CEis8(n,D+i)*sigma(i,3)*24/(12-n)
    return sum

def CTB121(n,l,D): #the Fourier coefficients C(D,r) of T(l)B_{12,n}^1 for any r in Z and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*CB121(n,D*l**2/i**2)*i**(11-(n+1)/2)
    return sum

def CTB122(n,l,D): #the Fourier coefficients C(D,r) of T(l)B_{12,n}^2 for any r in Z and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*CB122(n,D*l**2/i**2)*i**(11-(n+1)/2)
    return sum

def T12(n,l): #the matrix of T(l) on B121(n) and B122(n)
    G = matrix([[CB121(n,-1),CB122(n,-1)],[CB121(n,-2),CB122(n,-2)]])
    VB121 = matrix([[CTB121(n,l,-1)],[CTB121(n,l,-2)]])
    GB121 = G.inverse()*VB121
    VB122 = matrix([[CTB122(n,l,-1)],[CTB122(n,l,-2)]])
    GB122 = G.inverse()*VB122
    return matrix([[GB121[0,0],GB121[1,0]],[GB122[0,0],GB122[1,0]]])

def Cpsi12(D): # Fourier coefficients C(D,0) of the first Hecke eigenform in J_{12,D_1}, computed by diagonalizing T(3)
    return -CB121(1,D)*6/319+CB122(1,D)*1440/2233

def CTpsi12(l,D): #the Fourier coefficients C(D,0) of T(l)psi_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+Cpsi12(D*l**2/i**2)*kronecker((-1)**((1-1)/2)*D*8,i)*i**(11-(1+1)/2)
    return sum

def lpsi12(l,D): #the eigenvalue of psi_{12} for T(l)
    return CTpsi12(l,D)/Cpsi12(D)

def Calp12(D): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{12,D_1}, computed by diagonalizing T(3)
    return CB121(1,D)*325/319-CB122(1,D)*1440/2233

def CTalp12(l,D): #the Fourier coefficients C(D,0) of T(l)alpha_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+Calp12(D*l**2/i**2)*kronecker((-1)**((1-1)/2)*D*8,i)*i**(11-(1+1)/2)
    return sum

def lalp12(l,D): #the eigenvalue of alpha_{12} for T(l)
    return CTalp12(l,D)/Calp12(D)

def CB1231(D): #the Fourier coefficients C(D,r) of the cusp form B_{12,3}^1 of weight 12 for D_3 for any r in Z^3
    sum=-D*CE4(3,D)
    for i in xrange(1,-D+1):
        sum=sum+(5*sigma(i,1)+100*sigma(i,7)+504*sigma(i,5)*(D+i+5/24))*CE4(3,D+i)
        for j in xrange(1,-D-i+1):
            sum=sum-2520*sigma(i,5)*sigma(j,1)*CE4(3,D+i+j)
    return sum

def CB1232(D): #the Fourier coefficients C(D,0) of the cusp form B_{12,3}^2 of weight 12 for D_3
    sum=C8(3,D)
    for i in xrange(1,-D):
        sum=sum+240*sigma(i,3)*C8(3,D+i)
    return sum

def CTB1231(l,D): #the Fourier coefficients C(D,r) of T(l)B_{12,3}^1 for any r in Z^3 and D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*(i**9)*CB1231(D*l**2/i**2)
    return sum

def CTB1232(l,D): #the Fourier coefficients C(D,0) of T(l)B_{12,3}^2 for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*CB1232(D*l**2/i**2)*i**9
    return sum

def T123(l): #the matrix of T(l) on F_1 and F_2
    G = matrix([[CB1231(-1),CB1232(-1)],[CB1231(-2),CB1232(-2)]])
    VB1231 = matrix([[CTB1231(l,-1)],[CTB1231(l,-2)]])
    GB1231 = G.inverse()*VB1231
    VB1232 = matrix([[CTB1232(l,-1)],[CTB1232(l,-2)]])
    GB1232 = G.inverse()*VB1232
    return matrix([[GB1231[0,0],GB1231[1,0]],[GB1232[0,0],GB1232[1,0]]])

def CPsi12(D): # Fourier coefficients C(D,0) of the first Hecke eigenform in J_{12,D_3}, computed by diagonalizing T(3)
    return 2/9*CB1231(D)+35/9*CB1232(D)

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

def CA12(D): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{12,D_3}, computed by diagonalizing T(3)
    return 7/9*CB1231(D)-35/9*CB1232(D)

def CTA12(l,D): #the Fourier coefficients C(D,0) of T(l)A_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*CA12(D*l**2/i**2)*i**9
    return sum

def lA12(l,D): #the eigenvalue of A_{12} for T(l)
    return CTA12(l,D)/CA12(D)

def Cphi12(D): # Fourier coefficients C(D,0) of the first Hecke eigenform in J_{12,D_5}, computed by diagonalizing T(3)
    return CB121(5,D)*13/63-CB122(5,D)*80/189

def CTphi12(l,D): #the Fourier coefficients C(D,0) of T(l)phi_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+Cphi12(D*l**2/i**2)*kronecker((-1)**((5-1)/2)*D*8,i)*i**(11-(5+1)/2)
    return sum

def lphi12(l,D): #the eigenvalue of phi_{12} for T(l)
    return CTphi12(l,D)/Cphi12(D)

def Ckap12(D): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{12,D_5}, computed by diagonalizing T(3)
    return CB121(5,D)*50/63+CB122(5,D)*80/189

def CTkap12(l,D): #the Fourier coefficients C(D,0) of T(l)kappa_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+Ckap12(D*l**2/i**2)*kronecker((-1)**((5-1)/2)*D*8,i)*i**(11-(5+1)/2)
    return sum

def lkap12(l,D): #the eigenvalue of kappa_{12} for T(l)
    return CTkap12(l,D)/Ckap12(D)

def CPhi12(D): # Fourier coefficients C(D,0) of the first Hecke eigenform in J_{12,D_7}, computed by diagonalizing T(3)
    return 3*CB121(7,D)-CB122(7,D)*432/5

def CTPhi12(l,D): #the Fourier coefficients C(D,0) of T(l)Phi_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+CPhi12(D*l**2/i**2)*kronecker((-1)**((7-1)/2)*D*8,i)*i**(11-(7+1)/2)
    return sum

def lPhi12(l,D): #the eigenvalue of Phi_{12} for T(l)
    return CTPhi12(l,D)/CPhi12(D)

def CK12(D): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{12,D_7}, computed by diagonalizing T(3)
    return -2*CB121(7,D)+CB122(7,D)*432/5

def CTK12(l,D): #the Fourier coefficients C(D,0) of T(l)K_{12} for any D coprime to l
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+CK12(D*l**2/i**2)*kronecker((-1)**((7-1)/2)*D*8,i)*i**(11-(7+1)/2)
    return sum

def lK12(l,D): #the eigenvalue of K_{12} for T(l)
    return CTK12(l,D)/CK12(D)

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
