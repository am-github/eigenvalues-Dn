d7=load("d7.sobj")
d5=load("d5.sobj")
d3=load("d3.sobj")
d1=load("d1.sobj")
d83=load("d83.sobj")
t=load("t.sobj")

def tau(n):
    return t[n]

def rep(n,x): #anna's code; #Fourier coefficients C(D,r) of E_{4,D_n} for every r in Z^n
    if n==1:
        return d7[x]
    if n==3:
        return d5[x]
    if n==5:
        return d3[x]
    if n==7:
        return d1[x]
    fnc = ''
    sm = ''
    bd = 'x'
    fnc += 's=0; x='+str(x)+'; n ='+str(n)+';'
    for idx in range(8-n):
        fnc += 'for(i%s = 0, sqrt('%idx + bd + '), '
        bd += ' - i%s^2'%idx
        sm += '+ i%s^2'%idx
    fnc += 'if('+sm[1:]+'== x,p=0;'
    for idx in range(8-n):
        fnc += 'if(i%s == 0, p++);'%idx
    fnc += 's += 2^(8-n-p);)'+')'*(8-n)
    gp(fnc)
    rk = gp('s')
    return ZZ(rk)

def rep1(n,D,r): #Fourier coefficients C(D,r) of E_{4,D_n} for every r in (1/2+Z)^n
    fnc = ''
    sm = ''
    bd = '-2*x-(8-n)/4'
    fnc += 's=0; x='+str(D)+'; n ='+str(n)+'; r='+str(r)+';'
    for idx in range(8-n):
        fnc += 'for(i%s = 0, sqrt('%idx + bd + '),'
        bd += ' - i%s^2 -i%s'%(idx,idx)
        fnc += 'if('+bd+'<0, break);'
        sm += '+ i%s^2 + i%s'%(idx,idx)
    fnc += 'if('+sm[1:]+'== -2*x-(8-n)/4,'
    fnc += 's += 2^(8-n);)'+')'*(8-n)
    gp(fnc)
    rk = gp('s')
    return ZZ(rk)

def repminodd(n,D,r): #Fourier coefficients C(D,r) of Jacobi forms of weight 4+n and index D_n for every r in (1/2+Z)^n
    if not (-D-(8-n)/8).is_integer():
        raise ValueError('Incorrect value of D!')  
    sm2 = ''
    bd = '-2*x-(8-n)/4'
    sm = ''
    pr = ''
    fnc = 's=0; x='+str(D)+'; n ='+str(n)+'; r='+str(r)+';'
    for idx in range(8-n):
        fnc += 'for(i%s = 0, sqrt('%idx + bd + '),'
        bd += ' - i%s^2 -i%s'%(idx,idx)
        fnc += 'if('+bd+'<0, break);'
        sm2 += '+ i%s^2 + i%s'%(idx,idx)
        sm += '+ i%s'%idx
        pr += '* (i%s+1/2)'%idx
    fnc += 'if('+sm2[1:]+'== -2*x-(8-n)/4,p=0;'
    fnc += 's += 2^(8-n)*(-1)^(r+(8-n)/2'+sm+')'+pr+';)'+')'*(8-n)
    gp(fnc)
    rk = gp('s')
    return ZZ(rk)

def smart_quotient(num,den):
    poly = den.minpoly()
    [p, s, _] = poly.coefficients()
    denconj = -s - den
    return (num*denconj/p).maxima_methods().rootscontract().simplify_full().expand()

def smarter_quotient(num,den):
    poly = den.minpoly()
    [p, s, _] = poly.coefficients()
    denconj = -s - den
    Q=(num*denconj/p).maxima_methods().rootscontract().simplify_full().maxima_methods().rootscontract().simplify_full()
    return smart_quotient(Q.numerator(),Q.denominator())

def scalar_product(l1,l2):
    return sum([l1[i]*l2[i] for i in range(len(l1))])

def CT(fnc,k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)fnc for any r in Z^n and D coprime to l for any fnc of even weight k and index D_n
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum += kronecker((-1)**((n-1)/2)*D*8,i)*fnc(n,D*l**2/i**2,r*l/i)*i**(k-1-(n+1)/2)
    return sum

def l(fnc,k,n,l,D,r): #the eigenvalue of fnc for T(l) computed using the pair (D,r) in the support of D_n for any fnc of even weight k and index D_n
    return CT(fnc,k,n,l,D,r)/fnc(n,D,r)

def CE4(n,D,r): #the Fourier coefficients C(D,r) of the Eisenstein series of weight 4 and index D_n
    if (8*D)%8 in [1,3,5,7]:
        return rep1(n,D,r)
    return rep(n,-2*D)

def CE6(n,D,r): #the Fourier coefficients C(D,r) of E_{6,D_n}
    sum = -(D+(8-n)/24)*CE4(n,D,r)
    for i in xrange(1,-D+1):
        sum += (8-n)*sigma(i,1)*CE4(n,D+i,r)
    return sum

def CEis8(n,D,r): #the Fourier coefficients C(D,r) of E_{8,D_n}
    sum =-(D+(12-n)/24)*CE6(n,D,r)
    for i in xrange(1,-D+1):
        sum=sum+(12-n)*sigma(i,1)*CE6(n,D+i,r)
    return sum

def CE(k,n,D,r): #the Fourier coefficients C(D,r) of E_{k,D_n}
    if k==8:
        return CEis8(n,D,r)
    if k==6:
        return CE6(n,D,r)
    if k==4:
        return CE4(n,D,r)

def CEE(t,k,n,D,r): #the Fourier coefficients C(D,r) of E_tE_{k,D_n}
    if k==8 and n==3:
        sum = C8(n,D,r)
        for i in xrange(1,-D+1):    
            sum += -2*t*sigma(i,t-1)*C8(n,D+i,r)/bernoulli(t)
        return sum
    sum = CE(k,n,D,r)
    for i in xrange(1,-D+1):    
        sum += -2*t*sigma(i,t-1)*CE(k,n,D+i,r)/bernoulli(t)
    return sum

def CDE(k,n,D,r): #the Fourier coefficients C(D,r) of DeltaE_{k,D_n} for any r in Z
    sum = 0
    for i in xrange(1,-D+1):
        sum=sum+tau(i)*CE(k,n,D+i,r)
    return sum

def C8(n,D,r):  #the Fourier coefficients C(D,r) of the cusp form of weight 8 and index D_n for any r in Z^n
    if n==3:
        if r==0 and (-D)<2402:
            return d83[-D]
        sum=0
        for i in xrange(0,sqrt(-2*D-1)+1):
            if -2*D-1-i**2-i>=0:
                for j in xrange(0,sqrt(-2*D-1-i**2-i)+1):
                    if -2*D-1-i**2-i-j**2-j>=0:
                        for k in xrange(0,sqrt(-2*D-1-i**2-i-j**2-j)+1):
                            if -2*D-1-i**2-i-j**2-j-k**2-k>=0:
                                for l in xrange(0,sqrt(-2*D-1-i**2-i-j**2-j-k**2-k)+1):
                                    if r!=0 and i**2+i+j**2+j+k**2+k+l**2+l == (-2*D-1):
                                        sum += 16*(-1)**(i+j+k+l)*(i+1/2)*(j+1/2)*(k+1/2)*(l+1/2)
                                    if -8*D-4-4*(i**2+i+j**2+j+k**2+k+l**2+l) > 0 and (-8*D-4-4*(i**2+i+j**2+j+k**2+k+l**2+l)).is_square():
                                        sum += 32*(-1)**(i+j+k+l+2*r)*(i+1/2)*(j+1/2)*(k+1/2)*(l+1/2)
        return sum
    return CEE(4,4,n,D,r)-CEis8(n,D,r)*24**2/((8-n)*(12-n))

def C10(n,D,r):  #the Fourier coefficients C(D,r) of the cusp form of weight 10 and index D_n for any r in Z^n
    return CEE(6,4,n,D,r)*(8-n)/24+CEE(4,6,n,D,r)

def CTk(fnc,k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)fnc for any r in Z^n and D coprime to l for any fnc of even weight k and index D_n
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum += kronecker((-1)**((n-1)/2)*D*8,i)*fnc(k,n,D*l**2/i**2,r*l/i)*i**(k-1-(n+1)/2)
    return sum

def lk(fnc,k,n,l,D,r): #the eigenvalue of fnc for T(l) computed using the pair (D,r) in the support of D_n for any fnc of even weight k and index D_n
    return CTk(fnc,k,n,l,D,r)/fnc(k,n,D,r)

def CB1(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form B_{k,n}^1 of even weight k>8 and index D_n for any r in Z^n (n!=3)
    return CEE(k-4,4,n,D,r)-CEE(k-8,8,n,D,r)*576/((12-n)*(8-n))

def CB2(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form B_{k,n}^2 of even weight k>8 and index D_n for any r in Z^n (n!=3)
    return CEE(k-6,6,n,D,r)+CEE(k-8,8,n,D,r)*24/(12-n)

def CB3(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form B_{k,n}^3 of even weight k>14 and index D_n for any r in Z^n (n!=3)
    if k in [16,18,20]:
        sum=0
        for i in xrange(1,-D+1):
            sum += CE(k-12,n,D+i,r)*tau(i)
        return sum
    if k == 22:
        F=CuspForms(1,18).basis()[0]
        sum = 0
        for i in xrange(1,-D+1):
            sum += CE(k-18,n,D+i,r)*F[i]
        return sum

def CB4(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form B_{k,n}^4 of even weight k>18 and index D_n for any r in Z^n (n!=3)
    F=CuspForms(1,16).basis()[0]
    sum=0
    for i in xrange(1,-D+1):
        sum += CE(k-16,n,D+i,r)*F[i]
    return sum

def TT(k,n,l): #the matrix of T(l) on B_{k,n}^1 and B_{k,n}^2 for even k<16
    G = matrix([[CB1(k,n,-1,0),CB2(k,n,-1,0)],[CB1(k,n,-2,0),CB2(k,n,-2,0)]])
    VB1 = matrix([[CTk(CB1,k,n,l,-1,0)],[CTk(CB1,k,n,l,-2,0)]])
    GB1 = G.inverse()*VB1
    VB2 = matrix([[CTk(CB2,k,n,l,-1,0)],[CTk(CB2,k,n,l,-2,0)]])
    GB2 = G.inverse()*VB2
    return matrix([[GB1[0,0],GB1[1,0]],[GB2[0,0],GB2[1,0]]])

def TTT(k,n,l): #the matrix of T(l) on B_{k,n}^1, B_{k,n}^2 and B_{k,n}^3 for even 20>k>14
    G = matrix([[CB1(k,n,-1,0),CB2(k,n,-1,0),CB3(k,n,-1,0)],[CB1(k,n,-2,0),CB2(k,n,-2,0),CB3(k,n,-2,0)],[CB1(k,n,-4,0),CB2(k,n,-4,0),CB3(k,n,-4,0)]])
    VB1 = matrix([[CTk(CB1,k,n,l,-1,0)],[CTk(CB1,k,n,l,-2,0)],[CTk(CB1,k,n,l,-4,0)]])
    GB1 = G.inverse()*VB1
    VB2 = matrix([[CTk(CB2,k,n,l,-1,0)],[CTk(CB2,k,n,l,-2,0)],[CTk(CB2,k,n,l,-4,0)]])
    GB2 = G.inverse()*VB2
    VB3 = matrix([[CTk(CB3,k,n,l,-1,0)],[CTk(CB3,k,n,l,-2,0)],[CTk(CB3,k,n,l,-4,0)]])
    GB3 = G.inverse()*VB3
    return matrix([[GB1[0,0],GB1[1,0],GB1[2,0]],[GB2[0,0],GB2[1,0],GB2[2,0]],[GB3[0,0],GB3[1,0],GB3[2,0]]])

def TTTT(k,n,l): #the matrix of T(l) on B_{k,n}^1 and B_{k,n}^2, B_{k,n}^3 and B_{k,n}^4 for even k>18
    G = matrix([[CB1(k,n,-1,0),CB2(k,n,-1,0),CB3(k,n,-1,0),CB4(k,n,-1,0)],[CB1(k,n,-2,0),CB2(k,n,-2,0),CB3(k,n,-2,0),CB4(k,n,-2,0)],[CB1(k,n,-4,0),CB2(k,n,-4,0),CB3(k,n,-4,0),CB4(k,n,-4,0)],[CB1(k,n,-5,0),CB2(k,n,-5,0),CB3(k,n,-5,0),CB4(k,n,-5,0)]])
    VB1 = matrix([[CTk(CB1,k,n,l,-1,0)],[CTk(CB1,k,n,l,-2,0)],[CTk(CB1,k,n,l,-4,0)],[CTk(CB1,k,n,l,-5,0)]])
    GB1 = G.inverse()*VB1
    VB2 = matrix([[CTk(CB2,k,n,l,-1,0)],[CTk(CB2,k,n,l,-2,0)],[CTk(CB2,k,n,l,-4,0)],[CTk(CB2,k,n,l,-5,0)]])
    GB2 = G.inverse()*VB2
    VB3 = matrix([[CTk(CB3,k,n,l,-1,0)],[CTk(CB3,k,n,l,-2,0)],[CTk(CB3,k,n,l,-4,0)],[CTk(CB3,k,n,l,-5,0)]])
    GB3 = G.inverse()*VB3
    VB4 = matrix([[CTk(CB4,k,n,l,-1,0)],[CTk(CB4,k,n,l,-2,0)],[CTk(CB4,k,n,l,-4,0)],[CTk(CB4,k,n,l,-5,0)]])
    GB4 = G.inverse()*VB4
    return matrix([[GB1[0,0],GB1[1,0],GB1[2,0],GB1[3,0]],[GB2[0,0],GB2[1,0],GB2[2,0],GB2[3,0]],[GB3[0,0],GB3[1,0],GB3[2,0],GB3[3,0]],[GB4[0,0],GB4[1,0],GB4[2,0],GB4[3,0]]])

def C1(k,n,D,r): # Fourier coefficients C(D,r) of the first Hecke eigenform in J_{k,D_n}, computed by diagonalizing T(3)
    if k in [12,14]:
        E,P=TT(k,n,3).eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[0,0]+CB2(k,n,D,r)*S[0,1]
    if k in [16,18]:
        M=matrix(SR,TTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[0,0]+CB2(k,n,D,r)*S[0,1]+CB3(k,n,D,r)*S[0,2]
    if k == 20:
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[0,0]+CB2(k,n,D,r)*S[0,1]+CB3(k,n,D,r)*S[0,2]+CB4(k,n,D,r)*S[0,3]
    if k == 22:
        if n in [1,7]:
            T = TTTT(22,n,3)
            R.<x> = QQ[]
            poly = R(T.characteristic_polynomial())
            K.<y> = poly.splitting_field()
            T = matrix(K,T)
            _, P = T.eigenmatrix_right()
            S=P.inverse()
            return CB1(k,n,D,r)*S[0,0]+CB2(k,n,D,r)*S[0,1]+CB3(k,n,D,r)*S[0,2]+CB4(k,n,D,r)*S[0,3]
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[0,0]+CB2(k,n,D,r)*S[0,1]+CB3(k,n,D,r)*S[0,2]+CB4(k,n,D,r)*S[0,3]

def C2(k,n,D,r): # Fourier coefficients C(D,0) of the second Hecke eigenform in J_{k,D_n}, computed by diagonalizing T(3)
    if k in [12,14]:
        E,P=TT(k,n,3).eigenmatrix_right()
        S=P.inverse()
        return S[1,0]*CB1(k,n,D,r)+S[1,1]*CB2(k,n,D,r)
    if k in [16,18]:
        M=matrix(SR,TTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[1,0]+CB2(k,n,D,r)*S[1,1]+CB3(k,n,D,r)*S[1,2]
    if k == 20:
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[1,0]+CB2(k,n,D,r)*S[1,1]+CB3(k,n,D,r)*S[1,2]+CB4(k,n,D,r)*S[1,3]
    if k == 22:
        if n in [1,7]:
            T = TTTT(22,n,3)
            R.<x> = QQ[]
            poly = R(T.characteristic_polynomial())
            K.<y> = poly.splitting_field()
            T = matrix(K,T)
            _, P = T.eigenmatrix_right()
            S=P.inverse()
            return CB1(k,n,D,r)*S[1,0]+CB2(k,n,D,r)*S[1,1]+CB3(k,n,D,r)*S[1,2]+CB4(k,n,D,r)*S[1,3]
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[1,0]+CB2(k,n,D,r)*S[1,1]+CB3(k,n,D,r)*S[1,2]+CB4(k,n,D,r)*S[1,3]

def C3(k,n,D,r): # Fourier coefficients C(D,0) of the third Hecke eigenform in J_{k,D_n}, computed by diagonalizing T(3)
    if k in [16,18]:
        T=matrix(SR,TTT(k,n,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[2,0]+CB2(k,n,D,r)*S[2,1]+CB3(k,n,D,r)*S[2,2]
    if k == 20:
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[2,0]+CB2(k,n,D,r)*S[2,1]+CB3(k,n,D,r)*S[2,2]+CB4(k,n,D,r)*S[2,3]
    if k == 22:
        if n in [1,7]:
            T = TTTT(22,n,3)
            R.<x> = QQ[]
            poly = R(T.characteristic_polynomial())
            K.<y> = poly.splitting_field()
            T = matrix(K,T)
            _, P = T.eigenmatrix_right()
            S=P.inverse()
            return CB1(k,n,D,r)*S[2,0]+CB2(k,n,D,r)*S[2,1]+CB3(k,n,D,r)*S[2,2]+CB4(k,n,D,r)*S[2,3]
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[2,0]+CB2(k,n,D,r)*S[2,1]+CB3(k,n,D,r)*S[2,2]+CB4(k,n,D,r)*S[2,3]

def C4(k,n,D,r): # Fourier coefficients C(D,0) of the fourth Hecke eigenform in J_{k,D_n}, computed by diagonalizing T(3)
    if k == 20:
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[3,0]+CB2(k,n,D,r)*S[3,1]+CB3(k,n,D,r)*S[3,2]+CB4(k,n,D,r)*S[3,3]
    if k == 22:
        if n in [1,7]:
            T = TTTT(22,n,3)
            R.<x> = QQ[]
            poly = R(T.characteristic_polynomial())
            K.<y> = poly.splitting_field()
            T = matrix(K,T)
            _, P = T.eigenmatrix_right()
            S=P.inverse()
            return CB1(k,n,D,r)*S[3,0]+CB2(k,n,D,r)*S[3,1]+CB3(k,n,D,r)*S[3,2]+CB4(k,n,D,r)*S[3,3]
        M=matrix(SR,TTTT(k,n,3))
        E,P=M.eigenmatrix_right()
        S=P.inverse()
        return CB1(k,n,D,r)*S[3,0]+CB2(k,n,D,r)*S[3,1]+CB3(k,n,D,r)*S[3,2]+CB4(k,n,D,r)*S[3,3]

def l1parallel(k,n,l,D,r): #the eigenvalue of psi_{k,n} for T(l)
    if k in [12,14]:
        return lk(C1,k,n,l,D,r)
    if k == 22 and n in [1,7]:
        return lk(C1,k,n,l,D,r)
    return smart_quotient(CTk(C1,k,n,l,D,r),C1(k,n,D,r))

def l2parallel(k,n,l,D,r): #the eigenvalue of alp_{k,n} for T(l)
    if k in [12,14]:
        return lk(C2,k,n,l,D,r)
    if k == 22 and n in [1,7]:
        T = TTTT(22,n,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        if n==7:
            poly = x^3 - 12422194*x - 2645665785
            v = R(poly).roots(K)[0][0]
            k, from_k = K.subfield(v)
            b1 = (24*v^2 + 44712*v - 198755104)/979
            b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        if n==1:
            poly = x^3-x^2-2784108376*x + 1945534874860
            v = R(poly).roots(K)[0][0]
            k, from_k = K.subfield(v)
            b1 = 48*v- 16
            b2 = 4*v^2 + 4212*v-7424290408
        V, from_V, to_V  = K.vector_space()
        R.<B1,B2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,beta_1,beta_2]])
        return scalar_product(W.coordinates(to_V(CTk(C2,22,n,l,D,r)/C2(22,n,D,r))),[1,B1,B2])
    return smart_quotient(CTk(C2,k,n,l,D,r),C2(k,n,D,r))

def l3parallel(k,n,l,D,r): #the eigenvalue of ffrak_{k,n} for T(l)
    if k in [16,18]:
        return lk(C3,k,n,l,D,r)
    if k == 22 and n in [1,7]:
        T = TTTT(22,n,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        if n==7:
            poly = x^3 - 12422194*x - 2645665785
            v = R(poly).roots(K)[2][0]
            k, from_k = K.subfield(v)
            b1 = (24*v^2 + 44712*v - 198755104)/979
            b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        if n==1:
            poly = x^3-x^2-2784108376*x + 1945534874860
            v = R(poly).roots(K)[2][0]
            k, from_k = K.subfield(v)
            b1 = 48*v- 16
            b2 = 4*v^2 + 4212*v-7424290408
        V, from_V, to_V  = K.vector_space()
        R.<B1,B2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,beta_1,beta_2]])
        return scalar_product(W.coordinates(to_V(CTk(C3,22,n,l,D,r)/C3(22,n,D,r))),[1,B1,B2])
    if k in [20,22] and l==1:
        return smart_quotient(CTk(C3,k,n,l,D,r),C3(k,n,D,r))
    return smarter_quotient(CTk(C3,k,n,l,D,r),C3(k,n,D,r))

def l4parallel(k,n,l,D,r): #the eigenvalue of hfrak_{k,n} for T(l)
    if k == 22 and n in [1,7]:
        T = TTTT(22,n,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        if n==7:
            poly = x^3 - 12422194*x - 2645665785
            v = R(poly).roots(K)[1][0]
            k, from_k = K.subfield(v)
            b1 = (24*v^2 + 44712*v - 198755104)/979
            b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        if n==1:
            poly = x^3-x^2-2784108376*x + 1945534874860
            v = R(poly).roots(K)[1][0]
            k, from_k = K.subfield(v)
            b1 = 48*v- 16
            b2 = 4*v^2 + 4212*v-7424290408
        V, from_V, to_V  = K.vector_space()
        R.<B1,B2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,beta_1,beta_2]])
        return scalar_product(W.coordinates(to_V(CTk(C4,22,n,l,D,r)/C4(22,n,D,r))),[1,B1,B2])
    if k in [20,22] and l==1:
        return smart_quotient(CTk(C4,k,n,l,D,r),C4(k,n,D,r))
    return smarter_quotient(CTk(C4,k,n,l,D,r),C4(k,n,D,r))

def CTn(fnc,k,l,D,r): #the Fourier coefficients C(D,r) of T(l)fnc for any D coprime to l and fnc in J_{k,D_3}
    if l!=1 and gcd(D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-D)*8,i)*fnc(k,D*l**2/i**2,r*l/i)*i**(k-3)
    return sum

def ln(fnc,k,l,D,r): #the eigenvalue of fnc for T(l)
    return CTn(fnc,k,l,D,r)/fnc(k,D,r)

def CB31(k,D,r): #the Fourier coefficients C(D,r) of the cusp form B_{k,3}^1 for any r in Z^3
    sum=CE4(3,D,r)*5/24+CE6(3,D,r)
    for i in xrange(1,-D+1):
        sum=sum-2*(k-4)*sigma(i,k-5)*CE4(3,D+i,r)*5/(24*bernoulli(k-4))-2*(k-6)*sigma(i,k-7)*CE6(3,D+i,r)/bernoulli(k-6)
    return sum

def CB32(k,D,r): #the Fourier coefficients C(D,0) of the cusp form B_{k,3}^2
    sum=C8(3,D,r)
    if r==0:
        for i in xrange(1,-D):
            sum=sum-2*(k-8)*sigma(i,k-9)*C8(3,D+i,r)/bernoulli(k-8)
    else:
        for i in xrange(1,-D+1):
            sum = sum-2*(k-8)*sigma(i,k-9)*C8(3,D+i,r)/bernoulli(k-8)
    return sum

def TT3(k,l): #the matrix of T(l) on B_{k,3}^1 and B_{k,3}^2 for even k<16
    G = matrix([[CB31(k,-1,0),CB32(k,-1,0)],[CB31(k,-2,0),CB32(k,-2,0)]])
    V1 = matrix([[CTn(CB31,k,l,-1,0)],[CTn(CB31,k,l,-2,0)]])
    G1 = G.inverse()*V1
    V2 = matrix([[CTn(CB32,k,l,-1,0)],[CTn(CB32,k,l,-2,0)]])
    G2 = G.inverse()*V2
    return matrix([[G1[0,0],G1[1,0]],[G2[0,0],G2[1,0]]])

def TTT3(k,l): #the matrix of T(l) on B_{k,3}^1, B_{k,3}^2 and B_{k,3}^3 for even k > 14
    G = matrix([[CB31(k,-1,0),CB32(k,-1,0),CB3(k,3,-1,0)],[CB31(k,-2,0),CB32(k,-2,0),CB3(k,3,-2,0)],[CB31(k,-4,0),CB32(k,-4,0),CB3(k,3,-4,0)]])
    V1 = matrix([[CTn(CB31,k,l,-1,0)],[CTn(CB31,k,l,-2,0)],[CTn(CB31,k,l,-4,0)]])
    G1 = G.inverse()*V1
    V2 = matrix([[CTn(CB32,k,l,-1,0)],[CTn(CB32,k,l,-2,0)],[CTn(CB32,k,l,-4,0)]])
    G2 = G.inverse()*V2
    V3 = matrix([[CTk(CB3,k,3,l,-1,0)],[CTk(CB3,k,3,l,-2,0)],[CTk(CB3,k,3,l,-4,0)]])
    G3 = G.inverse()*V3
    return matrix([[G1[0,0],G1[1,0],G1[2,0]],[G2[0,0],G2[1,0],G2[2,0]],[G3[0,0],G3[1,0],G3[2,0]]])

def TTTT3(k,l): #the matrix of T(l) on B_{k,3}^1, B_{k,3}^2, B_{k,3}^3 and B_{k,3}^4 for even k > 18
    G = matrix([[CB31(k,-1,0),CB32(k,-1,0),CB3(k,3,-1,0),CB4(k,3,-1,0)],[CB31(k,-2,0),CB32(k,-2,0),CB3(k,3,-2,0),CB4(k,3,-2,0)],[CB31(k,-4,0),CB32(k,-4,0),CB3(k,3,-4,0),CB4(k,3,-4,0)],[CB31(k,-5,0),CB32(k,-5,0),CB3(k,3,-5,0),CB4(k,3,-5,0)]])
    V1 = matrix([[CTn(CB31,k,l,-1,0)],[CTn(CB31,k,l,-2,0)],[CTn(CB31,k,l,-4,0)],[CTn(CB31,k,l,-5,0)]])
    G1 = G.inverse()*V1
    V2 = matrix([[CTn(CB32,k,l,-1,0)],[CTn(CB32,k,l,-2,0)],[CTn(CB32,k,l,-4,0)],[CTn(CB32,k,l,-5,0)]])
    G2 = G.inverse()*V2
    V3 = matrix([[CTk(CB3,k,3,l,-1,0)],[CTk(CB3,k,3,l,-2,0)],[CTk(CB3,k,3,l,-4,0)],[CTk(CB3,k,3,l,-5,0)]])
    G3 = G.inverse()*V3
    V4 = matrix([[CTk(CB4,k,3,l,-1,0)],[CTk(CB4,k,3,l,-2,0)],[CTk(CB4,k,3,l,-4,0)],[CTk(CB4,k,3,l,-5,0)]])
    G4 = G.inverse()*V4
    return matrix([[G1[0,0],G1[1,0],G1[2,0],G1[3,0]],[G2[0,0],G2[1,0],G2[2,0],G2[3,0]],[G3[0,0],G3[1,0],G3[2,0],G3[3,0]],[G4[0,0],G4[1,0],G4[2,0],G4[3,0]]])

def CpsiD3(k,D,r): # Fourier coefficients C(D,r) of the first Hecke eigenform in J_{k,D_3}
    if k in [12,14]:
        T = matrix(SR,TT3(k,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[0,0]+CB32(k,D,r)*S[0,1]
    if k in [16,18]:
        T = matrix(SR,TTT3(k,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[0,0]+CB32(k,D,r)*S[0,1]+CB3(k,3,D,r)*S[0,2]
    if k in [20,22]:
        T = TTTT3(k,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        T = matrix(K,T)
        _, P = T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[0,0]+CB32(k,D,r)*S[0,1]+CB3(k,3,D,r)*S[0,2]+CB4(k,3,D,r)*S[0,3]

def CphiD3(k,D,r): # Fourier coefficients C(D,r) of the second Hecke eigenform in J_{k,D_3}
    if k in [12,14]:
        T = matrix(SR,TT3(k,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[1,0]+CB32(k,D,r)*S[1,1]
    if k in [16,18]:
        T = matrix(SR,TTT3(k,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[1,0]+CB32(k,D,r)*S[1,1]+CB3(k,3,D,r)*S[1,2]
    if k in [20,22]:
        T = TTTT3(k,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        T = matrix(K,T)
        _, P = T.eigenmatrix_right()
        S=P.inverse()   
        return CB31(k,D,r)*S[1,0]+CB32(k,D,r)*S[1,1]+CB3(k,3,D,r)*S[1,2]+CB4(k,3,D,r)*S[1,3]

def CdeltaD3(k,D,r): # Fourier coefficients C(D,r) of the third Hecke eigenform in J_{k,D_3}
    if k in [16,18]:
        T = matrix(SR,TTT3(k,3))
        E,P=T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[2,0]+CB32(k,D,r)*S[2,1]+CB3(k,3,D,r)*S[2,2]
    if k in [20,22]:
        T = TTTT3(k,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        T = matrix(K,T)
        _, P = T.eigenmatrix_right()
        S=P.inverse()
        return CB31(k,D,r)*S[2,0]+CB32(k,D,r)*S[2,1]+CB3(k,3,D,r)*S[2,2]+CB4(k,3,D,r)*S[2,3]

def CkappaD3(k,D,r): # Fourier coefficients C(D,r) of the fourth Hecke eigenform in J_{k,D_3}
    R.<x> = QQ[]
    poly = R(T.characteristic_polynomial())
    K.<y> = poly.splitting_field()
    T = matrix(K,T)
    _, P = T.eigenmatrix_right()
    S=P.inverse()
    return CB31(k,D,r)*S[3,0]+CB32(k,D,r)*S[3,1]+CB3(k,3,D,r)*S[3,2]+CB4(k,3,D,r)*S[3,3]

def lpsiD3parallel(k,l,D,r): #the eigenvalue of psi_{k,D_3} for T(l)
    if k in [12,20,22]:
        return ln(CpsiD3,k,l,D,r)
    return smart_quotient(CTn(CpsiD3,k,l,D,r),CpsiD3(k,D,r))

def lphiD3parallel(k,l,D,r): #the eigenvalue of phi_{k,D_3} for T(l)
    if k == 12:
        return ln(CphiD3,k,l,D,r)
    if k == 20:
        T = TTTT3(20,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        poly = x^3 - 12422194*x - 2645665785
        v = R(poly).roots(K)[0][0]
        k, from_k = K.subfield(v)
        b1 = (24*v^2 + 44712*v - 198755104)/979
        b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CphiD3,20,l,D,r)/CphiD3(20,D,r))),[1,beta_1,beta_2])
    if k == 22:
        T = TTTT3(22,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()    
        poly = x^3-175630027*x-142249227846
        v = R(poly).roots(K)[2][0]
        k, from_k = K.subfield(v)
        b1 = 72*v
        b2 = (216*v^2-262224*v-25290723888)/7
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CphiD3,22,l,D,r)/CphiD3(22,D,r))),[1,beta_1,beta_2])
    return smart_quotient(CTn(CphiD3,k,l,D,r),CphiD3(k,D,r))

def ldeltaD3parallel(k,l,D,r): #the eigenvalue of delta_{k,D_3} for T(l) and k>14
    if k == 20:
        T = TTTT3(20,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        poly = x^3 - 12422194*x - 2645665785
        v = R(poly).roots(K)[2][0]
        k, from_k = K.subfield(v)
        b1 = (24*v^2 + 44712*v - 198755104)/979
        b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CdeltaD3,20,l,D,r)/CdeltaD3(20,D,r))),[1,beta_1,beta_2])
    if k == 22:
        T = TTTT3(22,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()    
        poly = x^3-175630027*x-142249227846
        v = R(poly).roots(K)[0][0]
        k, from_k = K.subfield(v)
        b1 = 72*v
        b2 = (216*v^2-262224*v-25290723888)/7
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CdeltaD3,22,l,D,r)/CdeltaD3(22,D,r))),[1,beta_1,beta_2])
    return CTn(CdeltaD3,k,l,D,r)/CdeltaD3(k,D,r)

def lkappaD3parallel(k,l,D,r): #the eigenvalue of kappa_{k,D_3} for T(l) and k>18
    if k == 20:
        T = TTTT3(20,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()
        poly = x^3 - 12422194*x - 2645665785
        v = R(poly).roots(K)[1][0]
        k, from_k = K.subfield(v)
        b1 = (24*v^2 + 44712*v - 198755104)/979
        b2 = (-39144*v^2 + 84967848*v + 324169574624)/979
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CkappaD3,20,l,D,r)/CkappaD3(20,D,r))),[1,beta_1,beta_2])
    if k == 22:
        T = TTTT3(22,3)
        R.<x> = QQ[]
        poly = R(T.characteristic_polynomial())
        K.<y> = poly.splitting_field()    
        poly = x^3-175630027*x-142249227846
        v = R(poly).roots(K)[1][0]
        k, from_k = K.subfield(v)
        b1 = 72*v
        b2 = (216*v^2-262224*v-25290723888)/7
        V, from_V, to_V  = K.vector_space()
        R.<beta_1,beta_2> = PolynomialRing(QQ,2)
        R._latex_names = ['\\beta_1', '\\beta_2'] 
        W = V.subspace_with_basis([to_V(el) for el in [1,b1,b2]])
        return scalar_product(W.coordinates(to_V(CTn(CkappaD3,22,l,D,r)/CkappaD3(22,D,r))),[1,beta_1,beta_2])

def Cminodd(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form E_kpsi_{12-n,D_n} for the lattice D_n for any r in (1/2+Z)^n
    sum = repminodd(n,D,r)
    if k==0:
        return sum
    for i in xrange(1,-D+1):
        sum = sum-2*k*repminodd(n,D+i,r)*sigma(i,k-1)/bernoulli(k)
    return sum 

def CTminodd(k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)E_kpsi_{12-n,D_n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and (8*D)%l==0:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd(k,n,D*l**2/i**2,r*l/i)*i**(12-n+k-1-(n+1)/2)
    return sum

def lminoddparallel(k,n,l,D,r): #the eigenvalue of E_kpsi_{12-n,D_n} for T(l) for k<12
    return CTminodd(k,n,l,D,r)/Cminodd(k,n,D,r)

def Cminodd2(k,n,D,r): #the Fourier coefficients C(D,r) of the cusp form psi_{12-n,D_n}F_k of weight 12+k-n for the lattice D_n for any r in (1/2+Z)^n
    F=CuspForms(1,k).basis()[0]
    sum=0
    for i in xrange(1,-D+1):
        sum = sum+Cminodd(0,n,D+i,r)*F[i]
    return sum

def CTminodd2(k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)psi_{12-n,D_n}F_k for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminodd2(k,n,D*l**2/i**2,r*l/i)*i**(12+k-n-1-(n+1)/2)
    return sum

def Tminodd(k,n,p): #the matrix of T(p) on psi_{12-n}E_k and psi_{12-n}F_k
    if n==1:
        G = matrix([[Cminodd(k,n,-7/8,n/2),Cminodd2(k,n,-7/8,n/2)],[Cminodd(k,n,-23/8,n/2),Cminodd2(k,n,-23/8,n/2)]])
        VminoddE = matrix([[CTminodd(k,n,p,-(8-n)/8,n/2)],[CTminodd(k,n,p,-23/8,n/2)]])
        GminoddE = G.inverse()*VminoddE
        Vminodd = matrix([[CTminodd2(k,n,p,-(8-n)/8,n/2)],[CTminodd2(k,n,p,-23/8,n/2)]])
        Gminodd = G.inverse()*Vminodd
    else:
        G = matrix([[Cminodd(k,n,-(8-n)/8,n/2),Cminodd2(k,n,-(8-n)/8,n/2)],[Cminodd(k,n,-(16-n)/8,n/2),Cminodd2(k,n,-(16-n)/8,n/2)]])
        VminoddE = matrix([[CTminodd(k,n,p,-(8-n)/8,n/2)],[CTminodd(k,n,p,-(16-n)/8,n/2)]])
        GminoddE = G.inverse()*VminoddE
        Vminodd = matrix([[CTminodd2(k,n,p,-(8-n)/8,n/2)],[CTminodd2(k,n,p,-(16-n)/8,n/2)]])
        Gminodd = G.inverse()*Vminodd
    return matrix([[GminoddE[0,0],GminoddE[1,0]],[Gminodd[0,0],Gminodd[1,0]]])

def Cminoddk1(k,n,D,r): #the Fourier coefficients C(D,r) of the first Hecke eigenform in J_{12+k-n,D_n}, computed by diagonalizing Tminodd(k,n,p)
    if n in [1,3]:
        T = matrix(SR,Tminodd(k,n,3))
    else: 
        T = matrix(SR,Tminodd(k,n,5))
    E,P = T.eigenmatrix_right()
    S=P.inverse()
    return Cminodd(k,n,D,r)*S[0,0]+Cminodd2(k,n,D,r)*S[0,1]

def CTminoddk1(k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)psi_{12+k-n,D_n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminoddk1(k,n,D*l**2/i**2,r*l/i)*i**(12+k-n-1-(n+1)/2)
    return sum

def lminoddk1parallel(k,n,l,D,r): #the Hecke eigenvalues of psi_{12+k-n,D_n} for any r in (Z+1/2)^n and 8D coprime to l
    return smart_quotient(CTminoddk1(k,n,l,D,r),Cminoddk1(k,n,D,r))

def Cminoddk2(k,n,D,r): #the Fourier coefficients C(D,r) of the second Hecke eigenform in J_{12+k-n,D_n}, computed by diagonalizing Tminodd(k,n,p)
    if n in [1,3]:
        T = matrix(SR,Tminodd(k,n,3))
    else: 
        T = matrix(SR,Tminodd(k,n,5))
    E,P = T.eigenmatrix_right()
    S=P.inverse()
    return Cminodd(k,n,D,r)*S[1,0]+Cminodd2(k,n,D,r)*S[1,1]

def CTminoddk2(k,n,l,D,r): #the Fourier coefficients C(D,r) of T(l)phi_{12+k-n,D_n} for any r in (Z+1/2)^n and 8D coprime to l
    if l!=1 and gcd(8*D,l)>1:
        raise ValueError('the variable is not coprime to l!')
    if l%2==0:
        raise ValueError('l is not coprime to the level!')
    sum=0
    for i in l.divisors():
        sum=sum+kronecker((-1)**((n-1)/2)*D*8,i)*Cminoddk2(k,n,D*l**2/i**2,r*l/i)*i**(12+k-n-1-(n+1)/2)
    return sum

def lminoddk2parallel(k,n,l,D,r): #the Hecke eigenvalues of phi_{12+k-n,D_n} for any r in (Z+1/2)^n and 8D coprime to l
    return smart_quotient(CTminoddk2(k,n,l,D,r),Cminoddk2(k,n,D,r))
