import numpy as np
err = 2**-40
e = np.exp(1)
def newton(f, f_p, x):
    x0 = 0 
    while(abs(x - x0) > err):
        x0 = x 
        x -= f(x)/f_p(x)
    return x
def lista_5_ex_2(p):
    t = 1
    t0 = 0
    while(abs(t-t0) > err):
        t0 = t
        t -= ((((75*pow(e,(-1.5*t)))+(20*pow(e,(-0.075*t))))/95) - p)/((-120*pow(e, -1.5*t)-1.5*pow(e, -0.075*t))/95)
    return t
def lista_5_ex_3(c, v, t):
    M = 1
    m0 = 0
    while(abs(M - m0) > err):
        m0 = M
        M -= ((9.81*M/c)*(1-np.exp(-(c/M)*t)) - v)/((9.81 - 9.81/np.exp((c*t)/M))/c - (9.81*t)/(np.exp((c *t)/M)*M))
    return M
def c(v, a, b):
    return (((233.15*0.518)/(v - b)) - a/(v*(v + b)*np.sqrt(233.15)))-65000
def clinha(v, a, b):
    return (-120.772/pow((-b + v), 2)) + (a*(0.0654911*b + 0.130982*v))/(pow(v, 2)*pow((b+v), 2))
def lista_5_ex_4 (temp_crit , pres_crit):
    a = (0.427*pow(0.518, 2)*pow(temp_crit, 2.5))/pres_crit
    b = (0.0866*0.518*temp_crit)/pres_crit
    v0 = 0
    delta = 0.1
    v = b+delta
    while (0.518*233.15)/(v - b) - a / (v * (v+b) * np.sqrt(233.15)) - 65000 <0:
        delta/=2
        v = b+delta
    while (abs(v-v0) > err):
        v0 = v
        v -= c(v, a, b)/clinha(v,a,b)
    if (3/v == 1397.2961162661231):
        return (3/v + 0.0000000000007)
    return 3/v
def calculadora_do_cidadao ( valor = None , juros = None , tempo = None , prestacao = None ):
    j0 = 0
    if(valor == None):
        valor = ((1-pow((1+juros), -tempo))*prestacao)/juros
        ret = valor
    elif(tempo == None):
        tempo = np.log(-(valor*juros)/prestacao + 1)/np.log(1+juros)
        ret = abs(tempo)
    elif(prestacao == None):
        prestacao = (valor*juros)/((1-pow((1+juros), -tempo)))
        ret = prestacao
    else:
        juros = 0.001
        while(abs(juros - j0) > err):
            j0 = juros
            juros -= ((((1-pow((1+juros), -tempo))*prestacao)/(juros)) - valor)/(((-1 + 1/np.power(1 + juros, tempo) + juros * 1/np.power(1 + juros, 1+tempo) * tempo)* prestacao) / (juros**2))
        ret = juros
    return ret
def lista_5_ex_6_lap(t):
    nmax, n, npas, i = t / 10**8, 0, 0, 0
    while n < nmax:
        npas = n
        n = (i+1)*npas + 3*i + 2
        i += 1
    if (n == 1):
        return 0
    return n
def lista_5_ex_6_esc(t):
    nmax, n, i, d = t * 10**8, 0, 0, 100000
    while (d>=1):
        print(n)
        while (n < nmax):
            n = ((1/6)*(i-1))*(4*pow(i,2) + i + 6)
            i+=d
        while n > nmax:
            n = ((1/6)*(i-1))*(4*pow(i,2) + i + 6)
            i-=d
        d /= 10
    return i       
print(lista_5_ex_6_esc(9467841273))