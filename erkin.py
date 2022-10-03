from matplotlib.pylab import *
from scipy.optimize import least_squares
from IPython.display import display, Math


def print_sci(x, e=0, n=0, stri=''):
    pot = int(log10(x))
    mant = x/10**pot
    if e !=0 or isinf(e):
        pote = int(log10(e))
        err = e/10**pot
        if pot == 0:
            return display(Math((f'{stri} = %.{pot-pote+1}f\ (\pm %.{pot-pote+2}f)'%(mant,err)).replace('.',',')))
        else:
            return display(Math((f'{stri} = %.{pot-pote+1}f\ (\pm %.{pot-pote+2}f)\\times 10^{{%i}}'%(mant,err,pot)).replace('.',',')))
    else:
        numero=f'{stri} = %.{n}f'
        if pot !=0 :
            numero += '\\times 10^{{%i}}'
            return display(Math(f'{numero}'%(mant,pot)))
        else:
            return display(Math(numero%mant))

def modelo2(par, X, Y):
    k, n, Co = par
    ts = 1/k/(n-1)*(Y**(1-n)-Co**(1-n))
    return X-ts

def ERArrhenius(par, X, Y):
    R = 8.31446261815324
    k, Ea = par
    T = mean(X)
    y = k/exp(-Ea/R/T)*exp(-Ea/R/X)
    return Y-y

def errFit(hess_inv, resVariance):
    return np.sqrt( np.diag( hess_inv * resVariance))


def cinetica(X,Y, conversao = False, grafico=False):
    '''
    X - lista de tempos
    Y - lista de concentrações
    
    se conversão = True, Y é uma lista de conversões
    
    Devolve k, n, Co e os erros associados
    
    Se conversao = True, devolve n e kCAo^{n-1}
    
    V0.9©FGF2022
    '''
    if len(X)<=2:
        print("\x1b[31m Número de dados insuficientes \x1b[0m")
    else:
        if conversao:
            if X[0]!=0:
                X.insert(0,0)
                Y=[1-i for i in Y]
                Y.insert(0,1)
            else:
                Y=[1-i for i in Y]
        par0 = [0.1,2.1, 1]
        result = least_squares(modelo2, par0 ,bounds=[(0,-3,0),(1,3,inf)], args=(X,Y))
        if len(X)>3:
            dFit = errFit(linalg.inv(dot(result.jac.T, result.jac)),
                  (modelo2(result.x, X, Y)**2).sum()/(len(X)-len(par0) ) )
        else:
            print("\x1b[31m Número de dados insuficientes para calcular erros. \x1b[0m")
            dFit=[0, 0, 0]

        Cs = linspace(Y[-1],Y[0])
        tmod = 1/result.x[0]/(result.x[1]-1)*(Cs**(1-result.x[1])-Y[0]**(1-result.x[1]))
        if grafico:
            plot(X,Y,'X')
            plot(tmod,Cs)
            xlabel('t')
            ylabel('C')
            if conversao:
                ylabel('X')
        if conversao:
            stri=['kCo^{(n-1)}','n']
            for i in range(2):
                print_sci(result.x[i], dFit[i], stri=stri[i])
        else:
            stri=['k','n','Co']
            for i in range(3):
                print_sci(result.x[i], dFit[i], stri=stri[i])
        return(result, dFit)

def ERA(X,Y, Kelvin = False, grafico = False):
    '''
    X - lista de temperaturas (ºC)
    Y - lista de constantes de velocidade

    Se Kelvin = True X está em K.
 
    Devolve k e Ea e os erros associados
    
    Se Kelvin != True, em primeiro lugar as temperaturas em ºC são transformadas em K
    
    V0.9©FGF2022
    '''
    R = 8.31446261815324
    
    if not Kelvin:
        X=[i+273.15 for i in X]
    
    if len(X)<=2:
        print("\x1b[31m Número de dados insuficientes para calcular erros. \x1b[0m")
        
    par0 = [0.1,1e3]
    result = least_squares(ERArrhenius, par0 ,bounds=[(0,1000),(inf,inf)], args=(X,Y))
    
    if len(X)>2:
        dFit = errFit(linalg.inv(dot(result.jac.T, result.jac)),
                  (ERArrhenius(result.x, X, Y)**2).sum()/(len(X)-len(par0) ) )
    else:
        dFit=[0]*len(par0)

    Ts = linspace(min(X),max(X))
    
    ks = result.x[0]/exp(-result.x[1]/R/mean(X))*exp(-result.x[1]/R/Ts)
    if grafico:
        plot(X,Y,'X')
        plot(Ts,ks)
        xlabel('T (K)')
        ylabel('k')

    stri=['ko','Ea\ \ (J/mol)']
    result.x[0]=result.x[0]/exp(-result.x[1]/R/mean(X))
    dFit[0]=dFit[0]/exp(-result.x[1]/R/mean(X))
    for i in range(2):
        print_sci(result.x[i], dFit[i], stri=stri[i])
    return(result, dFit)


if __name__ == '__main__':
    X=[0, 10, 20, 40]
    Y=[1000, 500, 333, 200]
    fit, erros = cinetica(X,Y, grafico=False)
    
    X=[20, 70, 50]
    Y=[1, 32, 8]
    fit, erros = ERA(X,Y, grafico=True)
    show()