import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from IPython.display import display, Math


def _format_sci(number, e=0, n=0, label=''):
    """Return a LaTeX-formatted string for a number with optional error.

    Designed for display(Math(...)). Works with scalar floats.
    """
    try:
        x = float(number)
    except Exception:
        x = number

    if e is None:
        e = 0

    # handle zero or non-finite
    if not np.isfinite(x):
        return f"{label} = {x}"

    if e and np.isfinite(e) and e != 0:
        # scientific notation
        pot = int(np.floor(np.log10(abs(x)))) if x != 0 else 0
        mant = x / 10 ** pot
        pote = int(np.floor(np.log10(abs(e)))) if e != 0 else pot
        errmant = e / 10 ** pot
        if pot == 0:
            return f"{label} = {mant:.{pot-pote+1}f} \ (\pm {errmant:.{pot-pote+2}f})".replace('.', ',')
        else:
            return f"{label} = {mant:.{pot-pote+1}f} \ (\pm {errmant:.{pot-pote+2}f})\\times 10^{{{pot}}}".replace('.', ',')
    else:
        # fixed number of decimals
        if n <= 0:
            return f"{label} = {x}"
        pot = int(np.floor(np.log10(abs(x)))) if x != 0 else 0
        if pot != 0:
            mant = x / 10 ** pot
            return f"{label} = {mant:.{n}f}\\times 10^{{{pot}}}".replace('.', ',')
        else:
            return f"{label} = {x:.{n}f}".replace('.', ',')


def print_sci(x, e=0, n=0, stri=''):
    """Display a number (with optional error) nicely in Jupyter using Math().

    Behaviour tries to preserve the original output style in the repo but uses
    numpy for numerical ops and avoids mutating inputs.
    """
    s = _format_sci(x, e=e, n=n, label=stri)
    try:
        display(Math(s))
    except Exception:
        print(s)


def modelo2(par, X, Y):
    """Model function for kinetics fitting.

    par = [k, n, Co]
    Returns residual vector: model_times - observed_times
    """
    k, n, Co = par
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    # avoid division by zero or invalid n
    if k == 0 or (n - 1) == 0:
        # produce large residuals to steer optimizer away
        return np.ones_like(X) * 1e6

    ts = 1.0 / k / (n - 1.0) * (Y ** (1.0 - n) - Co ** (1.0 - n))
    return X - ts


def ERArrhenius(par, X, Y):
    """Residual for Arrhenius fit.

    par = [k0, Ea]
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    R = 8.31446261815324

    k0, Ea = par
    # model: k(T) = k0 * exp(-Ea/(R*T))
    y = k0 * np.exp(-Ea / (R * X))
    return Y - y


def errFit(jacobian, resVariance):
    """Estimate parameter standard errors from jacobian and residual variance.

    jacobian: result.jac (m x p)
    resVariance: scalar estimate of residual variance (RSS/(m-p))
    Returns vector of 1-sigma uncertainties for parameters.
    """
    # Use pseudo-inverse for stability
    try:
        cov = np.linalg.pinv(jacobian.T.dot(jacobian)) * resVariance
        return np.sqrt(np.abs(np.diag(cov)))
    except Exception:
        return np.zeros(jacobian.shape[1])


def cinetica(X, Y, conversao=False, grafico=False):
    """Fit kinetic parameters k, n, Co from batch data.

    X - list/array of times
    Y - list/array of concentrations (or conversion if conversao=True)

    Returns: (result, dFit) where result is the least_squares OptimizeResult and
    dFit is an array of parameter standard errors.
    """
    # Work on copies and convert to numpy arrays
    X = list(X) if isinstance(X, (list, tuple)) else np.array(X, dtype=float)
    Y = list(Y) if isinstance(Y, (list, tuple)) else np.array(Y, dtype=float)

    if conversao:
        # If provided Y as conversions, transform to concentrations in [0, Co]
        # Assume first concentration Co corresponds to conversion 0
        # Do not mutate original inputs
        Y_arr = np.asarray(Y, dtype=float)
        # If conversions are given as fractions of initial conc.
        # We need an estimate for Co; use 1.0 as default scale if not enough info
        Co_guess = 1.0
        # Convert conversions X_conv to concentrations C = Co * (1 - X_conv)
        Yc = Co_guess * (1.0 - Y_arr)
        # if X was given without 0 at start, insert t=0 and conc=Co
        X_arr = np.asarray(X, dtype=float)
        if X_arr[0] != 0:
            X_arr = np.concatenate(([0.0], X_arr))
            Yc = np.concatenate(([Co_guess], Yc))
    else:
        X_arr = np.asarray(X, dtype=float)
        Yc = np.asarray(Y, dtype=float)

    if X_arr.size < 3 or Yc.size < 3:
        raise ValueError("Número de dados insuficientes (precisa de >= 3 pontos)")

    # initial guesses
    par0 = np.array([0.1, 1.5, float(Yc[0]) if Yc[0] != 0 else 1.0])
    lb = np.array([0.0, -5.0, 0.0])
    ub = np.array([np.inf, 5.0, np.inf])

    result = least_squares(modelo2, par0, bounds=(lb, ub), args=(X_arr, Yc))

    # residual variance estimate
    m = X_arr.size
    p = par0.size
    rss = np.sum(modelo2(result.x, X_arr, Yc) ** 2)
    resVar = rss / max(1, (m - p))

    # parameter uncertainties
    dFit = errFit(result.jac, resVar) if m > p else np.zeros_like(result.x)

    # produce model curve if requested
    if grafico:
        Cs = np.linspace(Yc[-1], Yc[0], 100)
        try:
            tmod = 1.0 / result.x[0] / (result.x[1] - 1.0) * (Cs ** (1.0 - result.x[1]) - Yc[0] ** (1.0 - result.x[1]))
            plt.plot(X_arr, Yc, 'X', label='data')
            plt.plot(tmod, Cs, label='fit')
            plt.xlabel('t')
            plt.ylabel('C' if not conversao else 'X')
            plt.legend()
            plt.show()
        except Exception:
            pass

    # print results (maintain similar string labels)
    if conversao:
        labels = ['kCo^{(n-1)}', 'n']
        for i in range(2):
            print_sci(result.x[i], dFit[i] if i < len(dFit) else 0, stri=labels[i])
    else:
        labels = ['k', 'n', 'Co']
        for i in range(3):
            print_sci(result.x[i], dFit[i] if i < len(dFit) else 0, stri=labels[i])

    return result, dFit


def ERA(X, Y, Kelvin=False, grafico=False):
    """Fit Arrhenius parameters k0 and Ea from k vs T data.

    X - temperatures (ºC or K if Kelvin=True)
    Y - rate constants
    """
    X_arr = np.asarray(X, dtype=float)
    Y_arr = np.asarray(Y, dtype=float)
    if not Kelvin:
        X_arr = X_arr + 273.15

    if X_arr.size < 2:
        raise ValueError("Precisa de pelo menos 2 pontos para ajuste de Arrhenius")

    par0 = np.array([1.0, 5e4])  # k0, Ea (J/mol)
    lb = np.array([0.0, 0.0])
    ub = np.array([np.inf, np.inf])

    result = least_squares(ERArrhenius, par0, bounds=(lb, ub), args=(X_arr, Y_arr))

    m = X_arr.size
    p = par0.size
    rss = np.sum(ERArrhenius(result.x, X_arr, Y_arr) ** 2)
    resVar = rss / max(1, (m - p))
    dFit = errFit(result.jac, resVar) if m > p else np.zeros_like(result.x)

    if grafico:
        Ts = np.linspace(np.min(X_arr), np.max(X_arr), 100)
        ks = result.x[0] * np.exp(-result.x[1] / (8.31446261815324 * Ts))
        plt.plot(X_arr, Y_arr, 'X', label='data')
        plt.plot(Ts, ks, label='fit')
        plt.xlabel('T (K)')
        plt.ylabel('k')
        plt.legend()
        plt.show()

    labels = ['ko', 'Ea\\ \\ (J/mol)']
    for i in range(2):
        print_sci(result.x[i], dFit[i] if i < len(dFit) else 0, stri=labels[i])

    return result, dFit


if __name__ == '__main__':
    # small smoke tests when run as a script
    X = [0, 10, 20, 40]
    Y = [1000, 500, 333, 200]
    fit, erros = cinetica(X, Y, grafico=False)

    X = [20, 70, 50]
    Y = [1, 32, 8]
    fit, erros = ERA(X, Y, grafico=True)
