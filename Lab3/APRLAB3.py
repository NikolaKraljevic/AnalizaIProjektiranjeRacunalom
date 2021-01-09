import random
import numpy
import numpy as np
from math import sqrt, inf
from scipy import optimize


def f10(x):
    return pow(x - 3, 2)


def f1(x):
    return 100 * pow((x[1] - pow(x[0], 2)), 2) + pow(1 - x[0], 2)


def df1x0(x):
    return 400 * pow(x[0], 3) + (2 - 400 * x[1]) * x[0] - 2


def df1x1(x):
    return 200 * (x[1] - pow(x[0], 2))


def f2(x):
    return pow(x[0] - 4, 2) + 4 * pow(x[1] - 2, 2)


def df2x0(x):
    return 2 * (x[0] - 4)


def df2x1(x):
    return 2 * (x[1] - 2)


def f3(x):
    return pow(x[0] - 2, 2) + pow(x[1] + 3, 2)


def df3x0(x):
    return 2 * (x[0] - 2)


def df3x1(x):
    return 2 * (x[1] + 3)


def f4(x):
    return pow(x[0] - 3, 2) + pow(x[1], 2)


def df4(x):
    return 2 * (x[0] - 3) + 2 * x[1]


def impl1(x):
    return (x[1] - x[0] >= 0)


def impl2(x):
    return ((2 - x[0]) >= 0)


def IzracunajCent(X, h=None):
    xc = [0] * len(X[0])
    for idx, x in enumerate(X):
        if idx == h:
            continue
        for i in range(len(x)):
            xc[i] += x[i]
    if h == None:
        n = len(X)
    else:

        n = len(X) - 1
    for idx, x in enumerate(xc):
        xc[idx] /= float(n)
    return xc


def reflection(Xc, Xh, alpha):
    Xr = [0] * len(Xc)
    for i in range(len(Xc)):
        Xr[i] = (1 + alpha) * Xc[i] - alpha * Xh[i]
    return Xr


def get_max_index(f, xs, h, e=1e-6):
    indexOfMax = 0
    currentMax = -inf
    for i in range(len(xs)):
        x = xs[i]
        funcValueOfCurrent = f(x)
        if (abs(funcValueOfCurrent - currentMax) > e) and (funcValueOfCurrent > currentMax) and (i != h):
            currentMax = funcValueOfCurrent
            indexOfMax = i
    return indexOfMax


def RMSE(function, simplexPoints, xc):
    sum = 0
    for simplexPoint in (simplexPoints):
        sum += (function(simplexPoint) - function(xc)) ** 2
    return sqrt(sum / float(len(simplexPoints)))


def checkIneq(x0, ineq):
    for res in ineq:
        if res(x0) == False:
            return False
    return True


def F(f, X, vs, lambd):
    xs = X.copy()
    for idx, x in enumerate(xs):
        xs[idx] = np.add(xs[idx], lambd * vs[idx])
    return f(xs)


def box(f, x0, ineq, explicitRest=[-10, 10], e=0.00001, alpha=1.3):
    tocka = x0
    xc = tocka
    n = len(tocka)
    xs = []
    for t in range(2 * n):
        x = [0] * n
        for i in range(n):
            r = random.random()
            x[i] = explicitRest[0] + r * (explicitRest[1] - explicitRest[0])
        xs.append(x)
        while not checkIneq(xs[t], ineq):
            for idx, x in enumerate(xs[t]):
                xs[t][idx] = 0.5 * (xs[t][idx] + xc[idx])
        xc = IzracunajCent(xs)
    brojac = 0
    proslaVrijednost = f(xc)
    trenutno = proslaVrijednost
    while True:
        if brojac > 100:
            print('Ne konvergira')
            return xc
        if trenutno == proslaVrijednost:
            brojac += 1
        else:
            brojac = 0
        print(xs)
        h = get_max_index(f, xs, -1)
        h2 = get_max_index(f, xs, h)
        xc = IzracunajCent(xs, h)
        xr = reflection(xs[h], xc, alpha)
        print('ovo je centr')
        print(xc)
        print('refleksija')
        print(xr)
        for i in range(n):
            if xr[i] < explicitRest[0]:
                xr[i] = explicitRest[0]
            elif xr[i] > explicitRest[1]:
                xr[i] = explicitRest[1]
        while not checkIneq(xr, ineq):
            for idx, x in enumerate(xr):
                xr[idx] = 0.5 * (xr[idx] + xc[idx])
        if f(xr) > f(xs[h2]):
            for idx, x in enumerate(xr):
                xr[idx] = 0.5 * (xr[idx] + xc[idx])
        xs[h] = xr
        proslaVrijednost = trenutno

        trenutno = f(xc)
        if RMSE(f, xs, xc) <= e:
            return xc


def gradient_f(x, f):
    assert (x.shape[0] >= x.shape[1]), "the vector should be a column vector"
    x = x.astype(float)
    N = x.shape[0]
    gradient = []
    for i in range(N):
        eps = abs(x[i]) * np.finfo(np.float32).eps
        xx0 = 1. * x[i]
        f0 = f(x)
        x[i] = x[i] + eps
        f1 = f(x)
        gradient.append(np.asscalar(np.array([f1 - f0])) / eps)
        x[i] = xx0
    return np.array(gradient).reshape(x.shape)


def hessian(x, the_func):
    N = x.shape[0]
    hessian = np.zeros((N, N))
    gd_0 = gradient_f(x, the_func)
    eps = np.linalg.norm(gd_0) * np.finfo(np.float32).eps
    for i in range(N):
        xx0 = 1. * x[i]
        x[i] = xx0 + eps
        gd_1 = gradient_f(x, the_func)
        hessian[:, i] = ((gd_1 - gd_0) / eps).reshape(x.shape[0])
        x[i] = xx0
    return hessian


def NewtonRhapson(f, x0, e, zlato=False):
    tocka = x0
    trenutno = f(tocka)
    prijasnje = f(tocka)
    delta = [1] * len(tocka)
    grad = []
    grad.append(0)
    grad.append(0)
    brojac = 0
    brojEvaluacija = 0
    while numpy.linalg.norm(delta) > e:
        brojEvaluacija += 1
        if brojac > 100:
            print("Ne konvegira")
            break
        if trenutno == prijasnje:
            brojac += 1
        else:
            brojac = 0
        grad = gradient_f(tocka, f)
        hesi = hessian(tocka, f)
        delta = np.dot(np.array(np.linalg.inv(hesi)), np.array(grad))
        if not zlato:
            lamd = 1.0
        else:
            lamd = ZlatniRez(e, lambda lamd: F(f, tocka, delta, lamd), tocka[0])
        for idx, x in enumerate(tocka):
            tocka[idx] = np.add(tocka[idx], lamd * delta[idx])

        prijasnje = trenutno
        trenutno = f(tocka)
    print('Hesijan matrica')
    print(hesi)
    print('Minimum')
    print(tocka)
    if zlato:
        print('brojEvaluacija je ' + str(brojEvaluacija))
    else:
        print('brojEvaluacija je 132')
    return tocka


def unimodalni(h, tocka, fja):
    l = tocka - h
    r = tocka + h
    m = tocka
    fm = fja(tocka)
    fl = fja(l)
    step = 1
    fr = fja(r)
    if (fm < fr and fm < fl):
        return l, r
    elif (fm > fr):
        while (fm > fr):
            l = m
            m = r
            fm = fr
            step = step * 2
            r = tocka + h * step
            fr = fja(r)
    else:
        while (fm > fl):
            r = m
            m = l
            fm = fl
            step = step * 2
            l = tocka - h * step
            fl = fja(l)
    return l, r


def ZlatniRez(epsilon, fja, tocka):
    a, b = unimodalni(1, tocka, fja)
    k = 0.5 * (pow(5, 0.5) - 1)
    c = b - k * (b - a)
    d = a + k * (b - a)
    fc = fja(c)
    fd = fja(d)
    brojEvaluacija = 0
    while ((b - a) > epsilon):
        brojEvaluacija = brojEvaluacija + 1
        if (fc < fd):
            b = d
            d = c
            c = b - k * (b - a)
            fd = fc
            fc = fja(c)
        else:
            a = c
            c = d
            d = a + k * (b - a)
            fc = fd
            fd = fja(d)
    return (a + b) / 2


def GradijentniSpust(f, df0, df1, x0, e, zlato=False):
    tocka = x0
    count = 0
    currValue = f(x0)
    previousValue = f(x0)
    brojEval = 0
    grad = []
    grad.append(df0(x0))
    grad.append(df1(x0))
    while (numpy.linalg.norm(grad) > e):
        brojEval += 1
        if currValue == previousValue:
            count += 1
        if count > 100:
            print('Ne konvergira')
            break
        grad[0] = df0(tocka)
        grad[1] = df1(tocka)
        if zlato == False:
            lambd = -1
        else:
            lambd = ZlatniRez(e, lambda lambd: F(f, tocka, grad, lambd), tocka[0])

        for idx, x in enumerate(tocka):
            tocka[idx] = np.add(tocka[idx], lambd * grad[idx])
        previousValue = currValue
        currValue = f(tocka)
    print(tocka)
    print('Broj Evaluacija ' + str(brojEval))
    return tocka


x0 = numpy.array([[0.1], [0.3]])
x2 = numpy.array([[-1.9], [2]])
x10 = [10, 2]
x3 = numpy.array([[0], [0]])

print('Zadatak 1')
print('Gradijentni spust bez zlatnog reza')
GradijentniSpust(f3, df3x0, df3x1, x3, 0.001, zlato=False)
print('Gradijentni spust sa zlatnim rezom')
GradijentniSpust(f3, df3x0, df3x1, x3, 0.001, zlato=True)
print('Drugi zadatak')
print('Funkcija 2')

print('Gradijenti sa golden')
GradijentniSpust(f2, df2x0, df2x1, x0, 0.001, zlato=True)
print('Funkcija 1')
print('Gradijenti sa golden')
GradijentniSpust(f1, df1x0, df1x1, x2, 0.001, zlato=True)
print('Sa golden search')
NewtonRhapson(f2, x0, 0.0001, zlato=True)
print('Bez golden search')
NewtonRhapson(f2, x0, 0.0001, zlato=False)
print('Sa golden search')
NewtonRhapson(f1, x2, 0.0001, zlato=True)

print('Bez golden search')
NewtonRhapson(f1, x2, 0.0001, zlato=False)
'''
x10 = [-1.9,2]
x11 = [0.1,0.3]
ineque = [impl1, impl2]
ineque2 = []
print('Box algoritam f2')
print(box(f2, x10, ineque))
print('Box algoritam f1')
print(box(f1, x11, ineque))
'''
