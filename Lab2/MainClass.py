import numpy as np
import copy
import random
from math import sin,sqrt
def f10(x):
    return pow(x-3,2)
def f1(x):
    return 100*pow((x[1]-pow(x[0],2)),2)+pow(1-x[0],2)
def f2(x):
    return pow(x[0]-4,2)+4*pow(x[1]-2,2)
def f3(x):
    zbroj = 0
    for i in range(len(x)):
        zbroj += pow(x[i]-(i+1),2)
    return zbroj
def f4(x):
    return abs(pow(x[0],2)-pow(x[1],2))+ pow((pow(x[0],2)+pow(x[1],2)),0.5)
def f5(X):
    suma = 0
    for x in X:
        suma += x ** 2
    return 0.5 + (sin(sqrt(suma)) ** 2 - 0.5) / (1 + 0.001 * suma) ** 2
def unimodalni(h,tocka,fja):
    l = tocka-h
    r= tocka+h
    m= tocka
    fm = fja(tocka)
    fl= fja(l)
    step = 1
    fr = fja(r)
    if (fm<fr and fm<fl ):
        return l,r
    elif(fm>fr):
        while (fm>fr):
            l=m
            m=r
            fm = fr
            step = step*2
            r= tocka+h*step
            fr= fja(r)
    else:
        while(fm>fl):
            r=m
            m=l
            fm=fl
            step = step*2
            l= tocka - h*step
            fl = fja(l)
    print (l,r)
    return l,r
def unimodalni2(h,tocka,fja,i):
    pomocnatocka = tocka
    l = tocka[i]-h
    r= tocka[i]+h
    m= tocka[i]
    fm = fja(tocka)
    pomocnatocka[i]= l
    fl= fja(pomocnatocka)
    step = 1
    pomocnatocka[i]=r
    fr = fja(pomocnatocka)
    if (fm<fr and fm<fl ):
        return l,r
    elif(fm>fr):
        while (fm>fr):
            l=m
            m=r
            fm = fr
            step = step*2
            r= tocka[i]+h*step
            pomocnatocka[i]=r
            fr= fja(pomocnatocka)
    else:
        while(fm>fl):
            r=m
            m=l
            fm=fl
            step = step*2
            l= tocka[i] - h*step
            pomocnatocka[i] = l
            fl = fja(pomocnatocka)
    return l,r
def ZlatniRez(epsilon,fja,tocka):
    a,b = unimodalni(1,tocka,fja)
    k= 0.5*(pow(5,0.5)-1)
    c = b-k*(b-a)
    d= a+k*(b-a)
    fc = fja(c)
    fd = fja(d)
    brojEvaluacija = 0
    while((b-a)>epsilon):
        brojEvaluacija=brojEvaluacija+1
        if(fc<fd):
            b=d
            d=c
            c = b-k*(b-a)
            fd=fc
            fc = fja(c)
        else:
            a=c
            c=d
            d=a+k*(b-a)
            fc=fd
            fd=fja(d)
    return brojEvaluacija,(a+b)/2
def ZlatniRez2(epsilon,fja,tocka,i):
    a,b = unimodalni2(1,tocka,fja,i)

    k= 0.5*(pow(5,0.5)-1)
    c = b-k*(b-a)
    d= a+k*(b-a)
    tocka[i]=c
    fc = fja(tocka)
    tocka[i]=d
    fd = fja(tocka)
    brojEvaluacija = 0
    while((b-a)>epsilon):
        brojEvaluacija=brojEvaluacija+1
        if(fc<fd):
            b=d
            d=c
            c = b-k*(b-a)
            fd=fc
            tocka[i]= c
            fc = fja(tocka)
        else:
            a=c
            c=d
            d=a+k*(b-a)
            fc=fd
            tocka[i]=d
            fd=fja(tocka)
    return a+b/2
def Koordinatno_trazenje(x0,epsilon,fja):
    x=x0
    xn= []
    brojEvaluacija = 0
    while(True):
        brojEvaluacija+=1
        xs = x
        for i in range(len(x0)):
            x[i]+= ZlatniRez2(epsilon,fja,xs,i)
            dif=0
            for j in range(len(x)):
                dif+=abs(x[j]-xs[j])

        if(dif<=epsilon):
                break
    xn=x
    return xn,brojEvaluacija
def istrazi(xp,Dx,n,fja):
    xi= xp
    for i in range(n):
        P = fja(xi)
        xi[i] += Dx
        N = fja(xi)
        if (N>P):
            xi[i]-=2*Dx
            N = fja(xi)
            if(N>P):
                xi[i]+=Dx
    return xi
def HookeJeeves(x0,fja,dx,epsilon):
    Dx = dx
    xp,xb=x0,x0
    n= len(x0)
    brojEvaluacija = 0
    while(True):
        xn = istrazi(xp,Dx,n,fja)
        if(fja(xn)<fja(xb)):
            xp=2*xn-xb
            xb = xn
        else:
            Dx=Dx/2
            xp = xb

        brojEvaluacija+=1
        if(epsilon>Dx):
            break
    return brojEvaluacija,xb
def NelderMead(x0,alpha,beta,gama,epsilon,sigma,fja):
    dim = len(x0)
    prijasni_najbolji  = fja(x0)
    res = [[x0,prijasni_najbolji]]

    for i in range(dim):
        x = copy.copy(x0)
        x[i] = x[i]+1
        score = fja(x)
        res.append([x, score])
    BrojEvaluacija = 0
    while(True):

        BrojEvaluacija+=1
        res.sort(key=lambda x: x[1])
        if(abs(res[0][1]-res[-1][1])<epsilon):
            break
        #centroid
        xc = [0.]*dim
        for tup in res[:-1]:
            for i,c in enumerate(tup[0]):
                xc[i]+=c/(len(res)-1)
        #refleksija
        xr = xc +alpha*(xc-res[-1][0])
        rscore = fja(xr)
        if (res[0][1]<=rscore and rscore<res[-1][1]):
            del res[-1]
            res.append([xr,rscore])
            continue
        #ekspanzija
        if rscore<res[0][1]:
            xe = xc + gama*(xc -res[-1][0])
            escore = fja(xe)
            if escore<rscore:
                del res[-1]
                res.append([xe,escore])
                continue
            else:
                del res[-1]
                res.append([xr,rscore])
                continue
        #kontrakcija
        xk = xc + beta*(xc-res[-1][0])
        kscore = fja(xk)
        if kscore<res[-1][1]:
            del res[-1]
            res.append(([xk,kscore]))
            continue

        x1 = res[0][0]
        nres = []
        for tup in res:
            redx =  x1 + sigma*(tup[0]-x1)
            score = fja(redx)
            nres.append([redx,score])
        res = nres

    return res[0][0],BrojEvaluacija

Dx,e =0.5,0.000001
print('Zadatak 1 pocetna tocka 10')
BrojEvaluacija,minimium = ZlatniRez(e,f10,10)
print(BrojEvaluacija,minimium)
print('Zadatak 1 pocetna tocka 15')
BrojEvaluacija,minimium = ZlatniRez(e,f10,15)
print(BrojEvaluacija,minimium)
print('Zadatak 1 pocetna tocka 25')
BrojEvaluacija,minimium = ZlatniRez(e,f10,25)
print(BrojEvaluacija,minimium)
x0 = np.array([1,1])
xvek =np.array([4,4,4,4,4])
print ('2. Zadatak')
print('Funkcija 1')
x0 = np.array([-1.9,2])
min,broj=NelderMead(x0,1,0.5,2,0.000001,0.5,f1)
print('Needler Mead')
print(min,broj)
broj10,min15= HookeJeeves(x0,f1,Dx,e)
print("Hooke Jeeves")
print(broj10,min15)
min5,broj5=Koordinatno_trazenje(x0,e,f1)
print ('Koordinatno trazenje')
print (min3,broj5)
print('Funkcija 2')
x0 = np.array([0.1,0.3])
min,broj=NelderMead(x0,1,0.5,2,0.000001,0.5,f2)
print('Needler Mead')
print(min,broj)
broj2,min3= HookeJeeves(x0,f2,Dx,e)
print("Hooke Jeeves")
print(broj2,min3)
min5,broj5=Koordinatno_trazenje(x0,e,f2)
print ('Koordinatno trazenje')
print (min3,broj5)
print('Funkcija 3')
x0 = np.array([1,2,3,4,5])
min,broj=NelderMead(x0,1,0.5,2,0.000001,0.5,f3)
print('Needler Mead')
print(min,broj)
broj2,min3= HookeJeeves(x0,f3,Dx,e)
print("Hooke Jeeves")
print(broj2,min3)
min5,broj5=Koordinatno_trazenje(x0,e,f3)
print ('Koordinatno trazenje')
print (min3,broj5)
print('Funkcija 4')
x0 = np.array([1,1])
min,broj=NelderMead(x0,1,0.5,2,0.000001,0.5,f4)
print('Needler Mead')
print(min,broj)
broj2,min3= HookeJeeves(x0,f4,Dx,e)
print("Hooke Jeeves")
print(broj2,min3)
min5,broj5=Koordinatno_trazenje(x0,e,f4)
print ('Koordinatno trazenje')
print (min3,broj5)
print("Zadatak 3")
print('Funkcija 4')
x0 = np.array([20,20])
min,broj=NelderMead(x0,1,0.5,2,0.000001,0.5,f4)
print('Needler Mead')
print(min,broj)
broj2,min3= HookeJeeves(x0,f4,Dx,e)
print("Hooke Jeeves")
print(broj2,min3)
print("Zadatak 4")
x0 = np.array([0.5,0.5])
min10,broj10 = NelderMead(x0,1,0.5,2,0.000001,0.5,f1)
print ("Sigma je 0.5")
print (min10,broj10)
min11,broj11 = NelderMead(x0,1,0.5,2,0.000001,0.75,f1)
print ("Sigma je 0.75")
print(min11,broj11)
min11,broj11 = NelderMead(x0,1,0.5,2,0.000001,0.9,f1)
print ("Sigma je 0.9")
print(min11,broj11)
print("Zadatak 5")
brojac = 0
for i in range(50):
    x0 = np.array([random.randint(-50,50),random.randint(-50,50)])
    min,broj = HookeJeeves(x0,f5,Dx,e)
    if(f5(broj)<0.001):
        brojac+=1

print(brojac)

