import numpy as np
import math
import matplotlib.pyplot as plt
def analitickoRj(x,T,maxT):
    x_k = np.array([[1.0],[1.0]])
    t=0
    xs = []
    for r in range(np.size(x,0)):
        xs.append([x[r,0]])
    while t<maxT:
        x_kl = x_k
        x_kl[0]= x_k[0]*math.cos(t)+x_k[1]*math.sin(t)
        x_kl[1]= x_k[1]*math.cos(t)-x_k[0]*math.sin(t)
        for idx,r in enumerate(xs):
            r.append(x_kl[idx,0])
        print('t= ',t , 'x_k= ',x_kl)
        t+=T
        x_k = x_kl
    return xs
def crtaj(xs,ts):
    for idx,row in enumerate(xs):
        plt.plot(ts, row, label = 'x' + str(idx+1))
    plt.legend(loc='best')
    plt.grid()
    plt.show()


def runge_kutta(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts =[]
    ts.append(t)
    while t<maxT:
        k1 = np.matmul(A,x_k)+np.matmul(B,r)
        meda = x_k+k1*(T/2)
        k2 = np.matmul(A,meda)+np.matmul(B,r)
        meda = x_k+k2*(T/2)
        k3 = np.matmul(A,meda)+np.matmul(B,r)
        meda = x_k+k3*(T)
        k4 = np.matmul(A,meda)+np.matmul(B,r)
        x_k1 = x_k+ (k1+k2*2+k3*2+k4)*(T/6)
        for idx,rr in enumerate(xs):
            rr.append(x_k1[idx,0])
        t+= T
        print('t=',t,"\nx_k = ",x_k)
        x_k = x_k1
        ts.append(t)
    crtaj(xs,ts)
    return x_k
def runge_kutta4(A,x,B,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts =[]
    ts.append(t)
    while t<maxT:
        r = np.array([[t],[t]])
        k1 = np.matmul(A,x_k)+np.matmul(B,r)
        meda = x_k+k1*(T/2)
        r = np.array([[t+T/2], [t+T/2]])
        k2 = np.matmul(A,meda)+np.matmul(B,r)
        meda = x_k+k2*(T/2)
        r = np.array([[t + T / 2], [t + T / 2]])
        k3 = np.matmul(A,meda)+np.matmul(B,r)
        meda = x_k+k3*(T)
        r = np.array([[t + T ], [t + T ]])
        k4 = np.matmul(A,meda)+np.matmul(B,r)
        x_k1 = x_k+ (k1+k2*2+k3*2+k4)*(T/6)
        for idx,rr in enumerate(xs):
            rr.append(x_k1[idx,0])
        t+= T
        print('t=',t,"\nx_k = ",x_k)
        x_k = x_k1
        ts.append(t)
    crtaj(xs,ts)
    return x_k
def trapezni(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    meda = U-A*(T/2)
    meda2 = U+A*(T/2)
    R = np.matmul(np.linalg.inv(meda),meda2)
    S = np.matmul(np.linalg.inv(meda),B*T/2)
    while t<maxT:
        x_k1 = np.matmul(R,x_k)+np.matmul(S,r+r)
        for idx,rr in enumerate(xs):
            rr.append(x_k1[idx,0])
        t+=T
        print('t=',t,"\nx_k = ",x_k1)
        x_k = x_k1
        ts.append(t)
    crtaj(xs,ts)
    return x_k
def trapezni4(A,x,B,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    meda = U-A*(T/2)
    meda2 = U+A*(T/2)
    R = np.matmul(np.linalg.inv(meda),meda2)
    S = np.matmul(np.linalg.inv(meda),B*T/2)
    while t<maxT:
        r = np.array([[2*t+T],[2*t+T]])
        x_k1 = np.matmul(R,x_k)+np.matmul(S,r)
        for idx,rr in enumerate(xs):
            rr.append(x_k1[idx,0])
        t+=T
        print('t=',t,"\nx_k = ",x_k1)
        x_k = x_k1
        ts.append(t)
    crtaj(xs,ts)
    return x_k

def euler(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    brojac=0
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    medu = U+ A*T
    N = T*B
    while t<maxT:
        x_kl = np.matmul(medu,x_k)+np.matmul(N,r)
        for idx,rr in enumerate(xs):
            rr.append(x_kl[idx,0])
        t += T
        print('t= ',t,"\nx_k = ",x_kl)
        x_k = x_kl
        ts.append(t)
    crtaj(xs,ts)

    return x_k
def euler4(A,x,B,T,maxT):
    x_k = x
    xs = []
    brojac=0
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    medu = U+ A*T
    N = T*B
    while t<maxT:
        r = np.array([[t],[t]])
        x_kl = np.matmul(medu,x_k)+np.matmul(N,r)
        for idx,rr in enumerate(xs):
            rr.append(x_kl[idx,0])
        t += T
        print('t= ',t,"\nx_k = ",x_kl)
        x_k = x_kl
        ts.append(t)
    crtaj(xs,ts)

    return x_k

def obrnuti_euler(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t=0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    L = U - A*(T)
    P = np.linalg.inv(L)
    dsads = P
    ql = np.matmul(dsads,T*B)

    while t < maxT:
        x_kl = np.matmul(P, x_k)+np.matmul(ql,r)
        for idx, rr in enumerate(xs):
            rr.append(x_kl[idx, 0])
        #print('t= ', t, "\nx_k = ", x_kl)
        x_k = x_kl
        t += T
        ts.append(t)
    crtaj(xs,ts)
    return x_k
def obrnuti_euler4(A,x,B,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t=0
    ts = []
    ts.append(t)
    U = np.identity(np.size(A, 0))
    L = U - A*(T)
    P = np.linalg.inv(L)
    dsads = P
    ql = np.matmul(dsads,T*B)

    while t < maxT:
        r = np.array([[t+T],[t+T]])
        x_kl = np.matmul(P, x_k)+np.matmul(ql,r)
        for idx, rr in enumerate(xs):
            rr.append(x_kl[idx, 0])
        #print('t= ', t, "\nx_k = ", x_kl)
        x_k = x_kl
        t += T
        ts.append(t)
    crtaj(xs,ts)
    return x_k
def PECECE(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    U = np.identity(np.size(A, 0))
    N = T*B
    while t<maxT:
        x_kl =x_k +T*(np.matmul(A,x_k)+np.matmul(B,r))
        for idx,rr in enumerate(xs):
            rr.append(x_kl[idx,0])
        #print('t= ',t,"\nx_k = ",x_kl)
        x_kl += T*obrnuti_euler2(A,x_kl,B,r,T,T/2)
        x_k = x_kl
        t+=T
    return x_k
def PECE(A,x,B,r,T,maxT):
    x_k = x
    xs = []
    for rr in range(np.size(x,0)):
        xs.append([x[rr,0]])
    t = 0
    U = np.identity(np.size(A, 0))
    N = T*B
    while t<maxT:
        x_kl = np.matmul((U+A*T),x_k)+np.matmul(N,r)

        for idx,rr in enumerate(xs):
            rr.append(x_kl[idx,0])
        #print('t= ',t,"\nx_k = ",x_kl)

        x_kl = x_kl +T*trapezni(A,x_kl,B,r,T,T/2)
        x_k = x_kl
        t+=T
    return x_k
A = np.array([[0,1],[-1,0]])
x = np.array([[1],[1]])

B = np.array([[0,0],[0,0]])
r = np.array([[0],[0]])
A = np.array([[0,1],[-1,0]])
x = np.array([[1],[1]])
#runge_kutta(A,x,B,r,0.01,10)
#trapezni(A,x,B,r,0.01,10)
#euler(A,x,B,r,0.01,10)
#obrnuti_euler(A,x,B,r,0.01,10)
#print(PECECE(A,x,B,r,0.01,10))
#print(PECE(A,x,B,r,0.01,10))
#DRUGI
A = np.array([[0,-2],[1,-3]])
x = np.array([[1],[3]])
#runge_kutta(A,x,B,r,0.01,10)
#trapezni(A,x,B,r,0.01,10)
#euler(A,x,B,r,0.01,10)
#obrnuti_euler(A,x,B,r,0.01,10)
#TRECI
A = np.array([[0,-2],[1,-3]])
x = np.array([[1],[3]])
B = np.array([[2,0],[0,3]])
r = np.array([[1],[1]])
#euler(A,x,B,r,0.01,10)
#obrnuti_euler(A,x,B,r,0.01,10)
#trapezni(A,x,B,r,0.01,10)
#runge_kutta(A,x,B,r,0.01,10)
#CETVRTI
A = np.array([[1,-5],[1,-7]])
x = np.array([[-1],[3]])
B = np.array([[5,0],[0,3]])
#euler4(A,x,B,0.01,10)
#obrnuti_euler4(A,x,B,0.01,10)
#trapezni4(A,x,B,0.01,10)
#runge_kutta4(A,x,B,0.01,10)