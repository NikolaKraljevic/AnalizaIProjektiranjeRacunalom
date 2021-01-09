import numpy
import random
import math
import matplotlib.pyplot as plt
def f1(x):
    return 100*pow((x[1]-pow(x[0],2)),2)+pow(1-x[0],2)
def f3(x):
    sum = 0
    for idx,xi in enumerate(x):
        sum+= (xi-idx+1)**2
    return sum
def f6(x):
    sum = 0
    for xi in x:
        sum +=pow(xi,2)
    n1 = math.sin(math.sqrt(sum))**2-0.5
    n2 = (1+0.001*sum)**2
    return 0.5+n1/n2
def f7(x):
    sum = 0
    for xi in x:
        sum +=xi**2
    n1 = sum**0.25
    sin_x = 50*sum**0.1
    n2 = 1+math.sin((sin_x))**2
    return n1*n2
class Gen:
    def __init__(self,matrica,VrijednostFje):
        self.matrica = matrica
        self.vrijednostFje = VrijednostFje
    def __repr__(self):
        return "%s %s"%(self.matrica,self.vrijednostFje)
def IzracunajPreciznost(dg,gg,p):
    broj = math.floor(1+(gg-dg)*pow(10,p))
    n = math.ceil(math.log(broj,10)/math.log(2,10))
    return n
def makeOne(duljina):
    matrica = []
    for i in range(duljina):
        x = random.random()
        if x>0.5:
            matrica.append(1)
        else:
            matrica.append(0)
    return matrica
def IzracunajVrijednostFje(f,matrica,dg,gg,varijable):
    vrijednost = [0 for i in range(len(varijable))]
    prije = 0
    poslije = 0
    for k in range(len(varijable)):
        poslije += varijable[k]
        razlika = poslije-prije
        for i in range(razlika):
            if matrica[prije+i]==1:
                vrijednost[k] += pow(2,razlika-i-1)
        vrijednost[k] = dg+vrijednost[k]/(pow(2,razlika))*(gg-dg)
        prije += varijable[k]
    return f(vrijednost)
def VrFj(elem):
    return elem.vrijednostFje
def Krizanje(mat1,mat2):
    dijete = []
    k = random.randint(0,len(mat1))
    for i in range(k):
        dijete.append(mat1[i])
    for i in range(k,len(mat1)):
        dijete.append(mat2[i])
    return dijete
def KrizanjePom(mat1,mat2):
    x = random.random()
    dijete = []
    for i in range(len(mat1)):
        rjesenje = x*mat1[i]+(1-x)*mat2[i]
        dijete.append(rjesenje)
    return dijete
def Mutiraj(dijete,mut):
    brojmut = 0
    for i in range(len(dijete)):
        k = random.random()
        if k<=mut:
            brojmut+=1
            if dijete[i]==1:
                dijete[i]=0
            else:
                dijete[i]=1

    return dijete
def MutirajPom(dijete,mut):
    for i in range(len(dijete)):
        k = random.random()
        if k<mut:
            val = random.uniform(-1,1)
        else:
            val = 0
        dijete[i]+=val
    return dijete
def GenetskiAlg(f,dg,gg,varijable,n,p,mut):
    duljinaBin = 0
    brojac = 0
    prijasnaNajbolja = 0
    novaNajbolja= 0
    jedinke = []
    varijable = []
    prec = IzracunajPreciznost(dg,gg,p)
    for i in range(n):
        varijable.append(prec)
        duljinaBin += prec
    for i in range(n):
        matrica = makeOne(duljinaBin)
        vrijednost = IzracunajVrijednostFje(f,matrica,dg,gg,varijable)
        jedinke.append(Gen(matrica,vrijednost))

    jedinke.sort(key=VrFj)
    while(jedinke[0].vrijednostFje>0.001):
        turnir = []
        prijasnaNajbolja = novaNajbolja
        novaNajbolja = jedinke[0].vrijednostFje
        if brojac==100000:
            print('Nije pronadjen min u 100000')
            print(jedinke[0])
            return jedinke[0].vrijednostFje
        uturniru = []
        for i in range(3):
            x = random.randint(0,n-1)
            while x in uturniru:
                x = random.randint(0,n-1)
            uturniru.append(x)
            turnir.append(jedinke[x])
        turnir.sort(key=VrFj)
        dijete = Krizanje(turnir[0].matrica,turnir[1].matrica)
        if(jedinke[0].vrijednostFje<1.5):
            mut2 = mut
        else:
            mut2 = mut
        NovoDijete = Mutiraj(dijete,mut2)
        Val = IzracunajVrijednostFje(f,NovoDijete,dg,gg,varijable)
        del jedinke[-1]
        nekaj = Gen(NovoDijete,Val)
        jedinke.append(nekaj)
        jedinke.sort(key=VrFj)
        brojac+=1
    print(brojac)
    print(jedinke[0])
    return jedinke[0].vrijednostFje


def GenetskiAlgPomicnaTocka(f,dg,gg,varijable,n,p,mut,brojX,velicinaTurnira):
    brojac =0
    jedinke =[]
    varijable = []
    matrica = []
    for k in range(n):
        matrica = []
        for i in range(brojX):
            matrica.append(random.uniform(dg,gg))
        vrijednost = f(matrica)
        jedinke.append(Gen(matrica,vrijednost))
    jedinke.sort(key=VrFj)
    while(jedinke[0].vrijednostFje>0.001):
        turnir = []
        brojac +=1
        if brojac == 100000:
            print('Nije pronadjen min nakon 100000 evaluacija')
            print(jedinke[0])
            return(jedinke[0].vrijednostFje)
        uturniru = []
        for i in range(3):
            x = random.randint(0, n - 1)
            while x in uturniru:
                x = random.randint(0, n - 1)
            uturniru.append(x)
            turnir.append(jedinke[x])
        turnir.sort(key=VrFj)
        dijete = KrizanjePom(turnir[0].matrica,turnir[1].matrica)
        dijete = MutirajPom(dijete,mut)
        val = f(dijete)
        del jedinke[-1]
        nekaj = Gen(dijete,val)
        jedinke.append(nekaj)
        jedinke.sort(key=VrFj)
    print(brojac)
    print(jedinke[0])
    return jedinke[0].vrijednostFje
def prvi():
    varijable = []
    print('Prvi zadatak')
    print('Prva fja')

    GenetskiAlg(f1,-50,150,varijable,5,3,0.1)
    GenetskiAlgPomicnaTocka(f1,-50,150,varijable,5,3,0.1,2,3)
    print('Treca fja')

    GenetskiAlg(f3, -50, 150, varijable, 5, 3, 0.1)
    GenetskiAlgPomicnaTocka(f3, -50, 150, varijable, 5, 3, 0.1, 5, 3)
    print('Sesta fja')

    GenetskiAlg(f6, -50, 150, varijable, 5, 3, 0.1)
    GenetskiAlgPomicnaTocka(f3, -50, 150, varijable, 5, 3, 0.1, 2, 3)
def drugi():
    varijable = []
    dim = [1,3,6,10]
    for i in dim:
        rjzai[z] = []
        for k in range(10):
            rj =GenetskiAlgPomicnaTocka(f6,-50,150,varijable,10,3,0.1,i,3)
            rjzai[z].append(rj)
        z+=1
    fig1,ax1 = plt.subplots()
    ax1.boxplot(rjzai)
    plt.show()

def treci():
    varijable = []
    rjzai = [[] for i in range(4)]
    rjzai[0] = []
    rjzai[1]= []
    for k in range(20):
        rj = GenetskiAlgPomicnaTocka(f6,-50,150,varijable,10,3,0.1,3,3)
        ry =  GenetskiAlgPomicnaTocka(f6,-50,150,varijable,10,3,0.1,3,3)
        rjzai[0].append(rj)
        rjzai[1].append(ry)
    rjzai[2]= []
    rjzai[3]=[]
    for k in range(20):
        rj = GenetskiAlgPomicnaTocka(f7,-50,150,varijable,10,3,0.1,3,3)
        ry =  GenetskiAlgPomicnaTocka(f7,-50,150,varijable,10,3,0.1,3,3)
        rjzai[2].append(rj)
        rjzai[3].append(ry)
    fig1,ax1 = plt.subplots()
    ax1.boxplot(rjzai)
    plt.show()




def cetvrti():
    varijable = []
    rjzai = [ [] for i in range(4)]
    z = 0
    populacija = [30,50,100,200]
    mutacija = [0.1,0.3,0.6,0.9]
    for i in populacija:
        rjzai[z] = []
        for k in range(10):
            rj =GenetskiAlgPomicnaTocka(f6,-50,150,varijable,i,3,0.1,2,3)
            rjzai[z].append(rj)
        z+=1
    fig1,ax1 = plt.subplots()
    ax1.boxplot(rjzai)
    plt.show()
    varijable = []
    rjzai = [ [] for i in range(4)]
    z = 0
    populacija = [30,50,100,200]
    mutacija = [0.1,0.3,0.6,0.9]
    for i in mutacija:
        rjzai[z] = []
        for k in range(10):
            rj =GenetskiAlgPomicnaTocka(f6,-50,150,varijable,100,3,i,2,3)
            rjzai[z].append(rj)
        z+=1
    fig1,ax1 = plt.subplots()
    ax1.boxplot(rjzai)
    plt.show()
def peti():
    velicinaTurnira = [3,5,7,9]
    varijable = []
    rjzai = [[] for i in range(4)]
    z = 0
    for i in velicinaTurnira:
        rjzai[z] = []
        for k in range(50):
            rj = GenetskiAlgPomicnaTocka(f6, -50, 150, varijable, 100, 3, 0.1, 2, i)
            rjzai[z].append(rj)
        z += 1
    fig1, ax1 = plt.subplots()
    ax1.boxplot(rjzai)
    plt.show()

#prvi()
#drugi()
#treci()
cetvrti()
#peti()