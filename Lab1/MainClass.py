
def SaveResult(fileName,matrica):
    f = open(fileName,"w")
    for i in range(len(matrica.matrica)):
        for j in range(len(matrica.matrica[0])):
            f.write(" " +str(matrica.matrica[i][j]))
    f.close
def MatrixWithScalar(mat,scalar):
    for i in range(mat.rows):
        for j in range(mat.columns):
            mat.matrica[i][j]*=scalar
def MatrixMul(mat1,mat2):
    result = [[0 for i in range(mat2.rows)] for k in range(mat1.columns)]
    for i in range(mat1.columns):
        for k in range(mat2.rows):
            for j in range(mat2.columns):
                    result[i][k]+=mat1.matrica[i][j]*mat2.matrica[j][k]
    return result
def SplitLU(matrica):
    LMatrica = [[0 for i in range(matrica.rows)] for k in range(matrica.columns) ]
    for i in range(matrica.rows):
        for j in range(i+1):
            if i == j:
                LMatrica[i][i] = 1
            else:
                LMatrica[i][j] = matrica.matrica[i][j]
    UMatrica =  [[0 for i in range(matrica.rows)] for k in range(matrica.columns) ]
    for i in range(matrica.rows):
        for j in range(i,matrica.rows):
            UMatrica[i][j] = matrica.matrica[i][j]
    return LMatrica,UMatrica
def LUDecomposition(matrica):
    for i in range(matrica.columns-1):
        if (matrica.matrica[i][i] <0.0000001):
            print ("Pivot je približno 0")
            break
        for j in range(i+1,matrica.rows):
            matrica.matrica[j][i] /= matrica.matrica[i][i]
            for k in range(i+1,matrica.rows):
                matrica.matrica[j][k] -=matrica.matrica[j][i]*matrica.matrica[i][k]

def LUPDecomposition(matrica):
    brojPermutacija = 0
    Pivot = [[0 for d in range(matrica.rows)] for k in range(matrica.columns)]
    for r in range(matrica.rows):
        Pivot[r][r]=1
    for i in range(matrica.rows-1):

        pamti = i
        for j in range(i+1,matrica.columns):
            if (abs(matrica.matrica[j][i])>abs(matrica.matrica[pamti][i])):
                pamti = j

        if(abs(matrica.matrica[pamti][i])<0.000001):
            return None
        if(pamti != i):
            matrica.matrica[i],matrica.matrica[pamti]=matrica.matrica[pamti],matrica.matrica[i]
            Pivot[i],Pivot[pamti]=Pivot[pamti],Pivot[i]
            brojPermutacija += 1
        for j in range(i+1,matrica.rows):
            matrica.matrica[j][i] /= matrica.matrica[i][i]
            for k in range(i+1,matrica.rows):
                matrica.matrica[j][k] -=matrica.matrica[j][i]*matrica.matrica[i][k]
    return brojPermutacija,Pivot




def SubstitutionForward(matrica,vektor):
    for i in range(matrica.columns-1):
        for j in range(i+1,matrica.columns):
            vektor.matrica[j][0] -= matrica.matrica[j][i]*vektor.matrica[i][0]

def SubstitionBackwards(matrica,vektor):

    for i in range(matrica.columns-1,-1,-1):
        vektor.matrica[i][0]/=matrica.matrica[i][i]
        nekibroj = i-1
        for j in range(nekibroj+1):

            vektor.matrica[j][0]-= matrica.matrica[j][i]*vektor.matrica[i][0]

class Matrix:

    def __init__(self,text):
        file = open(text)
        file2 = open(text)
        j=0
        for line in file:
            elements = line.split()
            self.rows = len(elements)
            j+=1
        self.columns = j
        self.matrica = [[0 for i in range(self.rows)] for k in range(self.columns)]
        j=0
        for line in file2:
            elements = line.split()
            for i in range(len(elements)):
                self.matrica[j][i]=float(elements[i])
            j+=1
class Pivots:
    def __init__(self,matrica):
        self.rows = len(matrica)
        self.columns = len(matrica[0])
        self.matrica = matrica

Matrica1 = Matrix('Example.txt')
Vektor = Matrix('ExampleVector.txt')
LUDecomposition(Matrica1)
SubstitutionForward(Matrica1,Vektor)
SubstitionBackwards(Matrica1,Vektor)
print("Zadatak 2")
Matrica2LU = Matrix('zad2Matrix.txt')
Matrica2LUP = Matrix('zad2Matrix.txt')
Vektor2LU = Matrix('zad2Vector.txt')
Vektor2LUP = Matrix('zad2Vector.txt')
LUDecomposition(Matrica2LU)
brojPermutacija2,P =LUPDecomposition(Matrica2LUP)
P = Pivots(P)
print(P.matrica)
print(Matrica2LUP.matrica)
Vektor2LUP = MatrixMul(P,Vektor2LUP)
Vektor2LUP = Pivots(Vektor2LUP)
'''
SubstitutionForward(Matrica2LU,Vektor2LU)
SubstitionBackwards(Matrica2LU,Vektor2LU)
'''
SubstitutionForward(Matrica2LUP,Vektor2LUP)
SubstitionBackwards(Matrica2LUP,Vektor2LUP)
'''
print ("LU Dekompozicija")
print (Matrica2LU.matrica)
print(Vektor2LU.matrica)
'''
print ("LUP Dekompozicija")
print (Matrica2LUP.matrica)
print(Vektor2LUP.matrica)
SaveResult("zad2RJ.txt",Vektor2LUP)
print("Zadatak 3")
Matrica3LU = Matrix('zad3Matrix')
Matrica3LUP = Matrix('zad3Matrix')
LUDecomposition(Matrica3LU)
print(Matrica3LU.matrica)
brojPermutacija3,P3 = LUPDecomposition(Matrica3LUP)
LMat3,UMat3 = SplitLU(Matrica3LUP)
print ("P Matrica ")
print (P3)
print ("L Matrica")
print (LMat3)
print ("U Matrica")
print (UMat3)
print("Zadatak 4")
Matrica4LU = Matrix('zad4Matrix.txt')
Matrica4LUP = Matrix('zad4Matrix.txt')
Vektor4LU = Matrix('zad4Vector.txt')
Vektor4LUP = Matrix('zad4Vector.txt')
LUDecomposition(Matrica4LU)
brojPermutacija4,P4 =LUPDecomposition(Matrica4LUP)
P4= Pivots(P4)
Vektor4LUP = MatrixMul(P4,Vektor4LUP)
Vektor4LUP = Pivots(Vektor4LUP)
'''
SubstitutionForward(Matrica4LU,Vektor4LU)
SubstitionBackwards(Matrica4LU,Vektor4LU)
'''
SubstitutionForward(Matrica4LUP,Vektor4LUP)
SubstitionBackwards(Matrica4LUP,Vektor4LUP)
'''
print ("LU Dekompozicija : LU Matrica")
print (Matrica4LU.matrica)
print ("Rješenje")
print (Vektor4LU.matrica)
'''
print ("LUP Dekompozicija : LUP Matrica")
print (Matrica4LUP.matrica)
print ("Rješenje")
print (Vektor4LUP.matrica)
print("Zadatak 5")
Matrica5LU = Matrix('zad5Matrix.txt')
Matrica5LUP = Matrix('zad5Matrix.txt')
Vektor5LU = Matrix('zad5Vector.txt')
Vektor5LUP = Matrix('zad5Vector.txt')
LUDecomposition(Matrica5LU)
brojPermutacija5,P5 =LUPDecomposition(Matrica5LUP)
P5= Pivots(P5)
Vektor5LUP = MatrixMul(P5,Vektor5LUP)
Vektor5LUP = Pivots(Vektor5LUP)
'''
SubstitutionForward(Matrica5LU,Vektor5LU)
SubstitionBackwards(Matrica5LU,Vektor5LU)
'''
SubstitutionForward(Matrica5LUP,Vektor5LUP)
SubstitionBackwards(Matrica5LUP,Vektor5LUP)
'''
print ("LU Dekompozicija : LU Matrica")
print (Matrica5LU.matrica)
print ("Rješenje")
print (Vektor5LU.matrica)
'''
print ("LUP Dekompozicija : LUP Matrica")
print (Matrica5LUP.matrica)
print ("Rješenje")
print (Vektor5LUP.matrica)
Matrica6LU = Matrix('zad6Matrix.txt')
Matrica6LUP = Matrix('zad6Matrix.txt')
Vektor6LU = Matrix('zad6Vector.txt')
Vektor6LUP = Matrix('zad6Vector.txt')
LUDecomposition(Matrica6LU)
brojPermutacija6,P6 =LUPDecomposition(Matrica6LUP)
P6= Pivots(P6)
Vektor6LUP = MatrixMul(P6,Vektor6LUP)
Vektor6LUP = Pivots(Vektor6LUP)
SubstitutionForward(Matrica6LU,Vektor6LU)
SubstitionBackwards(Matrica6LU,Vektor6LU)
SubstitutionForward(Matrica6LUP,Vektor6LUP)
SubstitionBackwards(Matrica6LUP,Vektor6LUP)
print ("Zadatak 6")
print ("LU Dekompozicija : LU Matrica")
print (Matrica6LU.matrica)
print ("Rješenje")
print (Vektor6LU.matrica)
print ("LUP Dekompozicija : LUP Matrica")
print (Matrica6LUP.matrica)

print ("Rješenje")
print (Vektor6LUP.matrica)
print ("Zadatak 7")
Matrica7LUP = Matrix('zad7Matrix.txt')
Vektor7LUP = Matrix('zad7Matrix.txt')

brojPermutacija7,P7 =LUPDecomposition(Matrica7LUP)
P7= Pivots(P7)
Vektor7LUP = MatrixMul(P7,Vektor7LUP)
Vektor7LUP = Pivots(Vektor7LUP)
novi = [[0 for i in range(1)] for k in range(3) ]
novi = Pivots(novi)
for i in range(Matrica7LUP.rows):
    for j in range(Matrica7LUP.rows):
        novi.matrica[j][0]=P7.matrica[j][i]
    SubstitutionForward(Matrica7LUP,novi)
    SubstitionBackwards(Matrica7LUP,novi)
    for j in range(Matrica7LUP.rows):
        P7.matrica[j][i] = novi.matrica[j][0]
print (P7.matrica)
print ("Zadatak 8")
Matrica8LUP = Matrix('zad8Matrix.txt')
brojPermutacija8,P8 =LUPDecomposition(Matrica8LUP)
P8= Pivots(P8)
novi = [[0 for i in range(1)] for k in range(3) ]
novi = Pivots(novi)
for i in range(Matrica8LUP.rows):
    for j in range(Matrica8LUP.rows):
        novi.matrica[j][0]=P8.matrica[j][i]
    SubstitutionForward(Matrica8LUP,novi)
    SubstitionBackwards(Matrica8LUP,novi)
    for j in range(Matrica8LUP.rows):
        P8.matrica[j][i] = novi.matrica[j][0]
print (P8.matrica)
print ("Zadatak 9")
Matrica9LUP = Matrix('zad9Matrix.txt')
brojPermutacija9,P9 = LUPDecomposition(Matrica9LUP)
L9Mat,U9Mat = SplitLU(Matrica9LUP)
DetU =1
for i in range(len(U9Mat)):
    DetU*=U9Mat[i][i]
print("Determinanta")
print(pow(-1,brojPermutacija9)*DetU)
print ("Zadatak 10")

Matrica10LUP = Matrix('zad10Matrix.txt')
brojPermutacija10,P10 = LUPDecomposition(Matrica10LUP)
L10Mat,U10Mat = SplitLU(Matrica10LUP)
DetU =1
for i in range(len(U10Mat)):
    DetU*=U10Mat[i][i]
print("Determinanta")
print(pow(-1,brojPermutacija9)*DetU)

