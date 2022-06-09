import numpy as np
import matplotlib.pyplot as plt
import pylab





def main():
    T,A,G,B,B1,shtrix,G_shtri_1,G_shtri_2,mat,matrix_U,matrix_Bu,matrix_x0= matrix_read()
    u=72
    w=-3.1
    q=0
    teta=-0.6
    H=420
    delta_e=0
    delta_T=0
    t=0
    w0=np.pi
    dt=0.1
    x=list()
    x.append(u)
    x.append(w)
    x.append(q)
    x.append(teta)
    x.append(H)
    x.append(delta_e)
    x.append(delta_T)
    print("x:")
    for i in x:
        print(i)
    

    A_matrix_shtriX=A_shtrix(A,B,shtrix)
    G_matrix_shtriX=G_shtrix(A,G,G_shtri_1)
    G_2_matrix_shtriX=G_2_shtrix(A,G,G_shtri_2)
    ksi_eto= multiply_T_x(T, x)
    D2udt2=d2udt2(x,A_matrix_shtriX,G,G_matrix_shtriX,t,w0)
    D3hdt3=d3hdt3(x,A_matrix_shtriX,G_matrix_shtriX,G_2_matrix_shtriX,t,w0)
    Ksi_eto_2=A_ksi_eto(ksi_eto,D3hdt3,D2udt2,mat)
    U_bol=multiply_B1_ksi_eto_2(B1, Ksi_eto_2,matrix_U)
    Bu=multiply_B_u(U_bol, matrix_Bu)
    Xr=multiply_x0(A,x,Bu,dt,matrix_x0,shtrix,G_matrix_shtriX,mat,matrix_U,T,G,w0,B1,matrix_Bu,G_2_matrix_shtriX)
    grafics(Xr,dt)
    print("ksi_eto:")
    for i in ksi_eto:
        print(i)
    print("A:")
    for i in A:
        print(i)    
    
    print("B")
    for i in B:
        print(i)    


    print("A_matrix_shtriX:")
    for i in A_matrix_shtriX:
        print(i)    

    
    print("G_matrix_shtriX:")
    for i in G_matrix_shtriX:
        print(i)    

    print("G_2_matrix_shtriX:")
    for i in G_2_matrix_shtriX:
        print(i)  
    
    print("D2udt2:")
    print(D2udt2)

    print("D3hdt3:")
    print(D3hdt3)

    print("Ksi_eto_2")
    for i in Ksi_eto_2:
        print(i)  

    print("U_bol")
    for i in U_bol:
        print(i)  
    
    print("Bu")
    for i in Bu:
        print(i)  


    print("Xr:")
    for i in Xr:
        print(i)  


def A_shtrix(A,B,shtrix):


    A_matrix_shtrix=np.copy(shtrix)



    A_matrix_shtrix[0][0]=A[0][0]*A[0][0]+A[0][1]*A[1][0]
    A_matrix_shtrix[0][1]=A[0][0]*A[0][1]+A[0][1]*A[1][1]
    A_matrix_shtrix[0][2]=A[0][1]*A[1][2]+A[0][3]
    A_matrix_shtrix[0][3]=A[0][0]*A[0][3]
    A_matrix_shtrix[0][5]=A[0][0]*B[0][0]+A[0][1]*B[1][0]-B[0][0]/0.3
    A_matrix_shtrix[0][6]=A[0][0]*B[0][1]+A[0][1]*B[1][0]-B[0][1]/2
    A_matrix_shtrix[4][0]=A[4][3]*A[2][0]+A[4][1]*(A[1][0]*A[0][0]+A[1][0]*A[0][1]+A[1][2]*A[2][0])
    A_matrix_shtrix[4][1]=A[4][3]*A[2][1]+A[4][1]*(A[1][0]*A[0][1]+A[1][1]*A[1][1]+A[1][2]*A[2][1])
    A_matrix_shtrix[4][2]=A[4][3]*A[2][2]
    A_matrix_shtrix[4][3]=A[4][1]*(A[1][1]*A[1][2]+A[1][1]*A[2][2]+A[1][0]*A[0][3])
    A_matrix_shtrix[4][5]=A[4][1]*(A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0]-B[1][0]/0.3)+A[4][3]*B[2][0]
    A_matrix_shtrix[4][6]=A[4][3]*B[2][1]+A[4][1]*(A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1]-B[1][1]/2)

    return A_matrix_shtrix

def G_shtrix(A,G,G_shtrix):

    G_matrix_shtrix=np.copy(G_shtrix)

    G_matrix_shtrix[0][0]=A[0][0]*G[0][0]+A[0][1]*G[1][0]
    G_matrix_shtrix[0][1]=A[0][0]*G[0][1]+A[0][1]*G[1][1]
    G_matrix_shtrix[4][0]=A[4][3]*G[2][0]+A[4][1]*(A[1][0]*G[0][0]+A[1][1]*G[1][0]+A[1][2]*G[2][0])
    G_matrix_shtrix[4][1]=A[4][3]*G[2][1]+A[4][1]*(A[1][0]*G[0][1]+A[1][1]*G[1][1]+A[1][2]*G[2][1])

    return G_matrix_shtrix

   
def G_2_shtrix(A,G,G_shtrix):

    G_2_matrix_shtrix=np.copy(G_shtrix)


    G_2_matrix_shtrix[4][0]=A[4][1]*G[1][0]
    G_2_matrix_shtrix[4][1]=A[4][1]*G[1][1]

    return G_2_matrix_shtrix


def d2udt2(x,A_matrix_shtriX,G,G_matrix_shtriX,t,w0):

    D2udt2=A_matrix_shtriX[0][0]*x[0]+A_matrix_shtriX[0][1]*x[1]+A_matrix_shtriX[0][2]*x[2]+A_matrix_shtriX[0][3]*x[3]+A_matrix_shtriX[0][5]*x[5]+A_matrix_shtriX[0][6]*x[6]+G_matrix_shtriX[0][0]*-(np.sin(w0*t))+G_matrix_shtriX[0][1]*(np.cos(w0*t)-1)+G[0][0]*w0*-(np.cos(w0*t))+G[0][1]*-(w0)*np.sin(w0*t)

    return D2udt2

def d3hdt3(x,A_matrix_shtriX,G_matrix_shtriX,G_2_matrix_shtriX,t,w0):

    D3hdt3=A_matrix_shtriX[4][0]*x[0]+A_matrix_shtriX[4][1]*x[1]+A_matrix_shtriX[4][2]*x[2]+A_matrix_shtriX[4][3]*x[3]+A_matrix_shtriX[4][5]*x[5]+A_matrix_shtriX[4][6]*x[6]+G_matrix_shtriX[4][0]*-(np.sin(w0*t))+G_matrix_shtriX[4][1]*(np.cos(w0*t)-1)+G_2_matrix_shtriX[4][0]*w0*-(np.cos(w0*t))+G_2_matrix_shtriX[4][1]*-(w0)*np.sin(w0*t)

    return D3hdt3

def A_ksi_eto(ksi_eto,D3hdt3,D2udt2,mat):

    ksi_eto_2=np.copy(mat)
    ksi_eto_2[2]=D3hdt3-(-25.11*ksi_eto[1]-50.11*ksi_eto[2]+12.72*ksi_eto[3]+49.41*ksi_eto[4])-(177.42*ksi_eto[5]+222.92*ksi_eto[6])
    ksi_eto_2[4]=D2udt2-(-1.508*ksi_eto[1]-2.88*ksi_eto[2]+0.69*ksi_eto[3]+1.117*ksi_eto[4])-(10.598*ksi_eto[5]+12.835*ksi_eto[6])

    return ksi_eto_2




   


def multiply_T_x(T, x):
        ksi_eto= []
        ksi_eto=np.dot(T,x)

        return ksi_eto

def multiply_B1_ksi_eto_2(B1, Ksi_eto_2,matrix_U):
        U= np.copy(matrix_U)
        U[0]=B1[0][2]*Ksi_eto_2[2]+B1[0][4]*Ksi_eto_2[4]
        U[1]=B1[1][2]*Ksi_eto_2[2]+B1[1][4]*Ksi_eto_2[4]


        return U

def multiply_B_u(U_bol, matrix_Bu):
        Bu= np.copy(matrix_Bu)
        Bu[5]=U_bol[0]/0.3
        Bu[6]=U_bol[1]/2

        return Bu

def multiply_x0(A,x,Bu,dt,matrix_x0,shtrix,G_shtriX,mat,matrix_U,T,G,w0,B1,matrix_Bu,G_2_matrix_shtriX):
    X=np.copy(matrix_x0)
    Xr=np.copy(x)


    X=np.reshape(X,(-1, 1))
    Xr=np.reshape(Xr,(1, -1))
    x=np.reshape(x,(-1, 1))



    A_matrix_shtriX=np.copy(shtrix)
    G_matrix_shtriX=np.copy(G_shtriX)
    G_2_matrix_shtrix=np.copy(G_2_matrix_shtriX)
    matrix_Bu2= np.copy(matrix_Bu)
    mat2=np.copy(mat)
    matrix_U2= np.copy(matrix_U)

    t=0
    for i in np.arange (0,0.4,dt):
        t=t+dt
        if i==0:
            x=np.reshape(x,(-1, 1))
            X=np.dot((np.dot(A,x)+Bu),dt)+x
            X=np.reshape(X,(1, -1))
            Xr=np.append(Xr,X, axis = 0)
            x=X
        if i >0 :
            
            x=np.reshape(x,(-1, 1))
            ksi_eto= multiply_T_x(T, x)
            print("T:")
            for i in T:
                print(i)   
            D2udt2=d2udt2(x,A_matrix_shtriX,G,G_matrix_shtriX,t,w0)
            G_matrix_shtriX=np.copy(G_shtriX)
            D3hdt3=d3hdt3(x,A_matrix_shtriX,G_matrix_shtriX,G_2_matrix_shtrix,t,w0)
            Ksi_eto_2=A_ksi_eto(ksi_eto,D3hdt3,D2udt2,mat2)
            U_bol=multiply_B1_ksi_eto_2(B1, Ksi_eto_2,matrix_U2)
            Bu=multiply_B_u(U_bol, matrix_Bu2)
            X=np.dot((np.dot(A,x)+Bu),dt)+x
            X=np.reshape(X,(1, -1))
            Xr=np.append(Xr,X, axis = 0)
            x=X
    Xr=np.transpose(Xr)
    return Xr

def grafics(Xr,dt):
    plt.figure(figsize=(8, 8))

    dT=np.array([])
    t=0
    for i in np.arange (0,0.5,dt):
        dT = np.append(dT, [t])
        t=t+dt

    plt.subplot(2, 3, 1)
    plt.plot(dT, Xr[0], color='r')
    plt.title('U')
    plt.grid()

    plt.subplot(2, 3, 2)
    plt.plot(dT, Xr[1], color='r')
    plt.title('W')
    plt.grid()

    plt.subplot(2, 3, 3)
    plt.plot(dT, Xr[2], color='r')
    plt.title('q')
    plt.grid()

    plt.subplot(2, 3, 4)
    plt.plot(dT, Xr[3], color='r')
    plt.title('teta')
    plt.grid()

    plt.subplot(2, 3, 5)
    plt.plot(dT, Xr[4], color='r')
    plt.title('Hight')
    plt.grid()

    plt.subplot(2, 3, 6)
    plt.plot(dT, Xr[5], color='r')
    plt.title('delta_e')
    plt.grid()
    plt.show()


    plt.subplot(2, 3, 1)
    plt.plot(dT, Xr[6], color='r')
    plt.title('delta_T')
    plt.grid()


    plt.show()

def matrix_read(): # считывание матриц
    T=[]
    A=[]
    G=[]
    B=[]
    B1=[]
    shtrix=[]
    G_shtrix_1=[]
    G_shtrix_2=[]
    mat=[]
    matrix_U=[]
    matrix_Bu=[] 
    matrix_x0=[] 
      
    



   
    with open('matrix_T.txt') as m:
        for line in m:
            T.append(list(map(float, line.split())))
    
    with open('matrix_A.txt') as m:
        for line in m:
            A.append(list(map(float, line.split())))

    with open('matrix_G.txt') as m:
        for line in m:
            G.append(list(map(float, line.split())))

    with open('matrix_B.txt') as m:
        for line in m:
            B.append(list(map(float, line.split())))
    
    with open('matrix_B1.txt') as m:
        for line in m:
            B1.append(list(map(float, line.split())))
    
    with open('shtrix.txt') as m:
        for line in m:
            shtrix.append(list(map(float, line.split())))
    
    with open('G_shtrix_1.txt') as m:
        for line in m:
            G_shtrix_1.append(list(map(float, line.split())))
    with open('G_shtrix_2.txt') as m:
        for line in m:
            G_shtrix_2.append(list(map(float, line.split())))
    with open('mat.txt') as m:
        for line in m:
            mat.append(list(map(float, line.split())))
    with open('matrix_U.txt') as m:
        for line in m:
            matrix_U.append(list(map(float, line.split())))
    with open('matrix_Bu.txt') as m:
        for line in m:
            matrix_Bu.append(list(map(float, line.split())))
    with open('matrix_x0.txt') as m:
        for line in m:
            matrix_x0.append(list(map(float, line.split())))
  
    return T,A,G,B,B1,shtrix,G_shtrix_1,G_shtrix_2,mat,matrix_U,matrix_Bu,matrix_x0


if __name__ == "__main__":
    main() 