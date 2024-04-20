import numpy as np
import sympy as sp
cos=np.cos
sin=np.sin
def Trotx(ang):
    Tx = np.array([[1, 0, 0, 0],
                   [0, cos(ang), -sin(ang), 0],
                   [0, sin(ang), cos(ang), 0],
                   [0, 0, 0, 1]])
    return Tx

def Troty(ang):
    Ty = np.array([[cos(ang), 0, sin(ang), 0],
                   [0, 1, 0, 0],
                   [-sin(ang), 0, cos(ang), 0],
                   [0, 0, 0, 1]])
    return Ty

def Trotz(ang):
    Tz = np.array([[cos(ang), -sin(ang), 0, 0],
                   [sin(ang), cos(ang), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])
    return Tz
pi=np.pi

def Cartesianas(x,y,z):
    #Cilindricas
    p=np.sqrt(x^2 + y^2)
    theta_cil=np.atan2(y,x)
    #Esfericas
    r=np.sqrt(x^2+y^2+z^2)
    theta_esf=np.arccos(z,r)
    phi=np.atan2(y,x)
    print(p,'  ',theta_cil,'\n',r,'  ',theta_esf,'  ',phi,'\n')

def Cilindricas(p,theta_c,z):
    #Cartesianas
    x=p*np.cos(theta_c)
    y = p * np.sin(theta_c)
    #Esfericas
    r=np.sqrt(p^2+z^2)
    theta_e=np.atan2(p,z)
    phi=theta_c
    print(x, '  ', y, '\n', r, '  ', theta_e, '  ', phi, '\n')

def Esfericas(r,theta,phi):
    #Cartesianas
    x=r*np.sin(theta)*np.cos(phi)
    y=r*np.sin(theta)*np.sin(phi)
    z=r*np.cos(theta)

    #Cilindricas
    p=r*np.sin(theta)
    theta_c=phi
    z=r*np.cos(theta)
    print(x, '  ', y, '  ', z, '\n', p, '  ', theta_c, '  ', z, '\n')

def Cuat(eje,ang):
    w=np.cos(ang/2)
    eje=eje/np.linalg.norm(eje)
    eps=eje*np.sin(ang/2)
    return w, eps

def CuatInv(w,eps):
    theta=2*np.arctan2(np.linalg.norm(eps),w)
    u=eps/np.linalg.norm(eps)
    return theta, u

def CuatR(w,eps):
    R=np.array([[2*(w**2+eps[0]**2)-1,2*(eps[0]*eps[1]-w*eps[2]),2*(eps[0]*eps[2]+w*eps[1])],
                [2*(eps[0]*eps[1]+w*eps[2]),2*(w**2+eps[1]**2)-1,2*(eps[1]*eps[2]-w*eps[0])],
                [2*(eps[0]*eps[2]-w*eps[1]),2*(eps[1]*eps[2]+w*eps[0]),2*(w**2+eps[2]**2)-1]])
    return R

def cuaternion(R):
        #cambiar signo aqui si quieres otra solucion
    w=0.5*np.sqrt(1+R[0][0]+R[1][1]+R[2][2])
    if (np.round(w)==0):
        print('W es igual a ')
    ex=1/(4*w)*(R[2][1]-R[1][2])
    ey=1/(4*w)*(R[0][2]-R[2][0])
    ez=1/(4*w)*(R[1][0]-R[0][1])
    eps=np.array([[ex],[ey],[ez]])
    return w,eps

#conjugada con np.conj

def RotX(theta):
    Rot_X=np.array([[1,0,0],
                    [0,np.cos(theta),-np.sin(theta)],
                    [0,np.sin(theta),np.cos(theta)]])
    return Rot_X

def RotY(theta):
    Rot_Y=np.array([[np.cos(theta),0,np.sin(theta)],
                    [0,1,0],
                    [-np.sin(theta),0,np.cos(theta)]])
    return Rot_Y

def RotZ(theta):
    Rot_Z=np.array([[np.cos(theta),-np.sin(theta),0],
                    [np.sin(theta),np.cos(theta),0],
                    [0,0,1]])
    return Rot_Z
def Trasl(x,y,z):
    T = np.array([[1, 0, 0, x],
                  [0, 1, 0, y],
                  [0, 0, 1, z],
                  [0, 0, 0, 1]])
    return T
#INVERSA DE HOMOGENEA ES np.linalg.inv

#--------------SIMBOLICO---------------
def TraslSym(x,y,z):
    T = sp.Matrix([[1, 0, 0, x],
                  [0, 1, 0, y],
                  [0, 0, 1, z],
                  [0, 0, 0, 1]])
    return T

def sTrotx(ang):
    T = sp.Matrix([[1, 0,0,0],
    [0, sp.cos(ang),-sp.sin(ang),0],
    [0, sp.sin(ang), sp.cos(ang),0],
    [0, 0, 0, 1]])
    return T
def sTroty(ang):
    T = sp.Matrix([[sp.cos(ang), 0, sp.sin(ang), 0],
                   [0, 1, 0, 0],
                   [-sp.sin(ang), 0, sp.cos(ang), 0],
                   [0, 0, 0, 1]])
    return T
def sTrotz(ang):
    T = sp.Matrix([[sp.cos(ang),-sp.sin(ang),0,0],
    [sp.sin(ang), sp.cos(ang),0,0],
    [0,0,1,0],
    [0,0,0,1]])
    return T
#--------FIN DE SIMBOLICO---------------

def EulerSym(Var,theta):
    theta=sp.symbols(theta)
    if Var=='X':
        R=sp.Matrix([[1,0,0],
                    [0,sp.cos(theta),-sp.sin(theta)],
                    [0,sp.sin(theta),sp.cos(theta)]])
    elif Var=='Y':
        R = sp.Matrix([[sp.cos(theta), 0, sp.sin(theta)],
                          [0, 1, 0],
                          [-sp.sin(theta), 0, sp.cos(theta)]])
    elif Var=='Z':
        R = sp.Matrix([[sp.cos(theta), -sp.sin(theta), 0],
                              [sp.sin(theta), sp.cos(theta), 0],
                              [0, 0, 1]])
    return R

def RPY(roll,pitch,yaw):
    #euler zyx
    R=RotZ(roll)@RotY(pitch)@RotX(yaw)
    print(EulerSym('Z','roll')@EulerSym('Y','pitch')@EulerSym('X','yaw'))
    return R

def R2RPY(R):                #agregar aqui negativo y todo saldra bien, el signo del pitch sera distinto
    pitch=np.arctan2(-R[2][0],np.sqrt(R[2][1]**2+R[2][2]**2))
    roll=np.arctan2(R[1][0]/(np.cos(pitch)),R[0][0]/np.cos(pitch))
    yaw=np.arctan2(R[2][1]/np.cos(pitch),R[2][2]/np.cos(pitch))
    print('Pitch ',pitch,'\n Yaw ',yaw, '\Roll ',roll)

def Rodrigues(eje,ang):
    I=np.array([[1,0,0],[0,1,0],[0,0,1]])
    us=np.array([[0,-eje[2],eje[1]],[eje[2],0,-eje[0]],[-eje[1],eje[0],0]])
    magnitud=np.linalg.norm(eje)
    us=us/magnitud
    us2=us@us

    R=I+us*np.sin(ang)+us2*(1-np.cos(ang))
    return R

def ejeang(R):
    c = (R[0,0]+R[1,1]+R[2,2]-1.0)/2.0
        #cambiar sino aqui en s para tener la otra rpta
    s = np.sqrt((R[1,0]-R[0,1])**2+(R[2,0]-R[0,2])**2+(R[2,1]-R[1,2])**2)/2.0
    th = np.arctan2(s,c)
    u = 1.0/(2.*np.sin(th))*np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])
    if(th==pi):
        print('Singularidad')
        u=np.array([[np.sqrt(1/2*(R[0][0]+1))],[np.sqrt(1/2*(R[1][1]+1))],[np.sqrt(1/2*(R[2][2]+1))]])
        #si se quiere la otra opción se toman los mismo valores pero se ponen negativos
    return th,u

def dh(d, theta, a, alpha):
    # Escriba aqui la matriz de transformacion homogenea en funcion de los valores de d, theta, a, alpha
    T = np.array([[cos(theta),-cos(alpha)*sin(theta),sin(alpha)*sin(theta),a*cos(theta)],
                  [sin(theta),cos(alpha)*cos(theta),-sin(alpha)*cos(theta),a*sin(theta)],
                  [0,sin(alpha),cos(alpha),d],
                  [0,0,0,1]])
    return T

def dhSym(d, theta, a, alpha):

    # Escriba aqui la matriz de transformacion homogenea en funcion de los valores de d, theta, a, alpha
    T = sp.Matrix([[sp.cos(theta),-sp.cos(alpha)*sp.sin(theta),sp.sin(alpha)*sp.sin(theta),a*sp.cos(theta)],
                  [sp.sin(theta),sp.cos(alpha)*sp.cos(theta),-sp.sin(alpha)*sp.cos(theta),a*sp.sin(theta)],
                  [0,sp.sin(alpha),sp.cos(alpha),d],
                  [0,0,0,1]])
    return T

def sTdh(d, th, a, alpha):
    cth = sp.cos(th); sth = sp.sin(th)
    ca = sp.cos(alpha); sa = sp.sin(alpha)
    Tdh = sp.Matrix([[cth, -ca*sth,  sa*sth, a*cth],
                     [sth,  ca*cth, -sa*cth, a*sth],
                     [0,        sa,     ca,      d],
                     [0,         0,      0,      1]])
    return Tdh

def pasosDH(q):
        print("Eje z0 a lo largo de q1")

        for (i) in range(q-1):
            print(f"Eje z{i+1} en articulacion q{i+2}")

        print("\n")
        for (i) in range(1,q+1):
            print(f"Origen {i} en interseccion de z{i} y z{i - 1}")

        print("\n")

        for (i) in range(1,q):
            print(f"Eje x{i} en dirección de z{i-1} x z{i}")


        print("\n")
        print(f"Efector final: x{q} ortogonal a z{q-1} e intersecarlo")
        print(f"z{q} en dirección de z{q-1} hacia afuera\n")
        return 0

def tablaDH(i):
        print(f"di: distancia de { {i-1} } a [interseccion de z{i-1} con x{i}] en z{i-1}\n")
        print(f"thi: angulo de x{i-1} a x{i} alrededor de z{i-1}\n")
        print(f"ai: distancia de [interseccion de z{i-1} con x{i} a { {i} } en x{i}\n")
        print(f"alphi: ángulo de z{i-1} a z{i} alrededor de x{i}\n")
        return 0

def fkine(q):
    x = np.cos(q[0]) + np.cos(q[0]+q[1]);
    y = np.sin(q[0]) + np.sin(q[0]+q[1]);
    return np.array([x,y])

def JacobianoGeo(tipo,Tinicial,Tfinal,eje):
    if (tipo=='r' or tipo=='R'):
        zi_1 = Tinicial[0:3,eje-1]
        pi_1 = Tinicial[0:3,3]
        pn=Tfinal[0:3,3]

        Jvi=np.cross(zi_1,pn-pi_1)
        Jwi=zi_1

        JGi=np.concatenate((Jvi,Jwi))

    elif (tipo==('p' or 'P')):
        zi_1=Tinicial[0:3,eje-1]
        Jvi=zi_1
        Jwi=np.zeros(3)
        JGi=np.concatenate((Jvi,Jwi))
    return JGi

def PolinomioCubico(t0,tf,q0,qf,dq0,dqf):
    A=np.array([[t0**3,t0**2,t0,1],[tf**3,tf**2,tf,1],[3*t0**2,2*t0,1,0],[3*tf**2,2*tf,1,0]])
    B=np.array([[q0],[qf],[dq0],[dqf]])
    sol=np.linalg.inv(A)@B
    return sol

def PolinomioQuintico(t0,tf,q0,qf,dq0,dqf,ddq0,ddqf):
    A=np.array([[t0**5,t0**4,t0**3,t0**2,t0,1],
                [tf**5,tf**4,tf**3,tf**2,tf,1],
                [5*t0**4,4*t0**3,3*t0**2,2*t0,1,0],
                [5*tf**4,4*tf**3,3*tf**2,2*tf,1,0],
                [20*t0**3,12*t0**2,6*t0,2,0,0],
                [20*tf**3,12*tf**2,6*tf,2,0,0]])
    B=np.array([[q0],[qf],[dq0],[dqf],[ddq0],[ddqf]])
    sol=np.linalg.inv(A)@B
    return sol

#Formula del profe sin probar
def vel_trapezoidal(q0, qf, dqmax, tf, tb, tt):
    q = np.zeros(tt.shape)
    dq = np.zeros(tt.shape)
    ddq = np.zeros(tt.shape)
    for i in range(len(tt)):
        t = tt[i]
        if(t <= tb):
            q[i] = q0 + 1./2.*dqmax/tb*t**2
            dq[i] = dqmax/tb*t
            ddq[i] = dqmax/tb
        elif(t <= tf-tb):
            q[i] = q0 - 1./2.*tb*dqmax+dqmax*t
            dq[i] = dqmax
            ddq[i] = 0.0
        else:
            q[i] = qf - 1./2.*dqmax*(t-tf)**2/tb
            dq[i] = -dqmax/tb*(t-tf)
            ddq[i] = -dqmax/tb
    return q, dq, ddq
