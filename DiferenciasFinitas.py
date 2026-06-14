import numpy as np

def DiferenciasFinitas(a ,b ,c ,d ,m ,n ,TOL ,MAX_iter):
    if (m < 3 or n < 3):
        print("Error!! m y n deben ser mayor o igual a 3.")
        return
    
    h = (b-a)/n
    k = (d-c)/m

    #Marcar las posiciones en el espacio de los nodos internos
    #Pasos 2 a 3
    x = []
    for i in range(1, n):
        x.append(a + i*h)
    y = []
    for j in range(1, m):
        y.append(c + j*k)
    
    #Inicializar y llenar la matriz de 0
    #Paso 4 a 5
    w = np.zeros((n-1,m-1))

    lambd= (h**2)/k**2
    mu = 2*(1 + lambd)
    l = 1

    while(l <= MAX_iter):
        #Paso 7
        z = (-(h**2)*f(x[0], y[m-2]) + g(a, y[m-2]) + lambd*g(x[0], d) + lambd*w[0][m-3] + w[1][m-2]) / mu
        NORM = abs(z - w[0][m-2])
        w[0][m-2] = z

        # Paso 8
        for i in range(1, n-2):
            z = (-(h**2)*f(x[i], y[m-2]) + lambd*g(x[i], d) + w[i-1][m-2] + w[i+1][m-2] + lambd*w[i][m-3]) / mu
            if abs(w[i][m-2] - z) > NORM:
                NORM = abs(w[i][m-2] - z)
            w[i][m-2] = z

        #Paso 9
        z = (-(h**2)*f(x[n-2], y[m-2]) + g(b, y[m-2]) + lambd*g(x[n-2], d) + w[n-3][m-2] + lambd*w[n-2][m-3]) / mu
        if abs(w[n-2][m-2] - z) > NORM:
            NORM = abs(w[n-2][m-2] - z)
            w[n-2][m-2] = z

        #Paso 10 al 13
        for j in range(m-2, 1, -1):
            z = (-(h**2)*f(x[0], y[j-1]) + g(a, y[j-1]) + lambd*w[0][j] + lambd*w[0][j-2] + w[1][j-1]) / mu
            if abs(w[0][j-1] - z) > NORM:
                NORM = abs(w[0][j-1] - z)
            w[0][j-1] = z

            for i in range(1, n-2):
                z = (-(h**2)*f(x[i], y[j-1]) + w[i-1][j-1] + lambd*w[i][j] + w[i+1][j-1] + lambd*w[i][j-2]) / mu
                if abs(w[i][j-1] - z) > NORM:
                    NORM = abs(w[i][j-1] - z)
                w[i][j-1] = z

            z = (-(h**2)*f(x[n-2], y[j-1]) + g(b, y[j-1]) + w[n-3][j-1] + lambd*w[n-2][j] + lambd*w[n-2][j-2]) / mu
            if abs(w[n-2][j-1] - z) > NORM:
                NORM = abs(w[n-2][j-1] - z)
            w[n-2][j-1] = z

        # Paso 14
        z = (-(h**2)*f(x[0], y[0]) + g(a, y[0]) + lambd*g(x[0], c) + lambd*w[0][1] + w[1][0]) / mu
        if abs(w[0][0] - z) > NORM:
            NORM = abs(w[0][0] - z)
        w[0][0] = z
    
        # Paso 15
        for i in range(1, n-2):
            z = (-(h**2)*f(x[i], y[0]) + lambd*g(x[i], c) + w[i-1][0] + lambd*w[i][1] + w[i+1][0]) / mu
            if abs(w[i][0] - z) > NORM:
                NORM = abs(w[i][0] - z)
            w[i][0] = z

        # Paso 16
        z = (-(h**2)*f(x[n-2], y[0]) + g(b, y[0]) + lambd*g(x[n-2], c) + w[n-3][0] + lambd*w[n-2][1]) / mu
        if abs(w[n-2][0] - z) > NORM:
            NORM = abs(w[n-2][0] - z)
        w[n-2][0] = z

        #Paso 17 al 19
        if NORM <= TOL:
            for i in range(n-1):
                for j in range(m-1):
                    print(f"x={x[i]:.4f}, y={y[j]:.4f}, w={w[i][j]:.6f}")
            return w
        
        l+= 1

    print("Numero maximo de iiteraciones excedido!!")

def f(x, y):
    return -(np.cos(x +y) + np.cos(x -y))

def g(x, y):
    if np.isclose(y, 0):
        return np.cos(x)
    elif np.isclose(x, 0):
        return np.cos(y)
    elif np.isclose(x, np.pi):
        return -np.cos(y)
    elif np.isclose(y, np.pi/2):
        return 0



DiferenciasFinitas(0,np.pi,0,(np.pi)/2,3,4, np.e,5)